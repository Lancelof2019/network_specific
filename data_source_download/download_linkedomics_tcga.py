#!/usr/bin/env python3
"""
Download selected TCGA gene-level datasets from LinkedOmics.

Default behavior:
- Start from https://www.linkedomics.org/login.php
- Keep rows whose Cohort Source is TCGA and Permission is Y
- Visit each row's Data Download page
- Download datasets matching:
  Clinical, Methylation gene-level, miRNA gene-level, Mutation gene-level,
  RPPA gene-level, and SCNV gene-level Thresholded
- If a cancer type does not have one selected data class, log it and continue
- SCNV is strict: only Gene level + Thresholded is selected
- Skip missing datasets silently, but record them in the manifest

Run examples:
  python3 download_linkedomics_tcga.py --dry-run
  python3 download_linkedomics_tcga.py --out linked_tcga
  python3 download_linkedomics_tcga.py --out linked_tcga --workers 4
  python3 download_linkedomics_tcga.py --out linked_tcga --only BRCA STAD ACC
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import unquote, urljoin, urlparse
from urllib.request import Request, urlopen


BASE_URL = "https://www.linkedomics.org/login.php"
EXPECTED_DATASET_REASONS = [
    "Clinical",
    "Methylation gene-level",
    "miRNA gene-level",
    "Mutation gene-level",
    "RPPA gene-level",
    "SCNV gene-level thresholded",
]
USER_AGENT = (
    "Mozilla/5.0 (X11; Linux x86_64) "
    "AppleWebKit/537.36 (KHTML, like Gecko) "
    "Chrome/120.0 Safari/537.36"
)


@dataclass
class CancerCohort:
    cancer_type: str
    cohort_source: str
    cancer_id: str
    samples: str
    permission: str
    download_page: str


@dataclass
class DatasetFile:
    cancer_type: str
    cancer_id: str
    omics_dataset: str
    description: str
    level: str
    unit: str
    samples: str
    attributes: str
    file_type: str
    url: str
    selected_reason: str
    local_path: str = ""
    status: str = "planned"


@dataclass
class CancerResult:
    order: int
    cancer_id: str
    cancer_type: str
    status: str
    selected_count: int
    downloaded_count: int
    existing_count: int
    failed_count: int
    manifest_path: str
    log_path: str
    missing_expected: list[str] | None = None
    error: str = ""
    files: list[DatasetFile] | None = None


@dataclass
class DatasetDiscovery:
    files: list[DatasetFile]
    missing_expected: list[str]
    skipped_rows: list[str]


class TableParser(HTMLParser):
    """Small HTML table parser with link capture; avoids third-party packages."""

    def __init__(self, page_url: str):
        super().__init__()
        self.page_url = page_url
        self.rows: list[list[dict[str, object]]] = []
        self._row: list[dict[str, object]] | None = None
        self._cell: dict[str, object] | None = None

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attrs_dict = {k.lower(): v for k, v in attrs}
        if tag.lower() == "tr":
            self._row = []
        elif tag.lower() in {"td", "th"} and self._row is not None:
            self._cell = {"text": "", "links": []}
        elif tag.lower() == "a" and self._cell is not None:
            href = attrs_dict.get("href")
            if href:
                links = self._cell["links"]
                assert isinstance(links, list)
                links.append(urljoin(self.page_url, href))

    def handle_data(self, data: str) -> None:
        if self._cell is not None:
            self._cell["text"] = str(self._cell["text"]) + data

    def handle_endtag(self, tag: str) -> None:
        tag = tag.lower()
        if tag in {"td", "th"} and self._row is not None and self._cell is not None:
            self._cell["text"] = clean_text(str(self._cell["text"]))
            self._row.append(self._cell)
            self._cell = None
        elif tag == "tr" and self._row is not None:
            if self._row:
                self.rows.append(self._row)
            self._row = None


def clean_text(value: str) -> str:
    return re.sub(r"\s+", " ", value.replace("\xa0", " ")).strip()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def write_log(path: Path, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(f"{utc_now()} {message}\n")


def fetch_text(url: str, timeout: int = 60, retries: int = 3, delay: float = 2.0) -> str:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            req = Request(url, headers={"User-Agent": USER_AGENT})
            with urlopen(req, timeout=timeout) as response:
                raw = response.read()
            return raw.decode("utf-8", errors="replace")
        except (HTTPError, URLError, TimeoutError) as exc:
            last_error = exc
            if attempt < retries:
                time.sleep(delay * attempt)
    raise RuntimeError(f"Failed to fetch {url}: {last_error}") from last_error


def parse_tables(url: str, timeout: int = 60, retries: int = 3, delay: float = 2.0) -> list[list[dict[str, object]]]:
    html = fetch_text(url, timeout=timeout, retries=retries, delay=delay)
    parser = TableParser(url)
    parser.feed(html)
    return parser.rows


def cell_text(row: list[dict[str, object]], idx: int) -> str:
    if idx >= len(row):
        return ""
    return str(row[idx].get("text", ""))


def cell_links(row: list[dict[str, object]], idx: int) -> list[str]:
    if idx >= len(row):
        return []
    links = row[idx].get("links", [])
    return list(links) if isinstance(links, list) else []


def header_index_map(header_row: list[dict[str, object]]) -> dict[str, int]:
    mapping: dict[str, int] = {}
    for i, cell in enumerate(header_row):
        key = clean_text(str(cell.get("text", ""))).lower()
        key = key.replace("atrributes", "attributes")
        mapping[key] = i
    return mapping


def find_table_with_headers(rows: list[list[dict[str, object]]], required: set[str]) -> tuple[dict[str, int], list[list[dict[str, object]]]]:
    for i, row in enumerate(rows):
        mapping = header_index_map(row)
        if required.issubset(mapping.keys()):
            return mapping, rows[i + 1 :]
    raise RuntimeError(f"Could not find table with headers: {sorted(required)}")


def discover_tcga_cohorts(timeout: int = 60, retries: int = 3, delay: float = 2.0) -> list[CancerCohort]:
    rows = parse_tables(BASE_URL, timeout=timeout, retries=retries, delay=delay)
    required = {"cancer type", "cohort source", "cancer id", "samples", "permission", "data download"}
    idx, data_rows = find_table_with_headers(rows, required)

    cohorts: list[CancerCohort] = []
    for row in data_rows:
        source = cell_text(row, idx["cohort source"])
        permission = cell_text(row, idx["permission"])
        download_links = cell_links(row, idx["data download"])
        if source.upper() != "TCGA":
            continue
        if permission.upper() != "Y":
            continue
        if not download_links:
            continue
        cohorts.append(
            CancerCohort(
                cancer_type=cell_text(row, idx["cancer type"]),
                cohort_source=source,
                cancer_id=cell_text(row, idx["cancer id"]),
                samples=cell_text(row, idx["samples"]),
                permission=permission,
                download_page=download_links[0],
            )
        )
    return cohorts


def reason_to_select(omics_dataset: str, level: str) -> str | None:
    name = omics_dataset.lower()
    level_l = level.lower()

    if "clinical" in name or "clinical" in level_l:
        return "Clinical"

    if "methylation" in name:
        if "gene level" in name or level_l == "gene":
            return "Methylation gene-level"
        return None

    if "mirna" in name:
        if "gene level" in name or level_l in {"mirgene", "gene"}:
            return "miRNA gene-level"
        return None

    if "mutation" in name:
        if "gene level" in name or level_l == "gene":
            return "Mutation gene-level"
        return None

    if "rppa" in name:
        if "gene level" in name or level_l == "gene":
            return "RPPA gene-level"
        return None

    if "scnv" in name:
        # The requested SCNV choice is gene-level + thresholded, not focal and not log-ratio.
        if ("gene level" in name or level_l == "gene") and "threshold" in name:
            return "SCNV gene-level thresholded"
        return None

    return None


def target_family(omics_dataset: str, level: str) -> str | None:
    name = omics_dataset.lower()
    level_l = level.lower()
    if "clinical" in name or "clinical" in level_l:
        return "Clinical"
    if "methylation" in name:
        return "Methylation"
    if "mirna" in name:
        return "miRNA"
    if "mutation" in name:
        return "Mutation"
    if "rppa" in name:
        return "RPPA"
    if "scnv" in name:
        return "SCNV"
    return None


def skip_reason(omics_dataset: str, level: str) -> str:
    name = omics_dataset.lower()
    level_l = level.lower()
    if "scnv" in name:
        if not ("gene level" in name or level_l == "gene"):
            return "SCNV is not gene-level"
        if "threshold" not in name:
            return "SCNV gene-level exists but is not Thresholded"
    if any(key in name for key in ["methylation", "mirna", "mutation", "rppa"]):
        if "gene level" not in name and level_l not in {"gene", "mirgene"}:
            return "target family exists but row is not gene-level"
    return "not part of requested dataset rules"


def discover_dataset_files_with_report(
    cohort: CancerCohort,
    timeout: int = 60,
    retries: int = 3,
    delay: float = 2.0,
) -> DatasetDiscovery:
    rows = parse_tables(cohort.download_page, timeout=timeout, retries=retries, delay=delay)
    required = {"omics dataset", "link", "description", "level", "unit", "samples", "attributes", "file type"}
    idx, data_rows = find_table_with_headers(rows, required)

    files: list[DatasetFile] = []
    skipped_rows: list[str] = []
    for row in data_rows:
        omics_dataset = cell_text(row, idx["omics dataset"])
        if not omics_dataset:
            continue
        level = cell_text(row, idx["level"])
        reason = reason_to_select(omics_dataset, level)
        if not reason:
            if target_family(omics_dataset, level):
                skipped_rows.append(f"{omics_dataset} | level={level} | {skip_reason(omics_dataset, level)}")
            continue
        links = cell_links(row, idx["link"])
        if not links:
            skipped_rows.append(f"{omics_dataset} | level={level} | selected but missing download link")
            continue
        files.append(
            DatasetFile(
                cancer_type=cohort.cancer_type,
                cancer_id=cohort.cancer_id,
                omics_dataset=omics_dataset,
                description=cell_text(row, idx["description"]),
                level=cell_text(row, idx["level"]),
                unit=cell_text(row, idx["unit"]),
                samples=cell_text(row, idx["samples"]),
                attributes=cell_text(row, idx["attributes"]),
                file_type=cell_text(row, idx["file type"]),
                url=links[0],
                selected_reason=reason,
            )
        )
    present = {item.selected_reason for item in files}
    missing_expected = [reason for reason in EXPECTED_DATASET_REASONS if reason not in present]
    return DatasetDiscovery(files=files, missing_expected=missing_expected, skipped_rows=skipped_rows)


def discover_dataset_files(
    cohort: CancerCohort,
    timeout: int = 60,
    retries: int = 3,
    delay: float = 2.0,
) -> list[DatasetFile]:
    return discover_dataset_files_with_report(cohort, timeout=timeout, retries=retries, delay=delay).files


def safe_name(value: str) -> str:
    value = value.strip().replace(" ", "_")
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    value = re.sub(r"_+", "_", value).strip("._")
    return value or "unknown"


def cancer_package_dir(out_dir: Path, cancer_id: str) -> Path:
    return out_dir / f"TCGA-{safe_name(cancer_id)}"


def filename_from_url(url: str) -> str:
    name = unquote(Path(urlparse(url).path).name)
    return safe_name(name)


def set_local_paths(files: list[DatasetFile], out_dir: Path) -> None:
    for item in files:
        cancer_dir = cancer_package_dir(out_dir, item.cancer_id) / "data"
        item.local_path = str(cancer_dir / filename_from_url(item.url))


def download_file(url: str, output_path: Path, overwrite: bool, timeout: int, retries: int, delay: float) -> str:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists() and output_path.stat().st_size > 0 and not overwrite:
        return "exists"

    part_path = output_path.with_suffix(output_path.suffix + ".part")
    last_error: Exception | None = None

    for attempt in range(1, retries + 1):
        try:
            req = Request(url, headers={"User-Agent": USER_AGENT})
            with urlopen(req, timeout=timeout) as response, part_path.open("wb") as handle:
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    handle.write(chunk)
            os.replace(part_path, output_path)
            return "downloaded"
        except (HTTPError, URLError, TimeoutError, OSError) as exc:
            last_error = exc
            if part_path.exists():
                part_path.unlink(missing_ok=True)
            if attempt < retries:
                time.sleep(delay * attempt)

    raise RuntimeError(f"Failed to download {url}: {last_error}") from last_error


def write_manifest(files: list[DatasetFile], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(files[0]).keys()) if files else [
            "cancer_type",
            "cancer_id",
            "omics_dataset",
            "description",
            "level",
            "unit",
            "samples",
            "attributes",
            "file_type",
            "url",
            "selected_reason",
            "local_path",
            "status",
        ])
        writer.writeheader()
        for item in files:
            writer.writerow(asdict(item))


def write_json_manifest(files: list[DatasetFile], cohorts: list[CancerCohort], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "source": BASE_URL,
        "cohort_count": len(cohorts),
        "file_count": len(files),
        "selection_rules": ["Cohort Source == TCGA", "Permission == Y"] + EXPECTED_DATASET_REASONS,
        "cohorts": [asdict(c) for c in cohorts],
        "files": [asdict(f) for f in files],
    }
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def read_manifest(path: Path) -> list[DatasetFile]:
    if not path.exists():
        return []
    fields = DatasetFile.__dataclass_fields__.keys()
    with path.open("r", newline="", encoding="utf-8") as handle:
        return [DatasetFile(**{field: row.get(field, "") for field in fields}) for row in csv.DictReader(handle)]


def write_status(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def process_cancer_worker(task: dict[str, object]) -> CancerResult:
    order = int(task["order"])
    cohort = CancerCohort(**task["cohort"])  # type: ignore[arg-type]
    options = task["options"]  # type: ignore[assignment]
    assert isinstance(options, dict)

    out_dir = Path(str(options["out_dir"]))
    dry_run = bool(options["dry_run"])
    overwrite = bool(options["overwrite"])
    redo_done = bool(options["redo_done"])
    timeout = int(options["timeout"])
    retries = int(options["retries"])
    delay = float(options["delay"])

    package_dir = cancer_package_dir(out_dir, cohort.cancer_id)
    package_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = package_dir / "manifest.csv"
    status_path = package_dir / "status.json"
    log_path = package_dir / "cancer.log"
    started_at = utc_now()

    try:
        if status_path.exists() and not overwrite and not redo_done:
            previous = json.loads(status_path.read_text(encoding="utf-8"))
            previous_status = str(previous.get("status", ""))
            skip_states = {"complete"} if not dry_run else {"complete", "dry_run"}
            if previous_status in skip_states:
                files = read_manifest(manifest_path)
                write_log(log_path, f"SKIP already processed status={previous_status} files={len(files)}")
                return CancerResult(
                    order=order,
                    cancer_id=cohort.cancer_id,
                    cancer_type=cohort.cancer_type,
                    status="skipped_done",
                    selected_count=len(files),
                    downloaded_count=sum(1 for item in files if item.status == "downloaded"),
                    existing_count=sum(1 for item in files if item.status == "exists"),
                    failed_count=sum(1 for item in files if item.status.startswith("failed")),
                    manifest_path=str(manifest_path),
                    log_path=str(log_path),
                    missing_expected=list(previous.get("missing_expected", [])),
                    files=files,
                )

        write_log(log_path, f"START cancer_id={cohort.cancer_id} cancer_type={cohort.cancer_type}")
        discovery = discover_dataset_files_with_report(cohort, timeout=timeout, retries=retries, delay=delay)
        files = discovery.files
        set_local_paths(files, out_dir)
        write_manifest(files, manifest_path)
        write_log(log_path, f"SELECTED files={len(files)}")
        for skipped in discovery.skipped_rows:
            write_log(log_path, f"SKIP_NOT_MATCHED {skipped}")
        for missing in discovery.missing_expected:
            write_log(log_path, f"MISSING_EXPECTED {missing}")

        if dry_run:
            finished_at = utc_now()
            write_status(
                status_path,
                {
                    "status": "dry_run",
                    "cancer_id": cohort.cancer_id,
                    "cancer_type": cohort.cancer_type,
                    "started_at": started_at,
                    "finished_at": finished_at,
                    "selected_count": len(files),
                    "downloaded_count": 0,
                    "existing_count": 0,
                    "failed_count": 0,
                    "missing_expected": discovery.missing_expected,
                    "manifest_path": str(manifest_path),
                    "log_path": str(log_path),
                },
            )
            write_log(log_path, "DONE dry_run")
            return CancerResult(
                order=order,
                cancer_id=cohort.cancer_id,
                cancer_type=cohort.cancer_type,
                status="dry_run",
                selected_count=len(files),
                downloaded_count=0,
                existing_count=0,
                failed_count=0,
                manifest_path=str(manifest_path),
                log_path=str(log_path),
                missing_expected=discovery.missing_expected,
                files=files,
            )

        for i, item in enumerate(files, start=1):
            write_log(log_path, f"DOWNLOAD_START {i}/{len(files)} dataset={item.omics_dataset} url={item.url}")
            try:
                item.status = download_file(
                    item.url,
                    Path(item.local_path),
                    overwrite=overwrite,
                    timeout=timeout,
                    retries=retries,
                    delay=delay,
                )
                write_log(log_path, f"DOWNLOAD_DONE {i}/{len(files)} status={item.status} path={item.local_path}")
            except Exception as exc:
                item.status = f"failed: {exc}"
                write_log(log_path, f"DOWNLOAD_FAILED {i}/{len(files)} error={exc}")
            write_manifest(files, manifest_path)
            time.sleep(delay)

        downloaded_count = sum(1 for item in files if item.status == "downloaded")
        existing_count = sum(1 for item in files if item.status == "exists")
        failed_count = sum(1 for item in files if item.status.startswith("failed"))
        status = "complete" if failed_count == 0 else "partial_failed"
        finished_at = utc_now()
        write_status(
            status_path,
            {
                "status": status,
                "cancer_id": cohort.cancer_id,
                "cancer_type": cohort.cancer_type,
                "started_at": started_at,
                "finished_at": finished_at,
                "selected_count": len(files),
                "downloaded_count": downloaded_count,
                "existing_count": existing_count,
                "failed_count": failed_count,
                "missing_expected": discovery.missing_expected,
                "manifest_path": str(manifest_path),
                "log_path": str(log_path),
            },
        )
        write_log(
            log_path,
            f"DONE status={status} selected={len(files)} downloaded={downloaded_count} existing={existing_count} failed={failed_count}",
        )
        return CancerResult(
            order=order,
            cancer_id=cohort.cancer_id,
            cancer_type=cohort.cancer_type,
            status=status,
            selected_count=len(files),
            downloaded_count=downloaded_count,
            existing_count=existing_count,
            failed_count=failed_count,
            manifest_path=str(manifest_path),
            log_path=str(log_path),
            missing_expected=discovery.missing_expected,
            files=files,
        )

    except Exception as exc:
        write_log(log_path, f"ERROR {exc}")
        write_status(
            status_path,
            {
                "status": "failed",
                "cancer_id": cohort.cancer_id,
                "cancer_type": cohort.cancer_type,
                "started_at": started_at,
                "finished_at": utc_now(),
                "selected_count": 0,
                "downloaded_count": 0,
                "existing_count": 0,
                "failed_count": 1,
                "missing_expected": [],
                "manifest_path": str(manifest_path),
                "log_path": str(log_path),
                "error": str(exc),
            },
        )
        return CancerResult(
            order=order,
            cancer_id=cohort.cancer_id,
            cancer_type=cohort.cancer_type,
            status="failed",
            selected_count=0,
            downloaded_count=0,
            existing_count=0,
            failed_count=1,
            manifest_path=str(manifest_path),
            log_path=str(log_path),
            missing_expected=[],
            error=str(exc),
            files=[],
        )


def resolve_worker_count(requested_workers: int, cohort_count: int) -> int:
    if cohort_count <= 0:
        return 1
    cpu_count = os.cpu_count() or 1
    if requested_workers == 0:
        return min(cpu_count, cohort_count)
    return max(1, min(requested_workers, cohort_count))


def main() -> int:
    parser = argparse.ArgumentParser(description="Download selected TCGA gene-level datasets from LinkedOmics.")
    parser.add_argument("--out", default="linkedomics_tcga_downloads", help="Output directory.")
    parser.add_argument("--dry-run", action="store_true", help="Only build manifests; do not download files.")
    parser.add_argument("--only", nargs="*", help="Optional cancer IDs to include, e.g. BRCA STAD ACC.")
    parser.add_argument("--max-cancers", type=int, default=None, help="Limit number of cancer cohorts for testing.")
    parser.add_argument(
        "--workers",
        type=int,
        default=min(4, os.cpu_count() or 1),
        help="Parallel cancer workers. Use 0 to use all CPUs. Default: min(4, CPU count).",
    )
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing files.")
    parser.add_argument("--redo-done", action="store_true", help="Reprocess cancers that already have status.json marked complete.")
    parser.add_argument("--delay", type=float, default=0.5, help="Delay between page/file requests in seconds.")
    parser.add_argument("--timeout", type=int, default=120, help="Network timeout in seconds.")
    parser.add_argument("--retries", type=int, default=3, help="Retries per page/file.")
    args = parser.parse_args()

    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    run_log = out_dir / "run.log"
    csv_manifest = out_dir / "download_manifest.csv"
    json_manifest = out_dir / "download_manifest.json"

    write_log(run_log, f"RUN_START dry_run={args.dry_run} out={out_dir}")
    print("Discovering TCGA cohorts...", file=sys.stderr)
    cohorts = discover_tcga_cohorts(timeout=args.timeout, retries=args.retries, delay=args.delay)
    if args.only:
        wanted = {x.upper() for x in args.only}
        cohorts = [c for c in cohorts if c.cancer_id.upper() in wanted]
    if args.max_cancers is not None:
        cohorts = cohorts[: args.max_cancers]

    worker_count = resolve_worker_count(args.workers, len(cohorts))
    print(f"Found {len(cohorts)} TCGA cohorts. Workers: {worker_count}", file=sys.stderr)
    write_log(run_log, f"COHORTS count={len(cohorts)} workers={worker_count}")

    options = {
        "out_dir": str(out_dir),
        "dry_run": args.dry_run,
        "overwrite": args.overwrite,
        "redo_done": args.redo_done,
        "timeout": args.timeout,
        "retries": args.retries,
        "delay": args.delay,
    }
    tasks = [{"order": i, "cohort": asdict(cohort), "options": options} for i, cohort in enumerate(cohorts)]

    results: list[CancerResult] = []
    if worker_count == 1:
        for i, task in enumerate(tasks, start=1):
            result = process_cancer_worker(task)
            results.append(result)
            print(
                f"[{i}/{len(tasks)}] {result.cancer_id} status={result.status} files={result.selected_count} failed={result.failed_count}",
                file=sys.stderr,
            )
    else:
        with ProcessPoolExecutor(max_workers=worker_count) as executor:
            future_to_cancer = {executor.submit(process_cancer_worker, task): task["cohort"] for task in tasks}
            for i, future in enumerate(as_completed(future_to_cancer), start=1):
                result = future.result()
                results.append(result)
                print(
                    f"[{i}/{len(tasks)}] {result.cancer_id} status={result.status} files={result.selected_count} failed={result.failed_count}",
                    file=sys.stderr,
                )

    results.sort(key=lambda item: item.order)
    files: list[DatasetFile] = []
    for result in results:
        if result.files:
            files.extend(result.files)

    write_manifest(files, csv_manifest)
    write_json_manifest(files, cohorts, json_manifest)
    summary_path = out_dir / "run_summary.json"
    summary = {
        "status": "dry_run" if args.dry_run else "complete",
        "started_from": BASE_URL,
        "out_dir": str(out_dir),
        "worker_count": worker_count,
        "cohort_count": len(cohorts),
        "file_count": len(files),
        "downloaded_count": sum(1 for item in files if item.status == "downloaded"),
        "existing_count": sum(1 for item in files if item.status == "exists"),
        "failed_count": sum(1 for item in files if item.status.startswith("failed")),
        "missing_expected_by_cancer": {
            result.cancer_id: result.missing_expected or []
            for result in results
            if result.missing_expected
        },
        "results": [asdict(result) | {"files": None} for result in results],
        "manifest_csv": str(csv_manifest),
        "manifest_json": str(json_manifest),
        "run_log": str(run_log),
    }
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    write_log(
        run_log,
        f"RUN_DONE cohorts={len(cohorts)} files={len(files)} failed={summary['failed_count']} summary={summary_path}",
    )

    print(f"Selected {len(files)} files.", file=sys.stderr)
    print(f"Manifest CSV: {csv_manifest}", file=sys.stderr)
    print(f"Manifest JSON: {json_manifest}", file=sys.stderr)
    print(f"Run summary: {summary_path}", file=sys.stderr)
    print("Dry run complete. No files downloaded." if args.dry_run else "Download finished.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
