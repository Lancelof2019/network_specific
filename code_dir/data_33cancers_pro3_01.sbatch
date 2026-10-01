#!/bin/bash



#SBATCH --job-name=tcga_multiomics

#SBATCH --account=project_2019755

#SBATCH --partition=small

#SBATCH --time=72:00:00

#SBATCH --nodes=1

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=32

#SBATCH --mem=50G

#SBATCH --array=0-32%13

#SBATCH --output=/users/zhanglia/project3_workspace/code_dir/code_data/slurm_%x_%A_%a.out

#SBATCH --error=/users/zhanglia/project3_workspace/code_dir/code_data/slurm_%x_%A_%a.err

#SBATCH --export=ALL



set -euo pipefail



PROJECT_ROOT="/users/zhanglia/project3_workspace/code_dir"

CODE_DIR="${PROJECT_ROOT}/code_data"

DATA_ROOT="${PROJECT_ROOT}/data_source/linkedomics_tcga_downloads"

SBATCH_ROOT="/scratch/project_2019755"

RESULT_ROOT="${SBATCH_ROOT}/result_data"



R_SCRIPT="${CODE_DIR}/patient_network_analysis.R"



CANCERS=(

TCGA-SKCM

TCGA-SARC

TCGA-PRAD

TCGA-STAD

TCGA-PCPG

TCGA-PAAD

TCGA-OV

TCGA-STES

TCGA-ACC

TCGA-THCA

TCGA-MESO

TCGA-THYM

TCGA-UCEC

TCGA-TGCT

TCGA-UCS

TCGA-ESCA

TCGA-COADREAD

TCGA-CHOL

TCGA-GBM

TCGA-CESC

TCGA-LUSC

TCGA-BLCA

TCGA-BRCA

TCGA-GBMLGG

TCGA-DLBC

TCGA-HNSC

TCGA-LUAD

TCGA-LIHC

TCGA-LGG

TCGA-KICH

TCGA-KIRP

TCGA-KIRC

TCGA-KIPAN

)



CANCER="${CANCERS[$SLURM_ARRAY_TASK_ID]}"



INPUT_DIR="${DATA_ROOT}/${CANCER}/data"

OUTPUT_DIR="${RESULT_ROOT}/${CANCER}"



REQUIRED_FILES=(

rnaseq.cct

methylation.cct

mutation.cbt

scnv.cgt

rppa.cct

clinical.tsi

)



echo "============================================================"

echo "TCGA MULTI-OMICS ANALYSIS"

echo "============================================================"

echo "Job ID            : ${SLURM_JOB_ID}"

echo "Array task ID     : ${SLURM_ARRAY_TASK_ID}"

echo "Cancer            : ${CANCER}"

echo "Node              : ${SLURMD_NODENAME:-$(hostname)}"

echo "Allocated CPUs    : ${SLURM_CPUS_PER_TASK:-unknown}"

echo "Input directory   : ${INPUT_DIR}"

echo "Output directory  : ${OUTPUT_DIR}"

echo "Start time        : $(date -Is)"

echo "============================================================"



if [[ ! -f "${R_SCRIPT}" ]]; then

    echo "ERROR: R script not found: ${R_SCRIPT}" >&2

    exit 1

fi



if [[ ! -d "${INPUT_DIR}" ]]; then

    echo "ERROR: input directory not found: ${INPUT_DIR}" >&2

    exit 1

fi



for FILE_NAME in "${REQUIRED_FILES[@]}"; do

    FILE_PATH="${INPUT_DIR}/${FILE_NAME}"



    if [[ ! -s "${FILE_PATH}" ]]; then

        echo "ERROR: required file missing or empty: ${FILE_PATH}" >&2

        exit 1

    fi



    echo "OK: ${FILE_PATH}"

done



mkdir -p "${OUTPUT_DIR}"



export ANALYSIS_CORES="${SLURM_CPUS_PER_TASK:-1}"



export OMP_NUM_THREADS=1

export MKL_NUM_THREADS=1

export OPENBLAS_NUM_THREADS=1

export BLIS_NUM_THREADS=1

export VECLIB_MAXIMUM_THREADS=1

export NUMEXPR_NUM_THREADS=1



module purge

module load r-env



if ! command -v Rscript >/dev/null 2>&1; then

    echo "ERROR: Rscript is not available." >&2

    exit 1

fi



echo "============================================================"

echo "RESOURCE INFORMATION"

echo "============================================================"

echo "Cancer          : ${CANCER}"

echo "Node            : ${SLURMD_NODENAME:-$(hostname)}"

echo "Analysis cores  : ${ANALYSIS_CORES}"

echo "Visible CPUs    : $(nproc)"

echo "Memory status:"

free -h || true

echo "============================================================"



srun --ntasks=1 --cpus-per-task="${SLURM_CPUS_PER_TASK}" Rscript --no-save "${R_SCRIPT}" "${CANCER}"



echo "============================================================"

echo "FINISHED SUCCESSFULLY"

echo "Cancer       : ${CANCER}"

echo "Finished     : $(date -Is)"

echo "Results      : ${OUTPUT_DIR}"

echo "============================================================"
