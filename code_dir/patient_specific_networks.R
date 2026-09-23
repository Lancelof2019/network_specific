setwd("/users/zhanglia/project3_workspace/code_dir/code_data")
# read datasets
rnaseq <- as.matrix(read.table("../raw_data/rnaseq.cct", header = TRUE, row.names = 1))
methylation <- as.matrix(read.table("../raw_data/methylation.cct", header = TRUE, row.names = 1))
mutation <- as.matrix(read.table("../raw_data/mutation.cbt", header = TRUE, row.names = 1))
scnv <- as.matrix(read.table("../raw_data/scnv.cgt", header = TRUE, row.names = 1))


# normalize rnaseq data gene-wise
rnaseq_z <- t(scale(t(rnaseq)))

# select genes that have top 100 variation in rnaseq data for heatmap visualization
vars <- apply(rnaseq_z, 1, var)
rnaseq_top100 <- rnaseq_z[order(vars, decreasing = TRUE)[1:100], ]
genes_to_select <- rownames(rnaseq_top100)

# extract these genes from other datsets
methylation_top100 <- methylation[rownames(methylation) %in% genes_to_select, , drop = FALSE]
mutation_top100 <- mutation[rownames(mutation) %in% genes_to_select, , drop = FALSE]
scnv_top100 <- scnv[rownames(scnv) %in% genes_to_select, , drop = FALSE]

# align column names in all datasets
sample_names <- colnames(rnaseq_top100)

align_cols <- function(mat, sample_names) {
  # if some columns are missing, fill them with NA
  missing <- setdiff(sample_names, colnames(mat))
  if (length(missing) > 0) {
    na_mat <- matrix(
      NA,
      nrow = nrow(mat),
      ncol = length(missing),
      dimnames = list(rownames(mat), missing)
    )
    mat <- cbind(mat, na_mat)
  }
  
  mat[, sample_names, drop = FALSE]
}

methylation_top100 <- align_cols(methylation_top100, sample_names)
mutation_top100  <- align_cols(mutation_top100,  sample_names)
scnv_top100  <- align_cols(scnv_top100,  sample_names)

# combine datasets based on patient ID
combined_omics_top100 <- rbind(
  rnaseq_top100,
  methylation_top100,
  mutation_top100,
  scnv_top100
)

# cluster patients based on combined omics data
clustered_col_order <- hclust(dist(t(combined_omics_top100)))


# create heatmaps and legends
# cluster genes in each set separately, and apply previous column order
library(circlize)
library(ComplexHeatmap)

rna_upper_clip <- quantile(rnaseq_top100, 0.95, na.rm = TRUE)
rna_lower_clip <- quantile(rnaseq_top100, 0.05, na.rm = TRUE)
rna_max <- max(rnaseq_top100, na.rm = TRUE)

rna_col_fun <- colorRamp2(
  c(rna_lower_clip, 0, rna_upper_clip),
  c("blue", "white", "red")
)


rnaseq_hm <- Heatmap(rnaseq_top100, name = "RNAseq", show_column_names = FALSE, col=rna_col_fun,
                     show_column_dend = TRUE, show_row_names = FALSE,
                     row_title_gp = gpar(fontsize = 15), show_row_dend = FALSE,
                     row_title = "RNA sequencing", row_title_side = "left",
                     show_heatmap_legend = FALSE,
                     cluster_rows = TRUE, cluster_columns = FALSE)

rna_legend <- Legend(
  title = "RNA seq",
  col_fun = rna_col_fun,
  at = c(rna_lower_clip, 0, rna_upper_clip, rna_max),
  labels = c(
    paste0(sprintf("%.1f", rna_lower_clip), " = 5%"), 
    0, 
    paste0(sprintf("%.1f", rna_upper_clip), " = 95%"), 
    paste0(sprintf("%.1f", rna_max), " = 5%")),
  legend_height = unit(7, "cm"),
  labels_gp = gpar(fontsize = 15),
  title_gp = gpar(fontsize = 20)
)


methylation_hm <- Heatmap(methylation_top100, name = "Methylation", 
                          na_col = "black", cluster_rows = TRUE, cluster_columns = FALSE,
                          row_title_gp = gpar(fontsize = 15), show_row_dend = FALSE,
                          row_title = "Methylation", row_title_side = "left",
                          show_heatmap_legend = FALSE,
                          show_column_names = FALSE, show_row_names = FALSE)

meth_legend <- Legend(
  at = c(-0.5, 0, 0.5),
  col_fun = colorRamp2(c(-0.5,0,0.5), c("blue", "white", "red")),
  title = "Methylation",
  legend_height = unit(7, "cm"),
  labels_gp = gpar(fontsize = 15),
  title_gp = gpar(fontsize = 20)
)

colors <- c("0" = "darkolivegreen2", "1" = "indianred1")

mutation_hm <- Heatmap(mutation_top100, name = "Mutation", col = colors,
                       na_col = "black", cluster_rows = TRUE, cluster_columns = FALSE,
                       row_title_gp = gpar(fontsize = 15), show_row_dend = FALSE,
                       show_row_names = FALSE, row_title = "Mutation", 
                       show_heatmap_legend = FALSE,
                       row_title_side = "left", show_column_names = FALSE)

mutation_legend <- Legend(
  at = c(0, 1),
  legend_gp = gpar(
    fill = c("#4DAF4A", "#E41A1C")
  ),
  title = "Mutation",
  legend_height = unit(7, "cm"),
  labels_gp = gpar(fontsize = 15),
  title_gp = gpar(fontsize = 20)
)

scnv_hm <- Heatmap(scnv_top100, name = "SCNV",
                   na_col = "black", cluster_rows = TRUE, cluster_columns = FALSE,
                   row_title_gp = gpar(fontsize = 15), show_row_dend = FALSE,
                   show_row_names = FALSE, row_title = "SCNV",
                   show_heatmap_legend = FALSE,
                   row_title_side = "left", show_column_names = FALSE)


scnv_legend <- Legend(
  at=c(-2, -1, 0, 1, 2),
  title = "SCNV",
  col_fun = colorRamp2(c(-2,0,2), c("blue", "white", "red")),
  legend_height = unit(7, "cm"),
  labels_gp = gpar(fontsize = 15),
  title_gp = gpar(fontsize = 20)
)

na_legend = Legend(
  labels = "NA",
  legend_gp = gpar(fill = "black"),
  title = "Missing",
  labels_gp = gpar(fontsize = 15),
  title_gp = gpar(fontsize = 20)
)


# create a heatmap list
heatmap_list = rnaseq_hm %v% methylation_hm %v% mutation_hm %v% scnv_hm


# draw a heatmap list with legends to pdf file
pdf("../figures/omics_heatmaps.pdf", width = 15, height = 15)
draw(heatmap_list, annotation_legend_list = list(rna_legend, meth_legend, 
                                                 mutation_legend, scnv_legend, 
                                                 na_legend), 
     column_title = "Samples", row_title = "Genes", 
     column_title_gp = gpar(fontsize = 30), row_title_gp = gpar(fontsize=30))
dev.off()

############################################################
# 0. WORKING DIRECTORY
############################################################

setwd("/users/zhanglia/project3_workspace/code_dir/code_data")

options(timeout = 600)


############################################################
# 1. READ DATASETS
############################################################

rnaseq <- as.matrix(
  read.table(
    "../raw_data/rnaseq.cct",
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
)

methylation <- as.matrix(
  read.table(
    "../raw_data/methylation.cct",
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
)

mutation <- as.matrix(
  read.table(
    "../raw_data/mutation.cbt",
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
)

scnv <- as.matrix(
  read.table(
    "../raw_data/scnv.cgt",
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
)


cat("RNA dimensions:", dim(rnaseq), "\n")
cat("Methylation dimensions:", dim(methylation), "\n")
cat("Mutation dimensions:", dim(mutation), "\n")
cat("SCNV dimensions:", dim(scnv), "\n")


############################################################
# 2. NORMALIZE RNA-SEQ GENE-WISE
############################################################

rnaseq_z <- t(
  scale(
    t(rnaseq)
  )
)


############################################################
# 3. EXTRACT ALL GENE IDS FROM ALL OMICS DATASETS
############################################################

gene_list <- union(
  row.names(methylation),
  union(
    row.names(mutation),
    union(
      row.names(scnv),
      row.names(rnaseq_z)
    )
  )
)

gene_list <- sort(gene_list)


patient_ids <- union(
  colnames(methylation),
  union(
    colnames(mutation),
    union(
      colnames(scnv),
      colnames(rnaseq_z)
    )
  )
)


cat(
  "Total unique genes:",
  length(gene_list),
  "\n"
)

cat(
  "Total unique patient/sample IDs:",
  length(patient_ids),
  "\n"
)


############################################################
# 4. ALIGN ROWS BASED ON GENE IDS
############################################################

align_rows <- function(mat, gene_list) {
  
  # Find genes missing from this omics matrix
  missing <- setdiff(
    gene_list,
    rownames(mat)
  )
  
  # Add missing genes as NA rows
  if (length(missing) > 0) {
    
    na_mat <- matrix(
      NA_real_,
      ncol = ncol(mat),
      nrow = length(missing),
      dimnames = list(
        missing,
        colnames(mat)
      )
    )
    
    mat <- rbind(
      mat,
      na_mat
    )
  }
  
  # Make sure gene order is identical
  mat[
    gene_list,
    ,
    drop = FALSE
  ]
}


methylation_aligned <- align_rows(
  methylation,
  gene_list
)

mutation_aligned <- align_rows(
  mutation,
  gene_list
)

scnv_aligned <- align_rows(
  scnv,
  gene_list
)

rnaseq_aligned <- align_rows(
  rnaseq_z,
  gene_list
)


############################################################
# 5. ADD OMICS TYPE TO SAMPLE COLUMN NAMES
############################################################

colnames(rnaseq_aligned) <- paste0(
  colnames(rnaseq_aligned),
  "_RNA"
)

colnames(methylation_aligned) <- paste0(
  colnames(methylation_aligned),
  "_METH"
)

colnames(mutation_aligned) <- paste0(
  colnames(mutation_aligned),
  "_MUT"
)

colnames(scnv_aligned) <- paste0(
  colnames(scnv_aligned),
  "_SCNV"
)


############################################################
# 6. COMBINE ALL OMICS HORIZONTALLY
############################################################

horiz_combined_omics <- cbind(
  rnaseq_aligned,
  methylation_aligned,
  mutation_aligned,
  scnv_aligned
)


horiz_combined_omics <- horiz_combined_omics[
  order(
    rownames(horiz_combined_omics)
  ),
  ,
  drop = FALSE
]


cat(
  "Combined omics dimensions:",
  nrow(horiz_combined_omics),
  "genes x",
  ncol(horiz_combined_omics),
  "measurements\n"
)


############################################################
# 7. DOWNLOAD OMNIPATH INTERACTIONS
############################################################

# omnipath_interactions() currently gives the ncbi_tax_id
# full_join error, therefore download directly from
# the OmniPath REST API.

omnipath_url <- paste0(
  "https://omnipathdb.org/interactions?",
  "genesymbols=yes&",
  "datasets=omnipath&",
  "organisms=9606&",
  "fields=sources,references,curation_effort&",
  "license=academic"
)


interactions <- read.delim(
  omnipath_url,
  sep = "\t",
  header = TRUE,
  quote = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


############################################################
# 8. CHECK OMNIPATH DATA
############################################################

required_columns <- c(
  "source_genesymbol",
  "target_genesymbol"
)


missing_columns <- setdiff(
  required_columns,
  colnames(interactions)
)


if (length(missing_columns) > 0) {
  
  stop(
    paste(
      "Required OmniPath columns are missing:",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}


# Remove interactions without source / target gene symbols

interactions <- interactions[
  !is.na(interactions$source_genesymbol) &
    !is.na(interactions$target_genesymbol) &
    interactions$source_genesymbol != "" &
    interactions$target_genesymbol != "",
  ,
  drop = FALSE
]


cat(
  "Downloaded OmniPath interactions:",
  nrow(interactions),
  "\n"
)


############################################################
# 9. CONNECT TO ENSEMBL BIOMART
############################################################

library(biomaRt)


# Ensembl main server can sometimes be unavailable.
# Automatically try several mirrors.

ensembl_mirrors <- c(
  "www",
  "useast",
  "asia"
)


mart <- NULL


for (mirror_name in ensembl_mirrors) {
  
  cat(
    "Trying Ensembl mirror:",
    mirror_name,
    "...\n"
  )
  
  
  mart_try <- try(
    
    useEnsembl(
      biomart = "genes",
      dataset = "hsapiens_gene_ensembl",
      mirror = mirror_name
    ),
    
    silent = TRUE
  )
  
  
  if (!inherits(mart_try, "try-error")) {
    
    mart <- mart_try
    
    cat(
      "Successfully connected to Ensembl mirror:",
      mirror_name,
      "\n"
    )
    
    break
  }
  
  
  cat(
    "Failed to connect to:",
    mirror_name,
    "\n"
  )
}


if (is.null(mart)) {
  
  stop(
    paste(
      "Could not connect to any Ensembl BioMart mirror.",
      "Ensembl may currently be unavailable."
    )
  )
}


############################################################
# 10. FIND VALID GENE SYMBOLS IN OMICS DATA
############################################################

found_symbols_data <- getBM(
  attributes = "external_gene_name",
  filters = "external_gene_name",
  values = gene_list,
  mart = mart
)


found_symbols_data <- unique(
  found_symbols_data$external_gene_name
)


missing_symbols_data <- setdiff(
  gene_list,
  found_symbols_data
)


cat(
  "Genes not directly recognized in omics data:",
  length(missing_symbols_data),
  "\n"
)


############################################################
# 11. FIND VALID GENE SYMBOLS IN OMNIPATH
############################################################

interaction_gene_symbols <- unique(
  c(
    interactions$source_genesymbol,
    interactions$target_genesymbol
  )
)


found_symbols_interact <- getBM(
  attributes = "external_gene_name",
  filters = "external_gene_name",
  values = interaction_gene_symbols,
  mart = mart
)


found_symbols_interact <- unique(
  found_symbols_interact$external_gene_name
)


missing_symbols_interact <- setdiff(
  interaction_gene_symbols,
  found_symbols_interact
)


cat(
  "Genes not directly recognized in OmniPath:",
  length(missing_symbols_interact),
  "\n"
)


############################################################
# 12. FUNCTION TO FIND NEW GENE NAMES USING SYNONYMS
############################################################

get_synonym_mapping <- function(
    missing_symbols,
    mart
) {
  
  if (length(missing_symbols) == 0) {
    
    return(
      setNames(
        character(0),
        character(0)
      )
    )
  }
  
  
  mapping_table <- getBM(
    attributes = c(
      "external_gene_name",
      "external_synonym"
    ),
    filters = "external_synonym",
    values = missing_symbols,
    mart = mart
  )
  
  
  # Remove empty mappings
  mapping_table <- mapping_table[
    !is.na(mapping_table$external_gene_name) &
      !is.na(mapping_table$external_synonym) &
      mapping_table$external_gene_name != "" &
      mapping_table$external_synonym != "",
    ,
    drop = FALSE
  ]
  
  
  if (nrow(mapping_table) == 0) {
    
    return(
      setNames(
        character(0),
        character(0)
      )
    )
  }
  
  
  ##########################################################
  # Some old synonyms can map to multiple current genes.
  # Keep only unambiguous mappings.
  ##########################################################
  
  synonym_groups <- split(
    mapping_table$external_gene_name,
    mapping_table$external_synonym
  )
  
  
  synonym_groups <- lapply(
    synonym_groups,
    unique
  )
  
  
  synonym_groups <- synonym_groups[
    lengths(synonym_groups) == 1
  ]
  
  
  if (length(synonym_groups) == 0) {
    
    return(
      setNames(
        character(0),
        character(0)
      )
    )
  }
  
  
  mapping <- vapply(
    synonym_groups,
    function(x) {
      x[1]
    },
    character(1)
  )
  
  
  return(mapping)
}


############################################################
# 13. FIND ALTERNATIVE GENE NAMES
############################################################

new_names_data <- get_synonym_mapping(
  missing_symbols_data,
  mart
)


new_names_interact <- get_synonym_mapping(
  missing_symbols_interact,
  mart
)


cat(
  "Omics gene synonyms successfully updated:",
  length(new_names_data),
  "\n"
)

cat(
  "OmniPath gene synonyms successfully updated:",
  length(new_names_interact),
  "\n"
)


############################################################
# 14. REPLACE OLD GENE SYMBOLS IN OMICS DATA
############################################################

omics_rownames <- rownames(
  horiz_combined_omics
)


need_updating <- omics_rownames %in%
  names(new_names_data)


if (any(need_updating)) {
  
  omics_rownames[
    need_updating
  ] <- unname(
    new_names_data[
      omics_rownames[
        need_updating
      ]
    ]
  )
}


rownames(
  horiz_combined_omics
) <- omics_rownames


############################################################
# 15. REMOVE DUPLICATE GENES CREATED BY SYMBOL UPDATE
############################################################

# Example:
#
# old_symbol -> TP53
#
# but TP53 might already exist.
#
# In this case, keep the row containing the largest
# number of non-NA measurements.

collapse_duplicate_rows <- function(mat) {
  
  if (!anyDuplicated(
    rownames(mat)
  )) {
    
    return(mat)
  }
  
  
  gene_groups <- split(
    seq_len(
      nrow(mat)
    ),
    rownames(mat)
  )
  
  
  rows_to_keep <- vapply(
    
    gene_groups,
    
    function(index) {
      
      if (length(index) == 1) {
        
        return(index)
      }
      
      
      non_missing_counts <- rowSums(
        !is.na(
          mat[
            index,
            ,
            drop = FALSE
          ]
        )
      )
      
      
      index[
        which.max(
          non_missing_counts
        )
      ]
    },
    
    integer(1)
  )
  
  
  result <- mat[
    rows_to_keep,
    ,
    drop = FALSE
  ]
  
  
  rownames(result) <- names(
    gene_groups
  )
  
  
  return(result)
}


horiz_combined_omics <- collapse_duplicate_rows(
  horiz_combined_omics
)


############################################################
# 16. REPLACE OLD OMNIPATH SOURCE GENE SYMBOLS
############################################################

source_IDs <- interactions$source_genesymbol


need_updating <- source_IDs %in%
  names(new_names_interact)


if (any(need_updating)) {
  
  source_IDs[
    need_updating
  ] <- unname(
    new_names_interact[
      source_IDs[
        need_updating
      ]
    ]
  )
}


interactions$source_genesymbol <- source_IDs


############################################################
# 17. REPLACE OLD OMNIPATH TARGET GENE SYMBOLS
############################################################

target_IDs <- interactions$target_genesymbol


need_updating <- target_IDs %in%
  names(new_names_interact)


if (any(need_updating)) {
  
  target_IDs[
    need_updating
  ] <- unname(
    new_names_interact[
      target_IDs[
        need_updating
      ]
    ]
  )
}


interactions$target_genesymbol <- target_IDs


############################################################
# 18. FILTER INTERACTIONS BASED ON OMICS GENE LIST
############################################################

interactions <- subset(
  interactions,
  
  source_genesymbol %in%
    rownames(horiz_combined_omics) &
    
    target_genesymbol %in%
    rownames(horiz_combined_omics)
)


cat(
  "Interactions overlapping with omics genes:",
  nrow(interactions),
  "\n"
)


if (nrow(interactions) == 0) {
  
  stop(
    "No OmniPath interactions overlap with the omics data."
  )
}


############################################################
# 19. REMOVE DUPLICATED SOURCE-TARGET INTERACTIONS
############################################################

interaction_pair <- paste(
  interactions$source_genesymbol,
  interactions$target_genesymbol,
  sep = "->"
)


interactions <- interactions[
  !duplicated(
    interaction_pair
  ),
  ,
  drop = FALSE
]


cat(
  "Unique source-target interactions:",
  nrow(interactions),
  "\n"
)


############################################################
# 20. FUNCTION FOR SPEARMAN CORRELATION
#
# IMPORTANT:
# This preserves your ORIGINAL METHOD:
#
# RNA + METH + MUT + SCNV
# are horizontally concatenated first.
#
# Gene-gene correlation is then calculated across
# this complete multi-omics vector.
############################################################

correlation <- function(
    gene1,
    gene2
) {
  
  gene1_data <- as.numeric(
    horiz_combined_omics[
      gene1,
      ,
      drop = TRUE
    ]
  )
  
  
  gene2_data <- as.numeric(
    horiz_combined_omics[
      gene2,
      ,
      drop = TRUE
    ]
  )
  
  
  ##########################################################
  # Keep positions where both genes have valid data
  ##########################################################
  
  valid <- complete.cases(
    gene1_data,
    gene2_data
  )
  
  
  ##########################################################
  # Need at least 3 complete observations
  ##########################################################
  
  if (sum(valid) < 3) {
    
    return(
      c(
        ref_correlation = NA_real_,
        ref_p_value = NA_real_
      )
    )
  }
  
  
  gene1_data <- gene1_data[
    valid
  ]
  
  
  gene2_data <- gene2_data[
    valid
  ]
  
  
  ##########################################################
  # If a gene has exactly the same value everywhere,
  # correlation cannot be calculated.
  ##########################################################
  
  if (
    length(
      unique(gene1_data)
    ) < 2 ||
    length(
      unique(gene2_data)
    ) < 2
  ) {
    
    return(
      c(
        ref_correlation = NA_real_,
        ref_p_value = NA_real_
      )
    )
  }
  
  
  ##########################################################
  # Spearman correlation
  ##########################################################
  
  result <- suppressWarnings(
    cor.test(
      gene1_data,
      gene2_data,
      method = "spearman",
      exact = FALSE
    )
  )
  
  
  return(
    c(
      ref_correlation =
        unname(
          result$estimate
        ),
      
      ref_p_value =
        result$p.value
    )
  )
}


############################################################
# 21. CALCULATE CORRELATION FOR ALL OMNIPATH INTERACTIONS
############################################################

cat(
  "Calculating correlations for",
  nrow(interactions),
  "interactions...\n"
)


ref_correlations <- t(
  vapply(
    
    seq_len(
      nrow(interactions)
    ),
    
    function(i) {
      
      correlation(
        interactions$source_genesymbol[i],
        interactions$target_genesymbol[i]
      )
    },
    
    FUN.VALUE = c(
      ref_correlation = 0,
      ref_p_value = 0
    )
  )
)


colnames(
  ref_correlations
) <- c(
  "ref_correlation",
  "ref_p_value"
)


rownames(
  ref_correlations
) <- paste0(
  interactions$source_genesymbol,
  ".",
  interactions$target_genesymbol
)


############################################################
# 22. ADD CORRELATIONS TO INTERACTIONS
############################################################

# IMPORTANT FIX:
#
# Original code:
#
# correlation = ref_correlations[1]
#
# is wrong because [1] extracts only ONE value.
#
# Correct version:

interactions$correlation <- ref_correlations[
  ,
  "ref_correlation"
]


interactions$p_value <- ref_correlations[
  ,
  "ref_p_value"
]


############################################################
# 23. CALCULATE FDR
############################################################

interactions$FDR <- p.adjust(
  interactions$p_value,
  method = "BH"
)


############################################################
# 24. REMOVE INTERACTIONS WHERE CORRELATION IS NA
############################################################

interactions <- interactions[
  !is.na(
    interactions$correlation
  ),
  ,
  drop = FALSE
]


cat(
  "Interactions with valid correlations:",
  nrow(interactions),
  "\n"
)


############################################################
# 25. FILTER EDGES BASED ON CORRELATION
############################################################

filter <- abs(
  interactions$correlation
) > 0.1


interactions_filtered <- interactions[
  filter,
  ,
  drop = FALSE
]


cat(
  "Interactions with |correlation| > 0.1:",
  nrow(interactions_filtered),
  "\n"
)


if (nrow(interactions_filtered) == 0) {
  
  stop(
    "No interactions passed |correlation| > 0.1."
  )
}


############################################################
# 26. CONSTRUCT REFERENCE NETWORK
############################################################

library(igraph)


reference_nw <- graph_from_data_frame(
  interactions_filtered[
    ,
    c(
      "source_genesymbol",
      "target_genesymbol"
    )
  ],
  directed = TRUE
)


############################################################
# 27. ASSIGN CORRELATIONS AS EDGE WEIGHTS
############################################################

E(reference_nw)$weight <- interactions_filtered$correlation


############################################################
# 28. ADD P-VALUE AND FDR TO EDGES
############################################################

E(reference_nw)$p_value <- interactions_filtered$p_value

E(reference_nw)$FDR <- interactions_filtered$FDR


############################################################
# 29. KEEP ONLY THE BIGGEST CONNECTED SUBGRAPH
############################################################

network_components <- components(
  reference_nw,
  mode = "weak"
)


biggest_subgraph <- which.max(
  network_components$csize
)


reference_nw <- induced_subgraph(
  reference_nw,
  
  which(
    network_components$membership ==
      biggest_subgraph
  )
)


############################################################
# 30. GIVE NAMES TO EDGES BASED ON SOURCE/TARGET GENE
############################################################

E(reference_nw)$name <- paste0(
  
  as_ids(
    tail_of(
      reference_nw,
      E(reference_nw)
    )
  ),
  
  ".",
  
  as_ids(
    head_of(
      reference_nw,
      E(reference_nw)
    )
  )
)


############################################################
# 31. FINAL NETWORK INFORMATION
############################################################

cat(
  "\n====================================\n"
)

cat(
  "FINAL REFERENCE NETWORK\n"
)

cat(
  "====================================\n"
)


cat(
  "Number of genes:",
  vcount(reference_nw),
  "\n"
)


cat(
  "Number of edges:",
  ecount(reference_nw),
  "\n"
)


cat(
  "Correlation range:",
  range(
    E(reference_nw)$weight,
    na.rm = TRUE
  ),
  "\n"
)


############################################################
# 32. EXPORT GRAPHML FOR CYTOSCAPE
############################################################

output_file <- "../reference_network.graphml"


write_graph(
  reference_nw,
  file = output_file,
  format = "graphml"
)


cat(
  "\nReference network successfully exported:\n"
)

cat(
  output_file,
  "\n"
)

# ============================================================
# 0. LOAD PACKAGES
# ============================================================

library(parallel)
library(igraph)
library(ggplot2)
library(tidyverse)


# ============================================================
# 1. PARALLEL SETTINGS
# ============================================================

# You requested 10 CPU cores from Slurm
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

# Do not allow each worker to create extra MKL/BLAS threads
Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1"
)

cat("Using", n_cores, "CPU cores\n")

cat(
  "SLURM_CPUS_PER_TASK =",
  Sys.getenv("SLURM_CPUS_PER_TASK"),
  "\n"
)

cat(
  "SLURM_CPUS_ON_NODE =",
  Sys.getenv("SLURM_CPUS_ON_NODE"),
  "\n"
)


# ============================================================
# 2. PREPARE FILTERED INTERACTIONS
# ============================================================

# Use the final filtered interaction set everywhere
edges <- as.data.frame(
  interactions_filtered,
  stringsAsFactors = FALSE
)

required_cols <- c(
  "source_genesymbol",
  "target_genesymbol"
)

if (!all(required_cols %in% colnames(edges))) {
  stop(
    "interactions_filtered must contain ",
    "'source_genesymbol' and 'target_genesymbol'."
  )
}

# Convert explicitly to character
edges$source_genesymbol <- as.character(
  edges$source_genesymbol
)

edges$target_genesymbol <- as.character(
  edges$target_genesymbol
)

# Permanent unique edge IDs
edges$edge_id <- sprintf(
  "edge_%07d",
  seq_len(nrow(edges))
)

source_genes <- edges$source_genesymbol
target_genes <- edges$target_genesymbol
edge_ids <- edges$edge_id

cat(
  "\nNumber of filtered interactions:",
  nrow(edges),
  "\n"
)


# ============================================================
# 3. RECALCULATE REFERENCE CORRELATIONS
#    USING EXACTLY THE SAME FILTERED INTERACTIONS
# ============================================================

# IMPORTANT:
# Do NOT use the old ref_correlations object here.
# The old object contains 32469 rows whereas the final
# network contains 24373 filtered interactions.
#
# We therefore recalculate reference correlations directly
# on interactions_filtered.

calculate_ref_correlation <- function(
    gene1,
    gene2
) {
  
  # Genes must exist
  if (
    !(gene1 %in% rownames(horiz_combined_omics)) ||
    !(gene2 %in% rownames(horiz_combined_omics))
  ) {
    
    return(
      c(
        ref_correlation = NA_real_,
        ref_p_value = NA_real_
      )
    )
  }
  
  
  gene1_data <- as.numeric(
    horiz_combined_omics[
      gene1,
      ,
      drop = TRUE
    ]
  )
  
  gene2_data <- as.numeric(
    horiz_combined_omics[
      gene2,
      ,
      drop = TRUE
    ]
  )
  
  
  valid <- (
    is.finite(gene1_data) &
      is.finite(gene2_data)
  )
  
  
  if (sum(valid) < 3) {
    
    return(
      c(
        ref_correlation = NA_real_,
        ref_p_value = NA_real_
      )
    )
  }
  
  
  x <- gene1_data[valid]
  y <- gene2_data[valid]
  
  
  # Correlation undefined for constant vectors
  if (
    length(unique(x)) < 2 ||
    length(unique(y)) < 2
  ) {
    
    return(
      c(
        ref_correlation = NA_real_,
        ref_p_value = NA_real_
      )
    )
  }
  
  
  result <- suppressWarnings(
    cor.test(
      x,
      y,
      method = "spearman",
      exact = FALSE
    )
  )
  
  
  c(
    ref_correlation =
      unname(result$estimate),
    
    ref_p_value =
      result$p.value
  )
}


cat(
  "\nCalculating reference correlations for",
  nrow(edges),
  "edges using",
  n_cores,
  "cores...\n"
)

ref_start_time <- Sys.time()


ref_list <- parallel::mclapply(
  
  X = seq_len(nrow(edges)),
  
  FUN = function(j) {
    
    calculate_ref_correlation(
      gene1 = source_genes[j],
      gene2 = target_genes[j]
    )
  },
  
  mc.cores = n_cores,
  mc.preschedule = TRUE
)


ref_correlations_use <- do.call(
  rbind,
  ref_list
)

rownames(ref_correlations_use) <- edge_ids

colnames(ref_correlations_use) <- c(
  "ref_correlation",
  "ref_p_value"
)


ref_end_time <- Sys.time()


cat(
  "\nReference correlations finished.\n"
)

cat(
  "Reference rows:",
  nrow(ref_correlations_use),
  "\n"
)

cat(
  "Filtered interaction rows:",
  nrow(edges),
  "\n"
)

cat(
  "Reference calculation time:",
  round(
    as.numeric(
      difftime(
        ref_end_time,
        ref_start_time,
        units = "mins"
      )
    ),
    2
  ),
  "minutes\n"
)


stopifnot(
  nrow(ref_correlations_use) ==
    nrow(edges)
)

cat(
  "Reference-edge alignment OK.\n"
)


# ============================================================
# 4. PRECOMPUTE PATIENT COLUMN MASKS
# ============================================================

# This avoids repeatedly constructing patient column names
# inside millions of correlation calculations.

patient_keep_masks <- lapply(
  
  patient_ids,
  
  function(patient_id) {
    
    patient_columns <- paste0(
      patient_id,
      c(
        "_RNA",
        "_METH",
        "_MUT",
        "_SCNV"
      )
    )
    
    !colnames(horiz_combined_omics) %in%
      patient_columns
  }
)


# ============================================================
# 5. LEAVE-ONE-OUT SPEARMAN CORRELATION
# ============================================================

loo_correlation <- function(
    gene1,
    gene2,
    leave_one_out
) {
  
  # Genes must exist
  if (
    !(gene1 %in% rownames(horiz_combined_omics)) ||
    !(gene2 %in% rownames(horiz_combined_omics))
  ) {
    
    return(
      c(
        correlation = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  
  # Remove four omics columns belonging to this patient
  keep <- patient_keep_masks[[leave_one_out]]
  
  
  gene1_data <- as.numeric(
    horiz_combined_omics[
      gene1,
      keep,
      drop = TRUE
    ]
  )
  
  gene2_data <- as.numeric(
    horiz_combined_omics[
      gene2,
      keep,
      drop = TRUE
    ]
  )
  
  
  valid <- (
    is.finite(gene1_data) &
      is.finite(gene2_data)
  )
  
  
  if (sum(valid) < 3) {
    
    return(
      c(
        correlation = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  
  x <- gene1_data[valid]
  y <- gene2_data[valid]
  
  
  if (
    length(unique(x)) < 2 ||
    length(unique(y)) < 2
  ) {
    
    return(
      c(
        correlation = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  
  result <- suppressWarnings(
    cor.test(
      x,
      y,
      method = "spearman",
      exact = FALSE
    )
  )
  
  
  c(
    correlation =
      unname(result$estimate),
    
    p_value =
      result$p.value
  )
}


# ============================================================
# 6. CALCULATE ONE PATIENT
# ============================================================

calculate_one_patient <- function(i) {
  
  result <- vapply(
    
    seq_len(nrow(edges)),
    
    FUN = function(j) {
      
      loo_correlation(
        
        gene1 =
          source_genes[j],
        
        gene2 =
          target_genes[j],
        
        leave_one_out =
          i
      )
    },
    
    FUN.VALUE = c(
      correlation = 0.0,
      p_value = 0.0
    )
  )
  
  
  # vapply result:
  # 2 x number_of_edges
  #
  # transpose:
  # number_of_edges x 2
  
  result <- t(result)
  
  
  rownames(result) <- edge_ids
  
  
  colnames(result) <- c(
    patient_ids[i],
    paste0(
      patient_ids[i],
      "_p"
    )
  )
  
  
  result
}


# ============================================================
# 7. PARALLEL LOO CALCULATION
#    10 CPU CORES
# ============================================================

cat(
  "\nStarting parallel leave-one-out calculation...\n"
)

cat(
  "Patients:",
  length(patient_ids),
  "\n"
)

cat(
  "Edges per patient:",
  nrow(edges),
  "\n"
)

cat(
  "CPU cores:",
  n_cores,
  "\n"
)


loo_start_time <- Sys.time()


patient_specific_correlations_list <-
  parallel::mclapply(
    
    X = seq_along(patient_ids),
    
    FUN = calculate_one_patient,
    
    mc.cores = n_cores,
    
    mc.preschedule = TRUE
  )


loo_end_time <- Sys.time()


cat(
  "\nLOO correlation calculation finished.\n"
)

cat(
  "Elapsed time:",
  round(
    as.numeric(
      difftime(
        loo_end_time,
        loo_start_time,
        units = "mins"
      )
    ),
    2
  ),
  "minutes\n"
)


# ============================================================
# 8. CHECK PARALLEL RESULTS
# ============================================================

failed_workers <- vapply(
  
  patient_specific_correlations_list,
  
  inherits,
  
  logical(1),
  
  what = "try-error"
)


if (any(failed_workers)) {
  
  stop(
    "LOO calculation failed for: ",
    paste(
      patient_ids[
        failed_workers
      ],
      collapse = ", "
    )
  )
}


# ============================================================
# 9. COMBINE LOO RESULTS
# ============================================================

patient_specific_correlations_loo <-
  do.call(
    cbind,
    patient_specific_correlations_list
  )


cat(
  "\nLOO matrix:\n"
)

print(
  dim(
    patient_specific_correlations_loo
  )
)


cat(
  "\nReference matrix:\n"
)

print(
  dim(
    ref_correlations_use
  )
)


# They MUST have exactly the same number of rows
stopifnot(
  nrow(ref_correlations_use) ==
    nrow(
      patient_specific_correlations_loo
    )
)


# ============================================================
# 10. COMBINE REFERENCE + LOO RESULTS
# ============================================================

patient_specific_correlations <-
  cbind(
    ref_correlations_use,
    patient_specific_correlations_loo
  )


cat(
  "\nFinal combined correlation matrix:\n"
)

print(
  dim(
    patient_specific_correlations
  )
)


# ============================================================
# 11. SEPARATE CORRELATIONS AND P-VALUES
# ============================================================

patient_correlation_matrix <-
  patient_specific_correlations[
    ,
    patient_ids,
    drop = FALSE
  ]


patient_pvalue_matrix <-
  patient_specific_correlations[
    ,
    paste0(
      patient_ids,
      "_p"
    ),
    drop = FALSE
  ]


rownames(
  patient_correlation_matrix
) <- edge_ids


rownames(
  patient_pvalue_matrix
) <- edge_ids


# ============================================================
# 12. SAVE EXPENSIVE CORRELATION RESULTS
# ============================================================

dir.create(
  "../result_data",
  showWarnings = FALSE,
  recursive = TRUE
)


saveRDS(
  ref_correlations_use,
  "../result_data/ref_correlations_filtered.rds"
)


saveRDS(
  patient_specific_correlations,
  "../result_data/patient_specific_correlations.rds"
)


saveRDS(
  patient_correlation_matrix,
  "../result_data/patient_correlation_matrix.rds"
)


saveRDS(
  patient_pvalue_matrix,
  "../result_data/patient_pvalue_matrix.rds"
)


cat(
  "\nCorrelation results saved successfully.\n"
)


# ============================================================
# 13. SELECT PATIENT FOR VISUALIZATION
# ============================================================

patient_id <- "TCGA.ZG.A9L5"


if (!(patient_id %in% patient_ids)) {
  
  stop(
    patient_id,
    " was not found in patient_ids."
  )
}


# ============================================================
# 14. CONSTRUCT PATIENT NETWORK
# ============================================================

network_edges <- data.frame(
  
  from =
    source_genes,
  
  to =
    target_genes,
  
  edge_id =
    edge_ids,
  
  stringsAsFactors =
    FALSE
)


patient_nw <-
  graph_from_data_frame(
    network_edges,
    directed = TRUE
  )


# Assign patient's LOO correlation
E(patient_nw)$correlation <-
  patient_correlation_matrix[
    E(patient_nw)$edge_id,
    patient_id
  ]


# ============================================================
# 15. LARGEST CONNECTED COMPONENT
# ============================================================

component_result <-
  components(
    patient_nw,
    mode = "weak"
  )


biggest <-
  which.max(
    component_result$csize
  )


patient_nw <-
  induced_subgraph(
    
    patient_nw,
    
    which(
      component_result$membership ==
        biggest
    )
  )


cat(
  "\nLargest network component:",
  vcount(patient_nw),
  "nodes and",
  ecount(patient_nw),
  "edges\n"
)


# ============================================================
# 16. HUMAN-READABLE EDGE NAMES
# ============================================================

E(patient_nw)$name <-
  paste0(
    
    as_ids(
      tail_of(
        patient_nw,
        E(patient_nw)
      )
    ),
    
    ".",
    
    as_ids(
      head_of(
        patient_nw,
        E(patient_nw)
      )
    )
  )


# ============================================================
# 17. EXPORT GRAPHML
# ============================================================

write_graph(
  
  patient_nw,
  
  file =
    paste0(
      "../",
      patient_id,
      "_patient_network.graphml"
    ),
  
  format =
    "graphml"
)


# ============================================================
# 18. PREPARE PATIENT OMICS
# ============================================================

patient_omics_columns <-
  paste0(
    
    patient_id,
    
    c(
      "_RNA",
      "_METH",
      "_MUT",
      "_SCNV"
    )
  )


missing_omics_columns <-
  setdiff(
    patient_omics_columns,
    colnames(
      horiz_combined_omics
    )
  )


if (
  length(
    missing_omics_columns
  ) > 0
) {
  
  stop(
    "Missing patient omics columns: ",
    paste(
      missing_omics_columns,
      collapse = ", "
    )
  )
}


patient_omics <-
  horiz_combined_omics[
    
    V(patient_nw)$name,
    
    patient_omics_columns,
    
    drop = FALSE
  ]


# ============================================================
# 19. RNA SCALING
# ============================================================

scale_rna <- function(x) {
  
  q <- quantile(
    
    x,
    
    probs = c(
      0.05,
      0.95
    ),
    
    na.rm = TRUE
  )
  
  
  # Avoid division by zero
  if (
    !all(
      is.finite(q)
    ) ||
    q[2] <= q[1]
  ) {
    
    result <- rep(
      0.5,
      length(x)
    )
    
    result[
      is.na(x)
    ] <- NA_real_
    
    return(result)
  }
  
  
  x_scaled <-
    (
      x - q[1]
    ) /
    (
      q[2] - q[1]
    )
  
  
  x_scaled[
    x_scaled < 0
  ] <- 0
  
  
  x_scaled[
    x_scaled > 1
  ] <- 1
  
  
  x_scaled
}


# ============================================================
# 20. SCALE PATIENT OMICS
# ============================================================

patient_omics_scaled <-
  patient_omics


# RNA
patient_omics_scaled[, 1] <-
  scale_rna(
    patient_omics[, 1]
  )


# Methylation
patient_omics_scaled[, 2] <-
  patient_omics[, 2] + 0.5


# Mutation stays as original
# Assumed to already be 0/1


# SCNV
patient_omics_scaled[, 4] <-
  (
    patient_omics[, 4] + 2
  ) / 4


# Clip all numeric values to [0,1]
patient_omics_scaled[
  patient_omics_scaled < 0
] <- 0


patient_omics_scaled[
  patient_omics_scaled > 1
] <- 1


# ============================================================
# 21. RADIAL GLYPH FUNCTION
# ============================================================

create_radial_glyph <- function(
    values,
    gene_name,
    na_radius = 0.3,
    max_radius = 1
) {
  
  df <- data.frame(
    
    omic = factor(
      
      c(
        "RNA",
        "METH",
        "MUT",
        "SCNV"
      ),
      
      levels = c(
        "RNA",
        "METH",
        "MUT",
        "SCNV"
      )
    ),
    
    value =
      as.numeric(
        values
      )
  )
  
  
  df$radius <-
    ifelse(
      is.na(df$value),
      na_radius,
      df$value
    )
  
  
  df$fill <-
    ifelse(
      is.na(df$value),
      "NA",
      as.character(
        df$omic
      )
    )
  
  
  ggplot(
    
    df,
    
    aes(
      x = omic,
      y = radius,
      fill = fill
    )
  ) +
    
    geom_col(
      width = 1,
      color = "black",
      linewidth = 0.2
    ) +
    
    coord_polar(
      start = 0
    ) +
    
    scale_y_continuous(
      limits = c(
        0,
        max_radius
      )
    ) +
    
    scale_fill_manual(
      values = c(
        RNA = "blue",
        METH = "red",
        MUT = "green",
        SCNV = "orange",
        "NA" = "grey80"
      )
    ) +
    
    theme_void() +
    
    theme(
      
      legend.position =
        "none",
      
      panel.background =
        element_rect(
          fill = "white",
          colour = NA
        ),
      
      plot.background =
        element_rect(
          fill = "transparent",
          colour = NA
        )
    )
}


# ============================================================
# 22. CREATE SVG GLYPHS
# ============================================================

svg_output_dir <-
  paste0(
    "../",
    patient_id,
    "_radial_glyph_svgs"
  )


dir.create(
  svg_output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


for (
  gene in rownames(
    patient_omics_scaled
  )
) {
  
  p <- create_radial_glyph(
    
    patient_omics_scaled[
      gene,
      ,
      drop = TRUE
    ],
    
    gene
  )
  
  
  ggsave(
    
    filename =
      file.path(
        svg_output_dir,
        paste0(
          gene,
          ".svg"
        )
      ),
    
    plot = p,
    
    width = 2,
    
    height = 2,
    
    units = "in"
  )
}


# ============================================================
# 23. CREATE CYTOSCAPE IMAGE PATHS
# ============================================================

svg_dir <-
  paste0(
    
    "file://",
    
    normalizePath(
      svg_output_dir
    )
  )


cyto_images <-
  data.frame(
    
    name =
      rownames(
        patient_omics_scaled
      ),
    
    radial_glyph =
      file.path(
        
        svg_dir,
        
        paste0(
          rownames(
            patient_omics_scaled
          ),
          ".svg"
        )
      ),
    
    stringsAsFactors = FALSE
  )


write.csv(
  
  cyto_images,
  
  paste0(
    "../",
    patient_id,
    "_node_graphic_paths.csv"
  ),
  
  row.names = FALSE
)


# ============================================================
# 24. PREPARE NETWORK CORRELATIONS FOR CENTRALITY
# ============================================================

component_edge_ids <-
  E(patient_nw)$edge_id


network_correlations <-
  patient_correlation_matrix[
    component_edge_ids,
    ,
    drop = FALSE
  ]


cat(
  "\nCentrality input:",
  nrow(network_correlations),
  "edges x",
  ncol(network_correlations),
  "patients\n"
)


# ============================================================
# 25. CENTRALITY FOR ONE PATIENT
# ============================================================

calculate_centrality_one_patient <-
  function(i) {
    
    # Worker gets its own graph copy
    g <- patient_nw
    
    
    correlations <-
      as.numeric(
        network_correlations[
          ,
          i
        ]
      )
    
    
    # Valid correlations only
    valid <-
      is.finite(
        correlations
      )
    
    
    # --------------------------------------------------------
    # Full output vectors
    # --------------------------------------------------------
    
    edge_betweenness_full <-
      setNames(
        
        rep(
          NA_real_,
          ecount(g)
        ),
        
        E(g)$edge_id
      )
    
    
    node_names <-
      V(g)$name
    
    
    betweenness_full <-
      setNames(
        
        rep(
          NA_real_,
          length(node_names)
        ),
        
        node_names
      )
    
    
    strength_full <-
      setNames(
        
        rep(
          NA_real_,
          length(node_names)
        ),
        
        node_names
      )
    
    
    closeness_full <-
      setNames(
        
        rep(
          NA_real_,
          length(node_names)
        ),
        
        node_names
      )
    
    
    # No valid edges
    if (sum(valid) == 0) {
      
      return(
        list(
          
          edge_betweenness =
            edge_betweenness_full,
          
          betweenness =
            betweenness_full,
          
          strength =
            strength_full,
          
          closeness =
            closeness_full
        )
      )
    }
    
    
    # --------------------------------------------------------
    # Remove NA edges
    # --------------------------------------------------------
    
    g_valid <-
      delete_edges(
        g,
        E(g)[!valid]
      )
    
    
    correlations_valid <-
      correlations[valid]
    
    
    # ========================================================
    # ASSOCIATION STRENGTH
    #
    # Strong correlation = strong connection
    # ========================================================
    
    association_strength <-
      abs(
        correlations_valid
      )
    
    
    # ========================================================
    # DISTANCE FOR SHORTEST PATH MEASURES
    #
    # strong correlation -> short distance
    # weak correlation   -> long distance
    # ========================================================
    
    epsilon <- 1e-6
    
    
    distance_weights <-
      1 /
      pmax(
        association_strength,
        epsilon
      )
    
    
    # ========================================================
    # EDGE BETWEENNESS
    # ========================================================
    
    eb <-
      igraph::edge_betweenness(
        
        g_valid,
        
        directed = TRUE,
        
        weights =
          distance_weights
      )
    
    
    edge_betweenness_full[
      E(g_valid)$edge_id
    ] <- eb
    
    
    # ========================================================
    # NODE BETWEENNESS
    # ========================================================
    
    bc <-
      igraph::betweenness(
        
        g_valid,
        
        directed = TRUE,
        
        weights =
          distance_weights,
        
        normalized = TRUE
      )
    
    
    betweenness_full[
      names(bc)
    ] <- bc
    
    
    # ========================================================
    # WEIGHTED DEGREE / STRENGTH
    #
    # IMPORTANT:
    # Use abs(correlation)
    # NOT 1/abs(correlation)
    # ========================================================
    
    sc <-
      igraph::strength(
        
        g_valid,
        
        mode = "all",
        
        weights =
          association_strength
      )
    
    
    strength_full[
      names(sc)
    ] <- sc
    
    
    # ========================================================
    # CLOSENESS
    # ========================================================
    
    E(g_valid)$weight <-
      distance_weights
    
    
    g_und <-
      as.undirected(
        
        g_valid,
        
        mode = "collapse",
        
        edge.attr.comb =
          list(
            weight = "min",
            "ignore"
          )
      )
    
    
    cc <-
      igraph::closeness(
        
        g_und,
        
        weights =
          E(g_und)$weight,
        
        normalized = TRUE
      )
    
    
    closeness_full[
      names(cc)
    ] <- cc
    
    
    # ========================================================
    # RETURN
    # ========================================================
    
    list(
      
      edge_betweenness =
        edge_betweenness_full,
      
      betweenness =
        betweenness_full,
      
      strength =
        strength_full,
      
      closeness =
        closeness_full
    )
  }


# ============================================================
# 26. PARALLEL CENTRALITY
#     10 CPU CORES
# ============================================================

cat(
  "\nStarting parallel centrality calculation...\n"
)


centrality_start <-
  Sys.time()


centrality_results <-
  parallel::mclapply(
    
    X =
      seq_along(
        patient_ids
      ),
    
    FUN =
      calculate_centrality_one_patient,
    
    mc.cores =
      n_cores,
    
    mc.preschedule =
      TRUE
  )


centrality_end <-
  Sys.time()


cat(
  "\nCentrality calculation finished.\n"
)


cat(
  "Elapsed time:",
  round(
    as.numeric(
      difftime(
        centrality_end,
        centrality_start,
        units = "mins"
      )
    ),
    2
  ),
  "minutes\n"
)


# ============================================================
# 27. CHECK CENTRALITY WORKERS
# ============================================================

centrality_failed <-
  vapply(
    
    centrality_results,
    
    inherits,
    
    logical(1),
    
    what = "try-error"
  )


if (any(centrality_failed)) {
  
  stop(
    "Centrality failed for patients: ",
    paste(
      patient_ids[
        centrality_failed
      ],
      collapse = ", "
    )
  )
}


# ============================================================
# 28. CONVERT CENTRALITIES TO MATRICES
# ============================================================

edge_betweenness_matrix <-
  do.call(
    
    cbind,
    
    lapply(
      centrality_results,
      `[[`,
      "edge_betweenness"
    )
  )


betw_centrality <-
  do.call(
    
    cbind,
    
    lapply(
      centrality_results,
      `[[`,
      "betweenness"
    )
  )


degree_centrality <-
  do.call(
    
    cbind,
    
    lapply(
      centrality_results,
      `[[`,
      "strength"
    )
  )


closeness_centrality <-
  do.call(
    
    cbind,
    
    lapply(
      centrality_results,
      `[[`,
      "closeness"
    )
  )


# ============================================================
# 29. ASSIGN PATIENT COLUMN NAMES
# ============================================================

colnames(
  edge_betweenness_matrix
) <- patient_ids


colnames(
  betw_centrality
) <- patient_ids


colnames(
  degree_centrality
) <- patient_ids


colnames(
  closeness_centrality
) <- patient_ids


# ============================================================
# 30. SAVE CENTRALITY RESULTS
# ============================================================

saveRDS(
  
  edge_betweenness_matrix,
  
  "../result_data/edge_betweenness_matrix.rds"
)


saveRDS(
  
  betw_centrality,
  
  "../result_data/betweenness_centrality.rds"
)


saveRDS(
  
  degree_centrality,
  
  "../result_data/degree_strength_centrality.rds"
)


saveRDS(
  
  closeness_centrality,
  
  "../result_data/closeness_centrality.rds"
)


saveRDS(
  
  list(
    
    reference_correlations =
      ref_correlations_use,
    
    patient_correlations =
      patient_correlation_matrix,
    
    patient_pvalues =
      patient_pvalue_matrix,
    
    edge_betweenness =
      edge_betweenness_matrix,
    
    betweenness =
      betw_centrality,
    
    strength =
      degree_centrality,
    
    closeness =
      closeness_centrality
    
  ),
  
  "../result_data/patient_network_analysis_all.rds"
)


# ============================================================
# 31. FINAL SUMMARY
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "Patient network analysis completed.\n"
)

cat(
  "CPU cores used:",
  n_cores,
  "\n"
)

cat(
  "Filtered edges:",
  nrow(edges),
  "\n"
)

cat(
  "Patients:",
  length(patient_ids),
  "\n"
)

cat(
  "Reference matrix:",
  nrow(ref_correlations_use),
  "x",
  ncol(ref_correlations_use),
  "\n"
)

cat(
  "Patient correlation matrix:",
  nrow(patient_correlation_matrix),
  "x",
  ncol(patient_correlation_matrix),
  "\n"
)

cat(
  "Results saved to ../result_data/\n"
)

cat(
  "========================================\n"
)
