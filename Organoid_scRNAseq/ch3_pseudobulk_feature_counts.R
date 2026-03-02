# =============================================================================
# ch3_pseudobulk_feature_counts.R
# Author: Manveer Chauhan
#
# Purpose: For each organoid replicate across the three atlases, pseudobulk
# raw counts (summed across all cells per replicate) and report:
#   - Number of genes with pseudobulk count >= 10  (RNA assay)
#   - Number of isoforms with pseudobulk count >= 5 (iso assay)
#
# Uses the same integrated Seurat objects as ch3_isoTWAS_DE_heatmaps.R.
# Loads one atlas at a time to minimise peak memory usage.
# Source this script directly into your console.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

obj_dir <- "/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space/output_files/integrated_objects"

atlas_files <- c(
  "1M_Org (Day 1)" = "1M_Org_integrated_harmony_consensus_with_isoform_assay.rds",
  "3M_Org (Day 3)" = "3M_Org_integrated_harmony_consensus_with_isoform_assay.rds",
  "6M_Org (Day 6)" = "6M_Org_integrated_harmony_consensus_with_isoform_assay.rds"
)

atlas_samples <- list(
  "1M_Org (Day 1)" = c("org_1A", "org_1B"),
  "3M_Org (Day 3)" = c("org_3A", "org_3B", "org_3C"),
  "6M_Org (Day 6)" = c("org_6A", "org_6B", "org_6C")
)

gene_threshold    <- 10
isoform_threshold <- 5

results <- data.frame(
  Atlas         = character(),
  Sample        = character(),
  N_Cells       = integer(),
  Genes_gte10   = integer(),
  Isoforms_gte5 = integer(),
  stringsAsFactors = FALSE
)

cat("\nLoading atlas objects and computing pseudobulks...\n\n")

for (atlas_name in names(atlas_files)) {

  rds_path <- file.path(obj_dir, atlas_files[[atlas_name]])
  cat(sprintf("[%s] Loading %s ...\n", atlas_name, basename(rds_path)))

  seu <- readRDS(rds_path)

  cat(sprintf("[%s] Running AggregateExpression (RNA + iso assays)...\n", atlas_name))

  # Pseudobulk raw counts per sample_id for both assays in one pass
  pb_gene <- AggregateExpression(seu, assays = "RNA", group.by = "sample_id",
                                 slot = "counts")$RNA
  pb_iso  <- AggregateExpression(seu, assays = "iso",  group.by = "sample_id",
                                 slot = "counts")$iso

  # Cell counts per sample from metadata
  cell_counts <- table(seu$sample_id)

  rm(seu); gc(verbose = FALSE)

  for (samp in atlas_samples[[atlas_name]]) {
    # AggregateExpression replaces underscores with dashes in column names
    samp_col   <- gsub("_", "-", samp)
    n_cells    <- as.integer(cell_counts[samp])
    n_genes    <- sum(pb_gene[, samp_col] >= gene_threshold)
    n_isoforms <- sum(pb_iso[,  samp_col] >= isoform_threshold)

    results <- rbind(results, data.frame(
      Atlas         = atlas_name,
      Sample        = samp,
      N_Cells       = n_cells,
      Genes_gte10   = n_genes,
      Isoforms_gte5 = n_isoforms,
      stringsAsFactors = FALSE
    ))
    cat(sprintf("  %-10s  %d cells | %d genes >=10 | %d isoforms >=5\n",
                samp, n_cells, n_genes, n_isoforms))
  }

  rm(pb_gene, pb_iso, cell_counts); gc(verbose = FALSE)
  cat("\n")
}

# =============================================================================
# Print formatted summary table
# =============================================================================
cat("=============================================================================\n")
cat(" Pseudobulk Feature Counts per Replicate\n")
cat(sprintf(" Gene threshold: >= %d counts | Isoform threshold: >= %d counts\n",
            gene_threshold, isoform_threshold))
cat("=============================================================================\n")
cat(sprintf("%-20s %-10s %10s %14s %14s\n",
            "Atlas", "Sample", "N Cells", "Genes >=10", "Isoforms >=5"))
cat(sprintf("%-20s %-10s %10s %14s %14s\n",
            "--------------------", "----------", "----------",
            "--------------", "--------------"))

current_atlas <- ""
for (i in seq_len(nrow(results))) {
  row        <- results[i, ]
  atlas_label <- if (row$Atlas != current_atlas) { current_atlas <- row$Atlas; row$Atlas } else ""
  cat(sprintf("%-20s %-10s %10s %14s %14s\n",
              atlas_label,
              row$Sample,
              format(row$N_Cells,       big.mark = ","),
              format(row$Genes_gte10,   big.mark = ","),
              format(row$Isoforms_gte5, big.mark = ",")))
}
cat("=============================================================================\n\n")
