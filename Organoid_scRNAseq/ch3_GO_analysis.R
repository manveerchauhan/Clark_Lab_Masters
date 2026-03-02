#!/usr/bin/env Rscript
# GO Enrichment Analysis for Consensus-Annotated Organoid Objects
# Author: Manveer Chauhan
# Description: Performs Gene Ontology enrichment analysis on consensus cell types
#              across three developmental timepoints (1M, 3M, 6M) using gene-level
#              expression data from integrated organoid objects.

# Load required libraries ----
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(tidyverse)
library(ggplot2)

# Set working directory
setwd("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space")

# Helper Functions ----

#' Extract clean ENSEMBL ID from Seurat rowname
#'
#' Converts rownames like "ENSG00000130844.19-ZNF331" to "ENSG00000130844"
#'
#' @param rowname Character vector of Seurat feature names
#' @return Character vector of clean ENSEMBL IDs (version suffix removed)
extract_ensembl_id <- function(rowname) {
  # Split on '-' to get ENSEMBL portion
  ensembl_with_version <- str_split_fixed(rowname, "-", 2)[, 1]

  # Remove version suffix (e.g., ".19" -> "")
  ensembl_clean <- str_replace(ensembl_with_version, "\\.\\d+$", "")

  return(ensembl_clean)
}

#' Convert ENSEMBL IDs to ENTREZ IDs with error handling
#'
#' @param ensembl_ids Character vector of ENSEMBL gene IDs
#' @param feature_names Character vector of original feature names (for tracking)
#' @return data.frame with columns: feature_name, ensembl_id, entrez_id, conversion_success
convert_ensembl_to_entrez <- function(ensembl_ids, feature_names) {
  # Create input dataframe
  input_df <- data.frame(
    feature_name = feature_names,
    ensembl_id = ensembl_ids,
    stringsAsFactors = FALSE
  )

  # Attempt conversion using bitr
  cat("  Converting", length(unique(ensembl_ids)), "unique ENSEMBL IDs to ENTREZ IDs...\n")

  conversion_result <- tryCatch({
    bitr(
      geneID = unique(ensembl_ids),
      fromType = "ENSEMBL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
  }, error = function(e) {
    warning("bitr conversion failed: ", e$message)
    return(data.frame(ENSEMBL = character(0), ENTREZID = character(0)))
  })

  # Rename columns for clarity
  colnames(conversion_result) <- c("ensembl_id", "entrez_id")

  # Keep only first ENTREZ ID per ENSEMBL ID to avoid duplicates
  # (Some ENSEMBL IDs map to multiple ENTREZ IDs)
  conversion_result <- conversion_result %>%
    dplyr::group_by(ensembl_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  # Merge back to original data
  output_df <- left_join(input_df, conversion_result,
                         by = "ensembl_id",
                         relationship = "many-to-one")

  # Calculate success rate
  success_rate <- sum(!is.na(output_df$entrez_id)) / nrow(output_df) * 100
  cat("  Conversion success rate:", round(success_rate, 2), "%\n")
  cat("  Successfully converted:", sum(!is.na(output_df$entrez_id)), "genes\n")
  cat("  Failed conversions:", sum(is.na(output_df$entrez_id)), "genes\n\n")

  return(output_df)
}

#' Perform GO enrichment analysis for a single timepoint
#'
#' @param timepoint Character, one of "1M_Org", "3M_Org", "6M_Org"
#' @param output_dir Path to output directory
#' @return List of GO results (BP, MF, CC) and plots
analyze_timepoint <- function(timepoint, output_dir) {

  cat("\n========================================\n")
  cat("Processing timepoint:", timepoint, "\n")
  cat("========================================\n\n")

  # 1. Load Seurat object ----
  obj_path <- paste0("./output_files/integrated_objects/",
                     timepoint, "_integrated_harmony_consensus_with_isoform_assay.rds")

  cat("Loading object:", obj_path, "\n")
  seu <- readRDS(obj_path)

  # Set default assay to RNA and identity to consensus_cell_type
  DefaultAssay(seu) <- "RNA"
  Idents(seu) <- "consensus_cell_type"

  cat("  Cells:", ncol(seu), "\n")
  cat("  Genes:", nrow(seu), "\n")
  cat("  Cell types:", length(unique(seu$consensus_cell_type)), "\n")
  cat("  Cell type distribution:\n")
  print(table(seu$consensus_cell_type))
  cat("\n")

  # Join layers for Seurat v5 compatibility
  cat("Joining data layers for FindAllMarkers...\n")
  seu <- JoinLayers(seu, assay = "RNA")
  cat("  Layers joined successfully\n\n")

  # 2. Find marker genes ----
  cat("Finding marker genes for each consensus cell type...\n")
  cat("  Parameters: only.pos = TRUE, Seurat defaults (min.pct=0.01, logfc.threshold=0.1)\n")
  cat("  Post-filtering: p_val_adj <= 0.05, avg_log2FC >= 0.5\n\n")

  markers <- FindAllMarkers(
    seu,
    assay = "RNA",
    only.pos = TRUE,
    verbose = TRUE
  )

  # Check if markers were found
  if (is.null(markers) || nrow(markers) == 0) {
    cat("\n  WARNING: No markers found for this timepoint!\n")
    cat("  Skipping this timepoint...\n\n")
    return(NULL)
  }

  # Filter for significance and logFC threshold
  markers_sig <- markers %>% dplyr::filter(p_val_adj <= 0.05 & avg_log2FC >= 0.5)

  cat("\n  Total significant markers:", nrow(markers_sig), "\n")
  cat("  Markers per cell type:\n")
  print(table(markers_sig$cluster))
  cat("\n")

  # 3. Convert gene IDs ----
  cat("Converting gene IDs to ENTREZ format...\n")

  # Extract ENSEMBL IDs from feature names
  markers_sig$ensembl_id <- extract_ensembl_id(markers_sig$gene)

  # Show sample of conversions
  cat("  Sample of ID extraction:\n")
  sample_markers <- head(markers_sig[, c("gene", "ensembl_id")], 5)
  print(sample_markers)
  cat("\n")

  # Convert to ENTREZ
  conversion_df <- convert_ensembl_to_entrez(
    markers_sig$ensembl_id,
    markers_sig$gene
  )

  # Merge back to markers
  markers_with_entrez <- markers_sig %>%
    left_join(
      conversion_df %>% dplyr::select(feature_name, entrez_id),
      by = c("gene" = "feature_name")
    )

  # Save full marker table
  marker_output_path <- file.path(output_dir, paste0(timepoint, "_markers.csv"))
  write.csv(markers_with_entrez, marker_output_path, row.names = FALSE)
  cat("Saved marker table to:", marker_output_path, "\n\n")

  # 4. Prepare gene lists for compareCluster ----
  cat("Preparing gene lists for GO enrichment...\n")

  # Remove markers without ENTREZ conversion
  markers_for_go <- markers_with_entrez %>% dplyr::filter(!is.na(entrez_id))

  cat("  Markers with valid ENTREZ IDs:", nrow(markers_for_go), "\n")
  cat("  Cell types with valid markers:", length(unique(markers_for_go$cluster)), "\n\n")

  # Split by cluster
  cluster_gene_list <- split(markers_for_go$entrez_id, markers_for_go$cluster)

  # Show gene counts per cluster
  cat("  Genes per cell type (for GO analysis):\n")
  gene_counts <- sapply(cluster_gene_list, length)
  print(gene_counts)
  cat("\n")

  # 5. GO Enrichment Analysis ----
  results <- list()

  # Biological Process
  cat("Running GO enrichment: Biological Process...\n")
  go_bp <- tryCatch({
    compareCluster(
      geneCluster = cluster_gene_list,
      fun = "enrichGO",
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )
  }, error = function(e) {
    warning("BP enrichment failed: ", e$message)
    return(NULL)
  })
  results[["BP"]] <- go_bp

  if (!is.null(go_bp)) {
    cat("  Enriched BP terms:", nrow(go_bp@compareClusterResult), "\n")
  }
  cat("\n")

  # Molecular Function
  cat("Running GO enrichment: Molecular Function...\n")
  go_mf <- tryCatch({
    compareCluster(
      geneCluster = cluster_gene_list,
      fun = "enrichGO",
      OrgDb = org.Hs.eg.db,
      ont = "MF",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )
  }, error = function(e) {
    warning("MF enrichment failed: ", e$message)
    return(NULL)
  })
  results[["MF"]] <- go_mf

  if (!is.null(go_mf)) {
    cat("  Enriched MF terms:", nrow(go_mf@compareClusterResult), "\n")
  }
  cat("\n")

  # Cellular Component
  cat("Running GO enrichment: Cellular Component...\n")
  go_cc <- tryCatch({
    compareCluster(
      geneCluster = cluster_gene_list,
      fun = "enrichGO",
      OrgDb = org.Hs.eg.db,
      ont = "CC",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )
  }, error = function(e) {
    warning("CC enrichment failed: ", e$message)
    return(NULL)
  })
  results[["CC"]] <- go_cc

  if (!is.null(go_cc)) {
    cat("  Enriched CC terms:", nrow(go_cc@compareClusterResult), "\n")
  }
  cat("\n")

  # 6. Generate visualizations ----
  cat("Generating dotplots...\n")

  ontologies <- c("BP", "MF", "CC")
  ontology_names <- c(
    "BP" = "Biological Process",
    "MF" = "Molecular Function",
    "CC" = "Cellular Component"
  )

  for (ont in ontologies) {
    go_result <- results[[ont]]

    if (is.null(go_result) || nrow(go_result@compareClusterResult) == 0) {
      cat("  Skipping", ont, "- no enrichment results\n")
      next
    }

    # Create dotplot
    plot_title <- paste0(timepoint, " - GO Enrichment (", ontology_names[ont], ")")

    p <- dotplot(go_result, showCategory = 10, title = plot_title) +
      theme(
        axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1, size = 9),
        axis.text.y = element_text(size = 5, angle = 10, hjust = 1)
      )

    # Save as PDF
    pdf_path <- file.path(output_dir, paste0(timepoint, "_", ont, "_dotplot.pdf"))
    pdf(pdf_path, width = 12, height = 10)
    print(p)
    dev.off()

    cat("  Saved", ont, "dotplot to:", pdf_path, "\n")

    # Save RDS object
    rds_path <- file.path(output_dir, paste0(timepoint, "_", ont, "_results.rds"))
    saveRDS(go_result, rds_path)
    cat("  Saved", ont, "results to:", rds_path, "\n")
  }

  cat("\nCompleted timepoint:", timepoint, "\n")

  return(results)
}

# Main Analysis ----

cat("\n")
cat("================================================================================\n")
cat("  GO Enrichment Analysis for Consensus-Annotated Organoid Objects\n")
cat("  Author: Manveer Chauhan\n")
cat("  Date:", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("================================================================================\n")

# Create output directory
output_dir <- "./output_files/GO_analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("\nCreated output directory:", output_dir, "\n")
}

# Define timepoints
timepoints <- c("1M_Org", "3M_Org", "6M_Org")

# Store all results
all_results <- list()

# Process each timepoint
for (tp in timepoints) {
  result <- analyze_timepoint(tp, output_dir)

  if (!is.null(result)) {
    all_results[[tp]] <- result
  } else {
    cat("WARNING: Skipping", tp, "due to no marker genes found\n\n")
    all_results[[tp]] <- list(BP = NULL, MF = NULL, CC = NULL)
  }
}

# Generate summary statistics ----
cat("\n")
cat("========================================\n")
cat("Summary Statistics\n")
cat("========================================\n\n")

summary_data <- list()

for (tp in timepoints) {
  cat("Timepoint:", tp, "\n")

  for (ont in c("BP", "MF", "CC")) {
    result <- all_results[[tp]][[ont]]

    if (!is.null(result) && nrow(result@compareClusterResult) > 0) {
      n_terms <- nrow(result@compareClusterResult)
      n_clusters <- length(unique(result@compareClusterResult$Cluster))

      cat("  ", ont, "- Enriched terms:", n_terms,
          "| Cell types with enrichment:", n_clusters, "\n")

      summary_data[[paste0(tp, "_", ont)]] <- list(
        timepoint = tp,
        ontology = ont,
        n_terms = n_terms,
        n_clusters = n_clusters
      )
    } else {
      cat("  ", ont, "- No enrichment\n")

      summary_data[[paste0(tp, "_", ont)]] <- list(
        timepoint = tp,
        ontology = ont,
        n_terms = 0,
        n_clusters = 0
      )
    }
  }
  cat("\n")
}

# Save summary table
summary_df <- do.call(rbind, lapply(summary_data, as.data.frame))
summary_path <- file.path(output_dir, "GO_analysis_summary.csv")
write.csv(summary_df, summary_path, row.names = FALSE)

cat("Saved summary statistics to:", summary_path, "\n")

cat("\n")
cat("================================================================================\n")
cat("  Analysis Complete!\n")
cat("  Output directory:", output_dir, "\n")
cat("================================================================================\n\n")

cat("Generated files:\n")
cat("  - [timepoint]_markers.csv: Marker genes with ENTREZ conversion\n")
cat("  - [timepoint]_[BP/MF/CC]_dotplot.pdf: GO enrichment dotplots\n")
cat("  - [timepoint]_[BP/MF/CC]_results.rds: Raw GO enrichment objects\n")
cat("  - GO_analysis_summary.csv: Summary statistics table\n\n")
