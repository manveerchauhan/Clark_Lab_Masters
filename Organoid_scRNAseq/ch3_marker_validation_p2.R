## Marker Gene Validation Script - Dotplot Version
# Author: Manveer Chauhan
# This script creates developmental stage-specific marker gene dotplots
# for consensus-annotated integrated organoid timepoints

# Load essential packages
library(Seurat)
library(tidyverse)

# Set working directory
setwd("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space")

# Create output directory
dir.create("./output_files/marker_validation", recursive = TRUE, showWarnings = FALSE)

cat("\n========================================\n")
cat("MARKER GENE DOTPLOT ANALYSIS\n")
cat("Developmental Stage-Specific Markers\n")
cat("========================================\n\n")

## DEFINE TIMEPOINT-SPECIFIC MARKER GENE LISTS ------
# Flat structure: one list per timepoint

markers_by_timepoint <- list(
  "1M_Org" = c("PAX6", "SOX2", "VIM", "RELN", "WNT3A", "WNT2B",
               "TBR1", "CTIP2", "EOMES", "HOPX"),

  "3M_Org" = c("PAX6", "SOX9", "VIM", "SATB2", "WNT3A", "WNT2B",
               "CTIP2", "FOG2", "NMDA", "AMPA", "VGLUT1", "EOMES", "HOPX"),

  "6M_Org" = c("SATB2", "NMDA", "AMPA", "VGLUT1", "EOMES", "VIM",
               "PV", "SST", "GAD1", "PDGFRA", "OLIG1", "OLIG2",
               "GFAP", "S100B", "HOPX")
)

# Define gene aliases for common synonyms
gene_aliases <- list(
  "VGLUT1" = "SLC17A7",
  "PSD95" = "DLG4",
  "CTIP2" = "BCL11B",
  "FOG2" = "ZFPM2",
  "NMDA" = "GRIN1",
  "AMPA" = "GRIA2"
)

cat("Defined timepoint-specific marker gene lists:\n")
for(tp in names(markers_by_timepoint)) {
  n_genes <- length(markers_by_timepoint[[tp]])
  cat("  ", tp, ": ", n_genes, " genes\n", sep = "")
}
cat("\n")

## FUNCTION 1: Resolve gene symbol to full feature name with aliases ------
resolve_gene_name_with_aliases <- function(gene_symbol, seurat_obj, verbose = FALSE) {

  # First try direct match
  pattern <- paste0('-', gene_symbol, '$')
  matches <- grep(pattern, rownames(seurat_obj), value = TRUE, ignore.case = TRUE)

  if(length(matches) > 0) {
    if(length(matches) == 1) {
      if(verbose) cat("    ", gene_symbol, " → ", matches[1], "\n", sep = "")
      return(matches[1])
    } else {
      warning("Multiple matches for ", gene_symbol, ": ", paste(matches, collapse = ", "))
      if(verbose) cat("    ", gene_symbol, " → ", matches[1], " (multiple matches)\n", sep = "")
      return(matches[1])
    }
  }

  # Try alias if no direct match
  if(gene_symbol %in% names(gene_aliases)) {
    alias <- gene_aliases[[gene_symbol]]
    if(verbose) cat("    Trying alias: ", gene_symbol, " → ", alias, "\n", sep = "")
    pattern_alias <- paste0('-', alias, '$')
    matches_alias <- grep(pattern_alias, rownames(seurat_obj), value = TRUE, ignore.case = TRUE)

    if(length(matches_alias) > 0) {
      if(verbose) cat("    ", gene_symbol, " (", alias, ") → ", matches_alias[1], "\n", sep = "")
      return(matches_alias[1])
    }
  }

  # No matches found
  if(verbose) cat("    ", gene_symbol, " → NOT FOUND\n", sep = "")
  return(NA)
}

## FUNCTION 2: Validate gene set against Seurat object ------
validate_gene_set <- function(gene_symbols, seurat_obj, timepoint_name, verbose = TRUE) {

  if(verbose) cat("\nValidating", length(gene_symbols), "genes for", timepoint_name, "\n")

  validation_results <- data.frame(
    gene_symbol = gene_symbols,
    full_name = NA_character_,
    found = FALSE,
    stringsAsFactors = FALSE
  )

  for(i in seq_along(gene_symbols)) {
    full_name <- resolve_gene_name_with_aliases(gene_symbols[i], seurat_obj, verbose = FALSE)
    validation_results$full_name[i] <- full_name
    validation_results$found[i] <- !is.na(full_name)
  }

  n_found <- sum(validation_results$found)
  n_missing <- sum(!validation_results$found)
  pct_found <- round(100 * n_found / length(gene_symbols), 1)

  if(verbose) {
    cat("  Found:", n_found, "/", length(gene_symbols), "(", pct_found, "%)\n")

    if(n_missing > 0) {
      missing_genes <- validation_results$gene_symbol[!validation_results$found]
      cat("  Missing:", paste(missing_genes, collapse = ", "), "\n")
    }
  }

  return(validation_results)
}

## FUNCTION 3: Create marker dotplot ------
create_marker_dotplot <- function(seurat_obj, marker_genes_symbols,
                                  validation_results, timepoint_name) {

  cat("\nCreating dotplot for", timepoint_name, "\n")

  # Get valid genes
  valid_genes <- validation_results$full_name[validation_results$found]
  valid_symbols <- validation_results$gene_symbol[validation_results$found]

  if(length(valid_genes) == 0) {
    cat("  WARNING: No valid genes found - creating empty plot\n")
    empty_plot <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste0(timepoint_name, "\n\nNo valid genes found"),
               size = 6, color = "gray50") +
      theme_void()
    return(empty_plot)
  }

  cat("  Plotting", length(valid_genes), "genes\n")

  # Create gene name mapping for cleaner Y-axis labels
  # Extract just the gene symbol part from ENSG-SYMBOL format
  gene_labels <- setNames(
    sapply(valid_genes, function(x) gsub(".*-", "", x)),
    valid_genes
  )

  # Create DotPlot
  dotplot <- DotPlot(seurat_obj,
                     features = valid_genes,
                     group.by = "consensus_cell_type",
                     dot.scale = 8) +
    coord_flip() +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 11, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9)
    ) +
    scale_color_gradient(low = "lightgrey", high = "#e31837", name = "Avg\nExpression") +
    scale_y_discrete(labels = gene_labels) +
    labs(
      title = paste0(timepoint_name, " - Developmental Marker Expression"),
      subtitle = paste0("Genes: ", paste(valid_symbols, collapse = ", ")),
      x = "Genes",
      y = "Consensus Cell Type"
    )

  return(dotplot)
}

## FUNCTION 4: Process one timepoint ------
process_timepoint <- function(timepoint_name, file_path, marker_genes) {

  cat("\n\n####################################################\n")
  cat("PROCESSING TIMEPOINT:", timepoint_name, "\n")
  cat("####################################################\n")

  # Load consensus-annotated Seurat object
  cat("\nLoading consensus-annotated object from:", file_path, "\n")
  seurat_obj <- readRDS(file_path)
  cat("  Loaded:", ncol(seurat_obj), "cells,", nrow(seurat_obj), "features\n")
  cat("  Consensus cell types:", length(unique(seurat_obj$consensus_cell_type)), "\n")

  # Join layers for Seurat v5 compatibility
  cat("  Joining layers for v5 compatibility...\n")
  seurat_obj <- JoinLayers(seurat_obj)
  cat("  Layers joined successfully\n")

  # Validate marker genes
  validation_results <- validate_gene_set(
    marker_genes,
    seurat_obj,
    timepoint_name,
    verbose = TRUE
  )

  # Create dotplot
  dotplot <- create_marker_dotplot(
    seurat_obj,
    marker_genes,
    validation_results,
    timepoint_name
  )

  # Save as single-page PDF
  pdf_file <- paste0("./output_files/marker_validation/",
                     timepoint_name, "_marker_dotplot.pdf")

  cat("\nSaving dotplot to:", pdf_file, "\n")
  pdf(pdf_file, width = 12, height = 8)
  print(dotplot)
  dev.off()

  cat("\n====================================\n")
  cat("COMPLETED:", timepoint_name, "\n")
  cat("====================================\n")

  return(list(
    validation_results = validation_results,
    dotplot = dotplot,
    pdf_file = pdf_file
  ))
}

## MAIN EXECUTION: Process all timepoints ------
cat("\n\n####################################################\n")
cat("STARTING MARKER DOTPLOT GENERATION\n")
cat("####################################################\n\n")

# Define timepoint file paths (consensus-annotated objects)
timepoint_files <- list(
  "1M_Org" = "./output_files/integrated_objects/1M_Org_integrated_harmony_consensus.rds",
  "3M_Org" = "./output_files/integrated_objects/3M_Org_integrated_harmony_consensus.rds",
  "6M_Org" = "./output_files/integrated_objects/6M_Org_integrated_harmony_consensus.rds"
)

# Process each timepoint
results <- list()

for(timepoint_name in names(timepoint_files)) {
  # Get timepoint-specific marker list
  marker_genes <- markers_by_timepoint[[timepoint_name]]

  results[[timepoint_name]] <- process_timepoint(
    timepoint_name,
    timepoint_files[[timepoint_name]],
    marker_genes
  )

  # Force garbage collection
  gc()
}

cat("\n\n####################################################\n")
cat("ALL DOTPLOTS COMPLETE\n")
cat("####################################################\n")
cat("\nOutputs saved in: ./output_files/marker_validation/\n")
cat("\nGenerated dotplot PDFs:\n")
for(tp in names(results)) {
  cat("  ", tp, ": ", basename(results[[tp]]$pdf_file), "\n", sep = "")

  # Print validation summary
  n_tested <- nrow(results[[tp]]$validation_results)
  n_found <- sum(results[[tp]]$validation_results$found)
  cat("    Genes: ", n_found, "/", n_tested, " found\n", sep = "")
}
cat("\nScript completed successfully!\n")
