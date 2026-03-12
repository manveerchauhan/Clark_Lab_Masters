#!/usr/bin/env Rscript
# NRXN/NLGN Isoform Dotplot Generator for Consensus Annotated Objects
# Author: Manveer Chauhan
# Description: Generate dotplots showing all isoforms of NRXN1/2/3 and NLGN1 genes
#              across consensus cell types for each developmental timepoint.
#              Uses consensus-annotated integrated objects with isoform assays.

library(Seurat)
library(tidyverse)
library(patchwork)
library(grid)
library(gridExtra)

setwd("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space")

# Source the isoform plotting functions for helper functions
source("../isoform_expression_plotter.R")

cat("\n=== Loading Consensus Annotated Seurat Objects ===\n")

# Load integrated consensus objects with isoform assays
one_month_org <- readRDS("./output_files/integrated_objects/1M_Org_integrated_harmony_consensus_with_isoform_assay.rds")
three_month_org <- readRDS("./output_files/integrated_objects/3M_Org_integrated_harmony_consensus_with_isoform_assay.rds")
six_month_org <- readRDS("./output_files/integrated_objects/6M_Org_integrated_harmony_consensus_with_isoform_assay.rds")

# Set identities to consensus cell types
Idents(one_month_org) <- "consensus_cell_type"
Idents(three_month_org) <- "consensus_cell_type"
Idents(six_month_org) <- "consensus_cell_type"

cat("Seurat objects loaded successfully\n")
cat("  1M_Org:", ncol(one_month_org), "cells,", nrow(one_month_org[["iso"]]), "isoforms\n")
cat("  3M_Org:", ncol(three_month_org), "cells,", nrow(three_month_org[["iso"]]), "isoforms\n")
cat("  6M_Org:", ncol(six_month_org), "cells,", nrow(six_month_org[["iso"]]), "isoforms\n")

# Create named list of Seurat objects
seurat_timepoints <- list(
  "1M_Org" = one_month_org,
  "3M_Org" = three_month_org,
  "6M_Org" = six_month_org
)

cat("\nConsensus cell types present:\n")
for(tp in names(seurat_timepoints)) {
  cell_types <- unique(seurat_timepoints[[tp]]$consensus_cell_type)
  cat("  ", tp, ":", length(cell_types), "cell types\n")
}

# Create output directory
output_dir <- "./output_files/nrxn_nlgn_isoform_dotplots"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ===== Function to generate dotplot for a gene across all timepoints =====
generate_gene_dotplot_report <- function(seurat_list, gene_symbol, output_filename) {

  cat("\n--- Generating dotplot report for", gene_symbol, "---\n")

  # Initialize PDF
  pdf(output_filename, width = 14, height = 10)

  # Loop through each timepoint
  for(tp_name in names(seurat_list)) {

    seu_obj <- seurat_list[[tp_name]]

    # Get all isoforms for this gene
    isoforms <- get_isoforms_for_gene(seu_obj, gene_symbol)

    if(length(isoforms) == 0) {
      cat("  ", tp_name, ": No isoforms found for", gene_symbol, "\n")

      # Create empty plot with message
      plot.new()
      text(0.5, 0.5, paste0("No isoforms found for ", gene_symbol, "\nin ", tp_name),
           cex = 1.5, col = "red")

      next
    }

    cat("  ", tp_name, ": Found", length(isoforms), "isoforms\n")

    # Get cell type information
    n_cells <- ncol(seu_obj)
    n_celltypes <- length(unique(seu_obj$consensus_cell_type))

    # Create dotplot
    tryCatch({

      # Create DotPlot with isoform assay
      p <- DotPlot(seu_obj,
                   features = isoforms,
                   assay = "iso",
                   cols = c("lightgrey", "red"),
                   dot.scale = 8) +
        coord_flip() +
        theme_classic() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 8),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 11, hjust = 0.5),
          legend.position = "right"
        ) +
        labs(
          title = paste0(gene_symbol, " Isoform Expression - ", tp_name),
          subtitle = paste0("n = ", length(isoforms), " isoforms | ",
                           n_cells, " cells | ",
                           n_celltypes, " cell types"),
          x = "Isoforms",
          y = "Consensus Cell Types"
        )

      # Print plot
      print(p)

    }, error = function(e) {
      cat("    Error creating dotplot for", gene_symbol, "in", tp_name, ":", e$message, "\n")

      # Create error plot
      plot.new()
      text(0.5, 0.5, paste0("Error creating dotplot for ", gene_symbol, "\n", e$message),
           cex = 1, col = "red")
    })
  }

  # Close PDF
  dev.off()

  cat("Saved dotplot report to:", output_filename, "\n")

  return(output_filename)
}

# ===== Generate dotplot reports for each gene =====
cat("\n=== Generating Dotplot Reports ===\n")

# NRXN1 (Neurexin 1)
cat("\n[1/4] Processing NRXN1...\n")
generate_gene_dotplot_report(
  seurat_list = seurat_timepoints,
  gene_symbol = "NRXN1",
  output_filename = file.path(output_dir, "NRXN1_isoform_dotplot.pdf")
)

# NRXN2 (Neurexin 2)
cat("\n[2/4] Processing NRXN2...\n")
generate_gene_dotplot_report(
  seurat_list = seurat_timepoints,
  gene_symbol = "NRXN2",
  output_filename = file.path(output_dir, "NRXN2_isoform_dotplot.pdf")
)

# NRXN3 (Neurexin 3)
cat("\n[3/4] Processing NRXN3...\n")
generate_gene_dotplot_report(
  seurat_list = seurat_timepoints,
  gene_symbol = "NRXN3",
  output_filename = file.path(output_dir, "NRXN3_isoform_dotplot.pdf")
)

# NLGN1 (Neuroligin 1)
cat("\n[4/4] Processing NLGN1...\n")
generate_gene_dotplot_report(
  seurat_list = seurat_timepoints,
  gene_symbol = "NLGN1",
  output_filename = file.path(output_dir, "NLGN1_isoform_dotplot.pdf")
)

# ===== Summary =====
cat("\n=== Dotplot Generation Complete ===\n")
cat("All dotplot reports have been generated successfully!\n")
cat("Reports saved to:", output_dir, "\n")
cat("\nGenerated dotplot PDFs:\n")
cat("  - NRXN1_isoform_dotplot.pdf (Neurexin 1)\n")
cat("  - NRXN2_isoform_dotplot.pdf (Neurexin 2)\n")
cat("  - NRXN3_isoform_dotplot.pdf (Neurexin 3)\n")
cat("  - NLGN1_isoform_dotplot.pdf (Neuroligin 1)\n")
cat("\nEach PDF contains:\n")
cat("  - 3 pages (one per timepoint: 1M, 3M, 6M)\n")
cat("  - Dotplot showing all isoforms (Y-axis) vs consensus cell types (X-axis)\n")
cat("  - Dot size = % cells expressing\n")
cat("  - Dot color = average expression level\n")
cat("\nPlots use:\n")
cat("  - Consensus cell type annotations (consensus_cell_type)\n")
cat("  - Isoform-level expression data (iso assay)\n")
