#!/usr/bin/env Rscript
# Isoform Expression Plotter for Consensus Annotated Objects
# Author: Manveer Chauhan
# Description: Generate isoform expression reports using the updated consensus-annotated
#              integrated objects with isoform assays. Plots isoform-level expression
#              on gene-level harmony UMAP with consensus cell type annotations.

library(Seurat)
library(tidyverse)
library(patchwork)

setwd("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space")

# Source the isoform plotting functions (now updated for consensus cell types)
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
dir.create("./output_files/isoform_expression_reports", recursive = TRUE, showWarnings = FALSE)

# ===== Example 1: Generate NRXN1 isoform report =====
cat("\n=== Example 1: Generating NRXN1 Isoform Report ===\n")

nrxn1_report <- generate_isoform_report(
  seurat_list = seurat_timepoints,
  gene_symbol = "NRXN1",
  output_filename = "./output_files/isoform_expression_reports/NRXN1_Consensus_Isoform_Report.pdf"
)

cat("\nNRXN1 report generated:", nrxn1_report, "\n")

# ===== Example 2: Generate reports for synaptic adhesion genes =====
cat("\n=== Example 2: Generating Reports for Synaptic Adhesion Genes ===\n")

# PTBP family (splicing factors)
generate_isoform_report(seurat_timepoints, "PTBP1",
                       "./output_files/isoform_expression_reports/PTBP1_Consensus_Isoform_Report.pdf")
generate_isoform_report(seurat_timepoints, "PTBP2",
                       "./output_files/isoform_expression_reports/PTBP2_Consensus_Isoform_Report.pdf")

# NRXN family (presynaptic adhesion)
generate_isoform_report(seurat_timepoints, "NRXN2",
                       "./output_files/isoform_expression_reports/NRXN2_Consensus_Isoform_Report.pdf")
generate_isoform_report(seurat_timepoints, "NRXN3",
                       "./output_files/isoform_expression_reports/NRXN3_Consensus_Isoform_Report.pdf")

# NLGN family (postsynaptic adhesion)
generate_isoform_report(seurat_timepoints, "NLGN1",
                       "./output_files/isoform_expression_reports/NLGN1_Consensus_Isoform_Report.pdf")
generate_isoform_report(seurat_timepoints, "NLGN2",
                       "./output_files/isoform_expression_reports/NLGN2_Consensus_Isoform_Report.pdf")
generate_isoform_report(seurat_timepoints, "NLGN3",
                       "./output_files/isoform_expression_reports/NLGN3_Consensus_Isoform_Report.pdf")
generate_isoform_report(seurat_timepoints, "NLGN4X",
                       "./output_files/isoform_expression_reports/NLGN4X_Consensus_Isoform_Report.pdf")
generate_isoform_report(seurat_timepoints, "NLGN4Y",
                       "./output_files/isoform_expression_reports/NLGN4Y_Consensus_Isoform_Report.pdf")


# ===== Example 3: Using individual plotting functions =====
cat("\n=== Example 3: Using Individual Plotting Functions ===\n")

# Get all NRXN1 isoforms from 6M organoids
nrxn1_isoforms_6M <- get_isoforms_for_gene(six_month_org, "NRXN1")
cat("Found", length(nrxn1_isoforms_6M), "NRXN1 isoforms in 6M_Org:\n")
cat(paste(nrxn1_isoforms_6M, collapse = "\n"), "\n\n")

# Create individual plots for the first isoform (if available)
if(length(nrxn1_isoforms_6M) > 0) {
  first_isoform <- nrxn1_isoforms_6M[1]
  cat("Creating example UMAP plot for:", first_isoform, "\n")

  # UMAP plot with custom gradient (plotted on gene-level harmony UMAP)
  p_umap <- plot_isoform_umap(six_month_org, first_isoform, "6M_Org")

  # Save plot
  ggsave("./output_files/isoform_expression_reports/Example_Individual_Isoform_UMAP_Consensus.png",
         p_umap, width = 8, height = 6, dpi = 300)

  cat("Saved example UMAP plot to: ./output_files/isoform_expression_reports/Example_Individual_Isoform_UMAP_Consensus.png\n")
}

cat("\n=== Consensus Isoform Plotter Complete ===\n")
cat("All isoform reports have been generated successfully!\n")
cat("Reports saved to: ./output_files/isoform_expression_reports/\n")
cat("\nReports generated:\n")
cat("  - PTBP1, PTBP2 (splicing factors)\n")
cat("  - NRXN1, NRXN2, NRXN3 (presynaptic adhesion)\n")
cat("  - NLGN1, NLGN2, NLGN3, NLGN4X, NLGN4Y (postsynaptic adhesion)\n")
cat("\nAll plots use:\n")
cat("  - Gene-level harmony UMAP (umap.harmony reduction)\n")
cat("  - Consensus cell type annotations\n")
cat("  - Isoform-level expression data (iso assay)\n")
