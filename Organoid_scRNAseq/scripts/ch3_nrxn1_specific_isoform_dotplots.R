#!/usr/bin/env Rscript
# NRXN1 Specific Isoform Dotplot Generator for Consensus Annotated Objects
# Author: Manveer Chauhan
# Description: Generate dotplots showing specific NRXN1 isoforms
#              across consensus cell types for each developmental timepoint.
#              Uses consensus-annotated integrated objects with isoform assays.

library(Seurat)
library(tidyverse)
library(patchwork)
library(grid)
library(gridExtra)

setwd("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space")

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

# ===== Define specific NRXN1 isoforms of interest =====
nrxn1_isoforms_of_interest <- c(
  "ENST00000626899.1-NRXN1",
  "ENST00000657235.1-NRXN1",
  "ENST00000405581.3-NRXN1",
  "ENST00000628515.2-NRXN1",
  "ENST00000626249.2-NRXN1",
  "ENST00000405472.7-NRXN1",
  "ENST00000404971.5-NRXN1",
  "ENST00000406316.6-NRXN1",
  "ENST00000625672.2-NRXN1",
  "ENST00000401669.7-NRXN1",
  "ENST00000474354.1-NRXN1",
  "ENST00000628364.2-NRXN1",
  "ENST00000342183.9-NRXN1"
)

cat("\n=== NRXN1 Isoforms of Interest ===\n")
cat("Specified", length(nrxn1_isoforms_of_interest), "isoforms:\n")
for(iso in nrxn1_isoforms_of_interest) {
  cat("  -", iso, "\n")
}

# ===== Function to find available isoforms in Seurat object =====
find_available_isoforms <- function(seu_obj, isoforms_list) {
  # Get all features in the isoform assay
  all_features <- rownames(seu_obj[["iso"]])

  # Find which specified isoforms are present
  available <- isoforms_list[isoforms_list %in% all_features]
  missing <- isoforms_list[!isoforms_list %in% all_features]

  return(list(available = available, missing = missing))
}

# ===== Function to generate dotplot for specific NRXN1 isoforms across all timepoints =====
generate_nrxn1_specific_dotplot_report <- function(seurat_list, isoforms_list, output_filename) {

  cat("\n--- Generating NRXN1 specific isoform dotplot report ---\n")

  # Initialize PDF
  pdf(output_filename, width = 14, height = 10)

  # Loop through each timepoint
  for(tp_name in names(seurat_list)) {

    seu_obj <- seurat_list[[tp_name]]

    # Find available isoforms
    isoform_check <- find_available_isoforms(seu_obj, isoforms_list)
    available_isoforms <- isoform_check$available
    missing_isoforms <- isoform_check$missing

    cat("  ", tp_name, ":\n")
    cat("    Available:", length(available_isoforms), "of", length(isoforms_list), "isoforms\n")

    if(length(missing_isoforms) > 0) {
      cat("    Missing isoforms:\n")
      for(m in missing_isoforms) {
        cat("      -", m, "\n")
      }
    }

    if(length(available_isoforms) == 0) {
      cat("    ERROR: No specified isoforms found in", tp_name, "\n")

      # Create empty plot with message
      plot.new()
      text(0.5, 0.5, paste0("No specified NRXN1 isoforms found\nin ", tp_name),
           cex = 1.5, col = "red")

      next
    }

    # Get cell type information
    n_cells <- ncol(seu_obj)
    n_celltypes <- length(unique(seu_obj$consensus_cell_type))

    # Create dotplot
    tryCatch({

      # Create DotPlot with isoform assay
      p <- DotPlot(seu_obj,
                   features = available_isoforms,
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
          title = paste0("NRXN1 Specific Isoform Expression - ", tp_name),
          subtitle = paste0("n = ", length(available_isoforms), "/", length(isoforms_list),
                           " isoforms present | ",
                           n_cells, " cells | ",
                           n_celltypes, " cell types"),
          x = "NRXN1 Isoforms",
          y = "Consensus Cell Types"
        )

      # Print plot
      print(p)

    }, error = function(e) {
      cat("    Error creating dotplot for", tp_name, ":", e$message, "\n")

      # Create error plot
      plot.new()
      text(0.5, 0.5, paste0("Error creating dotplot\n", e$message),
           cex = 1, col = "red")
    })
  }

  # Close PDF
  dev.off()

  cat("Saved dotplot report to:", output_filename, "\n")

  return(output_filename)
}

# ===== Generate the NRXN1 specific isoform dotplot report =====
cat("\n=== Generating NRXN1 Specific Isoform Dotplot Report ===\n")

generate_nrxn1_specific_dotplot_report(
  seurat_list = seurat_timepoints,
  isoforms_list = nrxn1_isoforms_of_interest,
  output_filename = file.path(output_dir, "NRXN1_specific_isoforms_dotplot.pdf")
)

# ===== Summary =====
cat("\n=== Dotplot Generation Complete ===\n")
cat("NRXN1 specific isoform dotplot report has been generated!\n")
cat("Report saved to:", file.path(output_dir, "NRXN1_specific_isoforms_dotplot.pdf"), "\n")
cat("\nGenerated PDF contains:\n")
cat("  - 3 pages (one per timepoint: 1M, 3M, 6M)\n")
cat("  - Dotplot showing specified NRXN1 isoforms (Y-axis) vs consensus cell types (X-axis)\n")
cat("  - Dot size = % cells expressing\n")
cat("  - Dot color = average expression level\n")
cat("\nSpecified NRXN1 isoforms (", length(nrxn1_isoforms_of_interest), " total):\n")
for(iso in nrxn1_isoforms_of_interest) {
  cat("  -", iso, "\n")
}
cat("\nPlots use:\n")
cat("  - Consensus cell type annotations (consensus_cell_type)\n")
cat("  - Isoform-level expression data (iso assay)\n")
