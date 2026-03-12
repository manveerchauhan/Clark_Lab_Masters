## Consensus Cell Type Annotation Script
# Author: Manveer Chauhan
# This script adds expert-curated consensus cell type annotations to integrated Seurat objects
# based on cluster validation results and marker expression analysis

# Load essential packages
library(Seurat)
library(tidyverse)
library(gridExtra)
library(grid)

# Set working directory
setwd("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space")

# Create output directory
dir.create("./output_files/consensus_annotation", recursive = TRUE, showWarnings = FALSE)

cat("\n========================================\n")
cat("CONSENSUS CELL TYPE ANNOTATION\n")
cat("========================================\n\n")

## CONSENSUS ANNOTATION MAPPINGS------
# Expert-curated cell type assignments based on:
# 1. Cluster validation against ground truth annotations
# 2. Marker gene expression analysis
# 3. Cell type proportion analysis across developmental timepoints

# 1M_Org consensus annotations (13 clusters)
consensus_1M <- c(
  "0" = "Apical Radial Glia (aRG)",
  "1" = "Newborn Deep Layer Projection Neurons",
  "2" = "Preplate Neurons",
  "3" = "Cajal-Retzius Cells",
  "4" = "Deep Layer Projection Neurons",
  "5" = "Radial Glia (Mixed)",
  "6" = "Apical Radial Glia (aRG) - Cycling",
  "7" = "Preplate Neurons",
  "8" = "Cortical Hem Cells",
  "9" = "Apical Radial Glia (aRG) - Cycling",
  "10" = "Subcortical Progenitors",
  "11" = "Radial Glia (Mixed)",
  "12" = "Preplate Neurons"
)

# 3M_Org consensus annotations (10 clusters)
consensus_3M <- c(
  "0" = "Corticothalamic Projection Neurons (Deep Layer)",
  "1" = "Callosal Projection Neurons (Upper Layer)",
  "2" = "Corticofugal Projection Neurons (Deep Layer)",
  "3" = "Callosal Projection Neurons (Upper Layer)",
  "4" = "Intermediate Progenitors (IP)",
  "5" = "Glutamatergic Neurons (Unspecified Layer)",
  "6" = "Outer Radial Glia (oRG)",
  "7" = "Apical Radial Glia (aRG) - Cycling",
  "8" = "Cortical Hem Cells",
  "9" = "Glutamatergic Neurons (Unspecified Layer)"
)

# 6M_Org consensus annotations (11 clusters)
consensus_6M <- c(
  "0" = "oRG to Astrocytes (Transitioning)",
  "1" = "Interneurons",
  "2" = "oRG to Oligodendrocyte Progenitors (Transitioning)",
  "3" = "Interneurons",
  "4" = "Callosal Projection Neurons (Upper Layer)",
  "5" = "Interneurons",
  "6" = "Glutamatergic Neurons",
  "7" = "Oligodendrocyte Progenitors",
  "8" = "Intermediate Progenitors (IP)",
  "9" = "Outer Radial Glia (oRG) - Cycling",
  "10" = "Intermediate Progenitors (IP)"
)

## FUNCTION 1: Apply consensus annotations------
apply_consensus_annotations <- function(seurat_obj, consensus_mapping, timepoint_name) {

  cat("\n=== Applying consensus annotations for", timepoint_name, "===\n")

  # Get cluster assignments
  clusters <- as.character(seurat_obj$seurat_clusters)

  # Map clusters to consensus cell types
  consensus_labels <- consensus_mapping[clusters]

  # Remove names to avoid Seurat metadata assignment issues
  consensus_labels <- unname(consensus_labels)

  # Add to metadata
  seurat_obj$consensus_cell_type <- consensus_labels

  # Print summary
  cat("  Total cells annotated:", length(consensus_labels), "\n")
  cat("  Number of unique cell types:", length(unique(consensus_labels)), "\n")

  # Print cell type distribution
  cell_type_counts <- table(consensus_labels)
  cat("\n  Cell type distribution:\n")
  for(ct in names(sort(cell_type_counts, decreasing = TRUE))) {
    cat("    ", ct, ": ", cell_type_counts[ct], " cells\n", sep = "")
  }

  return(seurat_obj)
}

## FUNCTION 2: Create visualization PDF------
create_consensus_visualization <- function(seurat_obj, timepoint_name, output_dir) {

  cat("\n=== Creating visualization for", timepoint_name, "===\n")

  # Prepare output filename
  pdf_file <- file.path(output_dir, paste0(timepoint_name, "_consensus_annotation.pdf"))

  # Open PDF device
  pdf(pdf_file, width = 14, height = 10)

  # Page 1: UMAP by cluster
  cat("  Generating Page 1: UMAP by cluster...\n")
  p1 <- DimPlot(seurat_obj,
                reduction = "umap.harmony",
                group.by = "seurat_clusters",
                label = TRUE,
                label.size = 5,
                repel = TRUE) +
    ggtitle(paste0(timepoint_name, " - Seurat Clusters")) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = "right"
    )
  print(p1)

  # Page 2: UMAP by consensus cell type
  cat("  Generating Page 2: UMAP by consensus cell type...\n")
  p2 <- DimPlot(seurat_obj,
                reduction = "umap.harmony",
                group.by = "consensus_cell_type",
                label = TRUE,
                label.size = 3.5,
                repel = TRUE) +
    ggtitle(paste0(timepoint_name, " - Consensus Cell Type Annotations")) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.text = element_text(size = 9)
    )
  print(p2)

  # Page 3: Summary table
  cat("  Generating Page 3: Summary table...\n")

  # Create summary data frame
  summary_df <- seurat_obj@meta.data %>%
    group_by(seurat_clusters, consensus_cell_type) %>%
    summarise(n_cells = n(), .groups = 'drop') %>%
    arrange(as.numeric(as.character(seurat_clusters)))

  # Calculate percentages
  summary_df <- summary_df %>%
    mutate(percentage = round(100 * n_cells / sum(n_cells), 2))

  # Rename columns for display
  summary_df <- summary_df %>%
    rename(
      "Cluster" = seurat_clusters,
      "Consensus Cell Type" = consensus_cell_type,
      "Cell Count" = n_cells,
      "Percentage" = percentage
    )

  # Create table plot
  table_grob <- tableGrob(
    summary_df,
    rows = NULL,
    theme = ttheme_default(
      base_size = 10,
      core = list(
        fg_params = list(hjust = 0, x = 0.05),
        bg_params = list(fill = c(rep(c("white", "grey95"), length.out = nrow(summary_df))))
      ),
      colhead = list(
        fg_params = list(fontface = "bold"),
        bg_params = list(fill = "steelblue", alpha = 0.7)
      )
    )
  )

  # Create title
  title <- textGrob(
    paste0(timepoint_name, " - Consensus Annotation Summary"),
    gp = gpar(fontsize = 16, fontface = "bold")
  )

  # Add subtitle with total cells
  subtitle <- textGrob(
    paste0("Total cells: ", sum(summary_df$`Cell Count`),
           " | Unique cell types: ", nrow(summary_df)),
    gp = gpar(fontsize = 12)
  )

  # Arrange table with titles
  grid.newpage()
  grid.draw(
    arrangeGrob(
      title,
      subtitle,
      table_grob,
      ncol = 1,
      heights = unit(c(0.5, 0.3, 9.2), "inches")
    )
  )

  # Close PDF device
  dev.off()

  cat("  PDF saved to:", pdf_file, "\n")

  return(pdf_file)
}

## FUNCTION 3: Save annotated object------
save_annotated_object <- function(seurat_obj, timepoint_name, output_dir) {

  cat("\n=== Saving annotated object for", timepoint_name, "===\n")

  # Prepare output filename
  rds_file <- file.path(output_dir, paste0(timepoint_name, "_integrated_harmony_consensus.rds"))

  # Save object
  saveRDS(seurat_obj, file = rds_file)

  cat("  RDS file saved to:", rds_file, "\n")
  cat("  File size:", round(file.size(rds_file) / 1e9, 2), "GB\n")

  return(rds_file)
}

####################################################
## MAIN EXECUTION
####################################################

cat("\n\n####################################################\n")
cat("STARTING CONSENSUS ANNOTATION\n")
cat("####################################################\n\n")

# Define timepoint mappings
timepoint_config <- list(
  "1M_Org" = list(
    input = "./output_files/integrated_objects/1M_Org_integrated_harmony.rds",
    consensus = consensus_1M
  ),
  "3M_Org" = list(
    input = "./output_files/integrated_objects/3M_Org_integrated_harmony.rds",
    consensus = consensus_3M
  ),
  "6M_Org" = list(
    input = "./output_files/integrated_objects/6M_Org_integrated_harmony.rds",
    consensus = consensus_6M
  )
)

# Process each timepoint
annotation_results <- list()

for(timepoint_name in names(timepoint_config)) {

  cat("\n\n####################################################\n")
  cat("PROCESSING TIMEPOINT:", timepoint_name, "\n")
  cat("####################################################\n")

  config <- timepoint_config[[timepoint_name]]

  # Step 1: Load integrated object
  cat("\nLoading integrated object...\n")
  seu_obj <- readRDS(config$input)
  cat("  Loaded:", ncol(seu_obj), "cells\n")
  cat("  Number of clusters:", length(unique(seu_obj$seurat_clusters)), "\n")

  # Step 2: Apply consensus annotations
  seu_obj <- apply_consensus_annotations(seu_obj, config$consensus, timepoint_name)

  # Step 3: Create visualization PDF
  pdf_file <- create_consensus_visualization(
    seu_obj,
    timepoint_name,
    "./output_files/consensus_annotation"
  )

  # Step 4: Save annotated object
  rds_file <- save_annotated_object(
    seu_obj,
    timepoint_name,
    "./output_files/integrated_objects"
  )

  # Store results
  annotation_results[[timepoint_name]] <- list(
    n_cells = ncol(seu_obj),
    n_clusters = length(unique(seu_obj$seurat_clusters)),
    n_cell_types = length(unique(seu_obj$consensus_cell_type)),
    pdf_file = pdf_file,
    rds_file = rds_file
  )

  cat("\n", timepoint_name, "processing complete!\n")
}

####################################################
## FINAL SUMMARY
####################################################

cat("\n\n####################################################\n")
cat("CONSENSUS ANNOTATION COMPLETE\n")
cat("####################################################\n\n")

cat("Summary of annotated objects:\n\n")
for(tp in names(annotation_results)) {
  res <- annotation_results[[tp]]
  cat(tp, ":\n")
  cat("  Cells:", res$n_cells, "\n")
  cat("  Clusters:", res$n_clusters, "\n")
  cat("  Cell types:", res$n_cell_types, "\n")
  cat("  PDF:", basename(res$pdf_file), "\n")
  cat("  RDS:", basename(res$rds_file), "\n\n")
}

cat("All outputs saved to:\n")
cat("  - ./output_files/consensus_annotation/\n")
cat("  - ./output_files/integrated_objects/\n\n")

cat("Script completed successfully!\n")
