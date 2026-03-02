## Cluster Validation Script
# Author: Manveer Chauhan
# This script compares optimal resolution clusters from new integrated Seurat objects
# against ground truth manual annotations from the original analysis

# Load essential packages
library(Seurat)
library(tidyverse)
library(gridExtra)
library(grid)

# Set working directory
setwd("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space")

# Create output directory
dir.create("./output_files/cluster_validation", recursive = TRUE, showWarnings = FALSE)

cat("\n========================================\n")
cat("CLUSTER VALIDATION ANALYSIS\n")
cat("========================================\n\n")

## CUSTOM ANNOTATION MAPPING------
# Map original detailed annotations to simplified, interpretable categories for thesis
custom_annotation_mapping <- list(
  # Radial Glia
  "apical radial glia (aRG)" = "Apical Radial Glia (aRG)",
  "aRG (cycling progenitors)" = "Apical Radial Glia (aRG) - Cycling",
  "aRG/Outer radial glia" = "Radial Glia (Mixed)",
  "Outer radial glia (oRG)" = "Outer Radial Glia (oRG)",
  "oRG (cycling progenitors)" = "Outer Radial Glia (oRG) - Cycling",

  # Transitioning/Hybrid Glia
  "oRG/Astroglia" = "oRG to Astrocytes (Transitioning)",
  "oRG/Oligo progenitors (cycling)" = "oRG to Oligodendrocyte Progenitors (Transitioning)",

  # Progenitors
  "Intermediate Progenitors (IP)" = "Intermediate Progenitors (IP)",
  "IP" = "Intermediate Progenitors (IP)",
  "Subcortical progenitors" = "Subcortical Progenitors",

  # Glial Cells
  "Astroglia" = "Astrocytes",
  "Oligo progenitors" = "Oligodendrocyte Progenitors",

  # Projection Neurons - Upper Layer
  "Callosal projection neurons (CPN)" = "Callosal Projection Neurons (Upper Layer)",
  "CPN" = "Callosal Projection Neurons (Upper Layer)",

  # Projection Neurons - Deep/Subcortical
  "CFuPN" = "Corticofugal Projection Neurons (Deep Layer)",
  "CFuPN (likely SCPN)" = "Subcortical Projection Neurons (Deep Layer)",
  "Corticofugal projection neurons (CFuPN) (likely CthPN)" = "Corticothalamic Projection Neurons (Deep Layer)",
  "Deep Layer projection neurons (DL PN)" = "Deep Layer Projection Neurons",
  "Newborn DL PN" = "Deep Layer Projection Neurons - Newborn",

  # Other Projection Neurons
  "Projection neurons (PN)" = "Glutamatergic Neurons (Unspecified Layer)",

  # Interneurons
  "Interneurons (IN)" = "Interneurons",

  # Early Developmental
  "Preplate/Subplate neurons" = "Preplate/Subplate Neurons",
  "Cajal Retzius Cells" = "Cajal-Retzius Cells",

  # Regional Identity
  "Cortical hem" = "Cortical Hem",
  "cortical hem" = "Cortical Hem",

  # Unknown
  "unknown" = "Unknown"
)

cat("Defined custom annotation mapping with", length(custom_annotation_mapping), "categories\n\n")

## FUNCTION 0: Apply custom annotation mapping------
apply_custom_annotations <- function(gt_obj, mapping) {
  # Apply custom annotation mapping to ground truth object

  cat("  Applying custom annotation mapping...\n")

  original_annotations <- gt_obj$cluster_annotations

  # Map each annotation
  custom_annotations <- sapply(original_annotations, function(ann) {
    if(ann %in% names(mapping)) {
      return(mapping[[ann]])
    } else {
      warning("Unmapped annotation found: '", ann, "' - keeping original")
      return(ann)  # Keep original if not in mapping
    }
  })

  # Add as new metadata column
  gt_obj$cluster_annotations_custom <- unname(custom_annotations)

  # Report mapping statistics
  n_original_unique <- length(unique(original_annotations))
  n_custom_unique <- length(unique(custom_annotations))

  cat("  Original annotations:", n_original_unique, "unique categories\n")
  cat("  Custom annotations:", n_custom_unique, "unique categories\n")

  return(gt_obj)
}

## FUNCTION 1: Visualize ground truth UMAPs------
visualize_ground_truth_umaps <- function(gt_obj, timepoint_name) {

  cat("\n=== Creating ground truth UMAP visualizations for", timepoint_name, "===\n")

  # Check for UMAP reduction
  if(!"umap.harm" %in% names(gt_obj@reductions)) {
    cat("  Warning: 'umap.harm' reduction not found in ground truth object\n")
    cat("  Available reductions:", paste(names(gt_obj@reductions), collapse = ", "), "\n")
    return(NULL)
  }

  # Create UMAP for original annotations
  cat("  Creating UMAP for original annotations...\n")
  p_original <- DimPlot(
    gt_obj,
    reduction = "umap.harm",
    group.by = "cluster_annotations",
    label = FALSE
  ) +
    ggtitle(paste0(timepoint_name, " - Ground Truth (Original Annotations)")) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    ) +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))

  # Create UMAP for custom annotations
  cat("  Creating UMAP for custom annotations...\n")
  p_custom <- DimPlot(
    gt_obj,
    reduction = "umap.harm",
    group.by = "cluster_annotations_custom",
    label = FALSE
  ) +
    ggtitle(paste0(timepoint_name, " - Ground Truth (Simplified Annotations)")) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    ) +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))

  # Save individual PDFs
  pdf_original <- paste0("./output_files/cluster_validation/",
                         timepoint_name, "_ground_truth_UMAP_original.pdf")
  ggsave(pdf_original, p_original, width = 12, height = 8)
  cat("  Saved original annotations UMAP:", pdf_original, "\n")

  pdf_custom <- paste0("./output_files/cluster_validation/",
                       timepoint_name, "_ground_truth_UMAP_custom.pdf")
  ggsave(pdf_custom, p_custom, width = 12, height = 8)
  cat("  Saved custom annotations UMAP:", pdf_custom, "\n")

  # Also create a combined PDF with both
  pdf_combined <- paste0("./output_files/cluster_validation/",
                         timepoint_name, "_ground_truth_UMAPs_combined.pdf")
  pdf(pdf_combined, width = 14, height = 10)

  grid.arrange(
    p_original, p_custom, nrow = 2,
    top = textGrob(
      paste0(timepoint_name, " - Ground Truth Annotations Comparison"),
      gp = gpar(fontsize = 14, fontface = "bold")
    )
  )

  dev.off()
  cat("  Saved combined UMAP PDF:", pdf_combined, "\n")

  return(list(original = p_original, custom = p_custom))
}

## FUNCTION 2: Analyze barcode overlap between new and ground truth objects------
analyze_barcode_overlap <- function(new_obj, gt_obj, timepoint_name) {

  cat("\n=== Barcode overlap analysis for", timepoint_name, "===\n")

  # Get normalized barcodes from both objects
  # NEW object: strip ALL leading underscores from cell names
  new_barcodes_raw <- colnames(new_obj)
  new_barcodes_clean <- sub("^_+", "", new_barcodes_raw)  # Strip one or more underscores

  # GROUND TRUTH: use Barcode metadata column (already clean)
  gt_barcodes_raw <- colnames(gt_obj)
  gt_barcodes_clean <- gt_obj$Barcode

  cat("  NEW object barcode format: ", new_barcodes_raw[1], " -> ", new_barcodes_clean[1], "\n", sep="")
  cat("  GT object barcode format: ", gt_barcodes_raw[1], " -> ", gt_barcodes_clean[1], "\n", sep="")

  # Find overlapping barcodes using cleaned versions
  overlap_barcodes_clean <- intersect(new_barcodes_clean, gt_barcodes_clean)

  # Create mapping back to original cell names
  # For NEW object: map clean barcode -> original cell name
  new_barcode_map <- setNames(new_barcodes_raw, new_barcodes_clean)
  # For GT object: map clean barcode -> original cell name
  gt_barcode_map <- setNames(gt_barcodes_raw, gt_barcodes_clean)

  # Get original cell names for overlapping barcodes
  overlap_cells_new <- new_barcode_map[overlap_barcodes_clean]
  overlap_cells_gt <- gt_barcode_map[overlap_barcodes_clean]

  # Calculate statistics
  n_new <- length(new_barcodes_raw)
  n_gt <- length(gt_barcodes_raw)
  n_overlap <- length(overlap_barcodes_clean)

  pct_new <- round(100 * n_overlap / n_new, 2)
  pct_gt <- round(100 * n_overlap / n_gt, 2)

  # Print statistics
  cat("  New integrated object:", n_new, "cells\n")
  cat("  Ground truth object:", n_gt, "cells\n")
  cat("  Overlapping cells:", n_overlap, "\n")
  cat("  Overlap %% of new:", pct_new, "%%\n")
  cat("  Overlap %% of ground truth:", pct_gt, "%%\n")

  # Return results
  return(list(
    overlap_cells_new = overlap_cells_new,
    overlap_cells_gt = overlap_cells_gt,
    overlap_barcodes_clean = overlap_barcodes_clean,
    stats = data.frame(
      timepoint = timepoint_name,
      n_new = n_new,
      n_gt = n_gt,
      n_overlap = n_overlap,
      pct_new = pct_new,
      pct_gt = pct_gt
    )
  ))
}

## FUNCTION 3: Create confusion matrix comparing clusters to annotations------
create_confusion_matrix <- function(new_obj, gt_obj, overlap_result, timepoint_name,
                                   annotation_column = "cluster_annotations") {

  cat("\n=== Creating confusion matrix for", timepoint_name, "(",annotation_column, ") ===\n")

  # Extract cluster assignments for overlapping cells using original cell names
  # Use [row, column] indexing to properly match cell names to metadata rows
  new_clusters <- new_obj@meta.data[overlap_result$overlap_cells_new, "seurat_clusters"]
  gt_annotations <- gt_obj@meta.data[overlap_result$overlap_cells_gt, annotation_column]

  # Check for missing data
  n_missing_new <- sum(is.na(new_clusters))
  n_missing_gt <- sum(is.na(gt_annotations))

  if(n_missing_new > 0) {
    cat("  Warning:", n_missing_new, "cells missing cluster assignments in new object\n")
  }
  if(n_missing_gt > 0) {
    cat("  Warning:", n_missing_gt, "cells missing annotations in ground truth\n")
  }

  # Remove cells with missing data
  valid_cells <- !is.na(new_clusters) & !is.na(gt_annotations)
  new_clusters <- new_clusters[valid_cells]
  gt_annotations <- gt_annotations[valid_cells]

  cat("  Valid cells for comparison:", length(new_clusters), "\n")
  cat("  Number of new clusters:", length(unique(new_clusters)), "\n")
  cat("  Number of ground truth annotations:", length(unique(gt_annotations)), "\n")

  # Create confusion matrix (clusters = rows, annotations = columns)
  confusion_mat <- table(
    Cluster = new_clusters,
    Annotation = gt_annotations
  )

  # Convert to data frame for easier manipulation
  confusion_df <- as.data.frame.matrix(confusion_mat)

  # Calculate row percentages (what % of each cluster maps to each annotation)
  confusion_pct <- confusion_df / rowSums(confusion_df) * 100

  # Print summary
  cat("  Confusion matrix dimensions:", nrow(confusion_df), "clusters x",
      ncol(confusion_df), "annotations\n")

  return(list(
    confusion_mat = confusion_mat,
    confusion_df = confusion_df,
    confusion_pct = confusion_pct,
    new_clusters = new_clusters,
    gt_annotations = gt_annotations
  ))
}

## FUNCTION 4: Generate distribution summary------
generate_distribution_summary <- function(confusion_result, timepoint_name) {

  cat("\n=== Distribution summary for", timepoint_name, "===\n")

  # Check if confusion matrix has any annotations
  if(ncol(confusion_result$confusion_mat) == 0) {
    cat("  No annotations found - skipping distribution summary\n")
    return(data.frame(
      Annotation = character(0),
      total_cells = integer(0),
      top_clusters = character(0)
    ))
  }

  # For each ground truth annotation, show top 3 clusters it maps to
  confusion_long <- as.data.frame(confusion_result$confusion_mat) %>%
    filter(Freq > 0)

  if(nrow(confusion_long) == 0) {
    cat("  No overlapping cells - skipping distribution summary\n")
    return(data.frame(
      Annotation = character(0),
      total_cells = integer(0),
      top_clusters = character(0)
    ))
  }

  summary_list <- confusion_long %>%
    group_by(Annotation) %>%
    slice_head(n = 3) %>%
    mutate(pct = round(100 * Freq / sum(confusion_result$confusion_mat[, as.character(Annotation[1])]), 1)) %>%
    summarize(
      total_cells = sum(confusion_result$confusion_mat[, as.character(Annotation[1])]),
      top_clusters = paste0(
        paste(Cluster, " (", pct, "%)", sep = ""),
        collapse = ", "
      ),
      .groups = 'drop'
    )

  cat("\nTop cluster assignments per annotation:\n")
  print(summary_list, n = Inf)

  return(summary_list)
}

## FUNCTION 5: Visualize confusion matrix------
visualize_confusion_matrix <- function(confusion_result, timepoint_name,
                                      annotation_type = "Original") {

  cat("\n=== Creating confusion matrix plot for", timepoint_name, "(", annotation_type, ") ===\n")

  # Check if confusion matrix is empty
  if(ncol(confusion_result$confusion_mat) == 0 || nrow(confusion_result$confusion_mat) == 0) {
    cat("  Empty confusion matrix - creating placeholder plots\n")

    # Create empty placeholder plots
    empty_plot <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = "No overlapping cells\nfor validation",
               size = 6, color = "gray50") +
      theme_void() +
      labs(title = paste0(timepoint_name, ": No Data Available"))

    return(list(count_plot = empty_plot, pct_plot = empty_plot))
  }

  # Convert confusion matrix to long format for ggplot
  confusion_long <- as.data.frame(confusion_result$confusion_mat) %>%
    mutate(
      Cluster = factor(Cluster),
      Annotation = factor(Annotation)
    )

  # Create heatmap
  p <- ggplot(confusion_long, aes(x = Annotation, y = Cluster, fill = Freq)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Freq), size = 3) +
    scale_fill_gradient(low = "white", high = "steelblue",
                        name = "Cell Count") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      panel.grid = element_blank()
    ) +
    labs(
      title = paste0(timepoint_name, ": Cluster vs Ground Truth Annotation"),
      subtitle = paste0("Annotation Type: ", annotation_type),
      x = "Ground Truth Annotation",
      y = "New Cluster (Optimal Resolution)"
    )

  # Also create a percentage-based heatmap
  confusion_pct_long <- as.data.frame(confusion_result$confusion_mat) %>%
    group_by(Cluster) %>%
    mutate(
      pct = round(100 * Freq / sum(Freq), 1),
      Cluster = factor(Cluster),
      Annotation = factor(Annotation)
    ) %>%
    ungroup()

  p_pct <- ggplot(confusion_pct_long, aes(x = Annotation, y = Cluster, fill = pct)) +
    geom_tile(color = "white") +
    geom_text(aes(label = paste0(pct, "%")), size = 3) +
    scale_fill_gradient(low = "white", high = "darkorange",
                        name = "% of Cluster") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      panel.grid = element_blank()
    ) +
    labs(
      title = paste0(timepoint_name, ": Row-Normalized (% of Each Cluster)"),
      subtitle = paste0("Annotation Type: ", annotation_type),
      x = "Ground Truth Annotation",
      y = "New Cluster (Optimal Resolution)"
    )

  return(list(count_plot = p, pct_plot = p_pct))
}

## FUNCTION 6: Save validation results (Enhanced for 2-page PDF)------
save_validation_results_enhanced <- function(overlap_result,
                                             confusion_result_original,
                                             confusion_result_custom,
                                             distribution_summary_original,
                                             distribution_summary_custom,
                                             plots_original,
                                             plots_custom,
                                             timepoint_name) {

  cat("\n=== Saving enhanced results for", timepoint_name, "===\n")

  # Save overlap statistics (once - same for both annotation types)
  overlap_file <- paste0("./output_files/cluster_validation/",
                         timepoint_name, "_overlap_stats.csv")
  write.csv(overlap_result$stats, overlap_file, row.names = FALSE)
  cat("  Saved overlap stats:", overlap_file, "\n")

  # Save confusion matrices - ORIGINAL annotations
  confusion_file_orig <- paste0("./output_files/cluster_validation/",
                                timepoint_name, "_confusion_matrix_original.csv")
  write.csv(confusion_result_original$confusion_df, confusion_file_orig, row.names = TRUE)
  cat("  Saved original confusion matrix:", confusion_file_orig, "\n")

  confusion_pct_file_orig <- paste0("./output_files/cluster_validation/",
                                    timepoint_name, "_confusion_matrix_original_pct.csv")
  write.csv(confusion_result_original$confusion_pct, confusion_pct_file_orig, row.names = TRUE)
  cat("  Saved original confusion matrix (%):", confusion_pct_file_orig, "\n")

  # Save confusion matrices - CUSTOM annotations
  confusion_file_custom <- paste0("./output_files/cluster_validation/",
                                  timepoint_name, "_confusion_matrix_custom.csv")
  write.csv(confusion_result_custom$confusion_df, confusion_file_custom, row.names = TRUE)
  cat("  Saved custom confusion matrix:", confusion_file_custom, "\n")

  confusion_pct_file_custom <- paste0("./output_files/cluster_validation/",
                                      timepoint_name, "_confusion_matrix_custom_pct.csv")
  write.csv(confusion_result_custom$confusion_pct, confusion_pct_file_custom, row.names = TRUE)
  cat("  Saved custom confusion matrix (%):", confusion_pct_file_custom, "\n")

  # Save distribution summaries - both types
  dist_file_orig <- paste0("./output_files/cluster_validation/",
                           timepoint_name, "_distribution_summary_original.csv")
  write.csv(distribution_summary_original, dist_file_orig, row.names = FALSE)
  cat("  Saved original distribution summary:", dist_file_orig, "\n")

  dist_file_custom <- paste0("./output_files/cluster_validation/",
                             timepoint_name, "_distribution_summary_custom.csv")
  write.csv(distribution_summary_custom, dist_file_custom, row.names = FALSE)
  cat("  Saved custom distribution summary:", dist_file_custom, "\n")

  # Create 2-page PDF
  plot_file <- paste0("./output_files/cluster_validation/",
                      timepoint_name, "_confusion_plots.pdf")
  pdf(plot_file, width = 12, height = 10)

  # Prepare shared subtitle
  stats <- overlap_result$stats
  subtitle_text <- sprintf(
    "New Object: %d cells  |  Ground Truth: %d cells  |  Overlapping: %d cells (%.1f%% / %.1f%%)",
    stats$n_new, stats$n_gt, stats$n_overlap, stats$pct_new, stats$pct_gt
  )

  # Page 1: Original annotations
  title_text_original <- paste(timepoint_name, "- Cluster Validation (Original Annotations)")

  grid.arrange(
    plots_original$count_plot, plots_original$pct_plot, nrow = 2,
    top = textGrob(
      paste(title_text_original, "\n", subtitle_text, sep = ""),
      gp = gpar(fontsize = 12, fontface = "bold")
    )
  )

  # Page 2: Custom simplified annotations
  title_text_custom <- paste(timepoint_name, "- Cluster Validation (Simplified Annotations)")

  grid.arrange(
    plots_custom$count_plot, plots_custom$pct_plot, nrow = 2,
    top = textGrob(
      paste(title_text_custom, "\n", subtitle_text, sep = ""),
      gp = gpar(fontsize = 12, fontface = "bold")
    )
  )

  dev.off()
  cat("  Saved 2-page PDF:", plot_file, "\n")
}

## MAIN EXECUTION: Process all timepoints------
cat("\n\n####################################################\n")
cat("STARTING CLUSTER VALIDATION\n")
cat("####################################################\n\n")

# Define timepoint mappings
timepoint_mapping <- list(
  "1M_Org" = list(
    new = "./output_files/integrated_objects/1M_Org_integrated_harmony.rds",
    gt = "/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/one_month_seurat.intergrated_harm.isofrom.rds"
  ),
  "3M_Org" = list(
    new = "./output_files/integrated_objects/3M_Org_integrated_harmony.rds",
    gt = "/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/three_month_seurat.intergrated_harm.isofrom.rds"
  ),
  "6M_Org" = list(
    new = "./output_files/integrated_objects/6M_Org_integrated_harmony.rds",
    gt = "/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/six_month_seurat.intergrated_harm.isofrom.rds"
  )
)

# Process each timepoint
validation_results <- list()

for(timepoint_name in names(timepoint_mapping)) {

  cat("\n\n####################################################\n")
  cat("PROCESSING TIMEPOINT:", timepoint_name, "\n")
  cat("####################################################\n")

  paths <- timepoint_mapping[[timepoint_name]]

  # Load objects
  cat("\nLoading objects...\n")
  new_obj <- readRDS(paths$new)
  cat("  Loaded new integrated object:", ncol(new_obj), "cells\n")

  gt_obj <- readRDS(paths$gt)
  cat("  Loaded ground truth object:", ncol(gt_obj), "cells\n")

  # Apply custom annotations to ground truth
  gt_obj <- apply_custom_annotations(gt_obj, custom_annotation_mapping)

  # Step 0: Visualize ground truth UMAPs
  gt_umap_plots <- visualize_ground_truth_umaps(gt_obj, timepoint_name)

  # Step 1: Analyze barcode overlap (once)
  overlap_result <- analyze_barcode_overlap(new_obj, gt_obj, timepoint_name)

  # Step 2: Create confusion matrices for BOTH annotation types
  cat("\n--- Processing ORIGINAL annotations ---\n")
  confusion_result_original <- create_confusion_matrix(
    new_obj, gt_obj, overlap_result, timepoint_name,
    annotation_column = "cluster_annotations"
  )

  cat("\n--- Processing CUSTOM simplified annotations ---\n")
  confusion_result_custom <- create_confusion_matrix(
    new_obj, gt_obj, overlap_result, timepoint_name,
    annotation_column = "cluster_annotations_custom"
  )

  # Step 3: Generate distribution summaries for both
  distribution_summary_original <- generate_distribution_summary(
    confusion_result_original, paste0(timepoint_name, " (Original)")
  )

  distribution_summary_custom <- generate_distribution_summary(
    confusion_result_custom, paste0(timepoint_name, " (Custom)")
  )

  # Step 4: Visualize both
  plots_original <- visualize_confusion_matrix(
    confusion_result_original, timepoint_name, "Original"
  )

  plots_custom <- visualize_confusion_matrix(
    confusion_result_custom, timepoint_name, "Simplified"
  )

  # Step 5: Save enhanced results (2-page PDF + multiple CSVs)
  save_validation_results_enhanced(
    overlap_result,
    confusion_result_original,
    confusion_result_custom,
    distribution_summary_original,
    distribution_summary_custom,
    plots_original,
    plots_custom,
    timepoint_name
  )

  # Store results
  validation_results[[timepoint_name]] <- list(
    gt_umap_plots = gt_umap_plots,
    overlap = overlap_result,
    confusion_original = confusion_result_original,
    confusion_custom = confusion_result_custom,
    distribution_original = distribution_summary_original,
    distribution_custom = distribution_summary_custom,
    plots_original = plots_original,
    plots_custom = plots_custom
  )

  cat("\n====================================\n")
  cat("VALIDATION COMPLETE:", timepoint_name, "\n")
  cat("====================================\n")

  # Force garbage collection
  gc()
}

cat("\n\n####################################################\n")
cat("ALL VALIDATIONS COMPLETE\n")
cat("####################################################\n")

# Print final summary
cat("\n=== VALIDATION SUMMARY ===\n")
for(timepoint_name in names(validation_results)) {
  result <- validation_results[[timepoint_name]]
  cat("\nTimepoint:", timepoint_name, "\n")
  cat("  Overlapping cells:", result$overlap$stats$n_overlap, "\n")
  cat("  Number of clusters:", nrow(result$confusion_original$confusion_df), "\n")
  cat("  Original annotations:", ncol(result$confusion_original$confusion_df), "categories\n")
  cat("  Custom annotations:", ncol(result$confusion_custom$confusion_df), "categories\n")
}

cat("\n=== ALL VALIDATION COMPLETE ===\n")
cat("Results saved in: ./output_files/cluster_validation/\n")
cat("  - Ground truth UMAP PDFs (Original, Custom, and Combined)\n")
cat("  - 2-page confusion matrix PDFs per timepoint (Original + Simplified annotations)\n")
cat("  - Separate CSVs for both annotation types\n")
