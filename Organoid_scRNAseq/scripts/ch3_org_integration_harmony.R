## Harmony Integration Script for Organoid Timepoints------
# Author: Manveer Chauhan
# This script integrates QC'd organoid samples by timepoint using Harmony

# Load essential packages
library(Seurat)
library(tidyverse)
library(harmony)
library(clustree)
library(gridExtra)
library(grid)

# Load cell cycle markers
load("/data/gpfs/projects/punim0646/manveer/cycle.rda")

# Load silhouette functions
source("/data/gpfs/projects/punim2251/scRNAseq_scripts/silhouette.R")

# Source QC functions (for quantitative_elbow and optimize_clustering)
source("./qc_functions.R")

# Set working directory
setwd("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space")

# Create output directories
dir.create("./output_files/integrated_objects", recursive = TRUE, showWarnings = FALSE)
dir.create("./output_files/integration_QC", recursive = TRUE, showWarnings = FALSE)

cat("\n========================================\n")
cat("HARMONY INTEGRATION PIPELINE\n")
cat("========================================\n\n")

## STEP 1: Load all QC'd Seurat objects------
cat("=== STEP 1: Loading QC'd Seurat objects ===\n")

# Define sample groupings by timepoint
sample_groups <- list(
  "1M_Org" = c("org_1A", "org_1B"),
  "3M_Org" = c("org_3A", "org_3B", "org_3C"),
  "6M_Org" = c("org_6A", "org_6B", "org_6C")
)

# Clustering resolution override (NA = use calculated optimal)
clustering_params <- list(
  "1M_Org" = 0.3,   # Custom override
  "3M_Org" = NA,    # Use calculated optimal
  "6M_Org" = 0.5    # Custom override
)

cat("\n=== Clustering Resolution Configuration ===\n")
cat("1M_Org: Resolution =", ifelse(is.na(clustering_params[["1M_Org"]]), "calculated", clustering_params[["1M_Org"]]), "\n")
cat("3M_Org: Resolution =", ifelse(is.na(clustering_params[["3M_Org"]]), "calculated", clustering_params[["3M_Org"]]), "\n")
cat("6M_Org: Resolution =", ifelse(is.na(clustering_params[["6M_Org"]]), "calculated", clustering_params[["6M_Org"]]), "\n\n")

# Load all Seurat objects
seurat_objects <- list()
for(timepoint in names(sample_groups)) {
  samples <- sample_groups[[timepoint]]
  cat("\nLoading samples for", timepoint, ":\n")

  for(sample_id in samples) {
    file_path <- paste0("./output_files/seu_objects/", sample_id, "_QCed_final.rds")
    seurat_objects[[sample_id]] <- readRDS(file_path)
    cat("  Loaded:", sample_id, "-", ncol(seurat_objects[[sample_id]]), "cells\n")
  }
}

cat("\n=== All objects loaded successfully ===\n")

## STEP 2: Merge samples by timepoint------
merge_samples_by_timepoint <- function(sample_ids, seurat_list, timepoint_name) {

  cat("\n=== Merging samples for", timepoint_name, "===\n")

  # Get first object
  merged_obj <- seurat_list[[sample_ids[1]]]

  # Add sample_id to metadata if not present
  if(!"sample_id" %in% colnames(merged_obj@meta.data)) {
    merged_obj$sample_id <- sample_ids[1]
  }

  # Merge with remaining samples
  if(length(sample_ids) > 1) {
    for(i in 2:length(sample_ids)) {
      obj <- seurat_list[[sample_ids[i]]]

      # Add sample_id to metadata
      if(!"sample_id" %in% colnames(obj@meta.data)) {
        obj$sample_id <- sample_ids[i]
      }

      merged_obj <- merge(merged_obj, y = obj,
                          add.cell.ids = c("", sample_ids[i]),
                          project = timepoint_name)
    }
  }

  # Add timepoint metadata
  merged_obj$timepoint <- timepoint_name

  cat("  Total cells after merge:", ncol(merged_obj), "\n")
  cat("  Samples:", paste(sample_ids, collapse = ", "), "\n")

  return(merged_obj)
}

## STEP 3: Pre-integration processing------
pre_integration_processing <- function(seurat_obj, timepoint_name) {

  cat("\n=== Pre-integration processing for", timepoint_name, "===\n")

  # Normalize data
  seurat_obj <- NormalizeData(seurat_obj)
  cat("  Normalization complete\n")

  # Find variable features
  seurat_obj <- FindVariableFeatures(seurat_obj,
                                     selection.method = 'vst',
                                     nfeatures = 2000)
  cat("  Found 2000 variable features\n")

  # Scale data (all features)
  all_features <- rownames(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, features = all_features)
  cat("  Data scaling complete\n")

  # Run PCA
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
  cat("  PCA complete\n")

  # Quantitative elbow analysis to determine optimal PCs (with +5 PC adjustment)
  elbow_result <- quantitative_elbow(seurat_obj, timepoint_name, max_pcs = 50, pc_adjustment = 5)
  optimal_pcs <- elbow_result$optimal_pcs_used  # Use the adjusted value

  return(list(
    seurat_obj = seurat_obj,
    optimal_pcs = optimal_pcs,
    elbow_plots = elbow_result$plots
  ))
}

## STEP 4: Harmony integration------
harmony_integration <- function(seurat_obj, timepoint_name, optimal_pcs) {

  cat("\n=== Running Harmony integration for", timepoint_name, "===\n")

  # Create before-integration UMAP for comparison
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:optimal_pcs,
                        reduction = "pca",
                        reduction.name = "umap.unintegrated")
  cat("  Created pre-integration UMAP\n")

  # Run Harmony integration
  seurat_obj <- IntegrateLayers(
    object = seurat_obj,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "harmony",
    group.by = "sample_id",
    verbose = TRUE
  )
  cat("  Harmony integration complete\n")

  # Run UMAP on harmony reduction
  seurat_obj <- RunUMAP(seurat_obj,
                        reduction = "harmony",
                        dims = 1:optimal_pcs,
                        reduction.name = "umap.harmony")
  cat("  Created post-integration UMAP\n")

  # Create comparison plots
  before_plot <- DimPlot(seurat_obj, reduction = "umap.unintegrated",
                         group.by = "sample_id") +
    ggtitle(paste0(timepoint_name, ": Before Integration")) +
    labs(color = "Sample")

  after_plot <- DimPlot(seurat_obj, reduction = "umap.harmony",
                        group.by = "sample_id") +
    ggtitle(paste0(timepoint_name, ": After Harmony Integration")) +
    labs(color = "Sample")

  return(list(
    seurat_obj = seurat_obj,
    plots = list(before = before_plot, after = after_plot)
  ))
}

## STEP 5: Post-integration clustering------
post_integration_clustering <- function(seurat_obj, timepoint_name, optimal_pcs, custom_resolution = NA) {

  cat("\n=== Post-integration clustering for", timepoint_name, "===\n")

  # Find neighbors using harmony reduction
  seurat_obj <- FindNeighbors(seurat_obj,
                               reduction = "harmony",
                               dims = 1:optimal_pcs)
  cat("  FindNeighbors complete (harmony reduction, dims 1:", optimal_pcs, ")\n")

  # Silhouette analysis for optimal resolution
  cat("  Running silhouette analysis...\n")
  sil_results <- optimize_silhouette(sobject = seurat_obj,
                                     test_res = seq(0.1, 1.2, by = 0.1),
                                     summary_plot = TRUE,
                                     reduction = "harmony")

  sil_df <- as.data.frame(sil_results) %>%
    dplyr::rename(avg_sil_vals = sil_vals) %>%
    group_by(num_clusters) %>%
    dplyr::slice(which.max(avg_sil_vals)) %>%
    ungroup() %>%
    arrange(desc(avg_sil_vals))

  # Get optimal resolution (calculated from silhouette)
  optimal_res <- sil_df$res_vals[1]
  cat("  Calculated optimal resolution:", optimal_res, "(silhouette:",
      round(sil_df$avg_sil_vals[1], 4), ")\n")

  # Check for custom resolution override
  if(!is.na(custom_resolution)) {
    final_res <- custom_resolution
    resolution_method <- "CUSTOM OVERRIDE"
    cat("  Using CUSTOM resolution:", final_res, "(calculated optimal was:", optimal_res, ")\n")
  } else {
    final_res <- optimal_res
    resolution_method <- "calculated"
    cat("  Using calculated optimal resolution:", final_res, "\n")
  }

  # Cluster with all test resolutions for clustree
  seurat_obj <- FindClusters(seurat_obj, resolution = seq(0.1, 1.2, by = 0.1))

  # Create clustree plot with resolution method subtitle
  clustree_plot <- clustree(seurat_obj) +
    ggtitle(paste0(timepoint_name, ": Clustree")) +
    labs(subtitle = paste0("Final resolution: ", final_res, " (", resolution_method, ")"))

  # Re-cluster with final resolution (either custom or calculated)
  seurat_obj <- FindClusters(seurat_obj, resolution = final_res)

  # Create final cluster UMAP with resolution method subtitle
  cluster_umap <- DimPlot(seurat_obj, reduction = "umap.harmony",
                          label = TRUE) +
    ggtitle(paste0(timepoint_name, ": Final Clusters")) +
    labs(subtitle = paste0("Resolution: ", final_res, " (", resolution_method, ")"),
         color = "Cluster")

  return(list(
    seurat_obj = seurat_obj,
    optimal_res = optimal_res,        # Calculated value
    final_res = final_res,            # Actually used value
    resolution_method = resolution_method,
    sil_results = sil_df,
    plots = list(clustree = clustree_plot, cluster_umap = cluster_umap)
  ))
}

## STEP 6: Generate visualization plots------
generate_qc_plots <- function(seurat_obj, timepoint_name) {

  cat("\n=== Generating QC plots for", timepoint_name, "===\n")

  # Cell cycle UMAP
  cellcycle_umap <- DimPlot(seurat_obj, reduction = "umap.harmony",
                            group.by = "Phase") +
    ggtitle(paste0(timepoint_name, ": Cell Cycle Phase")) +
    labs(color = "Phase")

  # nFeature UMAP
  nfeature_umap <- FeaturePlot(seurat_obj, reduction = "umap.harmony",
                               features = "nFeature_RNA") +
    ggtitle(paste0(timepoint_name, ": nFeature_RNA"))

  # nCount UMAP
  ncount_umap <- FeaturePlot(seurat_obj, reduction = "umap.harmony",
                             features = "nCount_RNA") +
    ggtitle(paste0(timepoint_name, ": nCount_RNA"))

  return(list(
    cellcycle = cellcycle_umap,
    nfeature = nfeature_umap,
    ncount = ncount_umap
  ))
}

## STEP 7: Generate comprehensive PDF report------
generate_integration_pdf <- function(integration_results, timepoint_name) {

  pdf_file <- paste0("./output_files/integration_QC/", timepoint_name,
                     "_integration_report.pdf")

  cat("\nGenerating PDF report:", pdf_file, "\n")

  pdf(pdf_file, width = 14, height = 10)

  # Page 1: Pre-integration Elbow Analysis
  grid.arrange(
    integration_results$pre_integration$elbow_plots$elbow_plot,
    integration_results$pre_integration$elbow_plots$quant_plot,
    nrow = 1, ncol = 2,
    top = textGrob(paste(timepoint_name, "- PCA Elbow Analysis (Pre-Integration)"),
                   gp = gpar(fontsize = 14, fontface = "bold"))
  )

  # Page 2: Before/After Integration
  grid.arrange(
    integration_results$harmony$plots$before,
    integration_results$harmony$plots$after,
    nrow = 1, ncol = 2,
    top = textGrob(paste(timepoint_name, "- Harmony Integration Comparison"),
                   gp = gpar(fontsize = 14, fontface = "bold"))
  )

  # Page 3: Clustree
  grid.arrange(
    integration_results$clustering$plots$clustree,
    nrow = 1, ncol = 1,
    top = textGrob(paste(timepoint_name, "- Clustering Resolution Tree"),
                   gp = gpar(fontsize = 14, fontface = "bold"))
  )

  # Page 4: Final Clusters & Silhouette Table
  grid.arrange(
    integration_results$clustering$plots$cluster_umap,
    tableGrob(head(integration_results$clustering$sil_results, 10)),
    nrow = 1, ncol = 2,
    top = textGrob(paste(timepoint_name, "- Final Clustering Results"),
                   gp = gpar(fontsize = 14, fontface = "bold"))
  )

  # Page 5: QC Metrics on Integrated UMAP
  grid.arrange(
    integration_results$qc_plots$cellcycle,
    integration_results$qc_plots$nfeature,
    integration_results$qc_plots$ncount,
    integration_results$clustering$plots$cluster_umap,
    nrow = 2, ncol = 2,
    top = textGrob(paste(timepoint_name, "- QC Metrics on Integrated Data"),
                   gp = gpar(fontsize = 14, fontface = "bold"))
  )

  dev.off()

  cat("PDF report saved:", pdf_file, "\n")
}

## STEP 8: Master integration function------
integrate_timepoint <- function(timepoint_name, sample_ids, seurat_list) {

  cat("\n\n####################################################\n")
  cat("INTEGRATING TIMEPOINT:", timepoint_name, "\n")
  cat("####################################################\n")

  # Step 2: Merge samples
  merged_obj <- merge_samples_by_timepoint(sample_ids, seurat_list, timepoint_name)

  # Step 3: Pre-integration processing
  pre_integration <- pre_integration_processing(merged_obj, timepoint_name)

  # Step 4: Harmony integration
  harmony_result <- harmony_integration(pre_integration$seurat_obj,
                                       timepoint_name,
                                       pre_integration$optimal_pcs)

  # Step 5: Post-integration clustering (with custom resolution if specified)
  clustering_result <- post_integration_clustering(harmony_result$seurat_obj,
                                                   timepoint_name,
                                                   pre_integration$optimal_pcs,
                                                   custom_resolution = clustering_params[[timepoint_name]])

  # Step 6: Generate QC plots
  qc_plots <- generate_qc_plots(clustering_result$seurat_obj, timepoint_name)

  # Combine all results
  results <- list(
    timepoint = timepoint_name,
    seurat_obj = clustering_result$seurat_obj,
    optimal_pcs = pre_integration$optimal_pcs,
    optimal_res = clustering_result$optimal_res,        # Calculated resolution
    final_res = clustering_result$final_res,            # Actually used resolution
    resolution_method = clustering_result$resolution_method,
    pre_integration = pre_integration,
    harmony = harmony_result,
    clustering = clustering_result,
    qc_plots = qc_plots
  )

  # Step 7: Generate PDF report
  generate_integration_pdf(results, timepoint_name)

  # Save integrated object
  output_file <- paste0("./output_files/integrated_objects/",
                        timepoint_name, "_integrated_harmony.rds")
  saveRDS(results$seurat_obj, file = output_file)
  cat("\nIntegrated object saved:", output_file, "\n")

  cat("\n====================================\n")
  cat("INTEGRATION COMPLETE:", timepoint_name, "\n")
  cat("====================================\n")

  return(results)
}

## STEP 9: Execute integration for all timepoints------
cat("\n\n####################################################\n")
cat("STARTING INTEGRATION FOR ALL TIMEPOINTS\n")
cat("####################################################\n\n")

integration_results <- list()

for(timepoint in names(sample_groups)) {
  integration_results[[timepoint]] <- integrate_timepoint(
    timepoint_name = timepoint,
    sample_ids = sample_groups[[timepoint]],
    seurat_list = seurat_objects
  )

  # Force garbage collection
  gc()
}

cat("\n\n####################################################\n")
cat("ALL TIMEPOINTS INTEGRATED SUCCESSFULLY\n")
cat("####################################################\n")

# Print summary
cat("\n=== INTEGRATION SUMMARY ===\n")
for(timepoint in names(integration_results)) {
  result <- integration_results[[timepoint]]
  cat("\nTimepoint:", timepoint, "\n")
  cat("  Samples:", paste(sample_groups[[timepoint]], collapse = ", "), "\n")
  cat("  Total cells:", ncol(result$seurat_obj), "\n")
  cat("  Optimal PCs:", result$optimal_pcs, "\n")
  cat("  Calculated optimal resolution:", result$optimal_res, "\n")
  cat("  Final resolution used:", result$final_res, "(", result$resolution_method, ")\n")
  cat("  Number of clusters:", length(unique(result$seurat_obj$seurat_clusters)), "\n")
}

cat("\n=== PER-TIMEPOINT INTEGRATION COMPLETE ===\n")
cat("Integrated objects saved in: ./output_files/integrated_objects/\n")
cat("PDF reports saved in: ./output_files/integration_QC/\n")

## STEP 10: All-Samples Cross-Timepoint Harmony Integration------
cat("\n\n####################################################\n")
cat("STEP 10: ALL-SAMPLES CROSS-TIMEPOINT INTEGRATION\n")
cat("####################################################\n\n")

# 10a: Define merge function for all samples across timepoints
merge_all_samples <- function(seurat_list, sample_groups) {

  cat("\n=== Merging all samples across timepoints ===\n")

  # Build list of objects with sample_id and timepoint metadata set per-object
  all_objects <- list()
  all_sample_ids <- c()

  for (timepoint in names(sample_groups)) {
    for (sample_id in sample_groups[[timepoint]]) {
      obj <- seurat_list[[sample_id]]
      obj$sample_id <- sample_id
      obj$timepoint <- timepoint
      all_objects[[sample_id]] <- obj
      all_sample_ids <- c(all_sample_ids, sample_id)
      cat("  Prepared:", sample_id, "(", timepoint, ") -", ncol(obj), "cells\n")
    }
  }

  # Merge: first object + list of remaining objects
  merged_obj <- merge(
    x = all_objects[[all_sample_ids[1]]],
    y = all_objects[all_sample_ids[-1]],
    add.cell.ids = all_sample_ids,
    project = "All_Samples"
  )

  cat("\n  Total cells after merge:", ncol(merged_obj), "\n")
  cat("  Cells per timepoint:\n")
  print(table(merged_obj$timepoint))
  cat("  Cells per sample:\n")
  print(table(merged_obj$sample_id))

  return(merged_obj)
}

# 10b: Define extended PDF function for all-samples integration (7 pages)
generate_all_samples_pdf <- function(integration_results, timepoint_name) {

  pdf_file <- paste0("./output_files/integration_QC/", timepoint_name,
                     "_integration_report.pdf")

  cat("\nGenerating PDF report:", pdf_file, "\n")

  pdf(pdf_file, width = 14, height = 10)

  # Page 1: Pre-integration Elbow Analysis
  grid.arrange(
    integration_results$pre_integration$elbow_plots$elbow_plot,
    integration_results$pre_integration$elbow_plots$quant_plot,
    nrow = 1, ncol = 2,
    top = textGrob(paste(timepoint_name, "- PCA Elbow Analysis (Pre-Integration)"),
                   gp = gpar(fontsize = 14, fontface = "bold"))
  )

  # Page 2: Before/After Integration
  grid.arrange(
    integration_results$harmony$plots$before,
    integration_results$harmony$plots$after,
    nrow = 1, ncol = 2,
    top = textGrob(paste(timepoint_name, "- Harmony Integration Comparison"),
                   gp = gpar(fontsize = 14, fontface = "bold"))
  )

  # Page 3: Clustree
  grid.arrange(
    integration_results$clustering$plots$clustree,
    nrow = 1, ncol = 1,
    top = textGrob(paste(timepoint_name, "- Clustering Resolution Tree"),
                   gp = gpar(fontsize = 14, fontface = "bold"))
  )

  # Page 4: Final Clusters & Silhouette Table
  grid.arrange(
    integration_results$clustering$plots$cluster_umap,
    tableGrob(head(integration_results$clustering$sil_results, 10)),
    nrow = 1, ncol = 2,
    top = textGrob(paste(timepoint_name, "- Final Clustering Results"),
                   gp = gpar(fontsize = 14, fontface = "bold"))
  )

  # Page 5: QC Metrics on Integrated UMAP
  grid.arrange(
    integration_results$qc_plots$cellcycle,
    integration_results$qc_plots$nfeature,
    integration_results$qc_plots$ncount,
    integration_results$clustering$plots$cluster_umap,
    nrow = 2, ncol = 2,
    top = textGrob(paste(timepoint_name, "- QC Metrics on Integrated Data"),
                   gp = gpar(fontsize = 14, fontface = "bold"))
  )

  # Page 6: PCA Colored by Batch (sample_id)
  grid.arrange(
    integration_results$pca_batch_plot,
    nrow = 1, ncol = 1,
    top = textGrob(paste(timepoint_name, "- PCA Colored by Batch (Post-Integration)"),
                   gp = gpar(fontsize = 14, fontface = "bold"))
  )

  # Page 7: Timepoint Distribution on UMAP
  grid.arrange(
    integration_results$timepoint_umap,
    nrow = 1, ncol = 1,
    top = textGrob(paste(timepoint_name, "- Timepoint Distribution on Harmony UMAP"),
                   gp = gpar(fontsize = 14, fontface = "bold"))
  )

  dev.off()

  cat("PDF report saved:", pdf_file, "\n")
}

# 10c: Free per-timepoint results from memory (already saved to disk)
cat("Freeing per-timepoint integration results from memory...\n")
rm(integration_results)
gc()

# 10d: Merge all 8 samples across timepoints
all_samples_merged <- merge_all_samples(seurat_objects, sample_groups)
rm(seurat_objects)
gc()

# 10e: Pre-integration processing (reuse existing function)
cat("\n=== Pre-integration processing for All_Samples ===\n")
all_pre_integration <- pre_integration_processing(all_samples_merged, "All_Samples")
rm(all_samples_merged)
gc()

# 10f: Harmony integration (reuse existing function, group.by = "sample_id")
all_harmony <- harmony_integration(
  all_pre_integration$seurat_obj,
  "All_Samples",
  all_pre_integration$optimal_pcs
)

# 10g: Post-integration clustering (auto resolution via silhouette)
all_clustering <- post_integration_clustering(
  all_harmony$seurat_obj,
  "All_Samples",
  all_pre_integration$optimal_pcs,
  custom_resolution = NA
)

# 10h: QC plots (reuse existing function)
all_qc_plots <- generate_qc_plots(all_clustering$seurat_obj, "All_Samples")

# 10i: PCA batch plot (PC1 vs PC2 colored by sample_id)
pca_batch_plot <- DimPlot(all_clustering$seurat_obj,
                          reduction = "pca",
                          dims = c(1, 2),
                          group.by = "sample_id") +
  ggtitle("All Samples: PCA Colored by Batch")

# 10j: Timepoint UMAP
timepoint_umap <- DimPlot(all_clustering$seurat_obj,
                           reduction = "umap.harmony",
                           group.by = "timepoint") +
  ggtitle("All Samples: Timepoint Distribution")

# 10k: Assemble results
all_samples_results <- list(
  timepoint = "All_Samples",
  seurat_obj = all_clustering$seurat_obj,
  optimal_pcs = all_pre_integration$optimal_pcs,
  optimal_res = all_clustering$optimal_res,
  final_res = all_clustering$final_res,
  resolution_method = all_clustering$resolution_method,
  pre_integration = all_pre_integration,
  harmony = all_harmony,
  clustering = all_clustering,
  qc_plots = all_qc_plots,
  pca_batch_plot = pca_batch_plot,
  timepoint_umap = timepoint_umap
)

# 10l: Generate 7-page PDF report
generate_all_samples_pdf(all_samples_results, "All_Samples")

# 10m: Save integrated object
all_samples_output <- "./output_files/integrated_objects/All_Samples_integrated_harmony.rds"
saveRDS(all_samples_results$seurat_obj, file = all_samples_output)
cat("\nAll-samples integrated object saved:", all_samples_output, "\n")

# 10n: Print summary
cat("\n====================================\n")
cat("ALL-SAMPLES INTEGRATION COMPLETE\n")
cat("====================================\n")
cat("  Total cells:", ncol(all_samples_results$seurat_obj), "\n")
cat("  Optimal PCs:", all_samples_results$optimal_pcs, "\n")
cat("  Calculated optimal resolution:", all_samples_results$optimal_res, "\n")
cat("  Final resolution used:", all_samples_results$final_res,
    "(", all_samples_results$resolution_method, ")\n")
cat("  Number of clusters:",
    length(unique(all_samples_results$seurat_obj$seurat_clusters)), "\n")
cat("\n  Cells per timepoint:\n")
print(table(all_samples_results$seurat_obj$timepoint))
cat("\n  Cells per sample:\n")
print(table(all_samples_results$seurat_obj$sample_id))

# 10o: Final cleanup
rm(all_pre_integration, all_harmony, all_clustering, all_qc_plots,
   pca_batch_plot, timepoint_umap, all_samples_results)
gc()

cat("\n=== ALL INTEGRATION COMPLETE (Per-Timepoint + All-Samples) ===\n")
cat("Integrated objects saved in: ./output_files/integrated_objects/\n")
cat("PDF reports saved in: ./output_files/integration_QC/\n")
