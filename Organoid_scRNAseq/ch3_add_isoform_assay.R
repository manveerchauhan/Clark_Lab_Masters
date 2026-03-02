## Add Isoform-Level Assay to Integrated Consensus Objects
# Author: Manveer Chauhan
# This script adds isoform-level transcriptome data as a separate assay to
# the integrated, consensus-annotated Seurat objects. Isoform data undergoes
# full QC processing including quantitative PC selection and silhouette-optimized
# clustering, enabling comparison with gene-level annotations.

# Load essential packages
library(Seurat)
library(tidyverse)
library(patchwork)
library(gridExtra)
library(grid)

# Set working directory
setwd("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space")

# Create output directory
dir.create("./output_files/isoform_integration", recursive = TRUE, showWarnings = FALSE)

# Source required functions
source("qc_functions.R")
source("/data/gpfs/projects/punim2251/scRNAseq_scripts/silhouette.R")

# Load cell cycle genes
load("/data/gpfs/projects/punim0646/manveer/cycle.rda")

cat("\n========================================\n")
cat("ISOFORM ASSAY INTEGRATION\n")
cat("========================================\n\n")

## FUNCTION 1: Load and merge isoform objects by timepoint------
load_and_merge_isoform_objects <- function(timepoint_name, sample_ids) {

  cat("\n=== Loading and merging isoform objects for", timepoint_name, "===\n")

  # Load individual isoform objects
  iso_objects <- list()

  for(sample_id in sample_ids) {
    file_path <- paste0("./output_files/seu_objects/", sample_id, "_isoform_seurat.rds")
    cat("  Loading", sample_id, "from:", file_path, "\n")

    iso_obj <- readRDS(file_path)
    cat("    Cells:", ncol(iso_obj), ", Features:", nrow(iso_obj), "\n")

    iso_objects[[sample_id]] <- iso_obj
  }

  # Merge objects with sample prefixes
  cat("\n  Merging", length(sample_ids), "isoform objects...\n")

  if(length(sample_ids) == 1) {
    # Single sample - just add prefix
    iso_merged <- iso_objects[[sample_ids[1]]]
    iso_merged <- RenameCells(iso_merged, add.cell.id = sample_ids[1])
  } else {
    # Multiple samples - merge with prefixes
    iso_merged <- merge(
      x = iso_objects[[sample_ids[1]]],
      y = iso_objects[sample_ids[-1]],
      add.cell.ids = sample_ids,
      project = paste0(timepoint_name, "_isoform")
    )
  }

  cat("    Merged cells:", ncol(iso_merged), "\n")
  cat("    Merged features:", nrow(iso_merged), "\n")

  # Add sample_id metadata
  cell_names <- colnames(iso_merged)
  sample_ids_extracted <- sapply(strsplit(cell_names, "_"), `[`, 1)
  iso_merged$sample_id <- sample_ids_extracted

  cat("    Sample distribution:\n")
  print(table(iso_merged$sample_id))

  return(iso_merged)
}

## FUNCTION 2: Filter isoform object to match integrated cells------
filter_isoform_to_integrated_cells <- function(iso_merged, integrated_obj) {

  cat("\n=== Filtering isoform object to match integrated cells ===\n")

  # Get integrated cell names and sample_id
  integrated_cells <- colnames(integrated_obj)
  integrated_sample_ids <- integrated_obj$sample_id
  cat("  Integrated object cells:", length(integrated_cells), "\n")
  cat("  Sample distribution:\n")
  print(table(integrated_sample_ids))

  # Get isoform cell names (with sample prefix: sample_BARCODE)
  iso_cells <- colnames(iso_merged)
  cat("\n  Isoform merged cells:", length(iso_cells), "\n")

  # Create expected isoform cell names based on integrated cells
  # Integrated cells have MANY inconsistent formats across samples:
  # Format 1: __BARCODE (double underscore)
  # Format 2: _BARCODE (single underscore)
  # Format 3: sample_id_BARCODE (clean sample prefix)
  # Format 4: _sample_id_BARCODE (underscore + sample prefix)
  # Expected isoform format: sample_id_BARCODE

  expected_iso_cells <- character(length(integrated_cells))

  for(i in seq_along(integrated_cells)) {
    cell_name <- integrated_cells[i]
    sample_id <- integrated_sample_ids[i]

    # Universal barcode extraction strategy:
    # Remove ALL prefix patterns, then reconstruct consistently
    barcode_clean <- cell_name

    # Step 1: Remove _sample_id_ pattern (e.g., _org_3B_BARCODE → BARCODE)
    barcode_clean <- gsub(paste0("^_", sample_id, "_"), "", barcode_clean)

    # Step 2: Remove sample_id_ pattern if present (e.g., org_3C_BARCODE → BARCODE)
    barcode_clean <- gsub(paste0("^", sample_id, "_"), "", barcode_clean)

    # Step 3: Remove any remaining leading underscores (single or multiple)
    barcode_clean <- gsub("^_+", "", barcode_clean)

    # Reconstruct as sample_id_BARCODE
    expected_iso_cells[i] <- paste0(sample_id, "_", barcode_clean)
  }

  cat("  Expected isoform cell names (first 5):\n")
  print(head(expected_iso_cells, 5))
  cat("  Actual isoform cell names (first 5):\n")
  print(head(iso_cells, 5))

  # Find overlapping cells
  overlapping_cells <- intersect(expected_iso_cells, iso_cells)
  cat("\n  Overlapping cells:", length(overlapping_cells), "\n")

  if(length(overlapping_cells) == 0) {
    stop("ERROR: No overlapping cells found! Check barcode formatting.")
  }

  # Calculate overlap statistics
  pct_integrated <- round(100 * length(overlapping_cells) / length(integrated_cells), 1)
  pct_isoform <- round(100 * length(overlapping_cells) / length(iso_cells), 1)

  cat("  Overlap statistics:\n")
  cat("    ", pct_integrated, "% of integrated cells found in isoform data\n", sep = "")
  cat("    ", pct_isoform, "% of isoform cells match integrated data\n", sep = "")

  # Subset isoform object to overlapping cells
  iso_filtered <- subset(iso_merged, cells = overlapping_cells)

  cat("  Filtered isoform object:\n")
  cat("    Cells:", ncol(iso_filtered), "\n")
  cat("    Features:", nrow(iso_filtered), "\n")

  # Create mapping between integrated cells and isoform cells
  cell_mapping <- data.frame(
    integrated_cell = integrated_cells[expected_iso_cells %in% overlapping_cells],
    iso_cell = expected_iso_cells[expected_iso_cells %in% overlapping_cells],
    stringsAsFactors = FALSE
  )

  return(list(
    iso_filtered = iso_filtered,
    overlapping_cells = overlapping_cells,
    cell_mapping = cell_mapping,
    stats = list(
      integrated_cells = length(integrated_cells),
      isoform_cells = length(iso_cells),
      overlapping = length(overlapping_cells),
      pct_integrated = pct_integrated,
      pct_isoform = pct_isoform
    )
  ))
}

## FUNCTION 3: Add isoform assay to integrated object------
add_isoform_assay_to_object <- function(integrated_obj, iso_filtered, cell_mapping) {

  cat("\n=== Adding isoform assay to integrated object ===\n")

  # Join layers for Seurat v5 compatibility
  cat("  Joining isoform layers...\n")
  iso_filtered <- JoinLayers(iso_filtered)

  # Extract counts matrix from isoform assay
  cat("  Extracting isoform counts matrix...\n")
  counts_table <- iso_filtered[["isoform"]]$counts

  cat("    Counts matrix dimensions:", nrow(counts_table), "x", ncol(counts_table), "\n")
  cat("    Cell mapping entries:", nrow(cell_mapping), "\n")

  # Reorder counts table to match integrated object using cell_mapping
  # cell_mapping has: integrated_cell (with _), iso_cell (with sample_id_)

  # First, subset integrated object to cells that have isoform data
  integrated_obj <- subset(integrated_obj, cells = cell_mapping$integrated_cell)
  cat("    Integrated object subsetted to", ncol(integrated_obj), "cells with isoform data\n")

  # Now reorder counts table to match integrated object order
  counts_table_ordered <- counts_table[, cell_mapping$iso_cell]

  # Rename columns to match integrated cell names
  colnames(counts_table_ordered) <- cell_mapping$integrated_cell

  cat("    Reordered counts matrix dimensions:", nrow(counts_table_ordered), "x",
      ncol(counts_table_ordered), "\n")

  # Verify column names match
  if(!all(colnames(counts_table_ordered) == colnames(integrated_obj))) {
    stop("ERROR: Cell names don't match after reordering!")
  }

  cat("    Cell names match verified!\n")

  # Create new assay
  cat("  Creating 'iso' assay...\n")
  integrated_obj[["iso"]] <- CreateAssay5Object(counts = counts_table_ordered)

  cat("  Isoform assay added successfully\n")
  cat("    Assays in object:", paste(names(integrated_obj@assays), collapse = ", "), "\n")

  return(integrated_obj)
}

## FUNCTION 4: Process isoform assay with full QC pipeline------
process_isoform_assay <- function(integrated_obj, timepoint_name) {

  cat("\n=== Processing isoform assay for", timepoint_name, "===\n")

  # Set default assay to iso
  DefaultAssay(integrated_obj) <- "iso"

  # Step 1: Normalize
  cat("\n  Step 1: Normalizing isoform data...\n")
  integrated_obj <- NormalizeData(integrated_obj, assay = "iso")

  # Step 2: Find variable features
  cat("  Step 2: Finding variable features...\n")
  integrated_obj <- FindVariableFeatures(integrated_obj,
                                         assay = "iso",
                                         selection.method = "vst",
                                         nfeatures = 2000)
  cat("    Found 2000 variable isoforms\n")

  # Get top 10 for plotting
  top10_iso <- head(VariableFeatures(integrated_obj), 10)

  # Create variable feature plot
  var_feature_plot <- VariableFeaturePlot(integrated_obj) +
    ggtitle(paste0(timepoint_name, ": Top 2000 Variable Isoforms")) +
    NoLegend()
  var_feature_plot <- LabelPoints(plot = var_feature_plot,
                                  points = top10_iso,
                                  labels = gsub(".*-", "", top10_iso),
                                  repel = TRUE)

  # Step 3: Scale data
  cat("  Step 3: Scaling isoform data...\n")
  integrated_obj <- ScaleData(integrated_obj, assay = "iso")

  # Step 4: Run PCA
  cat("  Step 4: Running PCA on isoform data...\n")
  integrated_obj <- RunPCA(integrated_obj,
                           assay = "iso",
                           features = VariableFeatures(integrated_obj),
                           reduction.name = "pca_iso",
                           verbose = FALSE)

  # Step 5: Quantitative elbow to determine optimal PCs
  cat("  Step 5: Determining optimal PCs using quantitative elbow...\n")

  # Get PCA standard deviations
  pca_stdev <- integrated_obj@reductions$pca_iso@stdev

  # Calculate variance explained
  var_explained <- pca_stdev^2 / sum(pca_stdev^2) * 100
  cumvar_explained <- cumsum(var_explained)

  # Method 1: Cumulative >90% and individual <5%
  opt_pcs_method1 <- which(cumvar_explained > 90 & var_explained < 5)[1]

  # Method 2: Consecutive change <0.1%
  var_change <- abs(diff(var_explained))
  opt_pcs_method2 <- which(var_change < 0.1)[1]

  # Take minimum and add 5
  optimal_pcs_calculated <- min(opt_pcs_method1, opt_pcs_method2, na.rm = TRUE)
  optimal_pcs_used <- optimal_pcs_calculated + 5

  cat("    Optimal PCs (calculated):", optimal_pcs_calculated, "\n")
  cat("    Optimal PCs (used with +5 adjustment):", optimal_pcs_used, "\n")

  # Create elbow plot
  elbow_data <- data.frame(
    PC = 1:length(var_explained),
    Variance = var_explained,
    Cumulative = cumvar_explained
  )

  elbow_plot <- ggplot(elbow_data, aes(x = PC, y = Variance)) +
    geom_line() +
    geom_point() +
    geom_vline(xintercept = optimal_pcs_calculated, linetype = "dashed", color = "blue") +
    geom_vline(xintercept = optimal_pcs_used, linetype = "dashed", color = "red") +
    annotate("text", x = optimal_pcs_calculated, y = max(var_explained) * 0.9,
             label = paste("Calculated:", optimal_pcs_calculated), color = "blue", hjust = -0.1) +
    annotate("text", x = optimal_pcs_used, y = max(var_explained) * 0.8,
             label = paste("Used:", optimal_pcs_used), color = "red", hjust = -0.1) +
    labs(title = paste0(timepoint_name, " - Isoform PCA Elbow Plot"),
         x = "Principal Component",
         y = "% Variance Explained") +
    theme_classic()

  # Step 6: Find neighbors
  cat("  Step 6: Finding neighbors on isoform PCA...\n")
  integrated_obj <- FindNeighbors(integrated_obj,
                                  reduction = "pca_iso",
                                  dims = 1:optimal_pcs_used,
                                  graph.name = c("iso_nn", "iso_snn"))

  # Step 7: Optimize clustering with silhouette analysis
  cat("  Step 7: Optimizing clustering resolution via silhouette analysis...\n")

  resolutions <- seq(0.1, 1.2, by = 0.1)
  sil_results <- data.frame(
    resolution = numeric(),
    avg_silhouette = numeric()
  )

  for(res in resolutions) {
    integrated_obj <- FindClusters(integrated_obj,
                                   resolution = res,
                                   graph.name = "iso_snn",
                                   verbose = FALSE)

    # Calculate silhouette using sourced functions
    clusters <- integrated_obj@meta.data[[paste0("iso_snn_res.", res)]]
    if(length(unique(clusters)) > 1) {
      sil <- silhouette_seurat(
        sobject = integrated_obj,
        cluster_col = paste0("iso_snn_res.", res),
        dims = 1:optimal_pcs_used,
        reduction = "pca_iso"
      )
      avg_sil <- silhouette_mean(sil)
    } else {
      avg_sil <- 0
    }

    sil_results <- rbind(sil_results, data.frame(resolution = res, avg_silhouette = avg_sil))
    cat("    Resolution", res, ": avg silhouette =", round(avg_sil, 4), "\n")
  }

  # Find optimal resolution
  optimal_res <- sil_results$resolution[which.max(sil_results$avg_silhouette)]
  cat("    Optimal resolution:", optimal_res, "\n")

  # Set optimal clustering
  integrated_obj <- FindClusters(integrated_obj,
                                 resolution = optimal_res,
                                 graph.name = "iso_snn",
                                 verbose = FALSE)

  # Store in dedicated metadata column
  integrated_obj$seurat_clusters_iso <- integrated_obj@meta.data[[paste0("iso_snn_res.", optimal_res)]]

  # Create silhouette plot
  sil_plot <- ggplot(sil_results, aes(x = resolution, y = avg_silhouette)) +
    geom_line() +
    geom_point() +
    geom_point(data = sil_results[sil_results$resolution == optimal_res, ],
               aes(x = resolution, y = avg_silhouette),
               color = "red", size = 3) +
    labs(title = paste0(timepoint_name, " - Isoform Silhouette Optimization"),
         x = "Resolution",
         y = "Average Silhouette Width") +
    theme_classic()

  # Step 8: Run UMAP
  cat("  Step 8: Running UMAP on isoform PCA...\n")
  integrated_obj <- RunUMAP(integrated_obj,
                            reduction = "pca_iso",
                            dims = 1:optimal_pcs_used,
                            reduction.name = "umap_iso",
                            verbose = FALSE)

  cat("  Isoform assay processing complete!\n")

  # Reset default assay
  DefaultAssay(integrated_obj) <- "RNA"

  return(list(
    seurat_obj = integrated_obj,
    plots = list(
      var_feature_plot = var_feature_plot,
      elbow_plot = elbow_plot,
      sil_plot = sil_plot
    ),
    optimal_pcs = optimal_pcs_used,
    optimal_res = optimal_res,
    top10_iso = top10_iso
  ))
}

## FUNCTION 5: Create comparison PDF------
create_comparison_pdf <- function(integrated_obj, timepoint_name, processing_results) {

  cat("\n=== Creating comparison PDF for", timepoint_name, "===\n")

  pdf_file <- paste0("./output_files/isoform_integration/",
                     timepoint_name, "_isoform_assay_report.pdf")

  pdf(pdf_file, width = 14, height = 10)

  # Page 1: Summary
  cat("  Generating Page 1: Summary...\n")

  summary_text <- paste0(
    "ISOFORM ASSAY INTEGRATION SUMMARY\n",
    "Timepoint: ", timepoint_name, "\n\n",
    "Dataset Information:\n",
    "  Total Cells: ", ncol(integrated_obj), "\n",
    "  Gene-level Features (RNA): ", nrow(integrated_obj[["RNA"]]), "\n",
    "  Isoform-level Features (iso): ", nrow(integrated_obj[["iso"]]), "\n",
    "  Consensus Cell Types: ", length(unique(integrated_obj$consensus_cell_type)), "\n\n",
    "Isoform Processing:\n",
    "  Optimal PCs: ", processing_results$optimal_pcs, "\n",
    "  Optimal Resolution: ", processing_results$optimal_res, "\n",
    "  Isoform Clusters: ", length(unique(integrated_obj$seurat_clusters_iso)), "\n",
    "  Gene-level Clusters: ", length(unique(integrated_obj$seurat_clusters)), "\n"
  )

  summary_plot <- ggplot() +
    annotate("text", x = 0, y = 1,
             label = summary_text,
             hjust = 0, vjust = 1,
             size = 5, family = "mono") +
    theme_void() +
    xlim(0, 1) + ylim(0, 1)

  grid.arrange(
    summary_plot,
    top = textGrob(paste(timepoint_name, "- Isoform Assay Integration Report"),
                   gp = gpar(fontsize = 16, fontface = "bold"))
  )

  # Page 2: Variable features + Elbow plot
  cat("  Generating Page 2: Variable features and elbow plot...\n")
  grid.arrange(
    processing_results$plots$var_feature_plot,
    processing_results$plots$elbow_plot,
    ncol = 2,
    top = textGrob("Isoform Data Processing - QC Plots",
                   gp = gpar(fontsize = 14, fontface = "bold"))
  )

  # Page 3: Silhouette plot
  cat("  Generating Page 3: Silhouette optimization...\n")
  print(processing_results$plots$sil_plot)

  # Page 4: Gene-level UMAP by consensus cell type
  cat("  Generating Page 4: Gene-level UMAP...\n")
  p_gene <- DimPlot(integrated_obj,
                    reduction = "umap.harmony",
                    group.by = "consensus_cell_type",
                    label = TRUE,
                    label.size = 3,
                    repel = TRUE) +
    ggtitle("Gene-level UMAP (Harmony Integration)") +
    theme(legend.text = element_text(size = 8))
  print(p_gene)

  # Page 5: Isoform-level UMAP by consensus cell type
  cat("  Generating Page 5: Isoform-level UMAP...\n")
  p_iso <- DimPlot(integrated_obj,
                   reduction = "umap_iso",
                   group.by = "consensus_cell_type",
                   label = TRUE,
                   label.size = 3,
                   repel = TRUE) +
    ggtitle("Isoform-level UMAP (No Integration)") +
    theme(legend.text = element_text(size = 8))
  print(p_iso)

  # Page 6: Side-by-side comparison - consensus cell types
  cat("  Generating Page 6: Side-by-side consensus comparison...\n")
  grid.arrange(
    p_gene + NoLegend(),
    p_iso + NoLegend(),
    ncol = 2,
    top = textGrob("Gene vs Isoform UMAP - Consensus Cell Types",
                   gp = gpar(fontsize = 14, fontface = "bold"))
  )

  # Page 7: Gene-level clusters
  cat("  Generating Page 7: Gene-level clusters...\n")
  p_gene_clust <- DimPlot(integrated_obj,
                          reduction = "umap.harmony",
                          group.by = "seurat_clusters",
                          label = TRUE,
                          label.size = 4) +
    ggtitle("Gene-level Clusters")
  print(p_gene_clust)

  # Page 8: Isoform-level clusters
  cat("  Generating Page 8: Isoform-level clusters...\n")
  p_iso_clust <- DimPlot(integrated_obj,
                         reduction = "umap_iso",
                         group.by = "seurat_clusters_iso",
                         label = TRUE,
                         label.size = 4) +
    ggtitle("Isoform-level Clusters")
  print(p_iso_clust)

  # Page 9: Side-by-side cluster comparison
  cat("  Generating Page 9: Side-by-side cluster comparison...\n")
  grid.arrange(
    p_gene_clust + NoLegend(),
    p_iso_clust + NoLegend(),
    ncol = 2,
    top = textGrob("Gene vs Isoform Clustering",
                   gp = gpar(fontsize = 14, fontface = "bold"))
  )

  dev.off()

  cat("  PDF saved:", pdf_file, "\n")

  return(pdf_file)
}

## FUNCTION 6: Process one timepoint------
process_timepoint <- function(timepoint_name, sample_ids, integrated_file) {

  cat("\n\n####################################################\n")
  cat("PROCESSING TIMEPOINT:", timepoint_name, "\n")
  cat("####################################################\n")

  # Step 1: Load integrated consensus object
  cat("\nStep 1: Loading integrated consensus object...\n")
  integrated_obj <- readRDS(integrated_file)
  cat("  Cells:", ncol(integrated_obj), "\n")
  cat("  Assays:", paste(names(integrated_obj@assays), collapse = ", "), "\n")

  # Step 2: Load and merge isoform objects
  iso_merged <- load_and_merge_isoform_objects(timepoint_name, sample_ids)

  # Step 3: Filter isoform object to match integrated cells
  filter_result <- filter_isoform_to_integrated_cells(iso_merged, integrated_obj)
  iso_filtered <- filter_result$iso_filtered
  cell_mapping <- filter_result$cell_mapping

  # Step 4: Add isoform assay to integrated object
  integrated_obj <- add_isoform_assay_to_object(integrated_obj, iso_filtered, cell_mapping)

  # Step 5: Process isoform assay
  processing_results <- process_isoform_assay(integrated_obj, timepoint_name)
  integrated_obj <- processing_results$seurat_obj

  # Step 6: Create comparison PDF
  pdf_file <- create_comparison_pdf(integrated_obj, timepoint_name, processing_results)

  # Step 7: Save updated object
  cat("\n=== Saving updated object ===\n")
  output_file <- paste0("./output_files/integrated_objects/",
                       timepoint_name, "_integrated_harmony_consensus_with_isoform_assay.rds")
  saveRDS(integrated_obj, file = output_file)
  cat("  Saved:", output_file, "\n")
  cat("  File size:", round(file.size(output_file) / 1e9, 2), "GB\n")

  cat("\n====================================\n")
  cat("PROCESSING COMPLETE:", timepoint_name, "\n")
  cat("====================================\n")

  return(list(
    n_cells = ncol(integrated_obj),
    n_gene_features = nrow(integrated_obj[["RNA"]]),
    n_iso_features = nrow(integrated_obj[["iso"]]),
    optimal_pcs = processing_results$optimal_pcs,
    optimal_res = processing_results$optimal_res,
    n_iso_clusters = length(unique(integrated_obj$seurat_clusters_iso)),
    overlap_stats = filter_result$stats,
    output_file = output_file,
    pdf_file = pdf_file
  ))
}

## MAIN EXECUTION: Process all timepoints------
cat("\n\n####################################################\n")
cat("STARTING ISOFORM ASSAY INTEGRATION\n")
cat("####################################################\n\n")

# Define timepoint configurations
timepoint_config <- list(
  "1M_Org" = list(
    samples = c("org_1A", "org_1B"),
    integrated_file = "./output_files/integrated_objects/1M_Org_integrated_harmony_consensus.rds"
  ),
  "3M_Org" = list(
    samples = c("org_3A", "org_3B", "org_3C"),
    integrated_file = "./output_files/integrated_objects/3M_Org_integrated_harmony_consensus.rds"
  ),
  "6M_Org" = list(
    samples = c("org_6A", "org_6B", "org_6C"),
    integrated_file = "./output_files/integrated_objects/6M_Org_integrated_harmony_consensus.rds"
  )
)

# Process each timepoint
results <- list()

for(timepoint_name in names(timepoint_config)) {
  config <- timepoint_config[[timepoint_name]]

  results[[timepoint_name]] <- process_timepoint(
    timepoint_name,
    config$samples,
    config$integrated_file
  )

  # Force garbage collection
  gc()
}

cat("\n\n####################################################\n")
cat("ALL TIMEPOINTS PROCESSED\n")
cat("####################################################\n")
cat("\nOutputs saved:\n")
for(tp in names(results)) {
  res <- results[[tp]]
  cat("\n", tp, ":\n", sep = "")
  cat("  Cells: ", res$n_cells, "\n", sep = "")
  cat("  Gene features: ", res$n_gene_features, "\n", sep = "")
  cat("  Isoform features: ", res$n_iso_features, "\n", sep = "")
  cat("  Isoform optimal PCs: ", res$optimal_pcs, "\n", sep = "")
  cat("  Isoform optimal resolution: ", res$optimal_res, "\n", sep = "")
  cat("  Isoform clusters: ", res$n_iso_clusters, "\n", sep = "")
  cat("  Cell overlap: ", res$overlap_stats$pct_integrated, "%\n", sep = "")
  cat("  RDS: ", basename(res$output_file), "\n", sep = "")
  cat("  PDF: ", basename(res$pdf_file), "\n", sep = "")
}

cat("\nScript completed successfully!\n")
