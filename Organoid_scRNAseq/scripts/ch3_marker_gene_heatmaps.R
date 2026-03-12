#!/usr/bin/env Rscript
# Marker Gene Expression Heatmaps for Consensus-Annotated Organoid Objects
# Author: Manveer Chauhan
# Description: Generates heatmaps showing top marker genes per consensus cell type
#              at both gene-level (RNA) and isoform-level (iso) across three
#              developmental timepoints (1M, 3M, 6M).
#
# Approach:
#   1. Find markers using same parameters as ch3_GO_analysis.R
#   2. Select top 10 markers per cell type (sorted by avg_log2FC)
#   3. Generate pseudobulked expression heatmaps using code from ch3_isoTWAS_DE_heatmaps.R

# Load required libraries ----
library(Seurat)
library(tidyverse)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(AnnotationHub)
library(ensembldb)

# Set working directory
setwd("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space")

# Helper Functions ----

#' Resolve transcript IDs to gene symbols using AnnotationHub
#'
#' Adapted from Aim1.1_getGeneTxMetadata.R
#' Uses AnnotationHub EnsDb to resolve -NA suffixes to proper gene symbols
#'
#' @param transcript_ids Character vector of transcript IDs (format: ENST######.#-GENE or ENST######.#-NA)
#' @return data.frame with tx_id_original, tx_id_clean, gene_symbol_resolved, tx_is_canonical
resolve_transcript_symbols <- function(transcript_ids) {

  cat("  Initializing AnnotationHub...\n")

  # Initialize AnnotationHub
  ah <- AnnotationHub()
  human_ens <- query(ah, c("Homo sapiens", "EnsDb"))

  # Use latest EnsDb version
  edb <- human_ens[[length(human_ens)]]
  cat("  Using EnsDb:", human_ens$title[length(human_ens)], "\n")

  # Extract transcript metadata
  cat("  Extracting transcript metadata from EnsDb...\n")
  tx_metadata <- transcripts(edb, return.type = "data.frame") %>%
    dplyr::select(tx_id, tx_external_name, tx_is_canonical, tx_biotype) %>%
    # Strip version numbers from tx_id for matching
    dplyr::mutate(tx_id_clean = str_remove(tx_id, "\\..*"))

  cat("  EnsDb contains", nrow(tx_metadata), "transcripts\n")

  # Parse input transcript IDs
  input_df <- data.frame(
    tx_id_original = transcript_ids,
    stringsAsFactors = FALSE
  ) %>%
    # Extract ENST ID (strip version and gene suffix)
    dplyr::mutate(
      tx_id_clean = str_remove(str_extract(tx_id_original, "^ENST[0-9]+"), "\\..*"),
      current_suffix = str_extract(tx_id_original, "-.*$")
    )

  # Join with AnnotationHub metadata
  cat("  Matching transcript IDs with EnsDb...\n")
  resolved <- input_df %>%
    left_join(tx_metadata, by = "tx_id_clean") %>%
    dplyr::mutate(
      # Use tx_external_name from AnnotationHub, fallback to current suffix if not found
      gene_symbol_resolved = coalesce(tx_external_name, str_remove(current_suffix, "^-"))
    ) %>%
    dplyr::select(tx_id_original, tx_id_clean, gene_symbol_resolved, tx_is_canonical, tx_biotype)

  # Report statistics
  n_input <- nrow(resolved)
  n_resolved <- sum(!is.na(resolved$gene_symbol_resolved))
  n_na_input <- sum(grepl("-NA$", resolved$tx_id_original))
  n_na_fixed <- sum(grepl("-NA$", resolved$tx_id_original) &
                      !is.na(resolved$gene_symbol_resolved) &
                      resolved$gene_symbol_resolved != "NA")

  cat("  Resolution statistics:\n")
  cat("    Input transcripts:", n_input, "\n")
  cat("    Successfully resolved:", n_resolved, "(", round(100*n_resolved/n_input, 1), "%)\n")
  cat("    Input -NA transcripts:", n_na_input, "\n")
  cat("    Fixed -NA transcripts:", n_na_fixed, "(", round(100*n_na_fixed/n_na_input, 1), "%)\n")

  return(resolved)
}

#' Fix feature names in marker dataframe by resolving -NA suffixes
#'
#' IMPORTANT: Adds 'gene_display' column with resolved names while preserving
#' original 'gene' column for feature lookup in Seurat object
#'
#' @param marker_df Marker dataframe from FindAllMarkers (must have 'gene' column)
#' @return Marker dataframe with additional 'gene_display' column
fix_feature_names <- function(marker_df) {

  if (is.null(marker_df) || nrow(marker_df) == 0) {
    cat("  No markers to fix\n")
    return(marker_df)
  }

  cat("Fixing -NA suffixes in feature names...\n")

  # Extract unique feature names
  unique_features <- unique(marker_df$gene)
  cat("  Unique features to resolve:", length(unique_features), "\n")

  # Count -NA features
  n_na_features <- sum(grepl("-NA$", unique_features))
  cat("  Features with -NA suffix:", n_na_features, "\n\n")

  # Resolve to proper gene symbols
  resolved <- resolve_transcript_symbols(unique_features)

  # Create display names: tx_id_clean.version-gene_symbol_resolved
  # Need to preserve the version number from original
  resolved <- resolved %>%
    dplyr::mutate(
      # Extract version from original (e.g., ".1" from "ENST00000123456.1-NA")
      version = str_extract(tx_id_original, "\\.[0-9]+"),
      # Construct display name: ENST######.#-GENE
      gene_display = paste0(tx_id_clean, version, "-", gene_symbol_resolved)
    )

  # Map back to marker dataframe
  # Keep ORIGINAL 'gene' for Seurat lookup, add 'gene_display' for heatmap labels
  marker_df_fixed <- marker_df %>%
    left_join(resolved %>% dplyr::select(tx_id_original, gene_display),
              by = c("gene" = "tx_id_original")) %>%
    dplyr::mutate(gene_display = coalesce(gene_display, gene))

  cat("\n  ✓ Feature name resolution complete\n")
  cat("  Note: 'gene' column preserved for Seurat lookup, 'gene_display' added for labels\n\n")

  return(marker_df_fixed)
}

#' Find marker genes for a single timepoint
#'
#' Uses same parameters as ch3_GO_analysis.R:
#'   - FindAllMarkers with only.pos = TRUE
#'   - Filter: p_val_adj <= 0.05 & avg_log2FC >= 0.5
#'
#' @param seu.obj Seurat object with consensus_cell_type metadata
#' @param assay_type Character, "RNA" or "iso"
#' @param timepoint_name Character, e.g., "1M_Org"
#' @return Filtered marker dataframe
find_markers_per_timepoint <- function(seu.obj, assay_type, timepoint_name) {

  cat("\n========================================\n")
  cat("Finding markers for:", timepoint_name, "-", assay_type, "\n")
  cat("========================================\n\n")

  # Set assay and identity
  DefaultAssay(seu.obj) <- assay_type
  Idents(seu.obj) <- "consensus_cell_type"

  cat("  Cells:", ncol(seu.obj), "\n")
  cat("  Features:", nrow(seu.obj), "\n")
  cat("  Cell types:", length(unique(seu.obj$consensus_cell_type)), "\n")
  cat("  Cell type distribution:\n")
  print(table(seu.obj$consensus_cell_type))
  cat("\n")

  # Join layers ONLY for RNA assay (iso already joined)
  if (assay_type == "RNA") {
    cat("Joining RNA data layers for FindAllMarkers...\n")
    seu.obj <- JoinLayers(seu.obj, assay = "RNA")
    cat("  Layers joined successfully\n\n")
  } else {
    cat("iso assay already has unified layers - no JoinLayers needed\n\n")
  }

  # Find marker genes
  cat("Running FindAllMarkers...\n")
  cat("  Parameters: only.pos = TRUE, Seurat defaults (min.pct=0.01, logfc.threshold=0.1)\n")
  cat("  Post-filtering: p_val_adj <= 0.05, avg_log2FC >= 0.5\n\n")

  markers <- FindAllMarkers(
    seu.obj,
    assay = assay_type,
    only.pos = TRUE,
    verbose = TRUE
  )

  # Check if markers were found
  if (is.null(markers) || nrow(markers) == 0) {
    cat("\n  WARNING: No markers found for this timepoint!\n")
    return(NULL)
  }

  cat("\n  Total markers before filtering:", nrow(markers), "\n")

  # Filter for significance and logFC threshold
  markers_sig <- markers %>%
    dplyr::filter(p_val_adj <= 0.05 & avg_log2FC >= 0.5)

  cat("  Significant markers after filtering:", nrow(markers_sig), "\n")
  cat("  Markers per cell type:\n")
  print(table(markers_sig$cluster))
  cat("\n")

  return(markers_sig)
}

#' Select top N markers per cell type
#'
#' Sorts by avg_log2FC descending and takes top N per cluster.
#' If fewer than N markers available, uses all available markers.
#'
#' @param markers_df Filtered marker dataframe
#' @param n_per_celltype Integer, number of markers to select per cell type (default 10)
#' @return Dataframe with top markers
select_top_markers <- function(markers_df, n_per_celltype = 10) {

  if (is.null(markers_df) || nrow(markers_df) == 0) {
    cat("  No markers to select from\n")
    return(NULL)
  }

  cat("Selecting top", n_per_celltype, "markers per cell type...\n")

  top_markers <- markers_df %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(desc(avg_log2FC)) %>%
    dplyr::slice_head(n = n_per_celltype) %>%
    dplyr::ungroup()

  cat("  Total markers selected:", nrow(top_markers), "\n")
  cat("  Markers per cell type:\n")

  marker_counts <- top_markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(n_markers = n(), .groups = "drop")

  print(marker_counts)
  cat("\n")

  # Warn about cell types with fewer than N markers
  low_marker_celltypes <- marker_counts %>%
    dplyr::filter(n_markers < n_per_celltype)

  if (nrow(low_marker_celltypes) > 0) {
    cat("  NOTE: The following cell types have fewer than", n_per_celltype, "markers:\n")
    print(low_marker_celltypes)
    cat("\n")
  }

  return(top_markers)
}

#' Generate pseudobulked expression heatmap
#'
#' Adapted from ch3_isoTWAS_DE_heatmaps.R:generate_DE_heatmap()
#'
#' @param seu.obj Seurat object
#' @param marker_df Dataframe with selected markers (must have 'gene' column, optionally 'gene_display')
#' @param timepoint_name Character, e.g., "1M_Org"
#' @param assay_type Character, "RNA" or "iso"
#' @return pheatmap object or NULL if failed
generate_marker_heatmap <- function(seu.obj, marker_df, timepoint_name, assay_type) {

  if (is.null(marker_df) || nrow(marker_df) == 0) {
    cat("  No markers to plot\n")
    return(NULL)
  }

  # Use ORIGINAL gene names for Seurat feature lookup
  featureList <- marker_df$gene

  cat("Generating", assay_type, "heatmap for", timepoint_name, "with",
      length(featureList), "features\n")

  # Extract pseudobulked expression values for each consensus cell type
  # Uses ORIGINAL feature names (e.g., ENST00000123456.1-NA)
  pseudobulk.Expression.list <- AggregateExpression(
    seu.obj,
    features = featureList,
    assays = assay_type,
    group.by = "consensus_cell_type",
    return.seurat = FALSE
  )

  # Extract the matrix for the specific assay
  pseudobulk.Expression.mat <- pseudobulk.Expression.list[[assay_type]] %>%
    as.data.frame()

  cat("  Pseudobulk matrix dimensions:", nrow(pseudobulk.Expression.mat), "features x",
      ncol(pseudobulk.Expression.mat), "cell types\n")

  # Check if features were dropped during aggregation
  n_dropped <- length(featureList) - nrow(pseudobulk.Expression.mat)
  if (n_dropped > 0) {
    cat("  ⚠ Warning:", n_dropped, "features dropped during aggregation (likely all zeros)\n")
  }

  # Map display names to rownames if available
  # This replaces -NA suffixes with resolved gene symbols in heatmap
  if ("gene_display" %in% colnames(marker_df)) {
    cat("  Applying resolved gene names to heatmap labels...\n")

    # Create mapping from original to display names
    name_map <- marker_df %>%
      dplyr::select(gene, gene_display) %>%
      dplyr::distinct()

    # Get current rownames (original feature names from Seurat)
    current_rownames <- rownames(pseudobulk.Expression.mat)

    # Map to display names
    display_names <- name_map$gene_display[match(current_rownames, name_map$gene)]

    # Update rownames (fallback to original if mapping fails)
    rownames(pseudobulk.Expression.mat) <- ifelse(is.na(display_names),
                                                   current_rownames,
                                                   display_names)

    cat("  ✓ Rownames updated with resolved gene symbols\n")
  }

  # VALIDATION: Check minimum feature count
  if (nrow(pseudobulk.Expression.mat) < 3) {
    warning("Too few features (", nrow(pseudobulk.Expression.mat),
            ") for meaningful heatmap. Minimum 3 required. Skipping plot.")
    return(NULL)
  }

  # VALIDATION: Check for variance
  row_vars <- apply(pseudobulk.Expression.mat, 1, var, na.rm = TRUE)
  n_variable <- sum(row_vars > 0, na.rm = TRUE)
  cat("  Features with variance:", n_variable, "/", nrow(pseudobulk.Expression.mat), "\n")

  if (n_variable < 2) {
    warning("Insufficient variance in features for clustering (",
            n_variable, " variable features). Skipping plot.")
    return(NULL)
  }

  cat("  Cell types (x-axis labels):", paste(colnames(pseudobulk.Expression.mat), collapse = ", "), "\n")

  # Create heatmap title
  assay_label <- ifelse(assay_type == "RNA", "Gene-Level", "Isoform-Level")
  plot_title <- paste0(timepoint_name, " - Top Marker ", assay_label, " Expression\n(Row-scaled Z-score)")

  # Wrap pheatmap in tryCatch for robust error handling
  marker_heatmap <- tryCatch({
    pheatmap::pheatmap(
      pseudobulk.Expression.mat,
      cluster_rows = TRUE,
      show_rownames = TRUE,
      color = viridis::inferno(100),
      breaks = seq(-3, 3, length.out = 101),
      border_color = NA,
      fontsize = 10,
      scale = "row",
      fontsize_row = 6,
      fontsize_col = 9,
      height = 20,
      angle_col = 45,
      main = plot_title
    )
  }, error = function(e) {
    warning("pheatmap failed for ", timepoint_name, ": ", e$message)
    cat("  ✗ Heatmap generation failed - error:", e$message, "\n")
    return(NULL)
  }, warning = function(w) {
    # Suppress warnings but continue
    suppressWarnings(
      pheatmap::pheatmap(
        pseudobulk.Expression.mat,
        cluster_rows = TRUE,
        show_rownames = TRUE,
        color = viridis::inferno(100),
        breaks = seq(-3, 3, length.out = 101),
        border_color = NA,
        fontsize = 10,
        scale = "row",
        fontsize_row = 6,
        fontsize_col = 9,
        height = 20,
        angle_col = 45,
        main = plot_title
      )
    )
  })

  if (is.null(marker_heatmap)) {
    cat("  ✗ Heatmap generation returned NULL\n")
  } else {
    cat("  ✓ Heatmap generated successfully\n")
  }

  return(marker_heatmap)
}

# Main Analysis ----

cat("\n")
cat("================================================================================\n")
cat("  Marker Gene Expression Heatmaps for Consensus-Annotated Organoid Objects\n")
cat("  Author: Manveer Chauhan\n")
cat("  Date:", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("================================================================================\n")

# Create output directory
output_dir <- "./output_files/marker_gene_heatmaps"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("\nCreated output directory:", output_dir, "\n")
}

# Load consensus objects
cat("\nLoading consensus-annotated objects with isoform assays...\n")
one_month_org <- readRDS("./output_files/integrated_objects/1M_Org_integrated_harmony_consensus_with_isoform_assay.rds")
three_month_org <- readRDS("./output_files/integrated_objects/3M_Org_integrated_harmony_consensus_with_isoform_assay.rds")
six_month_org <- readRDS("./output_files/integrated_objects/6M_Org_integrated_harmony_consensus_with_isoform_assay.rds")
cat("  ✓ All objects loaded\n")

# Define timepoints and objects
timepoints <- c("1M_Org", "3M_Org", "6M_Org")
seu_objects <- list(
  "1M_Org" = one_month_org,
  "3M_Org" = three_month_org,
  "6M_Org" = six_month_org
)

# Define assays to analyze
assays <- c("RNA", "iso")
assay_labels <- c("RNA" = "Gene-Level", "iso" = "Isoform-Level")

# Store all results for summary
all_results <- list()

# Process each assay type separately
for (assay_type in assays) {

  cat("\n")
  cat("################################################################################\n")
  cat("##  ANALYZING ASSAY:", assay_type, "(", assay_labels[assay_type], ")  ##\n")
  cat("################################################################################\n")

  # Determine output PDF filename
  assay_suffix <- ifelse(assay_type == "RNA", "RNA", "ISO")
  pdf_filename <- file.path(output_dir, paste0("marker_genes_", assay_suffix, "_heatmaps.pdf"))

  # Open PDF device (each timepoint will be on a separate page)
  pdf(pdf_filename, width = 14, height = 10)

  page_count <- 0

  # Process each timepoint
  for (tp in timepoints) {

    seu.obj <- seu_objects[[tp]]

    # Step 1: Find markers
    markers <- find_markers_per_timepoint(seu.obj, assay_type, tp)

    if (is.null(markers) || nrow(markers) == 0) {
      cat("  Skipping", tp, "- no markers found\n\n")
      all_results[[paste0(assay_type, "_", tp)]] <- list(
        markers = NULL,
        top_markers = NULL,
        heatmap_success = FALSE
      )
      next
    }

    # Step 2: Select top 10 markers per cell type
    top_markers <- select_top_markers(markers, n_per_celltype = 10)

    if (is.null(top_markers) || nrow(top_markers) == 0) {
      cat("  Skipping", tp, "- no top markers selected\n\n")
      all_results[[paste0(assay_type, "_", tp)]] <- list(
        markers = markers,
        top_markers = NULL,
        heatmap_success = FALSE
      )
      next
    }

    # Step 2.5: Fix -NA suffixes for iso assay using AnnotationHub
    if (assay_type == "iso") {
      cat("\n")
      top_markers <- fix_feature_names(top_markers)
    }

    # Step 3: Save top markers to CSV
    csv_filename <- file.path(output_dir, paste0(tp, "_", assay_suffix, "_top_markers.csv"))
    write.csv(top_markers, csv_filename, row.names = FALSE)
    cat("  ✓ Saved top markers to:", csv_filename, "\n\n")

    # Step 4: Generate heatmap
    cat("Plotting heatmap (page", page_count + 1, ")...\n")
    heatmap_result <- generate_marker_heatmap(seu.obj, top_markers, tp, assay_type)

    if (!is.null(heatmap_result)) {
      page_count <- page_count + 1
      all_results[[paste0(assay_type, "_", tp)]] <- list(
        markers = markers,
        top_markers = top_markers,
        heatmap_success = TRUE
      )
      cat("  ✓ Heatmap added to PDF (page", page_count, ")\n\n")
    } else {
      cat("  ✗ Heatmap generation failed for", tp, "\n\n")
      all_results[[paste0(assay_type, "_", tp)]] <- list(
        markers = markers,
        top_markers = top_markers,
        heatmap_success = FALSE
      )
    }
  }

  # Close PDF device
  dev.off()

  if (page_count > 0) {
    cat("\n✓ Saved", page_count, "heatmap page(s) to:", pdf_filename, "\n")
  } else {
    # Remove empty PDF if no heatmaps were created
    if (file.exists(pdf_filename)) {
      file.remove(pdf_filename)
    }
    cat("\n✗ No valid heatmaps generated for", assay_type, "- PDF removed\n")
  }
}

# Generate summary statistics ----
cat("\n")
cat("================================================================================\n")
cat("  Summary Statistics\n")
cat("================================================================================\n\n")

summary_data <- data.frame(
  Assay = character(),
  Timepoint = character(),
  Total_Markers = integer(),
  Significant_Markers = integer(),
  Top_Markers_Selected = integer(),
  Heatmap_Success = logical(),
  stringsAsFactors = FALSE
)

for (assay_type in assays) {
  cat("Assay:", assay_type, "\n")

  for (tp in timepoints) {
    result_key <- paste0(assay_type, "_", tp)
    result <- all_results[[result_key]]

    if (!is.null(result)) {
      total_markers <- ifelse(is.null(result$markers), 0, nrow(result$markers))
      top_markers <- ifelse(is.null(result$top_markers), 0, nrow(result$top_markers))

      cat("  ", tp, "- Total significant markers:", total_markers,
          "| Top markers selected:", top_markers,
          "| Heatmap:", ifelse(result$heatmap_success, "✓", "✗"), "\n")

      summary_data <- rbind(summary_data, data.frame(
        Assay = assay_type,
        Timepoint = tp,
        Total_Markers = total_markers,
        Significant_Markers = total_markers,
        Top_Markers_Selected = top_markers,
        Heatmap_Success = result$heatmap_success,
        stringsAsFactors = FALSE
      ))
    } else {
      cat("  ", tp, "- No data\n")

      summary_data <- rbind(summary_data, data.frame(
        Assay = assay_type,
        Timepoint = tp,
        Total_Markers = 0,
        Significant_Markers = 0,
        Top_Markers_Selected = 0,
        Heatmap_Success = FALSE,
        stringsAsFactors = FALSE
      ))
    }
  }
  cat("\n")
}

# Save summary table
summary_path <- file.path(output_dir, "marker_heatmap_summary.csv")
write.csv(summary_data, summary_path, row.names = FALSE)

cat("Saved summary statistics to:", summary_path, "\n")

cat("\n")
cat("================================================================================\n")
cat("  Analysis Complete!\n")
cat("  Output directory:", output_dir, "\n")
cat("================================================================================\n\n")

cat("Generated files:\n")
cat("  - marker_genes_RNA_heatmaps.pdf: Gene-level marker heatmaps (3 pages)\n")
cat("  - marker_genes_ISO_heatmaps.pdf: Isoform-level marker heatmaps (3 pages)\n")
cat("  - [timepoint]_[RNA/ISO]_top_markers.csv: Top marker lists\n")
cat("  - marker_heatmap_summary.csv: Summary statistics table\n\n")
