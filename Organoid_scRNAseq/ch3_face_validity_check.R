# ============================================================================
# Single-Cell Isoform Expression Analysis
# Neurodevelopmental Models Analysis
# ============================================================================
# Author: [Manveer Chauhan]
# Date: [20.06.25]
# Description: Analysis of isoform-level expression patterns in brain organoids
#              across developmental timepoints (1, 3, and 6 months)
# ============================================================================

# REQUIRED LIBRARIES --------------------------------------------------------
library(Seurat)
library(tidyverse)
library(gridExtra)
library(grid)
library(patchwork)  # For plot_annotation function

# GLOBAL PARAMETERS ----------------------------------------------------------
# Directory paths
BASE_DIR <- "/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space"
setwd(BASE_DIR)

# Input file paths
ONE_MONTH_FILE <- file.path(BASE_DIR, "output_files/integrated_objects/1M_Org_integrated_harmony_consensus_with_isoform_assay.rds")
THREE_MONTH_FILE <- file.path(BASE_DIR, "output_files/integrated_objects/3M_Org_integrated_harmony_consensus_with_isoform_assay.rds")
SIX_MONTH_FILE <- file.path(BASE_DIR, "output_files/integrated_objects/6M_Org_integrated_harmony_consensus_with_isoform_assay.rds")

# Output directory
OUTPUT_DIR <- "./output_files/face_validity"

# Analysis parameters
IDENTITY_COLUMN <- "consensus_cell_type"
DEFAULT_ASSAY <- "RNA"  # Using gene-level (RNA) not isoform-level (iso)
DEFAULT_REDUCTION <- "umap.harmony"
DEFAULT_SLOT <- "counts"
DE_TEST_METHOD <- "wilcox"

# Features of interest for splicing analysis
SPLICING_FACTORS <- c("PTBP1", "PTBP2", "RBFOX1", "RBFOX2", "RBFOX3", "NOVA1", "NOVA2")
SYNAPTIC_GENES <- c("DLG4", "NRXN1", "NRXN2", "NRXN3", 
                    "NLGN1", "NLGN2", "NLGN3", "NLGN4X",
                    "NLGN4Y")
ALL_FEATURES <- c(SPLICING_FACTORS, SYNAPTIC_GENES)

# Plotting parameters
DEFAULT_NCOL <- 2
DEFAULT_POINT_SIZE <- 0.3
CHUNK_SIZE <- 4
PT_SIZE <- 1.2  # Point size for co-expression UMAPs
COEXPRESSION_COLOURS <- c("grey20", "cyan", "magenta")
#COEXPRESSION_COLOURS = c("grey90", "#0072B2", "orange")  # Original colors
#COEXPRESSION_COLOURS = c("grey90", "#7570b3", "#1b9e77")  # Alternative colors

# ============================================================================
# FUNCTION DEFINITIONS
# ============================================================================

#' Check available metadata in Seurat objects
#'
#' @param seu_obj Seurat object
#' @param obj_name Name of the object for reporting
check_metadata <- function(seu_obj, obj_name) {
  cat("Metadata columns in", obj_name, ":\n")
  cat(paste(colnames(seu_obj@meta.data), collapse = ", "), "\n\n")
  
  # Check for common pseudotime-related columns
  pseudotime_cols <- grep("pseudotime|monocle|trajectory|time", 
                          colnames(seu_obj@meta.data), 
                          value = TRUE, ignore.case = TRUE)
  
  if (length(pseudotime_cols) > 0) {
    cat("Found potential pseudotime columns:", paste(pseudotime_cols, collapse = ", "), "\n")
  } else {
    cat("No pseudotime-related columns found\n")
  }
  cat("---\n")
}

#' Gene symbol aliases for common gene names
#'
#' @description Maps commonly used gene names to their official symbols
gene_aliases <- list(
  "VGLUT1" = "SLC17A7",
  "PSD95" = "DLG4",
  "CTIP2" = "BCL11B",
  "FOG2" = "ZFPM2"
)

#' Resolve gene symbol to full feature name (ENSG-SYMBOL format)
#'
#' @param gene_symbol Short gene symbol (e.g., "PTBP1")
#' @param seurat_obj Seurat object to search in
#' @param verbose Print diagnostic messages
#' @return Full feature name (e.g., "ENSG00000011304.15-PTBP1") or NA if not found
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

#' Load and prepare Seurat objects
#'
#' @description Loads Seurat objects and sets the default identity
#' @return List containing the three Seurat objects
load_seurat_objects <- function() {
  cat("Loading Seurat objects...\n")
  
  # Load objects
  one_month_org <- readRDS(ONE_MONTH_FILE)
  three_month_org <- readRDS(THREE_MONTH_FILE)
  six_month_org <- readRDS(SIX_MONTH_FILE)
  
  # Set identities
  Idents(one_month_org) <- IDENTITY_COLUMN
  Idents(three_month_org) <- IDENTITY_COLUMN
  Idents(six_month_org) <- IDENTITY_COLUMN
  
  cat("Successfully loaded all objects\n")
  
  return(list(
    one_month = one_month_org,
    three_month = three_month_org,
    six_month = six_month_org
  ))
}

#' Find differentially expressed transcripts/isoforms
#'
#' @param seu.obj Seurat object
#' @param assay.type Assay to use for analysis (default: "iso")
#' @return DataFrame with differential expression results
returnTWAS_DE_Features <- function(seu.obj, assay.type = DEFAULT_ASSAY) {
  cat("Running FindAllMarkers for isoform-level analysis...\n")
  
  cellTypeMarkers_sc <- FindAllMarkers(
    seu.obj,
    assay = assay.type,
    slot = DEFAULT_SLOT,
    test.use = DE_TEST_METHOD
  )
  
  cat("Found", nrow(cellTypeMarkers_sc), "total markers\n")
  
  return(cellTypeMarkers_sc)
}

#' Create co-expression plot between two features
#'
#' @param seu.obj Seurat object
#' @param feature1 First gene/feature name
#' @param feature2 Second gene/feature name
#' @param reduction Dimensional reduction to use
#' @param title_prefix Prefix for plot title (e.g., "1 Month")
#' @return ggplot object
makeCoExpressionPlot <- function(seu.obj,
                                 feature1 = "PTBP1",
                                 feature2 = "PTBP2",
                                 reduction = DEFAULT_REDUCTION,
                                 title_prefix = "") {

  # Resolve gene symbols to full feature names (ENSG-SYMBOL format)
  feature1_full <- resolve_gene_name_with_aliases(feature1, seu.obj, verbose = FALSE)
  feature2_full <- resolve_gene_name_with_aliases(feature2, seu.obj, verbose = FALSE)

  # Check if features were found
  if(is.na(feature1_full)) {
    warning("Feature not found: ", feature1)
    return(NULL)
  }
  if(is.na(feature2_full)) {
    warning("Feature not found: ", feature2)
    return(NULL)
  }

  plt.output <- FeaturePlot(
    seu.obj,
    features = c(feature1_full, feature2_full),
    blend = TRUE,
    reduction = reduction,
    pt.size = PT_SIZE,
    cols = COEXPRESSION_COLOURS,
    order = TRUE
  )

  # Add title if prefix provided (use original gene symbols, not full names)
  if (title_prefix != "") {
    plot_title <- paste(title_prefix, ":", feature1, "vs", feature2, "Co-expression")
    plt.output <- plt.output +
      plot_annotation(title = plot_title, theme = theme(plot.title = element_text(hjust = 0.5, size = 10)))
  }

  # Add dark theme for better color contrast (& applies to all patchwork panels)
  plt.output <- plt.output &
    theme(
      panel.background = element_rect(fill = "black"),
      plot.background = element_rect(fill = "black"),
      panel.grid.major = element_line(color = "grey30", size = 0.2),
      panel.grid.minor = element_line(color = "grey20", size = 0.1),
      axis.text = element_text(color = "white"),
      axis.title = element_text(color = "white"),
      plot.title = element_text(color = "white", size = 8),
      legend.background = element_rect(fill = "black"),
      legend.text = element_text(color = "white"),
      legend.title = element_text(color = "white")
    )

  return(plt.output)
}

#' Plot features in chunks using violin plots
#'
#' @param seurat_obj Seurat object
#' @param features Vector of feature names
#' @param chunk_size Number of features per chunk
#' @param ncol Number of columns in plot
#' @param title_prefix Prefix for plot titles
plot_feature_chunks <- function(seurat_obj,
                                features,
                                chunk_size = CHUNK_SIZE,
                                ncol = DEFAULT_NCOL,
                                title_prefix = "") {

  # Resolve gene symbols to full feature names (ENSG-SYMBOL format)
  features_full <- sapply(features, function(f) {
    resolve_gene_name_with_aliases(f, seurat_obj, verbose = FALSE)
  })

  # Remove any features that weren't found
  features_full <- features_full[!is.na(features_full)]

  if(length(features_full) == 0) {
    warning("No features found for violin plots")
    return(NULL)
  }

  # Split features into chunks
  feature_chunks <- split(features_full, ceiling(seq_along(features_full) / chunk_size))

  # Loop over chunks and plot
  for (i in seq_along(feature_chunks)) {
    print(
      VlnPlot(seurat_obj, features = feature_chunks[[i]], ncol = ncol)
    )
  }
}

#' Create pseudotime comparison plots
#'
#' @param seu_obj Seurat object
#' @param title_prefix Prefix for plot titles
#' @param reduction Dimensional reduction to use
#' @return Combined ggplot object or NULL if pseudotime not available
plot_pseudotime_comparison <- function(seu_obj, 
                                       title_prefix = "", 
                                       reduction = DEFAULT_REDUCTION) {
  
  p1 <- DimPlot(
    seu_obj, 
    reduction = reduction, 
    group.by = IDENTITY_COLUMN,
    pt.size = DEFAULT_POINT_SIZE
  ) +
    ggtitle(paste(title_prefix, "- Cell types"))
  
  # Check if pseudotime exists in the metadata
  if ("pseudotime" %in% colnames(seu_obj@meta.data)) {
    p2 <- FeaturePlot(
      seu_obj, 
      features = "pseudotime", 
      reduction = reduction,
      pt.size = DEFAULT_POINT_SIZE
    ) +
      scale_color_viridis_c(name = "Pseudotime") +
      ggtitle(paste(title_prefix, "- Pseudotime"))
    
    return(p1 | p2)
  } else {
    cat("Warning: pseudotime not found in", title_prefix, "metadata. Skipping pseudotime plot.\n")
    return(p1)
  }
}

#' Create plots with and without cell type annotation labels shown
#'
#' @param seu_obj Seurat object
#' @param title_prefix Prefix for plot titles
#' @param reduction Dimensional reduction to use
#' @return Combined ggplot object or NULL if pseudotime not available
plot_annotation_comparison <- function(seu_obj,
                                       title_prefix = "",
                                       reduction = DEFAULT_REDUCTION) {

  p1 <- DimPlot(
    seu_obj,
    reduction = reduction,
    group.by = IDENTITY_COLUMN,
    pt.size = PT_SIZE
  ) +
    ggtitle(paste(title_prefix, "- Cell types")) +
    theme(
      panel.background = element_rect(fill = "black"),
      plot.background = element_rect(fill = "black"),
      panel.grid.major = element_line(color = "grey30", size = 0.2),
      panel.grid.minor = element_line(color = "grey20", size = 0.1),
      axis.text = element_text(color = "white"),
      axis.title = element_text(color = "white"),
      plot.title = element_text(color = "white", size = 10),
      legend.background = element_rect(fill = "black"),
      legend.text = element_text(color = "white"),
      legend.title = element_text(color = "white")
    )

  return(p1)
}

#' Generate ridge plots for multiple features
#'
#' @param seu_obj Seurat object
#' @param features Vector of feature names
#' @param ncol Number of columns in plot
#' @return ggplot object
plot_ridge_features <- function(seu_obj,
                                features = ALL_FEATURES,
                                ncol = DEFAULT_NCOL) {

  # Resolve gene symbols to full feature names (ENSG-SYMBOL format)
  features_full <- sapply(features, function(f) {
    resolve_gene_name_with_aliases(f, seu_obj, verbose = FALSE)
  })

  # Remove any features that weren't found
  features_full <- features_full[!is.na(features_full)]

  if(length(features_full) == 0) {
    warning("No features found for ridge plot")
    return(NULL)
  }

  theme_set(theme_minimal())
  return(RidgePlot(seu_obj, features = features_full, ncol = ncol))
}

#' Plot ridge features in chunks to avoid squashing
#'
#' @param seurat_obj Seurat object
#' @param features Vector of feature names
#' @param chunk_size Number of features per chunk (default 4-6)
#' @param ncol Number of columns in plot
#' @return NULL (prints plots directly)
plot_ridge_feature_chunks <- function(seurat_obj,
                                     features = ALL_FEATURES,
                                     chunk_size = 5,
                                     ncol = DEFAULT_NCOL) {

  # Resolve gene symbols to full feature names (ENSG-SYMBOL format)
  features_full <- sapply(features, function(f) {
    resolve_gene_name_with_aliases(f, seurat_obj, verbose = FALSE)
  })

  # Remove any features that weren't found
  features_full <- features_full[!is.na(features_full)]

  if(length(features_full) == 0) {
    warning("No features found for ridge plots")
    return(NULL)
  }

  # Split features into chunks
  feature_chunks <- split(features_full, ceiling(seq_along(features_full) / chunk_size))

  # Loop over chunks and plot
  theme_set(theme_minimal())
  for (i in seq_along(feature_chunks)) {
    print(
      RidgePlot(seurat_obj, features = feature_chunks[[i]], ncol = ncol)
    )
  }
}

#' Generate co-expression report for a single timepoint
#'
#' @param seu_obj Seurat object
#' @param timepoint_name Name of timepoint (e.g., "1M_Org", "3M_Org", "6M_Org")
#' @param timepoint_label Display label (e.g., "1 Month", "3 Month", "6 Month")
#' @description Generates co-expression plots with UMAP first page
#'              Designed to be called within a PDF device context
generate_coexpression_report <- function(seu_obj, timepoint_name, timepoint_label) {
  cat("Generating co-expression report for", timepoint_label, "...\n")

  # Page 1: Cell type annotations (with/without labels)
  cat("  - Cell type annotation UMAP\n")
  p_annotations <- plot_annotation_comparison(seu_obj, timepoint_label)
  print(p_annotations)

  # Pages 2-7: Co-expression analysis (6 gene pairs)
  cat("  - Co-expression plots\n")

  # PTBP1 vs PTBP2
  print(makeCoExpressionPlot(seu_obj, "PTBP1", "PTBP2", title_prefix = timepoint_label))

  # PTBP2 vs DLG4
  print(makeCoExpressionPlot(seu_obj, "PTBP2", "DLG4", title_prefix = timepoint_label))

  # PTBP1 vs DLG4
  print(makeCoExpressionPlot(seu_obj, "PTBP1", "DLG4", title_prefix = timepoint_label))

  # RBFOX1 vs RBFOX2
  print(makeCoExpressionPlot(seu_obj, "RBFOX1", "RBFOX2", title_prefix = timepoint_label))

  # RBFOX2 vs RBFOX3
  print(makeCoExpressionPlot(seu_obj, "RBFOX2", "RBFOX3", title_prefix = timepoint_label))

  # RBFOX1 vs RBFOX3
  print(makeCoExpressionPlot(seu_obj, "RBFOX1", "RBFOX3", title_prefix = timepoint_label))

  cat("Completed co-expression report for", timepoint_label, "\n\n")
}

#' Generate ridge plot report for a single timepoint
#'
#' @param seu_obj Seurat object
#' @param timepoint_name Name of timepoint (e.g., "1M_Org", "3M_Org", "6M_Org")
#' @param timepoint_label Display label (e.g., "1 Month", "3 Month", "6 Month")
#' @description Generates chunked ridge plots with UMAP first page
#'              Designed to be called within a PDF device context
generate_ridge_report <- function(seu_obj, timepoint_name, timepoint_label) {
  cat("Generating ridge plots report for", timepoint_label, "...\n")

  # Page 1: Cell type annotations (with/without labels)
  cat("  - Cell type annotation UMAP\n")
  p_annotations <- plot_annotation_comparison(seu_obj, timepoint_label)
  print(p_annotations)

  # Pages 2+: Ridge plots (chunked to avoid squashing)
  cat("  - Ridge plots (chunked)\n")
  plot_ridge_feature_chunks(seu_obj, ALL_FEATURES, chunk_size = 5, ncol = 2)

  cat("Completed ridge plots report for", timepoint_label, "\n\n")
}

#' Generate violin plot report for a single timepoint
#'
#' @param seu_obj Seurat object
#' @param timepoint_name Name of timepoint (e.g., "1M_Org", "3M_Org", "6M_Org")
#' @param timepoint_label Display label (e.g., "1 Month", "3 Month", "6 Month")
#' @description Generates chunked violin plots with UMAP first page
#'              Designed to be called within a PDF device context
generate_violin_report <- function(seu_obj, timepoint_name, timepoint_label) {
  cat("Generating violin plots report for", timepoint_label, "...\n")

  # Page 1: Cell type annotations (with/without labels)
  cat("  - Cell type annotation UMAP\n")
  p_annotations <- plot_annotation_comparison(seu_obj, timepoint_label)
  print(p_annotations)

  # Pages 2+: Violin plots by chunks
  cat("  - Violin plots (chunked)\n")
  plot_feature_chunks(seu_obj, ALL_FEATURES, chunk_size = 4, ncol = 2)

  cat("Completed violin plots report for", timepoint_label, "\n\n")
}

# ============================================================================
# MAIN ANALYSIS
# ============================================================================

# Create output directory ----------------------------------------------------
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
  cat("Created output directory:", OUTPUT_DIR, "\n")
}

# Load data ------------------------------------------------------------------
cat("\n=== LOADING DATA ===\n")
seurat_objects <- load_seurat_objects()
one_month_org <- seurat_objects$one_month
three_month_org <- seurat_objects$three_month
six_month_org <- seurat_objects$six_month

# Optional: Check metadata for debugging -------------------------------------
# Uncomment to inspect metadata structure
# cat("\n=== METADATA CHECK ===\n")
# check_metadata(one_month_org, "1 Month")
# check_metadata(three_month_org, "3 Month")
# check_metadata(six_month_org, "6 Month")

# Apply trace for FeaturePlot (adjusts color granularity) -------------------
trace("FeaturePlot", quote(brewer.gran <- 2), at = 22L)

# Generate face validity reports for each timepoint --------------------------
cat("\n=== GENERATING FACE VALIDITY REPORTS ===\n")

# Define timepoint configurations
timepoints <- list(
  list(obj = one_month_org, name = "1M_Org", label = "1 Month"),
  list(obj = three_month_org, name = "3M_Org", label = "3 Month"),
  list(obj = six_month_org, name = "6M_Org", label = "6 Month")
)

# Generate 3 separate PDF reports per timepoint (9 total PDFs)

# 1. Co-expression reports
cat("\n--- Generating Co-expression Reports ---\n")
for (tp in timepoints) {
  pdf_file <- file.path(OUTPUT_DIR, paste0(tp$name, "_coexpression_report.pdf"))

  cat("\nGenerating PDF:", pdf_file, "\n")
  pdf(file = pdf_file, width = 14, height = 10)

  generate_coexpression_report(tp$obj, tp$name, tp$label)

  dev.off()
  cat("Saved:", pdf_file, "\n")
}

# 2. Ridge plot reports
cat("\n--- Generating Ridge Plot Reports ---\n")
for (tp in timepoints) {
  pdf_file <- file.path(OUTPUT_DIR, paste0(tp$name, "_ridge_plots_report.pdf"))

  cat("\nGenerating PDF:", pdf_file, "\n")
  pdf(file = pdf_file, width = 14, height = 10)

  generate_ridge_report(tp$obj, tp$name, tp$label)

  dev.off()
  cat("Saved:", pdf_file, "\n")
}

# 3. Violin plot reports
cat("\n--- Generating Violin Plot Reports ---\n")
for (tp in timepoints) {
  pdf_file <- file.path(OUTPUT_DIR, paste0(tp$name, "_violin_plots_report.pdf"))

  cat("\nGenerating PDF:", pdf_file, "\n")
  pdf(file = pdf_file, width = 14, height = 10)

  generate_violin_report(tp$obj, tp$name, tp$label)

  dev.off()
  cat("Saved:", pdf_file, "\n")
}

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("All reports saved to:", OUTPUT_DIR, "\n")
cat("\nGenerated files (9 total):\n")
cat("\nCo-expression reports:\n")
cat("  - 1M_Org_coexpression_report.pdf\n")
cat("  - 3M_Org_coexpression_report.pdf\n")
cat("  - 6M_Org_coexpression_report.pdf\n")
cat("\nRidge plot reports:\n")
cat("  - 1M_Org_ridge_plots_report.pdf\n")
cat("  - 3M_Org_ridge_plots_report.pdf\n")
cat("  - 6M_Org_ridge_plots_report.pdf\n")
cat("\nViolin plot reports:\n")
cat("  - 1M_Org_violin_plots_report.pdf\n")
cat("  - 3M_Org_violin_plots_report.pdf\n")
cat("  - 6M_Org_violin_plots_report.pdf\n")