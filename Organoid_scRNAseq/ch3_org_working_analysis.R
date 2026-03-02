## QC Analysis Script for Gene-Level Seurat Objects------
# Load essential packages
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(Matrix)
library(gridExtra)
library(grid)
library(DoubletFinder)
library(clustree)

# Load cell cycle markers
load("/data/gpfs/projects/punim0646/manveer/cycle.rda")

# Load silhouette functions
source("/data/gpfs/projects/punim2251/scRNAseq_scripts/silhouette.R")

# Set working directory
setwd("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space")
dir.create("./output_files/QC", recursive = TRUE, showWarnings = FALSE)

# Sample information
sample_ids <- c("org_1A", "org_1B", "org_3A", "org_3B", "org_3C", "org_6A", "org_6B", "org_6C")
timepoints <- c("1M_Org", "1M_Org", "3M_Org", "3M_Org", "3M_Org", "6M_Org", "6M_Org", "6M_Org")

# FLAMES output directory
FLAMES_OUTPUT_DIR <- "/data/gpfs/projects/punim2251/Neurodevelopmental_Models/Cortical_Org_Diff/flame_outs_standard_ref"

## STEP 1: Create gene-level count matrices with gene symbols------

# Load isoform-gene dictionary (contains gene_id to gene_symbol mapping)
cat("=== Loading isoform-gene dictionary ===\n")
isoform_gene_dict <- read.csv("./output_files/ref_files/isoform_gene_dict.csv")

# Extract unique gene_id to gene_symbol mapping
gene_symbol_dict <- isoform_gene_dict %>%
  select(gene_id, gene_symbol) %>%
  distinct()

cat("Loaded", nrow(gene_symbol_dict), "unique gene mappings\n")

# Function to convert gene IDs to symbols in gene-level count matrices
convert_gene_counts_to_symbols <- function(gene_count_file, sample_id, gene_dict) {

  cat("\nProcessing:", sample_id, "\n")

  # Read gene count matrix
  cat("  Reading count matrix...\n")
  gene_counts <- read.csv(gene_count_file, row.names = 1, check.names = FALSE)

  # Convert NAs to zeros (FLAMES uses empty cells for zero counts)
  gene_counts[is.na(gene_counts)] <- 0
  cat("  Converted NA values to zeros\n")

  cat("  Original dimensions:", nrow(gene_counts), "genes x", ncol(gene_counts), "cells\n")

  # Create data frame with gene IDs
  gene_counts_df <- gene_counts %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene_id")

  # Remove version numbers from gene IDs (e.g., ENSG00000123.1 -> ENSG00000123)
  gene_counts_df$gene_id_clean <- str_replace(gene_counts_df$gene_id, "\\.\\d+$", "")

  # Also clean the dictionary gene IDs
  gene_dict_clean <- gene_dict %>%
    mutate(gene_id_clean = str_replace(gene_id, "\\.\\d+$", ""))

  # Merge with gene symbols
  cat("  Merging with gene symbols...\n")
  gene_counts_with_symbols <- gene_counts_df %>%
    left_join(gene_dict_clean, by = "gene_id_clean") %>%
    select(-gene_id_clean)

  # Check how many genes got symbols
  n_with_symbols <- sum(!is.na(gene_counts_with_symbols$gene_symbol))
  n_without_symbols <- sum(is.na(gene_counts_with_symbols$gene_symbol))
  cat("  Genes with symbols:", n_with_symbols, "\n")
  cat("  Genes without symbols:", n_without_symbols, "(will use gene_id)\n")

  # For genes without symbols, use the gene_id
  gene_counts_with_symbols <- gene_counts_with_symbols %>%
    mutate(gene_symbol = ifelse(is.na(gene_symbol), gene_id.x, gene_symbol))

  # Create combined gene_id_gene_symbol format for rownames (ensures uniqueness)
  gene_counts_with_symbols <- gene_counts_with_symbols %>%
    mutate(gene_name = paste0(gene_id.x, "_", gene_symbol))

  # Set combined names as rownames
  rownames(gene_counts_with_symbols) <- gene_counts_with_symbols$gene_name

  # Remove gene_id and gene_symbol columns, keeping only count data
  gene_counts_final <- gene_counts_with_symbols %>%
    select(-gene_id.x, -gene_id.y, -gene_symbol, -gene_name)

  cat("  Final dimensions:", nrow(gene_counts_final), "genes x", ncol(gene_counts_final), "cells\n")

  return(gene_counts_final)
}

# Process all samples and create gene symbol count matrices
cat("\n=== Converting gene IDs to symbols for all samples ===\n")

gene_symbol_count_matrices <- list()

for(i in seq_along(sample_ids)) {
  sample_id <- sample_ids[i]
  gene_count_file <- file.path(FLAMES_OUTPUT_DIR, paste0(sample_id, "_gene_count.csv"))

  gene_symbol_count_matrices[[sample_id]] <- convert_gene_counts_to_symbols(
    gene_count_file = gene_count_file,
    sample_id = sample_id,
    gene_dict = gene_symbol_dict
  )
}

cat("\n=== All count matrices converted successfully ===\n")

## STEP 2: Create Seurat objects from gene symbol count matrices------

# Function to create Seurat objects (adapted from ch3_qc_flames_v220_p2.R)
create_seurat_from_counts <- function(count_matrix_df, sample_id,
                                      assay_type = "RNA",
                                      MIN.CELL = 3,
                                      MIN.FEATURES = 20) {

  cat("\n=== Creating Seurat object for:", sample_id, "===\n")

  # Convert to matrix
  count_matrix <- as.matrix(count_matrix_df)

  # Check for and handle duplicate rownames
  if(any(duplicated(rownames(count_matrix)))) {
    cat("Warning: Found", sum(duplicated(rownames(count_matrix))), "duplicate gene names in", sample_id, "\n")
    cat("Making rownames unique...\n")
    rownames(count_matrix) <- make.unique(rownames(count_matrix))
  }

  # PRE-FILTER before converting to sparse matrix
  cat("Filtering genes: keeping genes detected in at least", MIN.CELL, "cells\n")
  genes_keep <- rowSums(count_matrix > 0, na.rm = TRUE) >= MIN.CELL
  count_matrix <- count_matrix[genes_keep, ]

  cat("Filtering cells: keeping cells with at least", MIN.FEATURES, "features\n")
  cells_keep <- colSums(count_matrix > 0, na.rm = TRUE) >= MIN.FEATURES
  count_matrix <- count_matrix[, cells_keep]

  cat("Filtered dimensions:", nrow(count_matrix), "genes x", ncol(count_matrix), "cells\n")

  # Convert to sparse matrix
  count_sparse <- Matrix::Matrix(count_matrix, sparse = TRUE)

  # Create Seurat object WITHOUT additional filtering (already filtered)
  seu_obj <- CreateSeuratObject(
    counts = count_sparse,
    project = sample_id,
    assay = assay_type
  )

  # Add sample metadata
  seu_obj$sample_id <- sample_id

  cat("Seurat object created successfully!\n")

  return(seu_obj)
}

# Create and save Seurat objects for all samples
cat("\n=== Creating Seurat objects for all samples ===\n")

seurat_objects <- list()

for(i in seq_along(sample_ids)) {
  sample_id <- sample_ids[i]

  # Create Seurat object
  seurat_objects[[sample_id]] <- create_seurat_from_counts(
    count_matrix_df = gene_symbol_count_matrices[[sample_id]],
    sample_id = sample_id,
    assay_type = "RNA",
    MIN.CELL = 3,
    MIN.FEATURES = 20
  )

  # Add timepoint metadata
  seurat_objects[[sample_id]]$timepoint <- timepoints[i]

  # Save immediately to disk
  output_file <- paste0("./output_files/seu_objects/", sample_id, "_gene_seurat_with_symbols.rds")
  saveRDS(seurat_objects[[sample_id]], file = output_file)
  cat("Saved:", output_file, "\n")

  # Print summary
  print(seurat_objects[[sample_id]])
}

cat("\n========================================\n")
cat("=== All Seurat objects created and saved successfully! ===\n")
cat("========================================\n")

## STEP 3: QC Analysis on Gene-Level Seurat Objects------

cat("\n=== Starting QC Analysis ===\n")

# Calculate mitochondrial and ribosomal percentages for all samples
for(sample_id in sample_ids) {
  cat("\nCalculating QC metrics for:", sample_id, "\n")

  # Mitochondrial genes (pattern matches MT- anywhere in feature name)
  # Feature names are formatted as: ENSG00000198888.2-MT-ND1
  seurat_objects[[sample_id]][["percent.mt"]] <- PercentageFeatureSet(
    seurat_objects[[sample_id]],
    pattern = "-MT-"
  )

  # Ribosomal genes (pattern matches RP[SL] anywhere in feature name)
  seurat_objects[[sample_id]][["percent.ribo"]] <- PercentageFeatureSet(
    seurat_objects[[sample_id]],
    pattern = "-RP[SL]"
  )

  cat("  Median percent.mt:", median(seurat_objects[[sample_id]]$percent.mt), "\n")
  cat("  Median percent.ribo:", median(seurat_objects[[sample_id]]$percent.ribo), "\n")
}

# Create QC violin plots for all samples
cat("\n=== Creating QC violin plots ===\n")

for(sample_id in sample_ids) {
  cat("Plotting:", sample_id, "\n")

  # Use VlnPlot with default labels
  # Default labels are clear: nFeature_RNA, nCount_RNA, percent.mt, percent.ribo
  p <- VlnPlot(
    seurat_objects[[sample_id]],
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
    ncol = 4,
    pt.size = 0.1
  )

  # Use png() device instead of ggsave() to avoid S4SXP error with Seurat v5
  output_file <- paste0("./output_files/QC/", sample_id, "_QC_violin.png")
  png(output_file, width = 16*300, height = 6*300, res = 300)
  print(p)
  dev.off()
  cat("  Saved:", output_file, "\n")
}

# Create scatter plots showing relationships between QC metrics
cat("\n=== Creating QC scatter plots ===\n")

for(sample_id in sample_ids) {
  cat("Plotting:", sample_id, "\n")

  # nCount vs nFeature - with intuitive labels
  p1 <- FeatureScatter(
    seurat_objects[[sample_id]],
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
  ) +
    xlab("Total UMI Counts (Library Size)") +
    ylab("Genes Detected per Cell") +
    ggtitle("Library Complexity")

  # nCount vs percent.mt - with intuitive labels
  p2 <- FeatureScatter(
    seurat_objects[[sample_id]],
    feature1 = "nCount_RNA",
    feature2 = "percent.mt"
  ) +
    xlab("Total UMI Counts (Library Size)") +
    ylab("Mitochondrial Content (%)") +
    ggtitle("Cell Viability")

  # Combine plots
  p <- p1 + p2 + plot_annotation(title = paste("QC Metrics -", sample_id))

  # Save plot
  output_file <- paste0("./output_files/QC/", sample_id, "_QC_scatter.png")
  ggsave(output_file, plot = p, width = 14, height = 6, dpi = 300)
  cat("  Saved:", output_file, "\n")
}

# Generate summary statistics table
cat("\n=== Generating summary statistics ===\n")

qc_summary <- data.frame()

for(i in seq_along(sample_ids)) {
  sample_id <- sample_ids[i]
  seu <- seurat_objects[[sample_id]]

  qc_summary <- rbind(qc_summary, data.frame(
    sample_id = sample_id,
    timepoint = timepoints[i],
    n_cells = ncol(seu),
    n_genes = nrow(seu),
    median_nCount = median(seu$nCount_RNA),
    median_nFeature = median(seu$nFeature_RNA),
    median_percent_mt = median(seu$percent.mt),
    median_percent_ribo = median(seu$percent.ribo),
    mean_nCount = mean(seu$nCount_RNA),
    mean_nFeature = mean(seu$nFeature_RNA),
    mean_percent_mt = mean(seu$percent.mt),
    mean_percent_ribo = mean(seu$percent.ribo)
  ))
}


# Print summary
cat("\n=== QC Summary Statistics ===\n")
print(qc_summary)

# Save summary table
write.csv(qc_summary, file = "./output_files/QC/qc_summary_statistics.csv", row.names = FALSE)
cat("\nSaved: ./output_files/QC/qc_summary_statistics.csv\n")

cat("\n========================================\n")
cat("=== QC Analysis Complete! ===\n")
cat("========================================\n")
cat("\nSeurat objects with QC metrics are stored in 'seurat_objects' list\n")
cat("QC plots saved in: ./output_files/QC/\n")


## STEP 4: Comprehensive Per-Sample QC Pipeline------
# Define QC Functions------

## Helper Function: Find Dominant Peak in Distribution
find_dominant_peak <- function(values, bw = "nrd0") {
  # Calculate density
  dens <- density(values, bw = bw, na.rm = TRUE)

  # Find local maxima (peaks) using derivative method
  # A peak is where the derivative changes from positive to negative
  peak_idx <- which(diff(sign(diff(dens$y))) == -2) + 1

  if(length(peak_idx) == 0) {
    # No peaks found, use mean as fallback
    return(list(
      center = mean(values, na.rm = TRUE),
      method = "mean",
      n_peaks_total = 0,
      n_peaks_top2 = 0,
      all_peaks = NULL,
      top2_peaks = NULL,
      selected_peak = mean(values, na.rm = TRUE),
      density_obj = dens
    ))
  }

  # Store all peak locations and their densities
  all_peak_values <- dens$x[peak_idx]
  all_peak_densities <- dens$y[peak_idx]

  # REFINEMENT: Select top 2 peaks by HEIGHT (density)
  # Then pick the RIGHTMOST one from those top 2
  if(length(peak_idx) == 1) {
    # Only one peak found, use it
    peak_value <- all_peak_values[1]
    top2_peak_values <- all_peak_values[1]
    n_top2 <- 1
  } else {
    # Get indices of top 2 peaks by density (height)
    top2_idx <- order(all_peak_densities, decreasing = TRUE)[1:min(2, length(all_peak_densities))]
    top2_peak_values <- all_peak_values[top2_idx]
    n_top2 <- length(top2_peak_values)

    # From top 2 peaks, select the RIGHTMOST one (highest x-value)
    peak_value <- max(top2_peak_values)
  }

  return(list(
    center = peak_value,
    method = "peak",
    n_peaks_total = length(peak_idx),
    n_peaks_top2 = n_top2,
    all_peaks = all_peak_values,
    top2_peaks = top2_peak_values,
    selected_peak = peak_value,
    density_obj = dens
  ))
}

## Helper Function: Plot Threshold Diagnostics
plot_threshold_diagnostics <- function(values, metric_name, peak_result,
                                       min_threshold, max_threshold, sample_id) {

  # Create data frame for ggplot
  plot_df <- data.frame(value = values)

  # Base histogram with density overlay
  p <- ggplot(plot_df, aes(x = value)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "lightblue", alpha = 0.6) +
    geom_density(color = "black", linewidth = 1) +
    theme_classic() +
    labs(
      title = paste0(sample_id, ": ", metric_name, " Distribution"),
      subtitle = paste0("Method: ", peak_result$method,
                       " | Selected peak: ", round(peak_result$center, 1),
                       " | Total peaks: ", peak_result$n_peaks_total,
                       " | Top 2 peaks used"),
      x = metric_name,
      y = "Density"
    )

  # Add all detected peaks as small gray points
  if(!is.null(peak_result$all_peaks) && length(peak_result$all_peaks) > 0) {
    # Get density values at all peak locations
    all_peak_densities <- approx(peak_result$density_obj$x,
                                 peak_result$density_obj$y,
                                 xout = peak_result$all_peaks)$y

    all_peak_df <- data.frame(
      x = peak_result$all_peaks,
      y = all_peak_densities
    )

    p <- p + geom_point(data = all_peak_df, aes(x = x, y = y),
                       color = "gray60", size = 2, shape = 19)
  }

  # Add top 2 peaks as larger orange points
  if(!is.null(peak_result$top2_peaks) && length(peak_result$top2_peaks) > 0) {
    top2_densities <- approx(peak_result$density_obj$x,
                             peak_result$density_obj$y,
                             xout = peak_result$top2_peaks)$y

    top2_df <- data.frame(
      x = peak_result$top2_peaks,
      y = top2_densities
    )

    p <- p + geom_point(data = top2_df, aes(x = x, y = y),
                       color = "orange", size = 4, shape = 19)
  }

  # Add vertical line for selected center (rightmost of top 2)
  p <- p + geom_vline(xintercept = peak_result$center,
                     color = "darkgreen", linetype = "dashed", linewidth = 1.2)

  # Add threshold lines
  p <- p + geom_vline(xintercept = min_threshold,
                     color = "red", linetype = "dotted", linewidth = 0.8) +
          geom_vline(xintercept = max_threshold,
                     color = "red", linetype = "dotted", linewidth = 0.8)

  # Add text annotations for thresholds
  y_max <- max(peak_result$density_obj$y) * 0.95
  p <- p + annotate("text", x = min_threshold, y = y_max,
                   label = paste0("MIN: ", round(min_threshold)),
                   angle = 90, vjust = -0.5, size = 3, color = "red") +
          annotate("text", x = max_threshold, y = y_max,
                   label = paste0("MAX: ", round(max_threshold)),
                   angle = 90, vjust = 1.5, size = 3, color = "red")

  return(p)
}

## Function 1: Cell/Gene Filtering with before/after plots
filter_cells_genes <- function(seurat_obj, sample_params) {

  sample_id <- sample_params$sample_id
  cat("\n=== STEP 1: Filtering cells for", sample_id, "===\n")

  # Store initial stats
  initial_cells <- ncol(seurat_obj)
  initial_genes <- nrow(seurat_obj)

  # Calculate MT percentage (already done in STEP 3, but ensure it exists)
  if(!"percent.mt" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "-MT-")
  }

  # Calculate adaptive thresholds or use custom - WITH PEAK DETECTION
  # For nFeature_RNA (genes per cell)
  if(is.na(sample_params$custom_min_genes) || is.na(sample_params$custom_max_genes)) {
    cat("  Detecting peaks in nFeature_RNA distribution...\n")
    peak_result_genes <- find_dominant_peak(seurat_obj$nFeature_RNA)
    center_genes <- peak_result_genes$center
    sd_genes <- sd(seurat_obj$nFeature_RNA)

    cat("  nFeature_RNA - Method:", peak_result_genes$method,
        "| Selected peak:", round(center_genes),
        "| SD:", round(sd_genes), "\n")
    cat("    Total peaks detected:", peak_result_genes$n_peaks_total,
        "| Top 2 peaks:", paste(round(peak_result_genes$top2_peaks), collapse = ", "), "\n")
  }

  if(is.na(sample_params$custom_min_genes)) {
    MIN_GENES <- round(center_genes - (1.5 * sd_genes))
    cat("  Using adaptive MIN_GENES:", MIN_GENES, "\n")
  } else {
    MIN_GENES <- sample_params$custom_min_genes
    peak_result_genes <- NULL  # No peak detection if custom value provided
    cat("  Using custom MIN_GENES:", MIN_GENES, "\n")
  }

  if(is.na(sample_params$custom_max_genes)) {
    MAX_GENES <- round(center_genes + (1.5 * sd_genes))
    cat("  Using adaptive MAX_GENES:", MAX_GENES, "\n")
  } else {
    MAX_GENES <- sample_params$custom_max_genes
    cat("  Using custom MAX_GENES:", MAX_GENES, "\n")
  }

  # For nCount_RNA (UMI counts per cell)
  if(is.na(sample_params$custom_max_counts)) {
    cat("  Detecting peaks in nCount_RNA distribution...\n")
    peak_result_counts <- find_dominant_peak(seurat_obj$nCount_RNA)
    center_counts <- peak_result_counts$center
    sd_counts <- sd(seurat_obj$nCount_RNA)

    cat("  nCount_RNA - Method:", peak_result_counts$method,
        "| Selected peak:", round(center_counts),
        "| SD:", round(sd_counts), "\n")
    cat("    Total peaks detected:", peak_result_counts$n_peaks_total,
        "| Top 2 peaks:", paste(round(peak_result_counts$top2_peaks), collapse = ", "), "\n")

    MAX_COUNTS <- round(center_counts + (1.5 * sd_counts))
    cat("  Using adaptive MAX_COUNTS:", MAX_COUNTS, "\n")
  } else {
    MAX_COUNTS <- sample_params$custom_max_counts
    peak_result_counts <- NULL  # No peak detection if custom value provided
    cat("  Using custom MAX_COUNTS:", MAX_COUNTS, "\n")
  }

  MIN_COUNTS <- sample_params$min_counts
  MT_THRESHOLD <- sample_params$mt_threshold

  # Create before-filtering plots
  vln_before_1 <- VlnPlot(seurat_obj, features = "nFeature_RNA", pt.size = 0.1) +
    ggtitle(paste0(sample_id, ": Genes Before Filtering")) + NoLegend()
  vln_before_2 <- VlnPlot(seurat_obj, features = "nCount_RNA", pt.size = 0.1) +
    ggtitle(paste0(sample_id, ": Counts Before Filtering")) + NoLegend()
  vln_before_3 <- VlnPlot(seurat_obj, features = "percent.mt", pt.size = 0.1) +
    ggtitle(paste0(sample_id, ": MT% Before Filtering")) + NoLegend()

  scatter_before <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_smooth(method = "lm") +
    ggtitle(paste0(sample_id, ": Before Filtering")) +
    NoLegend()

  # Apply filtering
  seurat_obj_filtered <- subset(seurat_obj,
                                subset = nFeature_RNA > MIN_GENES & nFeature_RNA < MAX_GENES &
                                         nCount_RNA > MIN_COUNTS & nCount_RNA < MAX_COUNTS &
                                         percent.mt < MT_THRESHOLD)

  # Calculate filtering stats
  final_cells <- ncol(seurat_obj_filtered)
  cells_removed <- initial_cells - final_cells
  pct_cells_kept <- round((final_cells / initial_cells) * 100, 2)

  cat("  Cells before:", initial_cells, "\n")
  cat("  Cells after:", final_cells, "\n")
  cat("  Cells removed:", cells_removed, "(", 100-pct_cells_kept, "%)\n")

  # Create after-filtering plots
  vln_after_1 <- VlnPlot(seurat_obj_filtered, features = "nFeature_RNA", pt.size = 0.1) +
    ggtitle(paste0(sample_id, ": Genes After Filtering")) + NoLegend()
  vln_after_2 <- VlnPlot(seurat_obj_filtered, features = "nCount_RNA", pt.size = 0.1) +
    ggtitle(paste0(sample_id, ": Counts After Filtering")) + NoLegend()
  vln_after_3 <- VlnPlot(seurat_obj_filtered, features = "percent.mt", pt.size = 0.1) +
    ggtitle(paste0(sample_id, ": MT% After Filtering")) + NoLegend()

  scatter_after <- FeatureScatter(seurat_obj_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_smooth(method = "lm") +
    ggtitle(paste0(sample_id, ": After Filtering")) +
    NoLegend()

  # Create summary table
  filter_summary <- data.frame(
    Metric = c("Cells", "Genes", "Median nCount", "Median nFeature", "Median MT%"),
    Before = c(initial_cells, initial_genes,
               median(seurat_obj$nCount_RNA),
               median(seurat_obj$nFeature_RNA),
               median(seurat_obj$percent.mt)),
    After = c(final_cells, nrow(seurat_obj_filtered),
              median(seurat_obj_filtered$nCount_RNA),
              median(seurat_obj_filtered$nFeature_RNA),
              median(seurat_obj_filtered$percent.mt))
  )

  # Create diagnostic plots for threshold determination
  diagnostic_plots <- list()

  if(!is.null(peak_result_genes)) {
    cat("  Generating nFeature_RNA diagnostic plot...\n")
    diagnostic_plots$genes_dist <- plot_threshold_diagnostics(
      values = seurat_obj$nFeature_RNA,
      metric_name = "nFeature_RNA (Genes per Cell)",
      peak_result = peak_result_genes,
      min_threshold = MIN_GENES,
      max_threshold = MAX_GENES,
      sample_id = sample_id
    )
  }

  if(!is.null(peak_result_counts)) {
    cat("  Generating nCount_RNA diagnostic plot...\n")
    diagnostic_plots$counts_dist <- plot_threshold_diagnostics(
      values = seurat_obj$nCount_RNA,
      metric_name = "nCount_RNA (UMI Counts per Cell)",
      peak_result = peak_result_counts,
      min_threshold = MIN_COUNTS,
      max_threshold = MAX_COUNTS,
      sample_id = sample_id
    )
  }

  # Create threshold summary table
  threshold_summary <- data.frame(
    Metric = c("nFeature_RNA", "nCount_RNA"),
    Method = c(
      ifelse(!is.null(peak_result_genes), peak_result_genes$method, "custom"),
      ifelse(!is.null(peak_result_counts), peak_result_counts$method, "custom")
    ),
    Selected_Peak = c(
      ifelse(!is.null(peak_result_genes), round(peak_result_genes$center), NA),
      ifelse(!is.null(peak_result_counts), round(peak_result_counts$center), NA)
    ),
    Total_Peaks = c(
      ifelse(!is.null(peak_result_genes), peak_result_genes$n_peaks_total, NA),
      ifelse(!is.null(peak_result_counts), peak_result_counts$n_peaks_total, NA)
    ),
    Top2_Peaks = c(
      ifelse(!is.null(peak_result_genes),
             paste(round(peak_result_genes$top2_peaks), collapse = ", "),
             NA),
      ifelse(!is.null(peak_result_counts),
             paste(round(peak_result_counts$top2_peaks), collapse = ", "),
             NA)
    ),
    SD = c(
      ifelse(!is.null(peak_result_genes), round(sd_genes), NA),
      ifelse(!is.null(peak_result_counts), round(sd_counts), NA)
    ),
    MIN_Threshold = c(MIN_GENES, MIN_COUNTS),
    MAX_Threshold = c(MAX_GENES, MAX_COUNTS)
  )

  # Return list with filtered object and plots
  return(list(
    seurat_obj = seurat_obj_filtered,
    plots = list(
      vln_before = list(vln_before_1, vln_before_2, vln_before_3),
      vln_after = list(vln_after_1, vln_after_2, vln_after_3),
      scatter_before = scatter_before,
      scatter_after = scatter_after,
      diagnostics = diagnostic_plots
    ),
    summary = filter_summary,
    threshold_summary = threshold_summary,
    thresholds = list(
      MIN_GENES = MIN_GENES,
      MAX_GENES = MAX_GENES,
      MIN_COUNTS = MIN_COUNTS,
      MAX_COUNTS = MAX_COUNTS,
      MT_THRESHOLD = MT_THRESHOLD
    )
  ))
}

# Source additional QC functions
source("qc_functions.R")


# QC Parameters Configuration (EASY TO MODIFY FOR ITERATION)------
qc_params <- tibble(
  sample_id = c("org_1A", "org_1B", "org_3A", "org_3B", "org_3C", "org_6A", "org_6B", "org_6C"),
  timepoint = c("1M_Org", "1M_Org", "3M_Org", "3M_Org", "3M_Org", "6M_Org", "6M_Org", "6M_Org"),

  # QC thresholds (NULL = use adaptive, or specify custom values)
  custom_min_genes = c(4000, 2500, NA, NA, NA, 2500, NA, 2500),
  custom_max_genes = c(NA, NA, NA, NA, NA, NA, NA, 10000),
  custom_max_counts = c(NA, NA, NA, NA, NA, 50000, 70000, NA),
  min_counts = c(500, 9500, 500, 500, 500, 500, 500, 500),
  mt_threshold = c(10, 10, 10, 10, 10, 10, 10, 10),
  doublet_rate = c(0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016)
)

cat("\n=== QC Parameters Configuration Loaded ===\n")
print(qc_params)


## Master Function: Process Complete QC Pipeline for One Sample
process_sample_qc <- function(sample_params, overwrite_intermediates = TRUE) {

  sample_id <- sample_params$sample_id
  cat("\n\n========================================\n")
  cat("PROCESSING SAMPLE:", sample_id, "\n")
  cat("========================================\n")

  # Create output directory for this sample
  sample_qc_dir <- paste0("./output_files/QC/", sample_id)
  dir.create(sample_qc_dir, recursive = TRUE, showWarnings = FALSE)

  # Define intermediate file paths
  filtered_path <- paste0("./output_files/seu_objects/", sample_id, "_filtered.rds")
  normalized_path <- paste0("./output_files/seu_objects/", sample_id, "_normalized.rds")
  scaled_pca_path <- paste0("./output_files/seu_objects/", sample_id, "_scaled_pca.rds")
  clustered_path <- paste0("./output_files/seu_objects/", sample_id, "_clustered.rds")
  final_path <- paste0("./output_files/seu_objects/", sample_id, "_QCed_final.rds")

  # STEP 1: Load initial Seurat object with gene symbols
  cat("\nLoading initial Seurat object...\n")
  input_path <- paste0("./output_files/seu_objects/", sample_id, "_gene_seurat_with_symbols.rds")
  seurat_obj <- readRDS(input_path)
  cat("Loaded:", ncol(seurat_obj), "cells,", nrow(seurat_obj), "genes\n")

  # STEP 2: Filter cells and genes
  if(!file.exists(filtered_path) || overwrite_intermediates) {
    filter_result <- filter_cells_genes(seurat_obj, sample_params)
    saveRDS(filter_result$seurat_obj, filtered_path)
    cat("Saved filtered object to:", filtered_path, "\n")
  } else {
    cat("Loading existing filtered object...\n")
    filter_result <- list(seurat_obj = readRDS(filtered_path))
  }

  # STEP 3: Normalization and variable features
  if(!file.exists(normalized_path) || overwrite_intermediates) {
    norm_result <- normalize_and_find_features(filter_result$seurat_obj, sample_id)
    saveRDS(norm_result$seurat_obj, normalized_path)
    cat("Saved normalized object to:", normalized_path, "\n")
  } else {
    cat("Loading existing normalized object...\n")
    norm_result <- list(seurat_obj = readRDS(normalized_path))
  }

  # STEP 4: Cell cycle scoring and scaling
  if(!file.exists(scaled_pca_path) || overwrite_intermediates) {
    cc_result <- cell_cycle_and_scale(norm_result$seurat_obj, sample_id)
    saveRDS(cc_result$seurat_obj, scaled_pca_path)
    cat("Saved scaled/PCA object to:", scaled_pca_path, "\n")
  } else {
    cat("Loading existing scaled/PCA object...\n")
    cc_result <- list(seurat_obj = readRDS(scaled_pca_path))
  }

  # STEP 5: Quantitative elbow analysis (with +5 PC adjustment)
  elbow_result <- quantitative_elbow(cc_result$seurat_obj, sample_id, max_pcs = 50, pc_adjustment = 5)
  optimal_pcs <- elbow_result$optimal_pcs_used  # Use the adjusted value

  # STEP 6: Clustering optimization
  if(!file.exists(clustered_path) || overwrite_intermediates) {
    cluster_result <- optimize_clustering(cc_result$seurat_obj, sample_id, optimal_pcs)
    saveRDS(cluster_result$seurat_obj, clustered_path)
    cat("Saved clustered object to:", clustered_path, "\n")
  } else {
    cat("Loading existing clustered object...\n")
    cluster_result <- list(seurat_obj = readRDS(clustered_path))
  }

  # STEP 7: UMAP and doublet removal
  if(!file.exists(final_path) || overwrite_intermediates) {
    doublet_result <- umap_and_doublets(cluster_result$seurat_obj,
                                        sample_id,
                                        optimal_pcs,
                                        sample_params$doublet_rate)
    saveRDS(doublet_result$seurat_obj, final_path)
    cat("Saved final QC'd object to:", final_path, "\n")
  } else {
    cat("Loading existing final object...\n")
    doublet_result <- list(seurat_obj = readRDS(final_path))
  }

  cat("\n========================================\n")
  cat("QC COMPLETE FOR:", sample_id, "\n")
  cat("========================================\n")

  # Return all results for PDF generation
  return(list(
    sample_id = sample_id,
    filter_result = filter_result,
    norm_result = norm_result,
    cc_result = cc_result,
    elbow_result = elbow_result,
    cluster_result = cluster_result,
    doublet_result = doublet_result,
    final_seurat = doublet_result$seurat_obj
  ))
}

## Function to Generate Comprehensive PDF Report
generate_qc_pdf <- function(qc_results) {

  sample_id <- qc_results$sample_id
  pdf_file <- paste0("./output_files/QC/", sample_id, "_complete_QC_report.pdf")

  cat("\nGenerating PDF report for", sample_id, "...\n")

  pdf(pdf_file, width = 14, height = 10)

  # Page 1: Threshold Diagnostics (NEW PAGE - FIRST!)
  if(length(qc_results$filter_result$plots$diagnostics) > 0) {
    # Build layout for diagnostic page
    diag_plots <- list()

    if(!is.null(qc_results$filter_result$plots$diagnostics$genes_dist)) {
      diag_plots[[length(diag_plots) + 1]] <- qc_results$filter_result$plots$diagnostics$genes_dist
    }

    if(!is.null(qc_results$filter_result$plots$diagnostics$counts_dist)) {
      diag_plots[[length(diag_plots) + 1]] <- qc_results$filter_result$plots$diagnostics$counts_dist
    }

    # Add threshold summary table
    diag_plots[[length(diag_plots) + 1]] <- tableGrob(
      qc_results$filter_result$threshold_summary,
      rows = NULL,
      theme = ttheme_default(base_size = 10)
    )

    # Arrange with 2 columns
    grid.arrange(
      grobs = diag_plots,
      ncol = 2,
      top = textGrob(paste(sample_id, "- Adaptive Threshold Determination"),
                     gp = gpar(fontsize = 14, fontface = "bold"))
    )
  }

  # Page 2: Filtering Summary
  grid.arrange(
    qc_results$filter_result$plots$vln_before[[1]],
    qc_results$filter_result$plots$vln_before[[2]],
    qc_results$filter_result$plots$vln_before[[3]],
    qc_results$filter_result$plots$vln_after[[1]],
    qc_results$filter_result$plots$vln_after[[2]],
    qc_results$filter_result$plots$vln_after[[3]],
    qc_results$filter_result$plots$scatter_before,
    qc_results$filter_result$plots$scatter_after,
    tableGrob(qc_results$filter_result$summary),
    nrow = 3, ncol = 3,
    top = textGrob(paste(sample_id, "- Filtering Summary"),
                   gp = gpar(fontsize = 14, fontface = "bold"))
  )

  # Page 3: Variable Features & Cell Cycle
  grid.arrange(
    qc_results$norm_result$plots$var_feature_plot,
    qc_results$cc_result$plots$cc_pc1_pc2,
    qc_results$cc_result$plots$cc_pc2_pc3,
    nrow = 2, ncol = 2,
    top = textGrob(paste(sample_id, "- Variable Features & Cell Cycle"),
                   gp = gpar(fontsize = 14, fontface = "bold"))
  )

  # Page 4: PCA Selection & Silhouette Analysis
  grid.arrange(
    qc_results$elbow_result$plots$elbow_plot,
    qc_results$elbow_result$plots$quant_plot,
    tableGrob(head(qc_results$cluster_result$sil_results, 10)),
    nrow = 2, ncol = 2,
    top = textGrob(paste(sample_id, "- PCA Selection & Silhouette Analysis"),
                   gp = gpar(fontsize = 14, fontface = "bold"))
  )

  # Page 5: Clustree
  grid.arrange(
    qc_results$cluster_result$plots$clustree,
    nrow = 1, ncol = 1,
    top = textGrob(paste(sample_id, "- Clustering Resolution Tree"),
                   gp = gpar(fontsize = 14, fontface = "bold"))
  )

  # Page 6: UMAP & Doublets
  grid.arrange(
    qc_results$doublet_result$plots$doublet_umap,
    qc_results$doublet_result$plots$final_umap,
    tableGrob(qc_results$doublet_result$doublet_stats),
    nrow = 2, ncol = 2,
    top = textGrob(paste(sample_id, "- UMAP & Doublet Removal"),
                   gp = gpar(fontsize = 14, fontface = "bold"))
  )

  # Page 7: QC Metrics & Cell Cycle on UMAP
  grid.arrange(
    qc_results$doublet_result$plots$nfeature_umap,
    qc_results$doublet_result$plots$ncount_umap,
    qc_results$doublet_result$plots$cellcycle_umap,
    qc_results$doublet_result$plots$final_umap,
    nrow = 2, ncol = 2,
    top = textGrob(paste(sample_id, "- QC Metrics & Cell Cycle on UMAP"),
                   gp = gpar(fontsize = 14, fontface = "bold"))
  )

  dev.off()

  cat("PDF report saved:", pdf_file, "\n")
}

## Execute QC pipeline for all samples------
cat("\n\n####################################################\n")
cat("STARTING COMPREHENSIVE QC PIPELINE FOR ALL SAMPLES\n")
cat("####################################################\n\n")

# Process each sample with QC pipeline
qc_results_list <- list()
for(i in 1:nrow(qc_params)) {

  cat("\n========== PROCESSING SAMPLE", i, "OF", nrow(qc_params), "==========\n")

  # Run full QC pipeline for this sample
  qc_results_list[[i]] <- process_sample_qc(qc_params[i,],
                                             overwrite_intermediates = TRUE)

  # Generate comprehensive PDF report
  generate_qc_pdf(qc_results_list[[i]])

  # Force garbage collection between samples
  gc()
}

cat("\n\n####################################################\n")
cat("ALL SAMPLES PROCESSED SUCCESSFULLY\n")
cat("####################################################\n")

# Print summary
cat("\n=== PROCESSING SUMMARY ===\n")
for(i in 1:length(qc_results_list)) {
  result <- qc_results_list[[i]]
  cat("\nSample:", result$sample_id, "\n")
  cat("  Cells after filtering:", ncol(result$seurat_obj), "\n")
  cat("  Optimal PCs:", result$optimal_pcs, "\n")
  cat("  Optimal resolution:", result$optimal_res, "\n")
  cat("  PDF report:", result$pdf_file, "\n")
}

cat("\n=== ALL QC ANALYSIS COMPLETE ===\n")
cat("Intermediate objects saved in: ./output_files/qc_intermediate/\n")
cat("PDF reports saved in: ./output_files/qc_reports/\n")

