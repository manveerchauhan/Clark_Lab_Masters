## Load supporting resources------
# Load cell cycle markers
load("/data/gpfs/projects/punim0646/manveer/cycle.rda")
# Import silhouette functions
source("/data/gpfs/projects/punim2251/scRNAseq_scripts/silhouette.R")

# FLAMES output directory for current project
FLAMES_OUTPUT_DIR <- "/data/gpfs/projects/punim2251/Neurodevelopmental_Models/Cortical_Org_Diff/flame_outs_standard_ref"
# Prefix for reading in formatted iso level count matrices from previous run
FORMATTED_ISO_MATRICES = "/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space/"

# Load ONLY essential packages
library(Seurat)      # For Seurat objects
library(Matrix)      # For sparse matrices
library(tidyverse)   # For data manipulation (tibble, purrr, dplyr)

# Set working directory and create folders for output files
setwd("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space")
dir.create("./output_files/seu_objects", recursive = TRUE, showWarnings = FALSE)

# Restore the sample metadata object
sample_metadata <- readRDS(file = "org_flames_v2.2.0_experiment_data.rds")

# Function to create Seurat objects------
create_seurat_from_counts <- function(count_matrix_df, sample_id, 
                                      assay_type = "RNA",
                                      MIN.CELL = 3,
                                      MIN.FEATURES = 1) {
  
  # Set the first column as rownames and remove it
  gene_names <- count_matrix_df[[1]]
  count_matrix_df <- count_matrix_df[, -1, drop = FALSE]
  
  # Check for and handle duplicate rownames
  if(any(duplicated(gene_names))) {
    cat("Warning: Found", sum(duplicated(gene_names)), "duplicate gene names in", sample_id, "\n")
    cat("Making rownames unique...\n")
    gene_names <- make.unique(gene_names)
  }
  
  rownames(count_matrix_df) <- gene_names
  
  # Convert to matrix
  count_matrix <- as.matrix(count_matrix_df)
  
  # PRE-FILTER before converting to sparse matrix
  # Filter genes: keep genes detected in at least MIN.CELL cells
  genes_keep <- rowSums(count_matrix > 0, na.rm = TRUE) >= MIN.CELL
  count_matrix <- count_matrix[genes_keep, ]
  
  # Filter cells: keep cells with at least MIN.FEATURES features
  cells_keep <- colSums(count_matrix > 0, na.rm = TRUE) >= MIN.FEATURES
  count_matrix <- count_matrix[, cells_keep]
  
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
  
  return(seu_obj)
}

# Process and save ONE sample at a time to avoid memory issues------
process_single_sample <- function(sample_idx, metadata) {
  
  sample_id <- metadata$sample_id[sample_idx]
  cat("\n=== Processing sample", sample_idx, "of", nrow(metadata), ":", sample_id, "===\n")
  
  # Process gene-level
  cat("Creating gene-level Seurat object...\n")
  gene_seu <- create_seurat_from_counts(
    metadata$gene_level_count_matrices[[sample_idx]], 
    sample_id,
    assay_type = "RNA",
    MIN.CELL = 3,
    MIN.FEATURES = 20
  )
  print(gene_seu)
  
  # Save immediately
  gene_file <- paste0("./output_files/seu_objects/", sample_id, "_gene_seurat.rds")
  saveRDS(gene_seu, file = gene_file)
  cat("Saved:", gene_file, "\n")
  
  # Clear from memory
  rm(gene_seu)
  gc()
  
  # Process isoform-level
  cat("Creating isoform-level Seurat object...\n")
  iso_seu <- create_seurat_from_counts(
    metadata$iso_level_count_matrices[[sample_idx]], 
    sample_id,
    assay_type = "isoform",
    MIN.CELL = 1,
    MIN.FEATURES = 1
  )
  print(iso_seu)
  
  # Save immediately
  iso_file <- paste0("./output_files/seu_objects/", sample_id, "_isoform_seurat.rds")
  saveRDS(iso_seu, file = iso_file)
  cat("Saved:", iso_file, "\n")
  
  # Clear from memory
  rm(iso_seu)
  gc()
  
  return(sample_id)
}

# Process all samples sequentially------
cat("\n========================================\n")
cat("Starting sequential processing of all samples...\n")
cat("========================================\n")

processed_samples <- c()
for(i in 1:nrow(sample_metadata)) {
  processed_sample <- process_single_sample(i, sample_metadata)
  processed_samples <- c(processed_samples, processed_sample)
}

cat("\n========================================\n")
cat("=== All samples processed successfully! ===\n")
cat("========================================\n")
cat("Processed samples:", paste(processed_samples, collapse = ", "), "\n")
cat("\nSeurat objects saved in: ./output_files/seu_objects/\n")
cat("\nTo load a specific object:\n")
cat("  gene_seu <- readRDS('./output_files/seu_objects/org_1A_gene_seurat.rds')\n")
cat("  iso_seu <- readRDS('./output_files/seu_objects/org_1A_isoform_seurat.rds')\n")