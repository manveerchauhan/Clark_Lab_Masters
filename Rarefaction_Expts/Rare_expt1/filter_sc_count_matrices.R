# Script to filter single-cell and single-nucleus count matrices
# Supports both ONT and PacBio data

library(tidyverse)
library(purrr)
library(readr)

# Set working directory to the pipeline location
setwd("/data/gpfs/projects/punim2251/Aim1_LongBench/longbench-analysis-pipeline")

# Define paths
SC_MATRICES_DIR <- "./sc_processed_matrices"
WHITELIST_DIR <- "/data/gpfs/projects/punim2251/LongBench_data/cell_line_bc_list"

# Load expanded sample sheet
sample_sheet <- read.csv("sample_sheet_expanded.csv")

# Function to load and process count matrices
process_count_matrices <- function(sample_sheet, data_type, sample_type, technology) {
  # Filter sample sheet for current data type and sample type
  current_samples <- sample_sheet %>%
    filter(data_type == !!data_type,
           sample_type == !!sample_type,
           technology == !!technology)
  
  # Load appropriate whitelist
  whitelist_file <- file.path(WHITELIST_DIR, 
                             paste0(sample_type, "/", 
                                   ifelse(technology == "ONT", "ont", "pacbio"),
                                   "_whitelist.txt"))
  whitelist <- readLines(whitelist_file)
  
  # Initialize lists to store matrices
  gene_pseudobulk_list <- list()
  iso_pseudobulk_list <- list()
  
  # Process each sample
  for (i in seq_len(nrow(current_samples))) {
    sample_info <- current_samples[i,]
    file_path <- sample_info$file_path
    
    # Read count matrix
    count_matrix <- read.csv(file_path, row.names = 1)
    
    # Filter for whitelisted barcodes
    count_matrix <- count_matrix[, colnames(count_matrix) %in% whitelist]
    
    # Create pseudobulk by summing across cells
    if (data_type == "gene") {
      gene_pseudobulk <- rowSums(count_matrix)
      gene_pseudobulk_list[[file_path]] <- data.frame(
        gene_id = names(gene_pseudobulk),
        counts = gene_pseudobulk,
        stringsAsFactors = FALSE
      )
    } else {
      iso_pseudobulk <- rowSums(count_matrix)
      iso_pseudobulk_list[[file_path]] <- data.frame(
        transcript_id = names(iso_pseudobulk),
        counts = iso_pseudobulk,
        stringsAsFactors = FALSE
      )
    }
  }
  
  return(list(
    gene_pseudobulk = gene_pseudobulk_list,
    iso_pseudobulk = iso_pseudobulk_list
  ))
}

# Function to filter features based on count thresholds
filter_features <- function(matrix_list, min_counts_per_cell_line, min_cell_lines) {
  filtered_list <- list()
  
  for (name in names(matrix_list)) {
    df <- matrix_list[[name]]
    
    # Count how many cell lines have counts >= min_counts_per_cell_line
    df$cell_lines_above_threshold <- ifelse(df$counts >= min_counts_per_cell_line, 1, 0)
    
    # Filter based on thresholds
    filtered_df <- df %>%
      filter(cell_lines_above_threshold >= min_cell_lines) %>%
      select(-cell_lines_above_threshold)
    
    filtered_list[[name]] <- filtered_df
  }
  
  return(filtered_list)
}

# Process and filter matrices for each combination
process_all_matrices <- function() {
  # Initialize lists to store results
  all_results <- list()
  
  # Process each technology and sample type combination
  for (tech in c("ONT", "PacBio")) {
    for (sample_type in c("sc", "sn")) {
      # Process gene matrices
      gene_results <- process_count_matrices(sample_sheet, "gene", sample_type, tech)
      filtered_gene <- filter_features(gene_results$gene_pseudobulk, 
                                     min_counts_per_cell_line = 5,  # Gene threshold
                                     min_cell_lines = 2)
      
      # Process transcript matrices
      tx_results <- process_count_matrices(sample_sheet, "transcript", sample_type, tech)
      filtered_tx <- filter_features(tx_results$iso_pseudobulk,
                                   min_counts_per_cell_line = 10,  # Transcript threshold
                                   min_cell_lines = 2)
      
      # Store results
      all_results[[paste0(tech, "_", sample_type, "_gene")]] <- filtered_gene
      all_results[[paste0(tech, "_", sample_type, "_tx")]] <- filtered_tx
    }
  }
  
  return(all_results)
}

# Run processing
results <- process_all_matrices()

# Save filtered matrices
for (name in names(results)) {
  saveRDS(results[[name]], 
          file.path(SC_MATRICES_DIR, 
                   paste0("filtered_", 
                         ifelse(grepl("gene$", name), "genePseudobulkDfs_", "isoPseudobulkDfs_"),
                         tolower(name), ".rds")))
}

# Save sample sheet used
write.csv(sample_sheet, 
          file.path(SC_MATRICES_DIR, "sc_sn_sample_sheet_used.csv"),
          row.names = FALSE) 