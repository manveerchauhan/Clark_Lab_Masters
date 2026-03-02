# Script to filter bulk count matrices
# Supports both ONT and PacBio data

library(tidyverse)
library(purrr)
library(readr)

# Set working directory to the pipeline location
setwd("/data/gpfs/projects/punim2251/Aim1_LongBench/longbench-analysis-pipeline")

# Define paths
BULK_MATRICES_DIR <- "./bulk_processed_matrices"

# Load expanded sample sheet
sample_sheet <- read.csv("sample_sheet_expanded.csv")

# Function to load and process bulk count matrices
process_bulk_matrices <- function(sample_sheet, technology) {
  # Filter sample sheet for current technology and bulk data
  current_samples <- sample_sheet %>%
    filter(technology == !!technology,
           sample_type == "bulk",
           data_type == "transcript")  # Bulk is transcript-level only
  
  # Initialize lists to store matrices
  bulk_iso_list <- list()
  
  # Process each sample
  for (i in seq_len(nrow(current_samples))) {
    sample_info <- current_samples[i,]
    file_path <- sample_info$file_path
    
    # Read count matrix
    count_matrix <- read.csv(file_path, row.names = 1)
  
    # Store in list
    bulk_iso_list[[file_path]] <- data.frame(
      transcript_id = rownames(count_matrix),
      counts = count_matrix$counts,
      stringsAsFactors = FALSE
    )
  }
  
  return(bulk_iso_list)
}

# Function to filter features based on count thresholds
filter_features <- function(matrix_list, min_total_counts, min_cell_lines) {
  filtered_list <- list()
  
  for (name in names(matrix_list)) {
    df <- matrix_list[[name]]
    
    # Calculate total counts and number of cell lines with counts > 0
    df$total_counts <- df$counts
    df$cell_lines_with_counts <- ifelse(df$counts > 0, 1, 0)
    
    # Filter based on thresholds
    filtered_df <- df %>%
      filter(total_counts >= min_total_counts,
             cell_lines_with_counts >= min_cell_lines) %>%
      select(-total_counts, -cell_lines_with_counts)
    
    filtered_list[[name]] <- filtered_df
  }
  
  return(filtered_list)
}

# Process and filter matrices for each technology
process_all_matrices <- function() {
  # Initialize lists to store results
  all_results <- list()
  
  # Process each technology
  for (tech in c("ONT", "PacBio")) {
    # Process transcript matrices
    tx_results <- process_bulk_matrices(sample_sheet, tech)
    filtered_tx <- filter_features(tx_results,
                                 min_total_counts = 10,  # Transcript threshold
                                 min_cell_lines = 2)
    
    # Store results
    all_results[[paste0(tech, "_tx")]] <- filtered_tx
  }
  
  return(all_results)
}

# Run processing
results <- process_all_matrices()

# Save filtered matrices
for (name in names(results)) {
  saveRDS(results[[name]], 
          file.path(BULK_MATRICES_DIR, 
                   paste0("bulk_isoDfs_", 
                         tolower(strsplit(name, "_")[[1]][1]), 
                         "_filtered.rds")))
}

# Save sample sheet used
write.csv(sample_sheet, 
          file.path(BULK_MATRICES_DIR, "bulk_sample_sheet_used.csv"),
          row.names = FALSE) 