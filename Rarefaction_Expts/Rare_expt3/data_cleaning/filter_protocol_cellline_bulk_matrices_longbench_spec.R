#!/usr/bin/env Rscript

# Protocol-CellLine Bulk Matrix Filtering Script - LongBench Specific
# Filters transcripts within each individual count matrix file for ≥5 counts

library(tidyverse)
library(purrr)
library(readr)

# Set working directory
setwd("/data/gpfs/projects/punim2251/LongBench_rarefaction_Yupei_converted")

# Define paths
PROPORTIONAL_SAMPLE_SHEET <- "bulk_sample_sheet_proportional.csv"
ABSOLUTE_SAMPLE_SHEET <- "bulk_sample_sheet_absolute.csv"
OUTPUT_DIR <- "filtered_matrices_longbench_spec"

# Create output directory
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Function to load and filter a single count matrix
load_and_filter_matrix <- function(file_path, min_counts = 5) {
  tryCatch({
    # Read CSV with transcript_id and counts columns
    count_data <- read.csv(file_path, stringsAsFactors = FALSE)
    
    # Check expected columns
    if ("transcript_id" %in% colnames(count_data) && "counts" %in% colnames(count_data)) {
      
      # Filter transcripts with ≥min_counts
      filtered_data <- count_data %>%
        filter(counts >= min_counts)
      
      # Convert to row names format
      filtered_matrix <- data.frame(
        counts = filtered_data$counts,
        row.names = filtered_data$transcript_id,
        stringsAsFactors = FALSE
      )
      
      # Return both filtered matrix and stats
      return(list(
        matrix = filtered_matrix,
        total_transcripts = nrow(count_data),
        filtered_transcripts = nrow(filtered_data),
        filtering_rate = nrow(filtered_data) / nrow(count_data)
      ))
      
    } else {
      stop("Expected columns 'transcript_id' and 'counts' not found")
    }
    
  }, error = function(e) {
    cat("Error loading", file_path, ":", as.character(e), "\n")
    return(NULL)
  })
}

# Function to process a protocol_cellline group (now much simpler)
process_protocol_cellline_group <- function(sample_group, min_counts = 5) {
  protocol_cellline <- unique(sample_group$protocol_cellline)
  cat("Processing", protocol_cellline, "with", nrow(sample_group), "depths\n")
  
  # Process each file independently
  filtered_matrices <- list()
  file_stats <- list()
  
  for (i in seq_len(nrow(sample_group))) {
    sample_info <- sample_group[i,]
    file_path <- sample_info$file_path
    depth <- as.character(sample_info$sampling_depth)  # Ensure it's character for list indexing
    
    # Filter this individual matrix
    result <- load_and_filter_matrix(file_path, min_counts)
    
    if (!is.null(result)) {
      filtered_matrices[[depth]] <- result$matrix
      file_stats[[depth]] <- data.frame(
        protocol_cellline = protocol_cellline,
        sampling_depth = depth,
        total_transcripts = result$total_transcripts,
        filtered_transcripts = result$filtered_transcripts,
        filtering_rate = result$filtering_rate,
        stringsAsFactors = FALSE
      )
      
      cat("  ", depth, ":", result$filtered_transcripts, "/", result$total_transcripts, 
          "(", round(result$filtering_rate * 100, 1), "%) transcripts retained\n")
    }
  }
  
  if (length(filtered_matrices) == 0) {
    cat("No valid matrices found for", protocol_cellline, "\n")
    return(NULL)
  }
  
  # Combine file statistics
  stats_df <- do.call(rbind, file_stats)
  
  # Create summary for this protocol_cellline
  summary_stats <- data.frame(
    protocol_cellline = protocol_cellline,
    depths_processed = nrow(stats_df),
    avg_total_transcripts = round(mean(stats_df$total_transcripts)),
    avg_filtered_transcripts = round(mean(stats_df$filtered_transcripts)),
    avg_filtering_rate = mean(stats_df$filtering_rate),
    min_counts_threshold = min_counts,
    stringsAsFactors = FALSE
  )
  
  return(list(
    matrices = filtered_matrices,
    file_stats = stats_df,
    summary = summary_stats
  ))
}

# Function to process a complete sample sheet
process_sample_sheet <- function(sample_sheet_path, output_suffix = "") {
  cat("=== Processing Sample Sheet:", sample_sheet_path, "===\n")
  
  # Load sample sheet
  sample_sheet <- read.csv(sample_sheet_path, stringsAsFactors = FALSE)
  cat("Loaded", nrow(sample_sheet), "samples\n")
  
  # Group by protocol_cellline
  protocol_groups <- split(sample_sheet, sample_sheet$protocol_cellline)
  cat("Found", length(protocol_groups), "unique protocol_cellline combinations\n\n")
  
  # Process each protocol_cellline group
  all_results <- list()
  summary_stats <- list()
  all_file_stats <- list()
  
  for (protocol_cellline in names(protocol_groups)) {
    group_data <- protocol_groups[[protocol_cellline]]
    
    # Process this group
    result <- process_protocol_cellline_group(group_data, min_counts = 5)
    
    if (!is.null(result)) {
      all_results[[protocol_cellline]] <- result
      summary_stats[[protocol_cellline]] <- result$summary
      all_file_stats[[protocol_cellline]] <- result$file_stats
      
      # Save filtered matrices for this protocol_cellline
      output_subdir <- file.path(OUTPUT_DIR, paste0(protocol_cellline, output_suffix))
      if (!dir.exists(output_subdir)) {
        dir.create(output_subdir, recursive = TRUE)
      }
      
      # Save each depth's filtered matrix
      for (depth in names(result$matrices)) {
        output_file <- file.path(output_subdir, paste0("filtered_", depth, ".csv"))
        
        # Convert matrix back to CSV format with transcript_id column
        matrix_data <- result$matrices[[depth]]
        csv_data <- data.frame(
          transcript_id = rownames(matrix_data),
          counts = matrix_data$counts,
          stringsAsFactors = FALSE
        )
        
        write.csv(csv_data, output_file, row.names = FALSE)
      }
      
      # Save detailed statistics as CSV files
      file_stats_file <- file.path(output_subdir, "file_stats.csv")
      summary_stats_file <- file.path(output_subdir, "summary_stats.csv")
      
      write.csv(result$file_stats, file_stats_file, row.names = FALSE)
      write.csv(result$summary, summary_stats_file, row.names = FALSE)
    }
  }
  
  # Combine and save summary statistics
  if (length(summary_stats) > 0) {
    summary_df <- do.call(rbind, summary_stats)
    summary_file <- file.path(OUTPUT_DIR, paste0("filtering_summary", output_suffix, ".csv"))
    write.csv(summary_df, summary_file, row.names = FALSE)
    
    # Combine detailed file statistics
    if (length(all_file_stats) > 0) {
      detailed_stats_df <- do.call(rbind, all_file_stats)
      detailed_file <- file.path(OUTPUT_DIR, paste0("detailed_filtering_stats", output_suffix, ".csv"))
      write.csv(detailed_stats_df, detailed_file, row.names = FALSE)
    }
    
    cat("\n=== Summary Statistics", output_suffix, "===\n")
    cat("Total protocol_cellline combinations processed:", nrow(summary_df), "\n")
    cat("Average transcripts per file:", round(mean(summary_df$avg_total_transcripts)), "\n")
    cat("Average filtered transcripts per file:", round(mean(summary_df$avg_filtered_transcripts)), "\n")
    cat("Average filtering rate:", round(mean(summary_df$avg_filtering_rate) * 100, 1), "%\n")
    cat("Summary saved to:", summary_file, "\n\n")
  }
  
  return(all_results)
}

# Main execution
main_filtering <- function() {
  cat("=== Protocol-CellLine Bulk Matrix Filtering - LongBench Specific ===\n")
  cat("Output directory:", OUTPUT_DIR, "\n")
  cat("Filtering criteria: ≥5 counts within each individual file\n\n")
  
  start_time <- Sys.time()
  
  # Process proportional sample sheet
  if (file.exists(PROPORTIONAL_SAMPLE_SHEET)) {
    proportional_results <- process_sample_sheet(PROPORTIONAL_SAMPLE_SHEET, "_proportional")
  }
  
  # Process absolute sample sheet  
  if (file.exists(ABSOLUTE_SAMPLE_SHEET)) {
    absolute_results <- process_sample_sheet(ABSOLUTE_SAMPLE_SHEET, "_absolute")
  }
  
  end_time <- Sys.time()
  total_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
  
  cat("=== Filtering Complete ===\n")
  cat("Total processing time:", round(total_time, 2), "minutes\n")
  cat("Results saved to:", OUTPUT_DIR, "\n")
  
  # Create overall summary
  overall_summary <- data.frame(
    processing_date = as.character(Sys.time()),
    proportional_sheet_processed = file.exists(PROPORTIONAL_SAMPLE_SHEET),
    absolute_sheet_processed = file.exists(ABSOLUTE_SAMPLE_SHEET),
    total_time_minutes = total_time,
    output_directory = OUTPUT_DIR,
    filtering_criteria = "≥5 counts per file",
    stringsAsFactors = FALSE
  )
  
  write.csv(overall_summary, file.path(OUTPUT_DIR, "overall_summary.csv"), row.names = FALSE)
  cat("Overall summary saved to:", file.path(OUTPUT_DIR, "overall_summary.csv"), "\n")
}

# Run the filtering
cat("Starting protocol_cellline individual file filtering...\n")
main_filtering()
cat("Done!\n") 