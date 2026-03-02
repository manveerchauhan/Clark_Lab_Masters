#!/usr/bin/env Rscript

# SC/SN Matrix Filtering Script
# Filters transcripts requiring ≥10 counts in ≥2 cell lines within each protocol+depth combination

library(tidyverse)
library(purrr)
library(readr)

# Set working directory
setwd("/data/gpfs/projects/punim2251/LongBench_rarefaction_Yupei_converted")

# Define paths
PROPORTIONAL_SAMPLE_SHEET <- "sc_sn_sample_sheet_proportional.csv"
ABSOLUTE_SAMPLE_SHEET <- "sc_sn_sample_sheet_absolute.csv"
OUTPUT_DIR <- "filtered_matrices_sc_sn"

# Create output directory
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Function to load a single count matrix
load_count_matrix <- function(file_path) {
  tryCatch({
    # Read CSV with transcript_id and counts columns
    count_data <- read.csv(file_path, stringsAsFactors = FALSE)
    
    # Check expected columns
    if ("transcript_id" %in% colnames(count_data) && "counts" %in% colnames(count_data)) {
      return(count_data)
    } else {
      stop("Expected columns 'transcript_id' and 'counts' not found")
    }
    
  }, error = function(e) {
    cat("Error loading", file_path, ":", as.character(e), "\n")
    return(NULL)
  })
}

# Function to filter transcripts across cell lines within a protocol+depth group
filter_protocol_depth_group <- function(sample_group, min_counts = 10, min_cell_lines = 2) {
  protocol <- unique(sample_group$protocol)
  depth <- unique(sample_group$sampling_depth)
  cat("Processing", protocol, "at depth", depth, "with", nrow(sample_group), "cell lines\n")
  
  # Load all matrices for this protocol+depth combination
  all_matrices <- list()
  
  for (i in seq_len(nrow(sample_group))) {
    sample_info <- sample_group[i,]
    file_path <- sample_info$file_path
    cell_line <- sample_info$cell_line
    
    # Load count matrix
    count_matrix <- load_count_matrix(file_path)
    
    if (!is.null(count_matrix)) {
      # Add cell line identifier
      count_matrix$cell_line <- cell_line
      all_matrices[[cell_line]] <- count_matrix
    }
  }
  
  if (length(all_matrices) == 0) {
    cat("No valid matrices found for", protocol, "at depth", depth, "\n")
    return(NULL)
  }
  
  # Combine all matrices
  combined_data <- do.call(rbind, all_matrices)
  
  # Calculate transcript statistics across cell lines
  transcript_stats <- combined_data %>%
    group_by(transcript_id) %>%
    summarise(
      total_counts = sum(counts, na.rm = TRUE),
      max_counts = max(counts, na.rm = TRUE),
      mean_counts = mean(counts, na.rm = TRUE),
      cell_lines_with_min_counts = sum(counts >= min_counts),
      cell_lines_present = n(),
      .groups = 'drop'
    )
  
  # Apply filtering criteria: ≥min_counts in ≥min_cell_lines
  filtered_transcripts <- transcript_stats %>%
    filter(cell_lines_with_min_counts >= min_cell_lines)
  
  # Create per-cell line filtered matrices
  filtered_matrices_per_cell_line <- list()
  
  for (cell_line in names(all_matrices)) {
    cell_line_data <- all_matrices[[cell_line]]
    
    # Keep only filtered transcripts
    filtered_data <- cell_line_data %>%
      filter(transcript_id %in% filtered_transcripts$transcript_id) %>%
      select(transcript_id, counts)
    
    filtered_matrices_per_cell_line[[cell_line]] <- filtered_data
  }
  
  # Create combined filtered transcript list
  combined_filtered <- filtered_transcripts %>%
    select(transcript_id, total_counts, max_counts, mean_counts, cell_lines_with_min_counts)
  
  # Return results
  return(list(
    combined_filtered = combined_filtered,
    per_cell_line_matrices = filtered_matrices_per_cell_line,
    stats = transcript_stats,
    summary = data.frame(
      protocol = protocol,
      sampling_depth = depth,
      total_transcripts = nrow(transcript_stats),
      filtered_transcripts = nrow(filtered_transcripts),
      filtering_rate = nrow(filtered_transcripts) / nrow(transcript_stats),
      cell_lines_processed = length(all_matrices),
      stringsAsFactors = FALSE
    )
  ))
}

# Function to process a complete sample sheet
process_sample_sheet <- function(sample_sheet_path, output_suffix = "") {
  cat("=== Processing Sample Sheet:", sample_sheet_path, "===\n")
  
  # Load sample sheet
  sample_sheet <- read.csv(sample_sheet_path, stringsAsFactors = FALSE)
  cat("Loaded", nrow(sample_sheet), "samples\n")
  
  # Group by protocol + sampling_depth
  sample_sheet$protocol_depth <- paste(sample_sheet$protocol, sample_sheet$sampling_depth, sep = "_")
  protocol_depth_groups <- split(sample_sheet, sample_sheet$protocol_depth)
  cat("Found", length(protocol_depth_groups), "unique protocol+depth combinations\n\n")
  
  # Process each protocol+depth group
  all_results <- list()
  summary_stats <- list()
  
  for (protocol_depth in names(protocol_depth_groups)) {
    group_data <- protocol_depth_groups[[protocol_depth]]
    
    # Filter this group
    result <- filter_protocol_depth_group(group_data, min_counts = 10, min_cell_lines = 2)
    
    if (!is.null(result)) {
      all_results[[protocol_depth]] <- result
      summary_stats[[protocol_depth]] <- result$summary
      
      # Save combined filtered transcript list
      output_file <- file.path(OUTPUT_DIR, paste0(protocol_depth, output_suffix, "_filtered.csv"))
      write.csv(result$combined_filtered, output_file, row.names = FALSE)
      
      # Save per-cell line filtered matrices in subdirectory
      output_subdir <- file.path(OUTPUT_DIR, paste0(protocol_depth, output_suffix, "_by_cell_line"))
      if (!dir.exists(output_subdir)) {
        dir.create(output_subdir, recursive = TRUE)
      }
      
      for (cell_line in names(result$per_cell_line_matrices)) {
        cell_line_file <- file.path(output_subdir, paste0(cell_line, "_filtered.csv"))
        write.csv(result$per_cell_line_matrices[[cell_line]], cell_line_file, row.names = FALSE)
      }
      
      # Save detailed statistics
      stats_file <- file.path(OUTPUT_DIR, paste0(protocol_depth, output_suffix, "_stats.csv"))
      write.csv(result$stats, stats_file, row.names = FALSE)
      
      summary_file <- file.path(OUTPUT_DIR, paste0(protocol_depth, output_suffix, "_summary.csv"))
      write.csv(result$summary, summary_file, row.names = FALSE)
    }
  }
  
  # Combine and save summary statistics
  if (length(summary_stats) > 0) {
    summary_df <- do.call(rbind, summary_stats)
    overall_summary_file <- file.path(OUTPUT_DIR, paste0("filtering_summary", output_suffix, ".csv"))
    write.csv(summary_df, overall_summary_file, row.names = FALSE)
    
    cat("\n=== Summary Statistics", output_suffix, "===\n")
    cat("Total protocol+depth combinations processed:", nrow(summary_df), "\n")
    cat("Average transcripts per combination:", round(mean(summary_df$total_transcripts)), "\n")
    cat("Average filtered transcripts per combination:", round(mean(summary_df$filtered_transcripts)), "\n")
    cat("Average filtering rate:", round(mean(summary_df$filtering_rate) * 100, 1), "%\n")
    cat("Summary saved to:", overall_summary_file, "\n\n")
  }
  
  return(all_results)
}

# Main execution
main_filtering <- function() {
  cat("=== SC/SN Matrix Filtering ===\n")
  cat("Output directory:", OUTPUT_DIR, "\n")
  cat("Filtering criteria: ≥10 counts in ≥2 cell lines within each protocol+depth\n\n")
  
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
    filtering_criteria = "≥10 counts in ≥2 cell lines",
    stringsAsFactors = FALSE
  )
  
  write.csv(overall_summary, file.path(OUTPUT_DIR, "overall_summary.csv"), row.names = FALSE)
  cat("Overall summary saved to:", file.path(OUTPUT_DIR, "overall_summary.csv"), "\n")
}

# Run the filtering
cat("Starting SC/SN cross-cell line filtering...\n")
main_filtering()
cat("Done!\n") 