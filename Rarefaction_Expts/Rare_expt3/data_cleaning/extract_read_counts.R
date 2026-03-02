#!/usr/bin/env Rscript
# Extract valid_best_aln counts from meta_info.json files across all directories
# and create a mapping dataset for merging with sample sheets.

library(jsonlite)
library(dplyr)

extract_read_counts <- function() {
  base_path <- "/data/gpfs/projects/punim2251/LongBench_rarefaction_Yupei"
  results <- data.frame()
  
  # Technology mappings for sample sheet consistency
  tech_mapping <- list(
    "pb_bulk" = "PacBio",
    "ont_bulk" = "ONT", 
    "dRNA_bulk" = "direct_RNA"
  )
  
  cat("Starting extraction...\n")
  
  # Process bulk data directories
  bulk_dirs <- c("pb_bulk", "ont_bulk", "dRNA_bulk")
  
  for (tech_dir in bulk_dirs) {
    tech_path <- file.path(base_path, tech_dir)
    
    if (!dir.exists(tech_path)) {
      cat("Warning:", tech_path, "does not exist\n")
      next
    }
    
    cat("Processing", tech_dir, "...\n")
    technology <- tech_mapping[[tech_dir]]
    if (is.null(technology)) technology <- tech_dir
    
    # Get all sampling depth directories
    depth_dirs <- list.dirs(tech_path, recursive = FALSE, full.names = FALSE)
    
    for (sampling_depth in depth_dirs) {
      depth_path <- file.path(tech_path, sampling_depth)
      cat("  Processing depth:", sampling_depth, "\n")
      
      # Get all cell line directories
      cell_line_dirs <- list.dirs(depth_path, recursive = FALSE, full.names = FALSE)
      
      for (cell_line in cell_line_dirs) {
        cell_line_path <- file.path(depth_path, cell_line)
        meta_file <- file.path(cell_line_path, "meta_info.json")
        
        if (file.exists(meta_file)) {
          tryCatch({
            meta_data <- fromJSON(meta_file)
            
            if (!is.null(meta_data$discard_table$valid_best_aln)) {
              valid_best_aln <- meta_data$discard_table$valid_best_aln
              
              # Create relative path
              rel_path <- file.path("LongBench_rarefaction_Yupei", tech_dir, sampling_depth, cell_line)
              rel_meta_path <- file.path(rel_path, "meta_info.json")
              
              new_row <- data.frame(
                technology = technology,
                data_type = "bulk",
                cell_line = cell_line,
                sampling_depth = sampling_depth,
                measurement_type = "absolute",
                valid_best_aln = valid_best_aln,
                directory_path = rel_path,
                meta_file_path = rel_meta_path,
                stringsAsFactors = FALSE
              )
              
              results <- rbind(results, new_row)
              cat("    ", cell_line, ":", format(valid_best_aln, big.mark = ","), "reads\n")
            } else {
              cat("    ", cell_line, ": No valid_best_aln found\n")
            }
          }, error = function(e) {
            cat("    ", cell_line, ": Error reading meta_info.json -", e$message, "\n")
          })
        } else {
          cat("    ", cell_line, ": meta_info.json not found\n")
        }
      }
    }
  }
  
  # Process sc_sn data directories
  sc_sn_path <- file.path(base_path, "sc_sn")
  
  if (dir.exists(sc_sn_path)) {
    cat("Processing sc_sn...\n")
    
    sc_sn_tech_mapping <- list(
      "pb_sc" = c("PacBio", "single_cell"),
      "pb_sn" = c("PacBio", "single_nucleus"),
      "ont_sc" = c("ONT", "single_cell"),
      "ont_sn" = c("ONT", "single_nucleus")
    )
    
    tech_data_dirs <- list.dirs(sc_sn_path, recursive = FALSE, full.names = FALSE)
    
    for (tech_data_name in tech_data_dirs) {
      if (!tech_data_name %in% names(sc_sn_tech_mapping)) next
      
      tech_data_path <- file.path(sc_sn_path, tech_data_name)
      technology <- sc_sn_tech_mapping[[tech_data_name]][1]
      data_type <- sc_sn_tech_mapping[[tech_data_name]][2]
      
      cat("  Processing", tech_data_name, "(", technology, data_type, ")...\n")
      
      # Get all sampling depth directories
      depth_dirs <- list.dirs(tech_data_path, recursive = FALSE, full.names = FALSE)
      
      for (sampling_depth in depth_dirs) {
        depth_path <- file.path(tech_data_path, sampling_depth)
        cat("    Processing depth:", sampling_depth, "\n")
        
        # Get all cell line directories
        cell_line_dirs <- list.dirs(depth_path, recursive = FALSE, full.names = FALSE)
        
        for (cell_line in cell_line_dirs) {
          cell_line_path <- file.path(depth_path, cell_line)
          meta_file <- file.path(cell_line_path, "meta_info.json")
          
          if (file.exists(meta_file)) {
            tryCatch({
              meta_data <- fromJSON(meta_file)
              
              if (!is.null(meta_data$discard_table$valid_best_aln)) {
                valid_best_aln <- meta_data$discard_table$valid_best_aln
                
                # Create relative path
                rel_path <- file.path("LongBench_rarefaction_Yupei", "sc_sn", tech_data_name, sampling_depth, cell_line)
                rel_meta_path <- file.path(rel_path, "meta_info.json")
                
                new_row <- data.frame(
                  technology = technology,
                  data_type = data_type,
                  cell_line = cell_line,
                  sampling_depth = sampling_depth,
                  measurement_type = "absolute",
                  valid_best_aln = valid_best_aln,
                  directory_path = rel_path,
                  meta_file_path = rel_meta_path,
                  stringsAsFactors = FALSE
                )
                
                results <- rbind(results, new_row)
                cat("      ", cell_line, ":", format(valid_best_aln, big.mark = ","), "reads\n")
              } else {
                cat("      ", cell_line, ": No valid_best_aln found\n")
              }
            }, error = function(e) {
              cat("      ", cell_line, ": Error reading meta_info.json -", e$message, "\n")
            })
          } else {
            cat("      ", cell_line, ": meta_info.json not found\n")
          }
        }
      }
    }
  }
  
  # Save results and create summary
  if (nrow(results) > 0) {
    # Sort for easier inspection
    results <- results %>%
      arrange(technology, data_type, cell_line, sampling_depth)
    
    # Save to CSV
    output_file <- "read_counts_metadata.csv"
    write.csv(results, output_file, row.names = FALSE)
    
    cat("\nExtraction complete!\n")
    cat("Total entries:", nrow(results), "\n")
    cat("Results saved to:", output_file, "\n")
    
    # Show summary statistics
    cat("\nSummary by technology and data type:\n")
    summary_stats <- results %>%
      group_by(technology, data_type) %>%
      summarise(
        count = n(),
        min_reads = min(valid_best_aln),
        max_reads = max(valid_best_aln),
        mean_reads = round(mean(valid_best_aln)),
        unique_cell_lines = n_distinct(cell_line),
        unique_depths = n_distinct(sampling_depth),
        .groups = 'drop'
      )
    
    print(summary_stats)
    
    return(results)
  } else {
    cat("No valid data found!\n")
    return(NULL)
  }
}

# Run the extraction
df <- extract_read_counts() 