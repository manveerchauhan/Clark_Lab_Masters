# Script that calculates the number of reads to pull from each bulk fastq file
# based on the proportion of reads belonging to each cell line in the single cell/nuclei data
# RUN THIS AFTER USING returnBarcodeReadInfo.py script***
#
# Usage:
#   Rscript createMatchedBulkFQ_helper_script.R <read_info_file> <seq_modality> <sample_id> <total_number_of_reads> <output_dir>
#
# Arguments:
#   read_info_file: path to CSV file containing read information
#   seq_modality: sequencing modality (e.g., "ont_sc")
#   sample_id: identifier for the sc/sn fastq being processed
#   total_number_of_reads: total number of reads (integer)
#   output_dir: directory for output files
#
# Example:
#   Rscript createMatchedBulkFQ_helper_script.R \
#     "/data/gpfs/projects/punim2251/LongBench_data/data_diagnostics/ont_sc_read_info.csv" \
#     "ont_sc" \
#     "test" \
#     91339154 \
#     "/data/gpfs/projects/punim2251/LongBench_data"

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

# Parse the arguments
read_info_file <- read.csv(args[1])
seq_modality <- args[2]
sample_id <- args[3]
total_number_of_reads <- as.numeric(args[4])
output_dir <- args[5]
setwd(output_dir)

# Check if the correct number of arguments is provided
if (length(args) != 5) {
  stop("Five arguments must be supplied:
       1. read_info_file path
       2. seq_modality
       3. sample_id
       4. total_number_of_reads
       5. output_dir")
}

# Parse with error checking
read_info_file <- read.csv(args[1])
if (!file.exists(args[1])) stop("Read info file does not exist")

seq_modality <- args[2]
sample_id <- args[3]

total_number_of_reads <- as.numeric(args[4])
if (is.na(total_number_of_reads)) stop("total_number_of_reads must be numeric")

output_dir <- args[5]
if (!dir.exists(output_dir)) stop("Output directory does not exist")
setwd(output_dir)

# Define whitelist containing barcodes and cell type labels----
ont_sc_whitelist <- read.csv("/data/gpfs/projects/punim2251/LongBench_data/sample_whitelists/ont_sc_whitelist_wLabels.csv")
pb_sc_whitelist <- read.csv("/data/gpfs/projects/punim2251/LongBench_data/sample_whitelists/pb_sc_whitelist_wLabels.csv")
pb_sn_whitelist <- read.csv("/data/gpfs/projects/punim2251/LongBench_data/sample_whitelists/pb_sn_whitelist_wLabels.csv")
ont_sn_whitelist <- read.csv("/data/gpfs/projects/punim2251/LongBench_data/sample_whitelists/ont_sn_whitelist_wLabels.csv")

if(seq_modality == "ont_sc"){
  barcode_cellType_reference <- ont_sc_whitelist
}else if(seq_modality == "pb_sc"){
  barcode_cellType_reference <- pb_sc_whitelist
}else if(seq_modality == "ont_sn"){
  barcode_cellType_reference <- ont_sn_whitelist
}else if(seq_modality == "pb_sn"){
  barcode_cellType_reference <- pb_sc_whitelist
} else{
  message("Invalid seq modality type... use 'ont_sc', 'pb_sc', 'ont_sn', or 'pb_sn'")
}

## Define functions being used-----
# Function to calculate counts and proportions of reads belonging to each cell type category
calculate_cell_category_counts <- function(df) {
  # Create a dataframe showing the relative proportion of reads to each annotation group
  count_prop_df <- df %>%
    # Group by category and sum counts
    group_by(merged_cellline_anno) %>%
    summarise(total_count = sum(count)) %>%
    # Sort by total count
    arrange(desc(total_count)) %>%
    # Calculate overall percentages
    mutate(
      overall_percent = round((total_count / sum(total_count) * 100), 2)
    ) %>%
    # Calculate percentages excluding doublets and unassigned
    mutate(
      clean_percent = case_when(
        merged_cellline_anno %in% c("doublet", "unassigned") ~ NA_real_,
        TRUE ~ (total_count / sum(total_count[!merged_cellline_anno %in% c("doublet", "unassigned")]))
      )
    )
  
  return(count_prop_df)
}

read_info_with_anno <- merge(
  read_info_file,
  barcode_cellType_reference[, c("barcode", "merged_cellline_anno")],
  by = "barcode",
  all.x = TRUE
)

category_proportions <- calculate_cell_category_counts(read_info_with_anno) %>% 
  na.omit() %>% 
  dplyr::select(merged_cellline_anno, clean_percent) %>% 
  dplyr::mutate(bulk_reads_required = round(clean_percent*total_number_of_reads))

if(sum(category_proportions$clean_percent) != 1){
  message("WARNING: Clean cell type proportions don't equal 1")
}

if(sum(category_proportions$bulk_reads_required) != total_number_of_reads){
  message(paste0("WARNING: Number of calculated bulk reads (", sum(category_proportions$bulk_reads_required),
                                                                   ") don't align with total # reads provided (",
                                                                   (total_number_of_reads)), ")")
  
  initial_difference = sum(category_proportions$bulk_reads_required) - total_number_of_reads
  
  absolute_difference = abs(initial_difference)
  
  message(paste("There is a", absolute_difference, "read discrepency"))
  
  message("To address the discrepency....")
  while(absolute_difference != 0){
    random_index <- sample(nrow(category_proportions), 1)
    
    if(initial_difference < 0){
      message(paste("Adding 1 read to:", category_proportions[[1]][random_index]))
      category_proportions$bulk_reads_required[random_index] <- category_proportions$bulk_reads_required[random_index] + 1
    }else{
      message(paste("Subtracting 1 read from:", category_proportions[[1]][random_index]))
      category_proportions$bulk_reads_required[random_index] <- category_proportions$bulk_reads_required[random_index] - 1
    }
    
    absolute_difference = absolute_difference - 1
  }
  
  message(paste0("Number of calculated bulk reads (", sum(category_proportions$bulk_reads_required),
                 ") match with total # reads provided now (",
                 (total_number_of_reads)), ")")
}

output_df <- category_proportions %>% 
  dplyr::select(-clean_percent)

write_csv(
  output_df,
  file.path(output_dir, paste0(seq_modality, "_", sample_id, "_bulkReadsNeeded.csv"))
)
