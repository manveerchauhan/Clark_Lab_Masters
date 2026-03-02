#!/usr/bin/env Rscript

# Script to extract all unique transcript biotype categories from count matrices
# Goes through both bulk and sc/sn sample sheets and examines all files

library(tidyverse)
library(readr)

# Set working directory
setwd("/data/gpfs/projects/punim2251/LongBench_rarefaction_Yupei_converted")

# Function to extract biotypes from a single file
extract_biotypes_from_file <- function(file_path) {
  tryCatch({
    cat("Processing:", basename(file_path), "\n")
    
    # Read the CSV file
    data <- read.csv(file_path, stringsAsFactors = FALSE)
    
    # Check if biotype column exists
    if ("biotype" %in% colnames(data)) {
      biotypes <- unique(data$biotype)
      # Remove any empty or NA values
      biotypes <- biotypes[!is.na(biotypes) & biotypes != ""]
      return(biotypes)
    } else {
      cat("  Warning: No 'biotype' column found in", basename(file_path), "\n")
      cat("  Available columns:", paste(colnames(data), collapse = ", "), "\n")
      return(character(0))
    }
    
  }, error = function(e) {
    cat("Error processing", file_path, ":", as.character(e), "\n")
    return(character(0))
  })
}

# Read both sample sheets
cat("=== READING SAMPLE SHEETS ===\n")

# Check if sample sheets exist
bulk_sheet_file <- "bulk_cleaned_sample_sheet_with_counts.csv"
sc_sn_sheet_file <- "sc_sn_cleaned_sample_sheet_with_counts.csv"

if (!file.exists(bulk_sheet_file)) {
  stop("Bulk sample sheet not found: ", bulk_sheet_file)
}

if (!file.exists(sc_sn_sheet_file)) {
  stop("SC/SN sample sheet not found: ", sc_sn_sheet_file)
}

bulk_sheet <- read.csv(bulk_sheet_file, stringsAsFactors = FALSE)
sc_sn_sheet <- read.csv(sc_sn_sheet_file, stringsAsFactors = FALSE)

cat("Loaded bulk sample sheet:", nrow(bulk_sheet), "files\n")
cat("Loaded SC/SN sample sheet:", nrow(sc_sn_sheet), "files\n")

# Combine all file paths
all_file_paths <- c(bulk_sheet$file_path, sc_sn_sheet$file_path)
cat("Total files to examine:", length(all_file_paths), "\n")

# Initialize containers for biotypes
all_biotypes <- character(0)
biotype_counts <- list()
files_processed <- 0
files_with_biotypes <- 0

cat("\n=== EXTRACTING BIOTYPES ===\n")

# Process files in batches to avoid overwhelming output
batch_size <- 50
total_files <- length(all_file_paths)

for (i in 1:total_files) {
  file_path <- all_file_paths[i]
  
  if (!file.exists(file_path)) {
    cat("File not found:", file_path, "\n")
    next
  }
  
  # Extract biotypes from this file
  file_biotypes <- extract_biotypes_from_file(file_path)
  
  files_processed <- files_processed + 1
  
  if (length(file_biotypes) > 0) {
    files_with_biotypes <- files_with_biotypes + 1
    
    # Add to master list
    all_biotypes <- unique(c(all_biotypes, file_biotypes))
    
    # Count occurrences
    for (biotype in file_biotypes) {
      if (biotype %in% names(biotype_counts)) {
        biotype_counts[[biotype]] <- biotype_counts[[biotype]] + 1
      } else {
        biotype_counts[[biotype]] <- 1
      }
    }
  }
  
  # Progress update every 50 files
  if (i %% batch_size == 0) {
    cat("Progress:", i, "/", total_files, "files processed\n")
    cat("Current unique biotypes found:", length(all_biotypes), "\n")
  }
}

# Convert biotype_counts to data frame for easier handling
biotype_df <- data.frame(
  biotype = names(biotype_counts),
  file_count = unlist(biotype_counts),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(file_count))

cat("\n=== SUMMARY RESULTS ===\n")
cat("Files processed:", files_processed, "/", total_files, "\n")
cat("Files with biotype data:", files_with_biotypes, "\n")
cat("Total unique biotypes found:", length(all_biotypes), "\n")

cat("\n=== UNIQUE BIOTYPE CATEGORIES ===\n")
# Sort biotypes alphabetically for easy reading
sorted_biotypes <- sort(all_biotypes)
for (i in 1:length(sorted_biotypes)) {
  cat(sprintf("%2d. %s\n", i, sorted_biotypes[i]))
}

cat("\n=== BIOTYPE FREQUENCY (Top 20) ===\n")
cat("Biotype\t\t\t\tFiles Present\n")
cat("-------\t\t\t\t-------------\n")
for (i in 1:min(20, nrow(biotype_df))) {
  cat(sprintf("%-30s\t%d\n", biotype_df$biotype[i], biotype_df$file_count[i]))
}

# Save detailed results
output_file <- "transcript_biotypes_summary.csv"
write.csv(biotype_df, output_file, row.names = FALSE)
cat("\nDetailed biotype frequency saved to:", output_file, "\n")

# Create a simple list file
biotype_list_file <- "unique_biotypes_list.txt"
writeLines(sorted_biotypes, biotype_list_file)
cat("Simple biotype list saved to:", biotype_list_file, "\n")

# Additional analysis: categorize biotypes
cat("\n=== BIOTYPE CATEGORIZATION ===\n")

protein_coding <- grep("protein_coding", sorted_biotypes, value = TRUE)
pseudogenes <- grep("pseudogene", sorted_biotypes, value = TRUE)
ncrna <- grep("RNA|rRNA|tRNA|snoRNA|snRNA|miRNA|lncRNA", sorted_biotypes, value = TRUE)
other <- setdiff(sorted_biotypes, c(protein_coding, pseudogenes, ncrna))

cat("Protein coding:", length(protein_coding), "\n")
if (length(protein_coding) > 0) cat("  -", paste(protein_coding, collapse = ", "), "\n")

cat("\nPseudogenes:", length(pseudogenes), "\n") 
if (length(pseudogenes) > 0) cat("  -", paste(pseudogenes, collapse = ", "), "\n")

cat("\nNon-coding RNAs:", length(ncrna), "\n")
if (length(ncrna) > 0) cat("  -", paste(ncrna, collapse = ", "), "\n")

cat("\nOther categories:", length(other), "\n")
if (length(other) > 0) cat("  -", paste(other, collapse = ", "), "\n")

cat("\n=== COMPLETED ===\n")
cat("Found", length(all_biotypes), "unique transcript biotype categories across", files_with_biotypes, "files\n") 