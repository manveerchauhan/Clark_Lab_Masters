# Script that is used to extract median reads per barcode for a given sc fastq file
# To use: Rscript findMedianReadsPerCell.R <output_directory> <fastq_Path> <sampleRate>

suppressPackageStartupMessages({
  library(tidyverse)
  library(ShortRead)
})

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if correct number of arguments are provided
if (length(args) != 4) {
  stop("You must provide three arguments: output_directory, fastq_path, sample_rate, log_filePrefix")
}

# Assign arguments to variables
outputDirectory <- args[1]
fastq_path <- args[2]
sample_rate <- args[3]
log_file_prefix <- args[4]

# Create the output directory if it doesn't exist
if (!dir.exists(outputDirectory)) {
  dir.create(outputDirectory, recursive = TRUE)
}

extractBarcodes <- function(fastq_file, sampleID) {
  fastq_data <- readFastq(fastq_file)
  ids <- id(fastq_data)
  
  # Convert the IDs to a character vector
  ids_vector <- as.character(ids)
  # Create a data frame from the character vector
  ids_df <- ids_vector %>%
    data.frame(Read_IDs = .) %>%
    mutate(Read_IDs = substr(Read_IDs, 1, 16)) %>% 
    group_by(Read_IDs) %>%
    summarise(Count = n()) %>%
    arrange(desc(Count)) %>% 
    mutate(SampleRate = sampleID)
  
  return(ids_df)
}
exportMedianValue <- function(sampleRate = sample_rate, 
                              medianReads = medianReadsPerCell){
  # Create a data frame with sample_rate and medianReadsPerCell
  output_data <- data.frame(SampleRate = sampleRate, 
                            MedianReadsPerCell = medianReads)
  
  # Specify the output file path
  output_file <- file.path(outputDirectory, 
                           paste0(log_file_prefix, "_MedianBarcodeReads_log.csv"))
  
  # Write the data to the CSV file. If the file already exists, append the new data.
  if (file.exists(output_file)) {
    write.table(output_data, file = output_file, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  } else {
    write.table(output_data, file = output_file, sep = ",", row.names = FALSE, col.names = TRUE)
  }
}

readsPerBarcode_df <- extractBarcodes(fastq_path, sample_rate)

medianReadsPerCell <- median(readsPerBarcode_df$Count)

exportMedianValue()
message("Added median reads per barcode into log")
