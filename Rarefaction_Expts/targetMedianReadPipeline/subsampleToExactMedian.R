# Script used to determine then use the sampling rate required to reach target median reads per cell using bbmap
# To use: Rscript subsampleToExactMedian.R <fastq_path> <target_median> <output_dir_path>

suppressPackageStartupMessages({
  library(tidyverse)
  library(ShortRead)
  library(fs)
})

### Initialise User-defined Variables from command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if correct number of arguments are provided
if (length(args) != 3) {
  stop("You must provide three arguments: fastq_path, target_median, outputDirPath")
}
fastq_file <- args[1]
target_Median <- as.numeric(args[2])
outputDir <- args[3]

### Initialise Global Default Variables
finegrain_increment <- 0.01
coarseMedianReadsFile <- paste0(outputDir, "/coarse_MedianBarcodeReads_log.csv")
fineMedianReadsFile <- paste0(outputDir, "/fine_MedianBarcodeReads_log.csv")
bash_subsampling_script <- "./2_fineSubsamplingScript.sh"

### Define Functions------
# Function to find sampling rates from coarse median log file to iterate between for fine-grained subsample search
calculateFineSampleRateBounds <- function(coarseMedianLog = coarseMedianReadsFile,
                                          desired_Median = target_Median){
  coarseData <- read.csv(coarseMedianLog)
  
  SampleRateBounds <- coarseData %>%
    mutate(diff = abs(MedianReadsPerCell - desired_Median)) %>% 
    arrange(diff) %>%
    dplyr::slice(1:2) %>% 
    arrange(SampleRate)
  
  sampleRate_LowerBound <- SampleRateBounds$SampleRate[1]
  sampleRate_UpperBound <- SampleRateBounds$SampleRate[2]
  
  return(c(sampleRate_LowerBound, sampleRate_UpperBound))
}

# Function to find the best sampling rate with median reads per barcode closest to target reads
determineFinalSamplingRate <- function(fineMedianLog = fineMedianReadsFile,
                                       desired_Median = target_Median){
  finalSampleRate <- read.csv(fineMedianLog) %>% 
    mutate(diff = abs(MedianReadsPerCell - desired_Median)) %>% 
    arrange(diff) %>% 
    dplyr::slice(1)
  
  return(finalSampleRate)
}

# Function that extracts occurences of each barcode in a fastq file
extractBarcodes <- function(fastqFile, 
                            sampleID = finalSampleRate) {
  fastq_data <- readFastq(fastqFile)
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
# Function that returns the file path of the most recently modified fastq file in output directory
returnMostRecentFQpath <- function(ouput_dir = outputDir){
  # List .fastq files in the output directory
  fastq_files <- dir_ls(ouput_dir, regexp = "\\.fastq$")
  # Get file modification times
  file_mtimes <- file.info(fastq_files)$mtime
  # Select the most recent file
  most_recent_FQ_file <- fastq_files[which.max(file_mtimes)]

  return(most_recent_FQ_file)
}

### Perform fine-grained optimal sampling rate search--------
sampleRateBounds <- calculateFineSampleRateBounds()
finegrain_sequence <- seq(from = sampleRateBounds[[1]], 
                          to = sampleRateBounds[[2]], 
                          by = finegrain_increment)

# call bash subsampling script to generate fine-grained log file displaying median reads per barcode
args <- c(fastq_file, "fine", "false", outputDir, as.character(finegrain_sequence))
system2("bash", args = c(bash_subsampling_script, args))


### Determine optimal sample rate for target median, make final fastq file & summary------
SamplingRateForTargetMedian <- determineFinalSamplingRate()
finalSampleRate <- SamplingRateForTargetMedian$SampleRate
finalSampleRateMedian <- SamplingRateForTargetMedian$MedianReadsPerCell
print(paste("Final Sample Rate:", finalSampleRate, "Final Median:", finalSampleRateMedian))

# Create final subsampled fastq with target median
args <- c(fastq_file, "final", "true", outputDir, finalSampleRate)
system2("bash", args = c(bash_subsampling_script, args))

# store the file path of the final rarified fastq file
finalFQ_file <- returnMostRecentFQpath()
print(paste0("THE FASTQ PATH BEING PASSED IN IS: ", finalFQ_file))

# Generate box plot
readsPerBarcode_df <- extractBarcodes(finalFQ_file, finalSampleRate)

## Make a boxplot and save it to show distribution of final subsampled fastq file
boxPlt <- ggplot(readsPerBarcode_df, aes(x = SampleRate, y = Count, fill = "red")) +
  geom_boxplot(alpha = 1) + 
  labs(title = "Distributions of Reads Per Cell",
       subtitle = paste("Target Median:", target_Median, "\nFinal Median:", finalSampleRateMedian),
       x = "Sample",
       y = "Number of Reads Per Barcode") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none") +
  theme_bw()

ggsave(filename = paste0(outputDir, "/", finalSampleRate, "_", 
                         finalSampleRateMedian, "_boxplot.png"), 
       plot = boxPlt, 
       width = 8, height = 6, dpi = 300)
