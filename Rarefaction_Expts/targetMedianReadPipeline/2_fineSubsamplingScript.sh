#!/bin/bash
# Script to create randomly subsampled fastq files from an original for a list of different sampling rates
# Usage- ./2_fineSubsamplingScript.sh <fast_q file> <true/false to keep rarified fastq files> <log_file_Prefix> <sample_rate1> <sample_rate2> <sample_rate3> ...
# **Make sure bbmap is installed in the current conda environment**

seed=42

# Define global variables from user-entered parameters
fastq_file="$1"
subsample_type="$2"
keepRarifiedFiles="$3"
output_dir="$4"

# Shift the first two arguments so that "$@" contains only the sampling rates
shift 4
# Define array of sampling rates from user input
sample_rates=("$@")

# Loop over the array of sample rates and perform subsampling, keeping track of total number of reads after
for rate in "${sample_rates[@]}"; do
    # Create an output file name based on the sampling rate
    output_file="${output_dir}/${rate}_percentRarified.fastq"

    # Subsample the FASTQ file using reformat.sh with the current sample rate
    reformat.sh in="$fastq_file" out="$output_file" samplerate="$rate" sampleseed="$seed" qin=33

    # Call the R script to process the subsampled file and extract the median reads per barcode
    Rscript findMedianReadsPerCell.R "$output_dir" "$output_file" "$rate" "$subsample_type"

    # Delete the temporary output file after processing if desired
        if [ "$keepRarifiedFiles" = "false" ]; then
        rm "$output_file"
        echo "Deleted $output_file"
    else
        echo "Kept $output_file"
    fi

    echo "Sampling rate $rate median reads per barcode stored"
done