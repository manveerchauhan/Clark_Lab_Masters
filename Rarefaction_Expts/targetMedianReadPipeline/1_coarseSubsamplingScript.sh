#!/bin/bash
# Script to create randomly subsampled fastq files from an original for a list of different sampling rates
# Usage- change sample_rates array to whatever you want to subsample to then run the following command in terminal:
# ./1_coarseSubsamplingScript.sh <fast_q file> <output_dir (eg. ./output)>
# **Make sure bbmap is installed in the current conda environment**

seed=42

# Define global variables from user-entered parameters
fastq_file="$1"
output_dir="$2"

# Define array of sampling rates
sample_rates=(1.0 0.90 0.80 0.70 0.60 0.50 0.40 0.30 0.20 0.10 0.01)

# Loop over the array of sample rates and perform subsampling, keeping track of total number of reads after
for rate in "${sample_rates[@]}"; do
    # Create an output file name based on the sampling rate
    output_file="${output_dir}/temp.fastq"

    # Subsample the FASTQ file using reformat.sh with the current sample rate
    reformat.sh in="$fastq_file" out="$output_file" samplerate="$rate" sampleseed="$seed" qin=33

    # Call the R script to process the subsampled file and extract the median reads per barcode
    Rscript findMedianReadsPerCell.R "$output_dir" "$output_file" "$rate" "coarse"

    # Delete the temporary output file after processing
    rm "$output_file"

    echo "Sampling rate $rate median reads per barcode stored"
done