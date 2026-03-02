#!/bin/bash
# Script to subsample a FASTQ file using predefined sampling rates
# Usage: ./subsample_fastq_with_rates.sh <input_fastq> <output_dir> <prefix> <sample_rate1> <sample_rate2> <sample_rate3> ...

# Check if minimum number of arguments provided
if [ $# -lt 4 ]; then
    echo "Usage: $0 <input_fastq> <output_dir> <prefix> <sample_rate1> [sample_rate2] [sample_rate3] ..."
    echo "Example: $0 input.fastq /path/to/output my_prefix 0.05 0.12 0.15 0.25 0.46"
    exit 1
fi

# Use the same seed as the original script for reproducibility
seed=42

# Define global variables from user-entered parameters
input_fastq="$1"
output_dir="$2"
prefix="$3"

# Shift the first three arguments so that "$@" contains only the sampling rates
shift 3
# Define array of sampling rates from user input
sample_rates=("$@")

# Check if input FASTQ file exists
if [ ! -f "$input_fastq" ]; then
    echo "Error: Input FASTQ file '$input_fastq' does not exist!"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Check if output directory is writable
if [ ! -w "$output_dir" ]; then
    echo "Error: Output directory '$output_dir' is not writable!"
    exit 1
fi

echo "Starting subsampling of: $input_fastq"
echo "Output directory: $output_dir"
echo "Prefix: $prefix"
echo "Sample rates: ${sample_rates[@]}"
echo "Using seed: $seed"
echo "----------------------------------------"

# Loop over the array of sample rates and perform subsampling
for rate in "${sample_rates[@]}"; do
    # Create an output file name based on the sampling rate and prefix
    output_file="${output_dir}/${rate}_${prefix}.fastq"
    
    echo "Processing sampling rate: $rate"
    echo "Output file: $output_file"
    
    # Subsample the FASTQ file using reformat.sh with the current sample rate
    reformat.sh in="$input_fastq" out="$output_file" samplerate="$rate" sampleseed="$seed" qin=33
    
    # Check if reformat.sh succeeded
    if [ $? -eq 0 ]; then
        echo "Successfully created: $output_file"
        
        # Optional: Display some basic stats about the output file
        total_reads=$(wc -l < "$output_file")
        read_count=$((total_reads / 4))
        echo "  -> Generated $read_count reads (sampling rate: $rate)"
    else
        echo "Error: Failed to create $output_file"
        exit 1
    fi
    
    echo "----------------------------------------"
done

echo "Subsampling completed successfully!"
echo "All output files saved to: $output_dir"