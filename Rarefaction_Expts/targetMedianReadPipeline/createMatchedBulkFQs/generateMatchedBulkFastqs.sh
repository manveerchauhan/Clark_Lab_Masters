#!/bin/bash
# Script to create randomly subsampled FASTQ files for a list of target cell lines with specific read counts
# Usage: ./generateMatchedBulkFastqs.sh <cell_line_info_csv> <input_fastq_directory>
# Example: ./generateMatchedBulkFastqs.sh cell_line_reads.csv /path/to/fastq/files
# 
# Description:
# This script reads a CSV file containing cell line names and their target read counts, then processes each
# cell line's corresponding FASTQ file to generate subsampled files with the specified number of reads.
# The script uses `bbmap`'s `reformat.sh` for subsampling and ensures reproducibility by setting a random seed.
# 
# Requirements:
# 1. Ensure `bbmap` is installed and available in the current conda environment.
# 2. Input FASTQ files must be named in the format `<cellline>.fastq` and located in the specified input directory.
# 3. The CSV file must have the following format (without headers):
#    <cell_line_name>,<target_read_count>
#
# Outputs:
# - Subsampled FASTQ files will be saved in a new directory called `subsampled_fastqs`.
# - A combined subsampled FASTQ file, containing all subsampled reads, will be created in the same directory.
#
# Note: Modify the `cellLine_read_info_file` and `input_dir` variables for testing purposes if needed.

# Set a random seed for reproducibility
seed=42

# Check if correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <cell_line_info_csv> <input_fastq_directory>"
    echo "Example: $0 cell_line_reads.csv /path/to/fastq/files"
    exit 1
fi

# Store command line arguments
cellLine_read_info_file="$1"
input_dir="$2"

# Declare the associative array
declare -A reads_required

# Read the CSV file into the associative array
while IFS=, read -r cellline reads; do
    if [[ $cellline != "merged_cellline_anno" ]]; then
        reads_required[$cellline]=$reads
    fi
done < "${cellLine_read_info_file}"

# Create output directory if it doesn't exist
mkdir -p subsampled_fastqs

# Process each cell line
for cellline in "${!reads_required[@]}"; do
    input_file="${input_dir}/${cellline}.fastq"
    output_file="subsampled_fastqs/${cellline}_subsampled.fastq"
    
    if [[ -f "$input_file" ]]; then
        echo "Subsampling ${cellline} to ${reads_required[$cellline]} reads..."
        reformat.sh in="$input_file" out="$output_file" samplereadstarget="${reads_required[$cellline]}" sampleseed="$seed"
    else
        echo "Warning: Input file not found for ${cellline}: ${input_file}"
    fi
done

# Merge all subsampled files into one combined file
echo "Merging all subsampled files into one combined file..."
cat subsampled_fastqs/*_subsampled.fastq > subsampled_fastqs/combined_subsampled.fastq

echo "Done! Combined file created at: subsampled_fastqs/combined_subsampled.fastq"