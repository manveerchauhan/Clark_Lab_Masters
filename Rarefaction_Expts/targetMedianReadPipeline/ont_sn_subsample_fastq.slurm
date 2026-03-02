#!/bin/bash
#SBATCH --partition="sapphire"
#SBATCH --nodes=1
#SBATCH --account="punim2251"
#SBATCH --ntasks=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=mschauhan@student.unimelb.edu.au
#SBATCH --cpus-per-task=4
#SBATCH --mem=70gb
#SBATCH --time=0-3:00:00
#SBATCH --job-name=fastq_subsample

# Load modules
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load R/4.4.0
module load Anaconda3/2024.02-1

# Activate conda env with python and bbmap
eval "$(conda shell.bash hook)"
conda activate bbmap_env
echo "Conda environment initialized"

# Change to the script directory
cd /data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/scripts/targetMedianReadPipeline

# Define input parameters
fastq="/data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/data/Longbench_All_run2/ont_sn_balanced_bcs/ont_sn_balanced_matched_reads.fastq"
outputDir="/data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/data/Longbench_All_run2/ont_sn_balanced_bcs"
prefix="ont_sn_subsampled"

# Define list of sample rates to use
sample_rates=(0.4)

echo "Starting subsampling job"
echo "Input FASTQ: $fastq"
echo "Output directory: $outputDir"
echo "Prefix: $prefix"
echo "Sample rates: ${sample_rates[@]}"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Working directory: $(pwd)"
echo "----------------------------------------"

# Check if the subsample script exists
if [ ! -f "./3_subsample_fastq_with_rates.sh" ]; then
    echo "Error: 3_subsample_fastq_with_rates.sh not found in current directory!"
    echo "Current directory: $(pwd)"
    echo "Files in directory:"
    ls -la
    exit 1
fi

# Make sure the script is executable
chmod +x ./3_subsample_fastq_with_rates.sh

# Get original file size for reporting
if [ -f "$fastq" ]; then
    original_size=$(wc -l < "$fastq")
    original_reads=$((original_size / 4))
    echo "Original file contains: $original_reads reads"
else
    echo "Warning: Could not access original file for size calculation"
fi

echo "----------------------------------------"
echo "Calling subsample script..."

# Record start time
start_time=$(date +%s)

# Call the existing subsample script with all the sample rates
./3_subsample_fastq_with_rates.sh "$fastq" "$outputDir" "$prefix" "${sample_rates[@]}"

# Check if the script succeeded
if [ $? -eq 0 ]; then
    # Calculate total processing time
    end_time=$(date +%s)
    total_duration=$((end_time - start_time))
    
    echo "----------------------------------------"
    echo "Subsampling completed successfully!"
    echo "Total processing time: ${total_duration}s"
    echo "Job completed at: $(date)"
    echo ""
    
    # Display generated files with their sizes
    echo "Generated files:"
    for rate in "${sample_rates[@]}"; do
        output_file="${outputDir}/${rate}_${prefix}.fastq"
        if [ -f "$output_file" ]; then
            file_size=$(du -h "$output_file" | cut -f1)
            output_size=$(wc -l < "$output_file")
            output_reads=$((output_size / 4))
            echo "  - $(basename $output_file): $output_reads reads ($file_size)"
        else
            echo "  - $(basename $output_file): FILE NOT FOUND!"
        fi
    done
    
    # Create a summary log file
    summary_file="${outputDir}/subsampling_summary_$(date +%Y%m%d_%H%M%S).log"
    {
        echo "Subsampling Summary"
        echo "==================="
        echo "Job ID: $SLURM_JOB_ID"
        echo "Date: $(date)"
        echo "Input file: $fastq"
        echo "Output directory: $outputDir"
        echo "Prefix: $prefix"
        echo "Total processing time: ${total_duration}s"
        if [ -n "$original_reads" ]; then
            echo "Original reads: $original_reads"
        fi
        echo ""
        echo "Sample rates processed: ${sample_rates[@]}"
        echo ""
        echo "Generated files:"
        for rate in "${sample_rates[@]}"; do
            output_file="${outputDir}/${rate}_${prefix}.fastq"
            if [ -f "$output_file" ]; then
                output_size=$(wc -l < "$output_file")
                output_reads=$((output_size / 4))
                file_size=$(du -h "$output_file" | cut -f1)
                if [ -n "$original_reads" ]; then
                    actual_rate=$(echo "scale=4; $output_reads / $original_reads" | bc -l)
                    echo "  $rate -> $(basename $output_file): $output_reads reads (actual rate: $actual_rate, size: $file_size)"
                else
                    echo "  $rate -> $(basename $output_file): $output_reads reads (size: $file_size)"
                fi
            else
                echo "  $rate -> $(basename $output_file): FILE NOT CREATED"
            fi
        done
    } > "$summary_file"
    
    echo ""
    echo "Summary log created: $summary_file"
    
else
    echo "Error: Subsampling script failed!"
    exit 1
fi