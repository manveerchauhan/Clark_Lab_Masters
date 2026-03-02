#!/bin/bash
#SBATCH --partition="cascade"
#SBATCH --nodes=1
#SBATCH --account="punim2251"
#SBATCH --ntasks=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=mschauhan@student.unimelb.edu.au
#SBATCH --cpus-per-task=4
#SBATCH --mem=200gb
#time in days-hours:mins:sec
#SBATCH --time=0-10:00:00

source /home/mschauhan/anaconda3/bin/activate /home/mschauhan/.conda/envs/bbmap_env

echo "Conda environment initialized"

# CSV file containing the number of bulk reads required from each cell line
cellLine_read_info_file="/data/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/data/LongBench_All/ont_sc/matched_ont_bulk/ont_sc_10percentRarified_bulkReadsNeeded.csv"
# Input directory containing original fastq files for each cell line
input_dir="/data/projects/punim2251/LongBench_data/transfer1/ont_bulk/fastq_files"
# Path to the generateMatchedBulkFastqs.sh script
generateMatchedBulkFastqsScript="/data/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/scripts/createMatchedBulkFQs/generateMatchedBulkFastqs.sh"

# Run script
${generateMatchedBulkFastqsScript} ${cellLine_read_info_file} ${input_dir}