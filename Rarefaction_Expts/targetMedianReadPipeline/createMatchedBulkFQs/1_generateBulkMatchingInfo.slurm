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
#SBATCH --time=0-04:00:00


# load modules needed
module load Anaconda3/2024.02-1
module load GCC/11.3.0  
module load OpenMPI/4.1.4 
module load R/4.4.0

echo "Modules loaded"

# activate conda environment
eval "$(conda shell.bash hook)"
conda activate PythonGenericEnv

echo "Conda environment activated"

# Define paths to python and R scripts
pythonScriptPath="/data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/scripts/createMatchedBulkFQs/returnBarcodeReadInfo.py"
RscriptPath="/data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/scripts/createMatchedBulkFQs/createMatchedBulkFQ_helper_script.R"

output_dir="/data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/data/LongBench_All/ont_sc/matched_ont_bulk" # make sure this is the correct output directory

input_fastq_path="/data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/data/LongBench_All/ont_sc/0.75_percentRarified.fastq" # make sure this is the correct rarified fastq file
prefix="75percentRarified"  # rarefaction prefix for the output csv files

# run the python script to get the number of reads for each barcode
output_csv_read_info_path="/data/gpfs/projects/punim2251/LongBench_data/data_diagnostics/${prefix}_ont_sc_read_info.csv" #
python ${pythonScriptPath} ${input_fastq_path} ${output_csv_read_info_path}
echo "Python script completed"

# run the R script to calculate the number of reads we need to pull from each sample's bulk fastq file to preserve the cell line proportions in our 'matched' bulk fastq files
# Note: total number of reads in the rarified fastq file = 68503055 ** MAKE SURE THIS IS THE CORRECT NUMBER **
Rscript ${RscriptPath} \
${output_csv_read_info_path} \
"ont_sc" \
${prefix} \
68503055 \
${output_dir}
echo "R script completed, output files saved in ${output_dir}"
