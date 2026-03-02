## Load supporting resources------
# Load cell cycle markers obtained from : https://hbctraining.github.io/scRNA-seq_online/lessons/06_SC_SCT_normalization.html
load("/data/gpfs/projects/punim0646/manveer/cycle.rda")
# import silhouette functions from silhouette R script (rrrSingleCellUtils github repo)
source("/data/gpfs/projects/punim2251/scRNAseq_scripts/silhouette.R")

# FLAMES output directory for current project
FLAMES_OUTPUT_DIR <- "/data/gpfs/projects/punim2251/Neurodevelopmental_Models/Cortical_Org_Diff/flame_outs_standard_ref"
# prefix for reading in formatted iso level count matrices from previous run
FORMATTED_ISO_MATRICES = "/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space/"

#install packages if required. Note some packages require installation via bioconductor. See installation instruction for each package to ensure installation is successful. 
library(rtracklayer)
library(Seurat)
library(DropletUtils)
library(gridExtra)
library(data.table)
library(BiocParallel)
library(celda)
library(SingleCellExperiment)
library(DoubletFinder)
library(stringr)
library(cowplot)
library(grid)
library(patchwork)
library(tidyverse)
library(Matrix)
library(tibble)
library(GenomicFeatures)
#library(Gviz)
#library(BSgenome.Hsapiens.UCSC.hg38)

# Set working directory and create folders for output files
setwd("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space")  # Set this to correct location


# Sample information
sample_metadata <- data.frame(
  sample_id = c("org_1A", "org_1B", "org_3A", "org_3B", "org_3C", "org_6A", "org_6B", "org_6C"),
  timepoint = c("1M_Org", "1M_Org", "3M_Org", "3M_Org", "3M_Org", "6M_Org", "6M_Org", "6M_Org"),
  timepoint_days = c(30, 30, 90, 90, 90, 180, 180, 180),
  replicate = c("A", "B", "A", "B", "C", "A", "B", "C"),
  stringsAsFactors = FALSE
)
dir.create("./output_files/ref_files", recursive = TRUE, showWarnings = FALSE)
dir.create("./output_files/counts", recursive = TRUE, showWarnings = FALSE)
dir.create("./output_files/seu_objects", recursive = TRUE, showWarnings = FALSE)
dir.create("./output_files/QC", recursive = TRUE, showWarnings = FALSE)
dir.create("./output_files/DE", recursive = TRUE, showWarnings = FALSE)


### Create naming dictionary for downstream isoform count matrix renaming-------
# Function to make csv naming resource 
make_isoform_gene_symbol_dict <- function(FLAMES_gtf, 
                                          reference_gtf, 
                                          output_file) {
  # Import the first GTF file (transcripts GTF)
  gtf1 <- import(FLAMES_gtf)
  gtf1_df <- as.data.frame(gtf1)
  
  # Select relevant columns from the first GTF
  selected_columns1 <- gtf1_df[, c("transcript_id", "gene_id")]
  unique_selected_cols <- unique(selected_columns1)
  
  # Import the second GTF file (reference GTF with gene symbols)
  gtf2 <- import(reference_gtf)
  gtf2_df <- as.data.frame(gtf2)
  
  # Select relevant columns from the second GTF
  selected_columns2 <- gtf2_df[, c("gene_name", "gene_id")]
  unique_gene_symbol <- unique(selected_columns2)
  
  # Merge the two data frames on 'gene_id'
  combined_data <- merge(unique_selected_cols, 
                         unique_gene_symbol, 
                         by = "gene_id", 
                         all.x = TRUE)
  
  # If 'gene_name' is missing, replace it with 'gene_id'
  combined_data$gene_symbol <- ifelse(is.na(combined_data$gene_name), 
                                      combined_data$gene_id, 
                                      combined_data$gene_name)
  
  # Select relevant columns
  final_combined_data <- combined_data[, c("transcript_id", "gene_id", "gene_symbol")]
  
  # Write to a CSV file
  write.csv(final_combined_data, file = file.path("output_files/ref_files", output_file), row.names = FALSE)
  
  
  return(final_combined_data)
}

# The FLAMES ref can be found in your selected output folder after running the Flames pipeline. 
FLAMES_gtf_file <- paste0(FLAMES_OUTPUT_DIR, "/isoform_annotated.gtf") #ensure file is unzipped
reference_gtf_file <- "/data/gpfs/projects/punim2251/Organoid_NDR_FLAMES_ref/modified_gencode_v47_primary_annotation.gtf" # ensure file is unzipped
output_file <- "isoform_gene_dict.csv"


# Call the helper function defined in code block above to create a dictionary containing corresponding gene information for each isoform
# This may take a few minutes 
isoform_gene_dict <- make_isoform_gene_symbol_dict(FLAMES_gtf_file,
                                                   reference_gtf_file,
                                                   output_file)


### Convert oarfish outputs into seurat compatible count matrix ./output_files/counts/-----
#This function reads in Oarfish count files and creates a csv file of count data. The function also appends the gene symbol to the ENSTID
process_oarfish_files_to_counts_matrix <- function(flames_output_folder, sample_name, output_dir) {
  
  # Read in the resource table (transcript_id, gene_id, gene_symbol)
  # Define the file paths based on the sample name
  count_matrix_path <- file.path(flames_output_folder, paste0(sample_name, ".count.mtx"))
  barcodes_path <- file.path(flames_output_folder, paste0(sample_name, ".barcodes.txt"))
  features_path <- file.path(flames_output_folder, paste0(sample_name, ".features.txt"))
  
  # Read the data
  counts <- readMM(count_matrix_path)
  barcodes <- readLines(barcodes_path)
  features <- read.delim(features_path, header = FALSE)
  
  # Transpose the matrix for Seurat compatibility
  counts <- t(counts)
  
  # Set row and column names
  rownames(counts) <- features$V1
  colnames(counts) <- barcodes
  
  # Convert to a data frame
  counts_df <- as.data.frame(counts)
  
  # Add transcript_id as the first column
  counts_df$transcript_id <- rownames(counts_df)
  counts_df <- counts_df[, c(ncol(counts_df), 1:(ncol(counts_df)-1))]
  
  # Merge with the resource table to add gene symbols
  df_genesymbol <- counts_df %>%
    left_join(isoform_gene_dict, by = "transcript_id")
  
  # Remove the gene_id column and reorder the columns
  df_genesymbol$gene_id <- NULL
  df_genesymbol <- df_genesymbol[, c(ncol(df_genesymbol), 1:(ncol(df_genesymbol)-1))]
  
  # Update row names to include gene symbol instead of transcript_id
  rownames(df_genesymbol) <- paste0(df_genesymbol$transcript_id, "_", df_genesymbol$gene_symbol)
  df_genesymbol$transcript_id <- NULL
  df_genesymbol$gene_symbol <- NULL
  
  # Write the output to a CSV file
  output_path <- file.path(output_dir, paste0("gene_symbol_", sample_name, "_counts.csv"))
  write.csv(df_genesymbol, output_path)
  
  cat("Processed sample:", sample_name, "\nOutput saved to:", output_path, "\n")
  
  return(df_genesymbol)
}


process_oarfish_files_to_counts_matrix(
  flames_output_folder = FLAMES_OUTPUT_DIR,
  sample_name = "org_1A",
  output_dir = "./output_files/counts/"
)

process_oarfish_files_to_counts_matrix(
  flames_output_folder = FLAMES_OUTPUT_DIR,
  sample_name = "org_1B",
  output_dir = "./output_files/counts/"
)

process_oarfish_files_to_counts_matrix(
  flames_output_folder = FLAMES_OUTPUT_DIR,
  sample_name = "org_3A",
  output_dir = "./output_files/counts/"
)

process_oarfish_files_to_counts_matrix(
  flames_output_folder = FLAMES_OUTPUT_DIR,
  sample_name = "org_3B",
  output_dir = "./output_files/counts/"
)

process_oarfish_files_to_counts_matrix(
  flames_output_folder = FLAMES_OUTPUT_DIR,
  sample_name = "org_3C",
  output_dir = "./output_files/counts/"
)

process_oarfish_files_to_counts_matrix(
  flames_output_folder = FLAMES_OUTPUT_DIR,
  sample_name = "org_6A",
  output_dir = "./output_files/counts/"
)

process_oarfish_files_to_counts_matrix(
  flames_output_folder = FLAMES_OUTPUT_DIR,
  sample_name = "org_6B",
  output_dir = "./output_files/counts/"
)

process_oarfish_files_to_counts_matrix(
  flames_output_folder = FLAMES_OUTPUT_DIR,
  sample_name = "org_6C",
  output_dir = "./output_files/counts/"
)
  




### Perform gene lvl qc


### RUN FROM HERE AFTER INITIAL RUN: Create and filter gene level seu objs-----
# Sample information
sample_metadata <- data.frame(
  sample_id = c("org_1A", "org_1B", "org_3A", "org_3B", "org_3C", "org_6A", "org_6B", "org_6C"),
  timepoint = c("1M_Org", "1M_Org", "3M_Org", "3M_Org", "3M_Org", "6M_Org", "6M_Org", "6M_Org"),
  timepoint_days = c(30, 30, 90, 90, 90, 180, 180, 180),
  replicate = c("A", "B", "A", "B", "C", "A", "B", "C"),
  
  # read in gene level count matrices to sample metadata
  gene_level_count_matrices = c(read.csv(paste0(FLAMES_OUTPUT_DIR, "/org_1A_gene_count.csv")),
                                read.csv(paste0(FLAMES_OUTPUT_DIR, "/org_1B_gene_count.csv")),
                                read.csv(paste0(FLAMES_OUTPUT_DIR, "/org_3A_gene_count.csv")),
                                read.csv(paste0(FLAMES_OUTPUT_DIR, "/org_3B_gene_count.csv")),
                                read.csv(paste0(FLAMES_OUTPUT_DIR, "/org_3C_gene_count.csv")),
                                read.csv(paste0(FLAMES_OUTPUT_DIR, "/org_6A_gene_count.csv")),
                                read.csv(paste0(FLAMES_OUTPUT_DIR, "/org_6B_gene_count.csv")),
                                read.csv(paste0(FLAMES_OUTPUT_DIR, "/org_6C_gene_count.csv"))),
  
  # read in isoform level count matrices
  iso_level_count_matrices = c(read.csv(paste0(FORMATTED_ISO_MATRICES, "output_files/counts/gene_symbol_org_1A_counts.csv")),
                               read.csv(paste0(FORMATTED_ISO_MATRICES, "output_files/counts/gene_symbol_org_1B_counts.csv")),
                               read.csv(paste0(FORMATTED_ISO_MATRICES, "output_files/counts/gene_symbol_org_3A_counts.csv")),
                               read.csv(paste0(FORMATTED_ISO_MATRICES, "output_files/counts/gene_symbol_org_3B_counts.csv")),
                               read.csv(paste0(FORMATTED_ISO_MATRICES, "output_files/counts/gene_symbol_org_3C_counts.csv")),
                               read.csv(paste0(FORMATTED_ISO_MATRICES, "output_files/counts/gene_symbol_org_6A_counts.csv")),
                               read.csv(paste0(FORMATTED_ISO_MATRICES, "output_files/counts/gene_symbol_org_6B_counts.csv")),
                               read.csv(paste0(FORMATTED_ISO_MATRICES, "output_files/counts/gene_symbol_org_6C_counts.csv"))),
  
  stringsAsFactors = FALSE
)


### test space----

sample_metadata <- tibble(
  sample_id      = c("org_1A","org_1B","org_3A","org_3B","org_3C","org_6A","org_6B","org_6C"),
  timepoint      = c("1M_Org","1M_Org","3M_Org","3M_Org","3M_Org","6M_Org","6M_Org","6M_Org"),
  timepoint_days = c(30, 30, 90, 90, 90, 180, 180, 180),
  replicate      = c("A","B","A","B","C","A","B","C"),
  
  gene_level_count_matrices = list(
    read.csv(file.path(FLAMES_OUTPUT_DIR, "org_1A_gene_count.csv"), check.names = FALSE),
    read.csv(file.path(FLAMES_OUTPUT_DIR, "org_1B_gene_count.csv"), check.names = FALSE),
    read.csv(file.path(FLAMES_OUTPUT_DIR, "org_3A_gene_count.csv"), check.names = FALSE),
    read.csv(file.path(FLAMES_OUTPUT_DIR, "org_3B_gene_count.csv"), check.names = FALSE),
    read.csv(file.path(FLAMES_OUTPUT_DIR, "org_3C_gene_count.csv"), check.names = FALSE),
    read.csv(file.path(FLAMES_OUTPUT_DIR, "org_6A_gene_count.csv"), check.names = FALSE),
    read.csv(file.path(FLAMES_OUTPUT_DIR, "org_6B_gene_count.csv"), check.names = FALSE),
    read.csv(file.path(FLAMES_OUTPUT_DIR, "org_6C_gene_count.csv"), check.names = FALSE)
  ),
  
  iso_level_count_matrices = list(
    read.csv(file.path(FORMATTED_ISO_MATRICES, "output_files/counts/gene_symbol_org_1A_counts.csv"), check.names = FALSE),
    read.csv(file.path(FORMATTED_ISO_MATRICES, "output_files/counts/gene_symbol_org_1B_counts.csv"), check.names = FALSE),
    read.csv(file.path(FORMATTED_ISO_MATRICES, "output_files/counts/gene_symbol_org_3A_counts.csv"), check.names = FALSE),
    read.csv(file.path(FORMATTED_ISO_MATRICES, "output_files/counts/gene_symbol_org_3B_counts.csv"), check.names = FALSE),
    read.csv(file.path(FORMATTED_ISO_MATRICES, "output_files/counts/gene_symbol_org_3C_counts.csv"), check.names = FALSE),
    read.csv(file.path(FORMATTED_ISO_MATRICES, "output_files/counts/gene_symbol_org_6A_counts.csv"), check.names = FALSE),
    read.csv(file.path(FORMATTED_ISO_MATRICES, "output_files/counts/gene_symbol_org_6B_counts.csv"), check.names = FALSE),
    read.csv(file.path(FORMATTED_ISO_MATRICES, "output_files/counts/gene_symbol_org_6C_counts.csv"), check.names = FALSE)
  )
)

# SAVE AND RELOAD ORG EXPERIMENT INPUT DATA OBJ-------
# Save sample metadata object to a file
saveRDS(sample_metadata, file = "org_flames_v2.2.0_experiment_data.rds")

# Restore the sample metadata object for future runs
sample_metadata = readRDS(file = "org_flames_v2.2.0_experiment_data.rds")

# create seurat objs------

# Updated function that handles duplicates
create_seurat_from_counts <- function(count_matrix_df, sample_id, assay_type = "RNA") {
  
  # Set the first column as rownames and remove it
  gene_names <- count_matrix_df[[1]]
  count_matrix_df <- count_matrix_df[, -1, drop = FALSE]
  
  # Check for and handle duplicate rownames
  if(any(duplicated(gene_names))) {
    cat("Warning: Found", sum(duplicated(gene_names)), "duplicate gene names in", sample_id, "\n")
    cat("Making rownames unique...\n")
    gene_names <- make.unique(gene_names)
  }
  
  rownames(count_matrix_df) <- gene_names
  
  # Convert to matrix
  count_matrix <- as.matrix(count_matrix_df)
  
  # PRE-FILTER before converting to sparse matrix
  # Filter genes: keep genes detected in at least 3 cells
  genes_keep <- rowSums(count_matrix > 0, na.rm = TRUE) >= 3
  count_matrix <- count_matrix[genes_keep, ]
  
  # Filter cells: keep cells with at least 1 feature
  cells_keep <- colSums(count_matrix > 0, na.rm = TRUE) >= 10
  count_matrix <- count_matrix[, cells_keep]
  
  # Now convert to sparse matrix
  count_sparse <- Matrix::Matrix(count_matrix, sparse = TRUE)
  
  # Create Seurat object WITHOUT additional filtering (already filtered)
  seu_obj <- CreateSeuratObject(
    counts = count_sparse,
    project = sample_id,
    assay = assay_type,
    min.cells = 0,      # Already filtered above
    min.features = 0    # Already filtered above
  )
  
  # Add sample metadata
  seu_obj$sample_id <- sample_id
  
  return(seu_obj)
}

### TEST SPACE 1: PROCESS AND SAVE SEU OBJS ONE AT A TIME-----
# Process and save ONE sample at a time to avoid memory overload
process_single_sample <- function(sample_idx, metadata) {
  
  sample_id <- metadata$sample_id[sample_idx]
  cat("\n=== Processing sample", sample_idx, "of", nrow(metadata), ":", sample_id, "===\n")
  
  # Process gene-level
  cat("Creating gene-level Seurat object...\n")
  gene_seu <- create_seurat_from_counts(
    metadata$gene_level_count_matrices[[sample_idx]], 
    sample_id,
    assay_type = "RNA"
  )
  
  # Save immediately
  gene_file <- paste0("./output_files/seu_objects/", sample_id, "_gene_seurat.rds")
  saveRDS(gene_seu, file = gene_file)
  cat("Saved:", gene_file, "\n")
  
  # Clear from memory
  rm(gene_seu)
  gc()
  
  # Process isoform-level
  cat("Creating isoform-level Seurat object...\n")
  iso_seu <- create_seurat_from_counts(
    metadata$iso_level_count_matrices[[sample_idx]], 
    sample_id,
    assay_type = "isoform"
  )
  
  # Save immediately
  iso_file <- paste0("./output_files/seu_objects/", sample_id, "_isoform_seurat.rds")
  saveRDS(iso_seu, file = iso_file)
  cat("Saved:", iso_file, "\n")
  
  # Clear from memory
  rm(iso_seu)
  gc()
  
  return(sample_id)
}

# Process all samples one by one
cat("Starting sequential processing of all samples...\n")
for(i in 1:nrow(sample_metadata)) {
  process_single_sample(i, sample_metadata)
}

cat("\n=== All samples processed successfully! ===\n")

# Later, when you need to work with them, load individually:
# gene_seu_1A <- readRDS("./output_files/seu_objects/org_1A_gene_seurat.rds")
# iso_seu_1A <- readRDS("./output_files/seu_objects/org_1A_isoform_seurat.rds")
## TEST SPACE 2:-----
# Extract the first sample's data
sample_1_id <- sample_metadata$sample_id[1]
gene_counts_1 <- sample_metadata$gene_level_count_matrices[[1]]
iso_counts_1 <- sample_metadata$iso_level_count_matrices[[1]]
length(sample_metadata$sample_id)

# Test this version
org1A_gene_seu_obj_raw <- create_seurat_from_counts(
  sample_metadata$gene_level_count_matrices[[1]], 
  sample_metadata$sample_id[1],
  assay_type = "RNA"
)

# Test on isoform too
test_iso <- create_seurat_from_counts(
  iso_counts_1,
  sample_1_id,
  assay_type = "isoform"
)

print(test_iso)


### CREATE AND SAVE GENE AND ISO LVL SEU OBJS----
# Apply function to create Seurat objects for all samples
sample_metadata <- sample_metadata %>%
  mutate(
    # Create gene-level Seurat objects
    gene_seurat = map2(
      gene_level_count_matrices, 
      sample_id,
      ~create_seurat_from_counts(.x, .y, assay_type = "RNA")
    ),
    
    # Create isoform-level Seurat objects
    isoform_seurat = map2(
      iso_level_count_matrices, 
      sample_id,
      ~create_seurat_from_counts(.x, .y, assay_type = "isoform")
    )
  )



# Optional: Save individual Seurat objects
for(i in 1:nrow(sample_metadata)) {
  # Save gene-level objects
  saveRDS(
    sample_metadata$gene_seurat[[i]], 
    file = paste0("./output_files/seu_objects/", 
                  sample_metadata$sample_id[i], "_gene_seurat.rds")
  )
  
  # Save isoform-level objects
  saveRDS(
    sample_metadata$isoform_seurat[[i]], 
    file = paste0("./output_files/seu_objects/", 
                  sample_metadata$sample_id[i], "_isoform_seurat.rds")
  )
}

# View the updated metadata structure
print(sample_metadata)
# Verify creation
cat("\n=== Summary of Created Seurat Objects ===\n")
for(i in 1:nrow(sample_metadata)) {
  cat("\nSample:", sample_metadata$sample_id[i], "\n")
  cat("  Gene-level:\n")
  print(sample_metadata$gene_seurat[[i]])
  cat("  Isoform-level:\n")
  print(sample_metadata$isoform_seurat[[i]])
}

# Save individual Seurat objects to disk
cat("\n=== Saving Seurat Objects ===\n")
for(i in 1:nrow(sample_metadata)) {
  # Save gene-level objects
  saveRDS(
    sample_metadata$gene_seurat[[i]], 
    file = paste0("./output_files/seu_objects/", 
                  sample_metadata$sample_id[i], "_gene_seurat.rds")
  )
  cat("Saved:", sample_metadata$sample_id[i], "gene Seurat object\n")
  
  # Save isoform-level objects
  saveRDS(
    sample_metadata$isoform_seurat[[i]], 
    file = paste0("./output_files/seu_objects/", 
                  sample_metadata$sample_id[i], "_isoform_seurat.rds")
  )
  cat("Saved:", sample_metadata$sample_id[i], "isoform Seurat object\n")
}

# Save the updated sample_metadata with Seurat objects
saveRDS(sample_metadata, file = "org_flames_v2.2.0_experiment_data_with_seurat.rds")
cat("\nSaved complete sample metadata object!\n")

# Quick summary
cat("\n=== Final Summary ===\n")
cat("Total samples processed:", nrow(sample_metadata), "\n")
cat("Gene Seurat objects created:", sum(!sapply(sample_metadata$gene_seurat, is.null)), "\n")
cat("Isoform Seurat objects created:", sum(!sapply(sample_metadata$isoform_seurat, is.null)), "\n")



## OLD CODE------

# go through and qc 1 organoid first manually, then make a report

### Create matching seu obj with iso counts----

# read in isoform level count matrices
org1A_unfiltered_iso_cm <- read.csv(paste0(FORMATTED_ISO_MATRICES, "output_files/counts/gene_symbol_org_1A_counts.csv"))
org1B_unfiltered_iso_cm <- read.csv(paste0(FORMATTED_ISO_MATRICES, "output_files/counts/gene_symbol_org_1B_counts.csv"))

org3A_unfiltered_iso_cm <- read.csv(paste0(FORMATTED_ISO_MATRICES, "output_files/counts/gene_symbol_org_3A_counts.csv"))
org3B_unfiltered_iso_cm <- read.csv(paste0(FORMATTED_ISO_MATRICES, "output_files/counts/gene_symbol_org_3B_counts.csv"))
org3C_unfiltered_iso_cm <- read.csv(paste0(FORMATTED_ISO_MATRICES, "output_files/counts/gene_symbol_org_3C_counts.csv"))

org6A_unfiltered_iso_cm <- read.csv(paste0(FORMATTED_ISO_MATRICES, "output_files/counts/gene_symbol_org_6A_counts.csv"))
org6B_unfiltered_iso_cm <- read.csv(paste0(FORMATTED_ISO_MATRICES, "output_files/counts/gene_symbol_org_6B_counts.csv"))
org6C_unfiltered_iso_cm <- read.csv(paste0(FORMATTED_ISO_MATRICES, "output_files/counts/gene_symbol_org_6C_counts.csv"))
