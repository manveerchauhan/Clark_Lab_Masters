# script to look at cell type proportions in seurat object
library(Seurat)
library(tidyverse)
library(gridExtra)
library(speckle)

theme_set(ggmin::theme_min())  # pre-set ggplot theme
##
#   Add cell-type proportion graphs using this: https://academic.oup.com/bioinformatics/article/38/20/4720/6675456
#   https://github.com/Oshlack/speckle
##

# import organoid seurat objects (consensus annotated with isoform assay)
one_month_org <- readRDS("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space/output_files/integrated_objects/1M_Org_integrated_harmony_consensus_with_isoform_assay.rds")
three_month_org <- readRDS("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space/output_files/integrated_objects/3M_Org_integrated_harmony_consensus_with_isoform_assay.rds")
six_month_org <- readRDS("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space/output_files/integrated_objects/6M_Org_integrated_harmony_consensus_with_isoform_assay.rds")

# Note: timepoint metadata already exists in these objects, no need to add

# function to create cell type proportion plot
MakeGlobalCellTypePropPlot <- function(seu.obj){
  plt <- plotCellTypeProps(clusters = seu.obj@meta.data$consensus_cell_type,
                    sample = seu.obj@meta.data$timepoint)

  return(plt)
}

# function to create cell type proportion plot with replicates
MakeReplicateCellTypePropPlot <- function(seu.obj){
  plt <- plotCellTypeProps(clusters = seu.obj@meta.data$consensus_cell_type,
                           sample = seu.obj@meta.data$orig.ident)

  return(plt)
}

Org_1M_global <- MakeGlobalCellTypePropPlot(one_month_org)
Org_3M_global <- MakeGlobalCellTypePropPlot(three_month_org)
Org_6M_global <- MakeGlobalCellTypePropPlot(six_month_org)

Org_1M_replicate <- MakeReplicateCellTypePropPlot(one_month_org)
Org_3M_replicate <- MakeReplicateCellTypePropPlot(three_month_org)
Org_6M_replicate <- MakeReplicateCellTypePropPlot(six_month_org)


org.global.prop.plts <- grid.arrange(Org_1M_global,
                                     Org_3M_global,
                                     Org_6M_global,
                                     ncol = 3)
org.replicate.prop.plts <- grid.arrange(Org_1M_replicate,
                                        Org_3M_replicate,
                                        Org_6M_replicate,
                                        ncol = 3)

# Create output directory if it doesn't exist
output_dir <- "./output_files/cell_type_proportions"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save plots to output directory
ggsave(file.path(output_dir, "org_global_cellTypeProp.pdf"), org.global.prop.plts, width = 30,
       height = 9, dpi = 300)
ggsave(file.path(output_dir, "org_replicate_cellTypeProp.pdf"), org.replicate.prop.plts, width = 30,
       height = 9, dpi = 300)

cat("Cell type proportion plots saved to:", output_dir, "\n")
