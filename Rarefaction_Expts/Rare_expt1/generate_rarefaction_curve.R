# Script to make rarefaction curve for matched bulk and sc/sn PB and ONT data using filtered matrices

library(tidyverse)
library(purrr)
library(readr)
library(tools) # for basename
library(stringr) # for str_match
library(patchwork)
library(gridExtra) # for publication-ready tables
library(grid) # for table styling

# Set working directory to the pipeline location
setwd("/data/gpfs/projects/punim2251/Aim1_LongBench/longbench-analysis-pipeline")

theme_set(theme_minimal())

# Define paths to filtered matrices
BULK_MATRICES_DIR <- "./bulk_processed_matrices"
SC_MATRICES_DIR <- "./sc_processed_matrices"

## Step 1: Load the filtered matrices -------
# Load filtered single-cell matrices
gene_pseudobulk_ont <- readRDS(file.path(SC_MATRICES_DIR, "filtered_genePseudobulkDfs_ont_sc_gene.rds"))
iso_pseudobulk_ont <- readRDS(file.path(SC_MATRICES_DIR, "filtered_isoPseudobulkDfs_ont_sc_tx.rds"))
gene_pseudobulk_pb <- readRDS(file.path(SC_MATRICES_DIR, "filtered_genePseudobulkDfs_pacbio_sc_gene.rds"))
iso_pseudobulk_pb <- readRDS(file.path(SC_MATRICES_DIR, "filtered_isoPseudobulkDfs_pacbio_sc_tx.rds"))

# Load filtered single-nucleus matrices
gene_pseudobulk_ont_sn <- readRDS(file.path(SC_MATRICES_DIR, "filtered_genePseudobulkDfs_ont_sn_gene.rds"))
iso_pseudobulk_ont_sn <- readRDS(file.path(SC_MATRICES_DIR, "filtered_isoPseudobulkDfs_ont_sn_tx.rds"))
gene_pseudobulk_pb_sn <- readRDS(file.path(SC_MATRICES_DIR, "filtered_genePseudobulkDfs_pacbio_sn_gene.rds"))
iso_pseudobulk_pb_sn <- readRDS(file.path(SC_MATRICES_DIR, "filtered_isoPseudobulkDfs_pacbio_sn_tx.rds"))

# Load filtered bulk matrices
bulk_iso_ont <- readRDS(file.path(BULK_MATRICES_DIR, "bulk_isoDfs_ont_filtered.rds"))
bulk_iso_pb <- readRDS(file.path(BULK_MATRICES_DIR, "bulk_isoDfs_pb_filtered.rds"))
# Load filtered bulk gene matrices
bulk_gene_ont <- readRDS(file.path(BULK_MATRICES_DIR, "bulk_geneDfs_ont_filtered.rds"))
bulk_gene_pb <- readRDS(file.path(BULK_MATRICES_DIR, "bulk_geneDfs_pb_filtered.rds"))

# Load expanded sample sheets for metadata
bulk_sample_sheet <- read.csv(file.path(BULK_MATRICES_DIR, "bulk_sample_sheet_used.csv"))
sc_sample_sheet <- read.csv(file.path(SC_MATRICES_DIR, "sc_sn_sample_sheet_used.csv"))

# Function to generate feature number summary from RDS data frame list
generateFeatureNumSummary <- function(df_list, tech_name, data_type, source_type, sample_sheet) {
  feature_nums <- data.frame()
  
  for (df_name in names(df_list)) {
    current_df <- df_list[[df_name]]
    # Extract rank from the name (e.g., ..._rank1)
    rank_match <- str_match(df_name, "rank(\\d+)")
    rank <- if (!is.na(rank_match[1,2])) as.integer(rank_match[1,2]) else NA
    # Try to match by rank_order if found, else fallback to basename
    if (!is.na(rank)) {
      sample_info <- sample_sheet %>% filter(rank_order == rank)
    } else {
      sample_info <- sample_sheet %>% filter(basename(file_path) == basename(df_name))
    }
    if (nrow(sample_info) == 0) {
      warning(sprintf("No sample sheet match for matrix name: %s", df_name))
      next
    }

    if (is.data.frame(current_df)) {
      feature_count <- nrow(current_df)
      seq_depth <- sample_info$sequencing_depth[1]
      # Extract total_reads if available in sample sheet
      total_reads <- if ("total_reads" %in% colnames(sample_info)) {
        sample_info$total_reads[1]
      } else {
        NA
      }
      sample_id <- if(source_type == "SC") {
        paste0(tech_name, " Single-Cell ", data_type)
      } else if(source_type == "SN") {
        paste0(tech_name, " Single-Nucleus ", data_type)
      } else {
        paste0(tech_name, " ", source_type, " ", data_type)
      }
      new_row <- data.frame(
        rankOrder = rank,
        sequencingDepth = seq_depth,
        total_reads = total_reads,
        featureNum = feature_count,
        sampleID = sample_id,
        technology = tech_name,
        sample_type = source_type,
        type = paste0(data_type, " Discovery")
      )
      feature_nums <- rbind(feature_nums, new_row)
    }
  }
  return(feature_nums)
}

## Step 2: Prepare Rarefaction Curve Dataframe Inputs -------
# Process single-cell data
gene_sc_featureNums_ont <- generateFeatureNumSummary(gene_pseudobulk_ont, "ONT", "Genes", "SC", sc_sample_sheet)
tx_sc_featureNums_ont <- generateFeatureNumSummary(iso_pseudobulk_ont, "ONT", "Isoforms", "SC", sc_sample_sheet)
gene_sc_featureNums_pb <- generateFeatureNumSummary(gene_pseudobulk_pb, "PacBio", "Genes", "SC", sc_sample_sheet)
tx_sc_featureNums_pb <- generateFeatureNumSummary(iso_pseudobulk_pb, "PacBio", "Isoforms", "SC", sc_sample_sheet)

# Process single-nucleus data
gene_sn_featureNums_ont <- generateFeatureNumSummary(gene_pseudobulk_ont_sn, "ONT", "Genes", "SN", sc_sample_sheet)
tx_sn_featureNums_ont <- generateFeatureNumSummary(iso_pseudobulk_ont_sn, "ONT", "Isoforms", "SN", sc_sample_sheet)
gene_sn_featureNums_pb <- generateFeatureNumSummary(gene_pseudobulk_pb_sn, "PacBio", "Genes", "SN", sc_sample_sheet)
tx_sn_featureNums_pb <- generateFeatureNumSummary(iso_pseudobulk_pb_sn, "PacBio", "Isoforms", "SN", sc_sample_sheet)

# Process bulk data
# Isoforms
tx_bulk_featureNums_ont <- generateFeatureNumSummary(bulk_iso_ont, "ONT", "Isoforms", "Bulk", bulk_sample_sheet)
tx_bulk_featureNums_pb <- generateFeatureNumSummary(bulk_iso_pb, "PacBio", "Isoforms", "Bulk", bulk_sample_sheet)
# Genes
gene_bulk_featureNums_ont <- generateFeatureNumSummary(bulk_gene_ont, "ONT", "Genes", "Bulk", bulk_sample_sheet)
gene_bulk_featureNums_pb <- generateFeatureNumSummary(bulk_gene_pb, "PacBio", "Genes", "Bulk", bulk_sample_sheet)

## Step 3: Plotting function ------
plotRarefactionCurve <- function(df_list_input,
                                 plt_title = "Rarefaction Curve",
                                 export_plt = FALSE,
                                 file_prefix = "unnamed_plt",
                                 txt_size = 14,
                                 file_width = 10,
                                 file_height = 6,
                                 removeLegend = FALSE,
                                 output_dir = "rarefaction_plots") {
  all_data <- do.call(rbind, df_list_input)
  gene_data <- all_data %>% filter(type == "Genes Discovery" | type == "Gene Discovery")
  iso_data  <- all_data %>% filter(type == "Isoforms Discovery" | type == "Isoform Discovery")

  # Color mapping for technologies
  color_map <- c(
    "PacBio" = "#DA1884",
    "ONT" = "#003F5C"
  )
  
  # Line type mapping for sample types
  linetype_map <- c(
    "SC" = "dashed",
    "SN" = "dotted", 
    "Bulk" = "solid"
  )
  
  # Point shape mapping for sample types
  shape_map <- c(
    "SC" = 16,    # solid circle
    "SN" = 17,    # solid triangle
    "Bulk" = 15   # solid square
  )

  gene_plt <- ggplot(gene_data, aes(x = factor(rankOrder), y = featureNum, 
                                    color = technology, linetype = sample_type, 
                                    shape = sample_type, group = sampleID)) +
    geom_line(size = 1, alpha = 0.7) +
    geom_point(size = 3, alpha = 0.8) +
    labs(x = "Sequencing Depth Rank",
         y = "Number of Features Detected",
         title = paste(plt_title, "- Genes"),
         color = "Technology",
         linetype = "Sample Type",
         shape = "Sample Type") +
    theme_minimal() +
    theme(text = element_text(size = txt_size),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right") +
    scale_x_discrete(labels = function(x) paste0("Rank ", x)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_color_manual(values = color_map) +
    scale_linetype_manual(values = linetype_map) +
    scale_shape_manual(values = shape_map)

  iso_plt <- ggplot(iso_data, aes(x = factor(rankOrder), y = featureNum, 
                                  color = technology, linetype = sample_type, 
                                  shape = sample_type, group = sampleID)) +
    geom_line(size = 1, alpha = 0.7) +
    geom_point(size = 3, alpha = 0.8) +
    labs(x = "Sequencing Depth Rank",
         y = "Number of Features Detected",
         title = paste(plt_title, "- Isoforms"),
         color = "Technology",
         linetype = "Sample Type",
         shape = "Sample Type") +
    theme_minimal() +
    theme(text = element_text(size = txt_size),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right") +
    scale_x_discrete(labels = function(x) paste0("Rank ", x)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
    scale_color_manual(values = color_map) +
    scale_linetype_manual(values = linetype_map) +
    scale_shape_manual(values = shape_map)

  combined_plt <- gene_plt + iso_plt + plot_layout(ncol = 2, guides = "collect")
  
  # Save plot if requested
  if (export_plt) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    # Save as PNG
    png_filename <- file.path(output_dir, paste0(file_prefix, ".png"))
    ggsave(png_filename, plot = combined_plt, width = file_width * 2, height = file_height, dpi = 300)
    
    # Save as PDF (vector format, good for publications)
    pdf_filename <- file.path(output_dir, paste0(file_prefix, ".pdf"))
    ggsave(pdf_filename, plot = combined_plt, width = file_width * 2, height = file_height, device = "pdf")
    
    # Save as SVG (vector format, good for editing)
    svg_filename <- file.path(output_dir, paste0(file_prefix, ".svg"))
    ggsave(svg_filename, plot = combined_plt, width = file_width * 2, height = file_height, device = "svg")
    
    # Save the data as RDS file
    rds_filename <- file.path(output_dir, paste0(file_prefix, "_data.rds"))
    saveRDS(all_data, rds_filename)
    
    # Create and save sequencing depth summary table
    depth_table <- all_data %>%
      select(technology, sample_type, rankOrder, sequencingDepth) %>%
      distinct() %>%
      arrange(technology, sample_type, rankOrder) %>%
      mutate(
        depth_metric = case_when(
          sample_type %in% c("SC", "SN") ~ "Median reads per cell",
          sample_type == "Bulk" ~ "Total reads",
          TRUE ~ "Unknown"
        )
      ) %>%
      select(technology, sample_type, rankOrder, sequencingDepth, depth_metric)
    
    # Save as CSV
    csv_filename <- file.path(output_dir, paste0(file_prefix, "_sequencing_depths.csv"))
    write.csv(depth_table, csv_filename, row.names = FALSE)
    
    # Save as RDS
    depth_rds_filename <- file.path(output_dir, paste0(file_prefix, "_sequencing_depths.rds"))
    saveRDS(depth_table, depth_rds_filename)
    
    cat("Plots and data saved as:\n")
    cat("- PNG:", png_filename, "\n")
    cat("- PDF:", pdf_filename, "\n") 
    cat("- SVG:", svg_filename, "\n")
    cat("- RDS data:", rds_filename, "\n")
    cat("- Sequencing depths CSV:", csv_filename, "\n")
    cat("- Sequencing depths RDS:", depth_rds_filename, "\n")
  }

  return(list(plot = combined_plt, data = all_data, depth_table = if(export_plt) depth_table else NULL))
}

plotRarefactionCurveContinuous <- function(df_list_input,
                                         plt_title = "Rarefaction Curve",
                                         export_plt = FALSE,
                                         file_prefix = "unnamed_plt",
                                         txt_size = 14,
                                         file_width = 10,
                                         file_height = 6,
                                         removeLegend = FALSE,
                                         output_dir = "rarefaction_plots",
                                         x_axis_metric = "sequencing_depth",  # or "total_reads"
                                         log_scale = FALSE) {
  all_data <- do.call(rbind, df_list_input)
  gene_data <- all_data %>% filter(type == "Genes Discovery" | type == "Gene Discovery")
  iso_data  <- all_data %>% filter(type == "Isoforms Discovery" | type == "Isoform Discovery")

  # Color mapping for technologies (same as original)
  color_map <- c(
    "PacBio" = "#DA1884",
    "ONT" = "#003F5C"
  )
  
  # Line type mapping for sample types (same as original)
  linetype_map <- c(
    "SC" = "dashed",
    "SN" = "dotted", 
    "Bulk" = "solid"
  )
  
  # Point shape mapping for sample types (same as original)
  shape_map <- c(
    "SC" = 16,    # solid circle
    "SN" = 17,    # solid triangle
    "Bulk" = 15   # solid square
  )
  
  # Determine x-axis variable and label
  if (x_axis_metric == "total_reads") {
    x_var <- "total_reads"
    x_label <- "Total Reads"
  } else {
    x_var <- "sequencingDepth"
    x_label <- "Sequencing Depth (Median reads per cell / Total reads for bulk)"
  }

  # Create gene plot
  gene_plt <- ggplot(gene_data, aes_string(x = x_var, y = "featureNum", 
                                          color = "technology", linetype = "sample_type", 
                                          shape = "sample_type", group = "sampleID")) +
    geom_line(size = 1, alpha = 0.7) +
    geom_point(size = 3, alpha = 0.8) +
    labs(x = x_label,
         y = "Number of Features Detected",
         title = paste(plt_title, "- Genes"),
         color = "Technology",
         linetype = "Sample Type",
         shape = "Sample Type") +
    theme_minimal() +
    theme(text = element_text(size = txt_size),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_color_manual(values = color_map) +
    scale_linetype_manual(values = linetype_map) +
    scale_shape_manual(values = shape_map)

  # Create isoform plot
  iso_plt <- ggplot(iso_data, aes_string(x = x_var, y = "featureNum", 
                                        color = "technology", linetype = "sample_type", 
                                        shape = "sample_type", group = "sampleID")) +
    geom_line(size = 1, alpha = 0.7) +
    geom_point(size = 3, alpha = 0.8) +
    labs(x = x_label,
         y = "Number of Features Detected",
         title = paste(plt_title, "- Isoforms"),
         color = "Technology",
         linetype = "Sample Type",
         shape = "Sample Type") +
    theme_minimal() +
    theme(text = element_text(size = txt_size),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
    scale_color_manual(values = color_map) +
    scale_linetype_manual(values = linetype_map) +
    scale_shape_manual(values = shape_map)

  # Apply log scale if requested
  if (log_scale) {
    gene_plt <- gene_plt + scale_x_log10(labels = scales::comma)
    iso_plt <- iso_plt + scale_x_log10(labels = scales::comma)
  } else {
    gene_plt <- gene_plt + scale_x_continuous(labels = scales::comma)
    iso_plt <- iso_plt + scale_x_continuous(labels = scales::comma)
  }

  combined_plt <- gene_plt + iso_plt + plot_layout(ncol = 2, guides = "collect")
  
  # Save plot if requested
  if (export_plt) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Create file suffix based on x-axis metric and log scale
    file_suffix <- paste0("_", x_axis_metric)
    if (log_scale) file_suffix <- paste0(file_suffix, "_log")
    
    # Save as PNG
    png_filename <- file.path(output_dir, paste0(file_prefix, file_suffix, ".png"))
    ggsave(png_filename, plot = combined_plt, width = file_width * 2, height = file_height, dpi = 300)
    
    # Save as PDF
    pdf_filename <- file.path(output_dir, paste0(file_prefix, file_suffix, ".pdf"))
    ggsave(pdf_filename, plot = combined_plt, width = file_width * 2, height = file_height, device = "pdf")
    
    # Save as SVG
    svg_filename <- file.path(output_dir, paste0(file_prefix, file_suffix, ".svg"))
    ggsave(svg_filename, plot = combined_plt, width = file_width * 2, height = file_height, device = "svg")
    
    # Save the data as RDS file
    rds_filename <- file.path(output_dir, paste0(file_prefix, file_suffix, "_data.rds"))
    saveRDS(all_data, rds_filename)
    
    cat("Continuous rarefaction plots and data saved as:\n")
    cat("- PNG:", png_filename, "\n")
    cat("- PDF:", pdf_filename, "\n") 
    cat("- SVG:", svg_filename, "\n")
    cat("- RDS data:", rds_filename, "\n")
  }

  return(list(plot = combined_plt, data = all_data))
}

# Function to create and display sequencing depth summary table
create_depth_summary_table <- function(all_data) {
  depth_table <- all_data %>%
    select(technology, sample_type, rankOrder, sequencingDepth) %>%
    distinct() %>%
    arrange(technology, sample_type, rankOrder) %>%
    mutate(
      depth_metric = case_when(
        sample_type %in% c("SC", "SN") ~ "Median reads per cell",
        sample_type == "Bulk" ~ "Total reads",
        TRUE ~ "Unknown"
      ),
      sequencingDepth_formatted = case_when(
        sample_type == "Bulk" ~ format(sequencingDepth, big.mark = ",", scientific = FALSE),
        TRUE ~ format(round(sequencingDepth, 0), big.mark = ",")
      )
    ) %>%
    select(technology, sample_type, rankOrder, sequencingDepth_formatted, depth_metric)
  
  return(depth_table)
}

# Function to create publication-ready table image using gridExtra
create_publication_table <- function(all_data, output_dir = "rarefaction_plots", file_prefix = "depth_table") {
  
  # Create the data for the table
  table_data <- all_data %>%
    select(technology, sample_type, rankOrder, sequencingDepth) %>%
    distinct() %>%
    arrange(technology, sample_type, rankOrder) %>%
    mutate(
      sample_type_full = case_when(
        sample_type == "SC" ~ "Single-Cell",
        sample_type == "SN" ~ "Single-Nucleus", 
        sample_type == "Bulk" ~ "Bulk",
        TRUE ~ sample_type
      ),
      sequencingDepth_formatted = case_when(
        sample_type == "Bulk" ~ format(sequencingDepth, big.mark = ",", scientific = FALSE),
        TRUE ~ format(round(sequencingDepth, 0), big.mark = ",")
      ),
      tech_sample = paste(technology, sample_type_full, sep = " ")
    ) %>%
    select(rankOrder, tech_sample, sequencingDepth_formatted) %>%
    pivot_wider(names_from = tech_sample, 
                values_from = sequencingDepth_formatted) %>%
    arrange(rankOrder) %>%
    # Ensure consistent column order
    select(rankOrder, 
           `ONT Single-Cell`, `ONT Single-Nucleus`, `ONT Bulk`,
           `PacBio Single-Cell`, `PacBio Single-Nucleus`, `PacBio Bulk`)
  
  # Debug: Check what columns we actually have
  cat("Available columns in table_data:", colnames(table_data), "\n")
  cat("First few rows of table_data:\n")
  print(head(table_data))
  
  # Prepare data for table display - handle missing columns gracefully
  display_table <- table_data
  
  # Get actual column names and create display names
  actual_cols <- colnames(display_table)
  display_names <- c("Rank")
  
  # Add display names for existing columns
  for (col in actual_cols[-1]) { # Skip first column (rankOrder)
    if (grepl("ONT Single-Cell", col)) display_names <- c(display_names, "ONT\nSingle-Cell")
    else if (grepl("ONT Single-Nucleus", col)) display_names <- c(display_names, "ONT\nSingle-Nucleus")
    else if (grepl("ONT Bulk", col)) display_names <- c(display_names, "ONT\nBulk")
    else if (grepl("PacBio Single-Cell", col)) display_names <- c(display_names, "PacBio\nSingle-Cell")
    else if (grepl("PacBio Single-Nucleus", col)) display_names <- c(display_names, "PacBio\nSingle-Nucleus")
    else if (grepl("PacBio Bulk", col)) display_names <- c(display_names, "PacBio\nBulk")
    else display_names <- c(display_names, col)
  }
  
  colnames(display_table) <- display_names
  
  # Create table grob with no row names
  table_grob <- tableGrob(display_table, 
                          rows = NULL,  # This removes row names
                          theme = ttheme_default(
                            core = list(
                              fg_params = list(cex = 0.9),
                              bg_params = list(fill = c("white", "#f8f9fa"), alpha = 0.8)
                            ),
                            colhead = list(
                              fg_params = list(cex = 1.0, fontface = "bold"),
                              bg_params = list(fill = "#e9ecef", alpha = 0.9)
                            )
                          ))
  
  # Create the final layout
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save as PNG
  png_filename <- file.path(output_dir, paste0(file_prefix, "_publication_table.png"))
  png(png_filename, width = 1200, height = 900, res = 150)
  
  grid.newpage()
  
  # Simple layout with just the table
  grid.draw(table_grob)
  
  dev.off()
  
  # Save as PDF
  pdf_filename <- file.path(output_dir, paste0(file_prefix, "_publication_table.pdf"))
  pdf(pdf_filename, width = 12, height = 9)
  
  grid.newpage()
  
  # Simple layout with just the table (same as PNG)
  grid.draw(table_grob)
  
  dev.off()
  
  cat("Publication table saved as:\n")
  cat("- PNG:", png_filename, "\n")
  cat("- PDF:", pdf_filename, "\n")
  
  return(table_grob)
}

# Function to prepare unified data with total reads for all sample types
prepare_unified_data <- function(df_list_input) {
  all_data <- do.call(rbind, df_list_input)
  
  # Create unified total_reads column
  all_data <- all_data %>%
    mutate(
      total_reads_unified = case_when(
        sample_type == "Bulk" ~ sequencingDepth,  # For bulk, sequencing_depth IS total reads
        sample_type %in% c("SC", "SN") & !is.na(total_reads) ~ total_reads,  # For SC/SN, use total_reads if available
        sample_type %in% c("SC", "SN") & is.na(total_reads) ~ sequencingDepth,  # Fallback to sequencingDepth if total_reads is NA
        TRUE ~ sequencingDepth  # Default fallback
      ),
      median_reads_per_cell = case_when(
        sample_type %in% c("SC", "SN") ~ sequencingDepth,  # For SC/SN, sequencing_depth IS median reads per cell
        TRUE ~ NA_real_  # Not applicable for bulk
      )
    )
  
  return(all_data)
}

plotRarefactionCurveUnified <- function(df_list_input,
                                       plt_title = "Unified Rarefaction Curve",
                                       export_plt = FALSE,
                                       file_prefix = "unified_rarefaction",
                                       txt_size = 14,
                                       file_width = 10,
                                       file_height = 6,
                                       output_dir = "rarefaction_plots",
                                       log_scale = FALSE,
                                       show_dual_axis = TRUE) {
  
  # Prepare unified data
  all_data <- prepare_unified_data(df_list_input)
  gene_data <- all_data %>% filter(type == "Genes Discovery" | type == "Gene Discovery")
  iso_data  <- all_data %>% filter(type == "Isoforms Discovery" | type == "Isoform Discovery")

  # Color mapping for technologies (same as original)
  color_map <- c(
    "PacBio" = "#DA1884",
    "ONT" = "#003F5C"
  )
  
  # Line type mapping for sample types (same as original)
  linetype_map <- c(
    "SC" = "dashed",
    "SN" = "dotted", 
    "Bulk" = "solid"
  )
  
  # Point shape mapping for sample types (same as original)
  shape_map <- c(
    "SC" = 16,    # solid circle
    "SN" = 17,    # solid triangle
    "Bulk" = 15   # solid square
  )

  # Create gene plot
  gene_plt <- ggplot(gene_data, aes(x = total_reads_unified, y = featureNum, 
                                   color = technology, linetype = sample_type, 
                                   shape = sample_type, group = sampleID)) +
    geom_line(size = 1, alpha = 0.7) +
    geom_point(size = 3, alpha = 0.8) +
    labs(x = "Total Reads",
         y = "Number of Features Detected",
         title = paste(plt_title, "- Genes"),
         color = "Technology",
         linetype = "Sample Type",
         shape = "Sample Type") +
    theme_minimal() +
    theme(text = element_text(size = txt_size),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right",
          plot.margin = margin(t = 20, r = 5, b = 5, l = 5, unit = "mm")) +  # Extra top margin for dual axis
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_color_manual(values = color_map) +
    scale_linetype_manual(values = linetype_map) +
    scale_shape_manual(values = shape_map)

  # Create isoform plot
  iso_plt <- ggplot(iso_data, aes(x = total_reads_unified, y = featureNum, 
                                 color = technology, linetype = sample_type, 
                                 shape = sample_type, group = sampleID)) +
    geom_line(size = 1, alpha = 0.7) +
    geom_point(size = 3, alpha = 0.8) +
    labs(x = "Total Reads",
         y = "Number of Features Detected",
         title = paste(plt_title, "- Isoforms"),
         color = "Technology",
         linetype = "Sample Type",
         shape = "Sample Type") +
    theme_minimal() +
    theme(text = element_text(size = txt_size),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right",
          plot.margin = margin(t = 20, r = 5, b = 5, l = 5, unit = "mm")) +  # Extra top margin for dual axis
    scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
    scale_color_manual(values = color_map) +
    scale_linetype_manual(values = linetype_map) +
    scale_shape_manual(values = shape_map)

  # Apply log scale if requested
  if (log_scale) {
    gene_plt <- gene_plt + scale_x_log10(
      labels = scales::comma,
      breaks = scales::log_breaks(n = 12)  # More breakpoints for log scale
    )
    iso_plt <- iso_plt + scale_x_log10(
      labels = scales::comma,
      breaks = scales::log_breaks(n = 12)  # More breakpoints for log scale
    )
  } else {
    gene_plt <- gene_plt + scale_x_continuous(
      labels = scales::comma,
      breaks = scales::pretty_breaks(n = 15)  # More breakpoints for linear scale
    )
    iso_plt <- iso_plt + scale_x_continuous(
      labels = scales::comma,
      breaks = scales::pretty_breaks(n = 15)  # More breakpoints for linear scale
    )
  }
  
  # Add secondary x-axis for median reads per cell if requested
  if (show_dual_axis) {
    # Calculate the range of total reads for SC/SN data to create secondary axis
    sc_sn_data <- all_data %>% filter(sample_type %in% c("SC", "SN"), !is.na(median_reads_per_cell))
    
    if (nrow(sc_sn_data) > 0) {
      # Create a mapping between total reads and median reads per cell for SC/SN data
      # We'll use a representative sample to establish the relationship
      sc_sn_summary <- sc_sn_data %>%
        group_by(technology, sample_type) %>%
        summarise(
          avg_ratio = mean(total_reads_unified / median_reads_per_cell, na.rm = TRUE),
          .groups = "drop"
        )
      
      # Use the median ratio across all SC/SN samples for the conversion
      avg_conversion_ratio <- median(sc_sn_summary$avg_ratio, na.rm = TRUE)
      
      # Create secondary axis transformation
      sec_trans <- function(x) x / avg_conversion_ratio
      sec_inv_trans <- function(x) x * avg_conversion_ratio
      
      # Add secondary axis to both plots
      if (log_scale) {
        gene_plt <- gene_plt + 
          scale_x_log10(
            labels = scales::comma,
            breaks = scales::log_breaks(n = 12),
            sec.axis = sec_axis(
              trans = sec_trans,
              name = "Median Reads per Cell (SC/SN equivalent)",
              labels = scales::comma,
              breaks = scales::log_breaks(n = 8)  # Fewer breaks on secondary axis to avoid crowding
            )
          )
        
        iso_plt <- iso_plt + 
          scale_x_log10(
            labels = scales::comma,
            breaks = scales::log_breaks(n = 12),
            sec.axis = sec_axis(
              trans = sec_trans,
              name = "Median Reads per Cell (SC/SN equivalent)",
              labels = scales::comma,
              breaks = scales::log_breaks(n = 8)  # Fewer breaks on secondary axis to avoid crowding
            )
          )
      } else {
        gene_plt <- gene_plt + 
          scale_x_continuous(
            labels = scales::comma,
            breaks = scales::pretty_breaks(n = 15),
            sec.axis = sec_axis(
              trans = sec_trans,
              name = "Median Reads per Cell (SC/SN equivalent)",
              labels = scales::comma,
              breaks = scales::pretty_breaks(n = 10)  # Fewer breaks on secondary axis to avoid crowding
            )
          )
        
        iso_plt <- iso_plt + 
          scale_x_continuous(
            labels = scales::comma,
            breaks = scales::pretty_breaks(n = 15),
            sec.axis = sec_axis(
              trans = sec_trans,
              name = "Median Reads per Cell (SC/SN equivalent)",
              labels = scales::comma,
              breaks = scales::pretty_breaks(n = 10)  # Fewer breaks on secondary axis to avoid crowding
            )
          )
      }
    }
  }

  combined_plt <- gene_plt + iso_plt + plot_layout(ncol = 2, guides = "collect")
  
  # Save plot if requested
  if (export_plt) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Create file suffix based on options
    file_suffix <- "_unified"
    if (show_dual_axis) file_suffix <- paste0(file_suffix, "_dual_axis")
    if (log_scale) file_suffix <- paste0(file_suffix, "_log")
    
    # Increase width for dual axis plots to provide more space
    plot_width <- if (show_dual_axis) file_width * 2.5 else file_width * 2
    plot_height <- if (show_dual_axis) file_height * 1.2 else file_height  # Slightly taller for dual axis
    
    # Save as PNG
    png_filename <- file.path(output_dir, paste0(file_prefix, file_suffix, ".png"))
    ggsave(png_filename, plot = combined_plt, width = plot_width, height = plot_height, dpi = 300)
    
    # Save as PDF with extra width for dual axis
    pdf_filename <- file.path(output_dir, paste0(file_prefix, file_suffix, ".pdf"))
    ggsave(pdf_filename, plot = combined_plt, width = plot_width, height = plot_height, device = "pdf")
    
    # Save as SVG
    svg_filename <- file.path(output_dir, paste0(file_prefix, file_suffix, ".svg"))
    ggsave(svg_filename, plot = combined_plt, width = plot_width, height = plot_height, device = "svg")
    
    # Save the unified data as RDS file
    rds_filename <- file.path(output_dir, paste0(file_prefix, file_suffix, "_data.rds"))
    saveRDS(all_data, rds_filename)
    
    cat("Unified rarefaction plots and data saved as:\n")
    cat("- PNG:", png_filename, "\n")
    cat("- PDF:", pdf_filename, "\n") 
    cat("- SVG:", svg_filename, "\n")
    cat("- RDS data:", rds_filename, "\n")
    cat("- Plot dimensions:", plot_width, "x", plot_height, "\n")
    
    # Print data summary
    cat("\nData summary:\n")
    summary_table <- all_data %>%
      group_by(technology, sample_type) %>%
      summarise(
        n_samples = n(),
        total_reads_range = paste0(scales::comma(min(total_reads_unified, na.rm = TRUE)), 
                                 " - ", scales::comma(max(total_reads_unified, na.rm = TRUE))),
        median_reads_per_cell_range = if(any(!is.na(median_reads_per_cell))) {
          paste0(scales::comma(min(median_reads_per_cell, na.rm = TRUE)), 
                 " - ", scales::comma(max(median_reads_per_cell, na.rm = TRUE)))
        } else {"N/A"},
        .groups = "drop"
      )
    print(summary_table)
  }

  return(list(plot = combined_plt, data = all_data))
}

# Create directories for plots
if (!dir.exists("rarefaction_plots")) {
  dir.create("rarefaction_plots")
}

## Step 4: Create and display rarefaction plots ------
# Combine all feature number data frames
all_feature_nums <- list(
  gene_sc_featureNums_ont,
  tx_sc_featureNums_ont,
  gene_sc_featureNums_pb,
  tx_sc_featureNums_pb,
  gene_sn_featureNums_ont,
  tx_sn_featureNums_ont,
  gene_sn_featureNums_pb,
  tx_sn_featureNums_pb,
  gene_bulk_featureNums_ont,
  gene_bulk_featureNums_pb,
  tx_bulk_featureNums_ont,
  tx_bulk_featureNums_pb
)

# Create combined rarefaction plot
rarefaction_result <- plotRarefactionCurve(
  df_list_input = all_feature_nums,
  plt_title = "Rarefaction Curve: Features Detected by Sequencing Depth Rank",
  export_plt = TRUE,
  file_prefix = "filtered_rarefaction_curve",
  txt_size = 16,
  file_height = 8,
  file_width = 12,
  output_dir = "rarefaction_plots"
)

# Create and display sequencing depth summary table
all_combined_data <- do.call(rbind, all_feature_nums)

# Debug: Check the structure of combined data
cat("\n=== DEBUGGING COMBINED DATA ===\n")
cat("Structure of all_combined_data:\n")
str(all_combined_data)
cat("\nUnique combinations of technology, sample_type, and rankOrder:\n")
unique_combos <- all_combined_data %>%
  select(technology, sample_type, rankOrder, sequencingDepth) %>%
  distinct() %>%
  arrange(technology, sample_type, rankOrder)
print(unique_combos)

depth_summary <- create_depth_summary_table(all_combined_data)

cat("\n=== SEQUENCING DEPTH SUMMARY BY RANK ===\n")
cat("Note: SC/SN depths are median reads per cell, Bulk depths are total reads\n\n")
print(depth_summary)

# Create publication-ready table image
cat("\n=== CREATING PUBLICATION TABLE ===\n")
publication_table <- create_publication_table(
  all_data = all_combined_data,
  output_dir = "rarefaction_plots",
  file_prefix = "filtered_rarefaction_curve"
)

## Step 5: Create unified rarefaction plots with dual x-axis ------
cat("\n=== CREATING UNIFIED RAREFACTION PLOTS WITH DUAL X-AXIS ===\n")

# Create unified rarefaction plot with dual x-axis (linear scale)
unified_result_linear <- plotRarefactionCurveUnified(
  df_list_input = all_feature_nums,
  plt_title = "Unified Rarefaction Curve: All Sample Types by Total Reads",
  export_plt = TRUE,
  file_prefix = "filtered_rarefaction_curve",
  txt_size = 16,
  file_height = 8,
  file_width = 12,
  output_dir = "rarefaction_plots",
  log_scale = FALSE,
  show_dual_axis = TRUE
)

# Create unified rarefaction plot with dual x-axis (log scale)
unified_result_log <- plotRarefactionCurveUnified(
  df_list_input = all_feature_nums,
  plt_title = "Unified Rarefaction Curve: All Sample Types by Total Reads (Log Scale)",
  export_plt = TRUE,
  file_prefix = "filtered_rarefaction_curve",
  txt_size = 16,
  file_height = 8,
  file_width = 12,
  output_dir = "rarefaction_plots",
  log_scale = TRUE,
  show_dual_axis = TRUE
)

# Create unified plot without dual axis for comparison
unified_result_single <- plotRarefactionCurveUnified(
  df_list_input = all_feature_nums,
  plt_title = "Unified Rarefaction Curve: All Sample Types by Total Reads (Single Axis)",
  export_plt = TRUE,
  file_prefix = "filtered_rarefaction_curve",
  txt_size = 16,
  file_height = 8,
  file_width = 12,
  output_dir = "rarefaction_plots",
  log_scale = FALSE,
  show_dual_axis = FALSE
)

cat("\nUnified rarefaction curves generated with:\n")
cat("- Bottom x-axis: Total reads (unified scale for all sample types)\n")
cat("- Top x-axis: Median reads per cell equivalent (for SC/SN reference)\n")
cat("- All sample types (Bulk, SC, SN) plotted on the same graph\n")
cat("- Both linear and log-scale versions\n")
cat("- Data normalization:\n")
cat("  * Bulk: total_reads = sequencing_depth\n")
cat("  * SC/SN: total_reads from CSV, median_reads_per_cell = sequencing_depth\n")
cat("\nAll plots saved in the 'rarefaction_plots' directory\n")

## Step 6: Export data for statistical analysis ------
cat("\n=== PREPARING DATA FOR STATISTICAL ANALYSIS ===\n")

# Create standardized data export for statistical analysis
statistical_data <- prepare_unified_data(all_feature_nums) %>%
  # Add additional metadata for statistical analysis
  mutate(
    # Create unique identifiers for each curve
    curve_id = paste(technology, sample_type, type, sep = "_"),
    # Copy sequencingDepth to total_reads where total_reads is NA (especially for bulk data)
    total_reads = ifelse(is.na(total_reads), sequencingDepth, total_reads),
    # Ensure total_reads_unified is never NA by using sequencingDepth as fallback
    total_reads_unified = ifelse(is.na(total_reads_unified), sequencingDepth, total_reads_unified),
    # Add log-transformed values for easier modeling
    log_total_reads = log10(total_reads_unified),
    log_features = log10(featureNum),
    # Add normalized feature counts (features per million reads)
    features_per_million = (featureNum / total_reads_unified) * 1e6
  ) %>%
  # Add sample size information in a separate step
  mutate(
    n_ranks = n(),
    .by = curve_id
  ) %>%
  # Ensure data is sorted by sequencing depth within each curve
  arrange(curve_id, total_reads_unified)

# Create metadata table for statistical analysis
curve_metadata <- statistical_data %>%
  group_by(curve_id, technology, sample_type, type) %>%
  summarise(
    n_points = n(),
    min_total_reads = min(total_reads_unified, na.rm = TRUE),
    max_total_reads = max(total_reads_unified, na.rm = TRUE),
    min_features = min(featureNum, na.rm = TRUE),
    max_features = max(featureNum, na.rm = TRUE),
    reads_range = max_total_reads - min_total_reads,
    features_range = max_features - min_features,
    .groups = "drop"
  )

# Save data for statistical analysis
stats_dir <- "statistical_analysis"
if (!dir.exists(stats_dir)) {
  dir.create(stats_dir, recursive = TRUE)
}

# Save the main dataset
statistical_data_file <- file.path(stats_dir, "rarefaction_data_for_stats.rds")
saveRDS(statistical_data, statistical_data_file)

# Save metadata
metadata_file <- file.path(stats_dir, "curve_metadata.rds")
saveRDS(curve_metadata, metadata_file)

# Save as CSV for easy inspection
csv_data_file <- file.path(stats_dir, "rarefaction_data_for_stats.csv")
write.csv(statistical_data, csv_data_file, row.names = FALSE)

csv_metadata_file <- file.path(stats_dir, "curve_metadata.csv")
write.csv(curve_metadata, csv_metadata_file, row.names = FALSE)

# Print summary for user
cat("Statistical analysis data exported:\n")
cat("- Main dataset (RDS):", statistical_data_file, "\n")
cat("- Main dataset (CSV):", csv_data_file, "\n")
cat("- Curve metadata (RDS):", metadata_file, "\n")
cat("- Curve metadata (CSV):", csv_metadata_file, "\n")

cat("\nDataset summary:\n")
cat("- Total data points:", nrow(statistical_data), "\n")
cat("- Unique curves:", length(unique(statistical_data$curve_id)), "\n")
cat("- Technologies:", paste(unique(statistical_data$technology), collapse = ", "), "\n")
cat("- Sample types:", paste(unique(statistical_data$sample_type), collapse = ", "), "\n")
cat("- Feature types:", paste(unique(statistical_data$type), collapse = ", "), "\n")

print(curve_metadata)

cat("\n=== RAREFACTION CURVE GENERATION COMPLETE ===\n")
cat("Next step: Run 'analyze_rarefaction_statistics.R' for statistical modeling\n") 
