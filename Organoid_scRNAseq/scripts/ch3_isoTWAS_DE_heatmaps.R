library(Seurat)
library(tidyverse)
library(pheatmap)
library(presto)
library(gridExtra)
library(grid)
library(patchwork)
library(RColorBrewer)

setwd("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space")

# Import consensus annotated seurat objects with isoform assay
one_month_org <- readRDS("./output_files/integrated_objects/1M_Org_integrated_harmony_consensus_with_isoform_assay.rds")
three_month_org <- readRDS("./output_files/integrated_objects/3M_Org_integrated_harmony_consensus_with_isoform_assay.rds")
six_month_org <- readRDS("./output_files/integrated_objects/6M_Org_integrated_harmony_consensus_with_isoform_assay.rds")

Idents(one_month_org) <- "consensus_cell_type"
Idents(three_month_org) <- "consensus_cell_type"
Idents(six_month_org) <- "consensus_cell_type"

# Define the trait-specific filtered database paths
trait_files <- list(
  #ASD = "../Updated_Mitch_Filtered_Database/ASD_filtered_5e-8_Sig_isoTWAS_developmental.csv",
  BP = "../Updated_Mitch_Filtered_Database/BP_filtered_5e-8_Sig_isoTWAS_developmental.csv",
  MDD = "../Updated_Mitch_Filtered_Database/MDD_filtered_5e-8_Sig_isoTWAS_developmental.csv",
  SCZ = "../Updated_Mitch_Filtered_Database/SCZ_filtered_5e-8_Sig_isoTWAS_developmental.csv"
)

# Read in all trait-specific TWAS data and combine
cat("Loading trait-specific GWAS-filtered TWAS databases...\n")
isoform_TWAS_list <- list()
for(trait_name in names(trait_files)) {
  file_path <- trait_files[[trait_name]]
  if(file.exists(file_path)) {
    trait_data <- read.csv(file_path, stringsAsFactors = FALSE)
    cat("Loaded", nrow(trait_data), "records for", trait_name, "\n")
    isoform_TWAS_list[[trait_name]] <- trait_data
  } else {
    cat("Warning: File not found:", file_path, "\n")
  }
}

# Combine all trait-specific data into one dataframe
isoform_TWAS_df <- bind_rows(isoform_TWAS_list)
cat("Combined dataset contains", nrow(isoform_TWAS_df), "total records across", 
    length(unique(isoform_TWAS_df$Trait)), "traits\n")

# Create the by-condition lookup
isoform_TWAS_byCondition <- isoform_TWAS_df %>% 
  dplyr::select(c("Transcript", "Trait"))

# Get unique traits
traits <- unique(isoform_TWAS_df$Trait)
cat("Traits to analyze:", paste(traits, collapse = ", "), "\n")

# Function that returns a dataframe with DETs overlapping with TWAS dataframe
returnTWAS_DE_Features <- function(seu.obj){
  # calculate DEGs
  cat("Running FindAllMarkers for isoform-level analysis...\n")
  cellTypeMarkers_sc <- FindAllMarkers(seu.obj,
                                       assay = "iso",
                                       layer = "data",  # Updated for Seurat v5 compatibility
                                       test.use = "wilcox")

  cat("Found", nrow(cellTypeMarkers_sc), "total markers\n")

  # Verbose filtering with step-by-step logging
  cat("\n--- Filtering Pipeline ---\n")

  # Step 1: p-value filter
  step1 <- cellTypeMarkers_sc %>% dplyr::filter(p_val_adj <= 0.05)
  cat("After p_val_adj <= 0.05 filter:", nrow(step1), "markers\n")

  # Step 2: log2FC filter REMOVED for TWAS integration
  # Justification: TWAS focuses on statistical associations, not effect sizes
  # Small fold changes can be biologically meaningful in psychiatric genetics
  step2 <- step1  # No filtering on log2FC
  cat("No log2FC filter applied - using all statistically significant markers:", nrow(step2), "markers\n")

  # Step 3: Bambu filter
  step3 <- step2 %>% dplyr::filter(!str_detect(gene, "Bambu"))
  cat("After Bambu exclusion:", nrow(step3), "markers\n")

  # Step 4: Add trimming and show diagnostic info
  step4 <- step3 %>% dplyr::mutate(trimmed_TxID = sub("\\..*", "", gene))

  cat("\n--- Diagnostic Info ---\n")
  cat("Sample feature names (first 5):\n")
  print(head(step4$gene, 5))
  cat("\nSample trimmed transcript IDs (first 5):\n")
  print(head(step4$trimmed_TxID, 5))
  cat("\nSample TWAS transcript IDs (first 5):\n")
  print(head(isoform_TWAS_df$Transcript, 5))

  # Step 5: Match with TWAS database
  DE_iso_filtered <- step4 %>%
    dplyr::filter(trimmed_TxID %in% isoform_TWAS_df$Transcript) %>%
    left_join(isoform_TWAS_byCondition, by = c("trimmed_TxID" = "Transcript"))

  cat("\n--- Final Results ---\n")
  cat("After TWAS matching:", nrow(DE_iso_filtered), "DETs overlap with TWAS transcripts\n")

  if(nrow(step4) > 0) {
    overlap_rate <- round(100 * nrow(DE_iso_filtered) / nrow(step4), 1)
    cat("Overlap rate:", overlap_rate, "%\n")
  }

  # Show trait distribution if matches found
  if(nrow(DE_iso_filtered) > 0) {
    cat("\nTrait distribution of matches:\n")
    print(table(DE_iso_filtered$Trait))
  }

  cat("-------------------------\n\n")

  return(DE_iso_filtered)
}

# Create dfs of differentially expressed isoforms that overlap with TWAS transcripts----
cat("\nProcessing 6-month organoids...\n")
TWAS_DE_6M_Org <- returnTWAS_DE_Features(six_month_org)

cat("\nProcessing 3-month organoids...\n")
TWAS_DE_3M_Org <- returnTWAS_DE_Features(three_month_org)

cat("\nProcessing 1-month organoids...\n")
TWAS_DE_1M_Org <- returnTWAS_DE_Features(one_month_org)

# split by TWAS trait
# FIXED: Use base R split() to avoid trait label mismatch between group_split() and setNames()
splitDataframeByTrait <- function(TWAS_DE_df){
  # split() automatically names list elements correctly by the grouping factor
  split(TWAS_DE_df, TWAS_DE_df$Trait)
}

TWAS_DE_6M_Org.split <- splitDataframeByTrait(TWAS_DE_6M_Org)
TWAS_DE_3M_Org.split <- splitDataframeByTrait(TWAS_DE_3M_Org)
TWAS_DE_1M_Org.split <- splitDataframeByTrait(TWAS_DE_1M_Org)

cat("Data split by trait. Available traits in results:\n")
cat("6M:", names(TWAS_DE_6M_Org.split), "\n")
cat("3M:", names(TWAS_DE_3M_Org.split), "\n")
cat("1M:", names(TWAS_DE_1M_Org.split), "\n")

# Function that collapses the TWAS DE dataframe to gene level list
CollapseIsoDEdf2Gene <- function(TWAS_DE_df){
  IsoGene_TWAS_dict <- isoform_TWAS_df %>%
    dplyr::select(c("Transcript", "HGNC"))

  TWAS_DE_df <- TWAS_DE_df %>%
    left_join(IsoGene_TWAS_dict, by = c("trimmed_TxID" = "Transcript"))

  gene_lvl_DE_list <- unique(TWAS_DE_df$HGNC)

  return(gene_lvl_DE_list)
}

# Function to prepare summary data for a specific trait across all timepoints
prepare_trait_summary_data <- function(trait_name,
                                       split_1M,
                                       split_3M,
                                       split_6M,
                                       total_gwas_records) {

  # Extract data for this trait from each timepoint
  data_1M <- if(trait_name %in% names(split_1M)) split_1M[[trait_name]] %>% mutate(Timepoint = "1M") else tibble()
  data_3M <- if(trait_name %in% names(split_3M)) split_3M[[trait_name]] %>% mutate(Timepoint = "3M") else tibble()
  data_6M <- if(trait_name %in% names(split_6M)) split_6M[[trait_name]] %>% mutate(Timepoint = "6M") else tibble()

  # Combine all timepoint data
  all_data <- bind_rows(data_1M, data_3M, data_6M)

  # Calculate timepoint totals
  timepoint_totals <- tibble(
    Timepoint = factor(c("1M", "3M", "6M"), levels = c("1M", "3M", "6M")),
    Hits = c(nrow(data_1M), nrow(data_3M), nrow(data_6M))
  )

  # Calculate cell type breakdown with explicit zeros for missing combinations
  if(nrow(all_data) > 0) {
    # Get all unique cell types across all timepoints
    all_celltypes <- unique(all_data$cluster)

    celltype_breakdown <- all_data %>%
      group_by(Timepoint, cluster) %>%
      summarise(Hits = n(), .groups = "drop") %>%
      # Complete to add explicit zeros for missing combinations
      complete(Timepoint = c("1M", "3M", "6M"),
               cluster = all_celltypes,
               fill = list(Hits = 0)) %>%
      rename(CellType = cluster) %>%
      mutate(Timepoint = factor(Timepoint, levels = c("1M", "3M", "6M")))
  } else {
    # No hits at all - create empty structure
    celltype_breakdown <- tibble(
      Timepoint = factor(character(), levels = c("1M", "3M", "6M")),
      CellType = character(),
      Hits = numeric()
    )
  }

  # Calculate summary statistics
  total_de_hits <- nrow(all_data)
  hit_rate_pct <- if(total_gwas_records > 0) {
    round(100 * total_de_hits / total_gwas_records, 1)
  } else {
    0
  }

  stats <- list(
    total_gwas_records = total_gwas_records,
    total_de_hits = total_de_hits,
    hit_rate_pct = hit_rate_pct,
    n_celltypes_1M = length(unique(data_1M$cluster)),
    n_celltypes_3M = length(unique(data_3M$cluster)),
    n_celltypes_6M = length(unique(data_6M$cluster))
  )

  return(list(
    timepoint_totals = timepoint_totals,
    celltype_breakdown = celltype_breakdown,
    stats = stats
  ))
}

# Function to create trait summary visualization page
create_trait_summary_page <- function(trait_name,
                                      split_1M,
                                      split_3M,
                                      split_6M,
                                      total_gwas_records) {

  cat("  Creating summary page for", trait_name, "...\n")

  # Prepare data
  summary_data <- prepare_trait_summary_data(trait_name, split_1M, split_3M,
                                              split_6M, total_gwas_records)

  cat("    ✓ Total hits:", summary_data$stats$total_de_hits, "\n")
  cat("    ✓ Hit rate:", paste0(summary_data$stats$hit_rate_pct, "%"), "\n")
  cat("    ✓ Unique cell types:",
      max(summary_data$stats$n_celltypes_1M,
          summary_data$stats$n_celltypes_3M,
          summary_data$stats$n_celltypes_6M), "\n")

  # Panel A: Timepoint progression bar chart
  panel_A <- ggplot(summary_data$timepoint_totals,
                    aes(x = Timepoint, y = Hits, fill = Timepoint)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = paste0("n=", Hits)),
              vjust = -0.5, size = 5, fontface = "bold") +
    scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Blues")) +
    labs(title = paste0(trait_name, " isoTWAS-DE Hits Across Developmental Timepoints"),
         x = "Developmental Timepoint",
         y = "Number of isoTWAS-DE Hits") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none",
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 10)) +
    ylim(0, max(summary_data$timepoint_totals$Hits) * 1.15)

  # Panel B: Summary statistics table
  stats_df <- data.frame(
    Metric = c("Total GWAS-filtered records",
               "Total DE hits (all timepoints)",
               "Hit rate",
               "",
               "1M hits",
               "3M hits",
               "6M hits",
               "",
               "Unique cell types (1M)",
               "Unique cell types (3M)",
               "Unique cell types (6M)"),
    Value = c(as.character(summary_data$stats$total_gwas_records),
              as.character(summary_data$stats$total_de_hits),
              paste0(summary_data$stats$hit_rate_pct, "%"),
              "",
              as.character(summary_data$timepoint_totals$Hits[1]),
              as.character(summary_data$timepoint_totals$Hits[2]),
              as.character(summary_data$timepoint_totals$Hits[3]),
              "",
              as.character(summary_data$stats$n_celltypes_1M),
              as.character(summary_data$stats$n_celltypes_3M),
              as.character(summary_data$stats$n_celltypes_6M))
  )

  panel_B_grob <- gridExtra::tableGrob(stats_df,
                                       rows = NULL,
                                       theme = gridExtra::ttheme_minimal(
                                         base_size = 10,
                                         padding = unit(c(4, 4), "mm")))

  # Panel C: Cell type breakdown heatmap
  if(nrow(summary_data$celltype_breakdown) > 0) {
    panel_C <- ggplot(summary_data$celltype_breakdown,
                      aes(x = Timepoint, y = CellType, fill = Hits)) +
      geom_tile(color = "white", size = 0.5) +
      geom_text(aes(label = Hits), size = 3, fontface = "bold") +
      scale_fill_viridis_c(option = "inferno",
                           direction = 1,
                           begin = 0.1,
                           end = 0.9) +
      labs(title = "Cell Type Distribution of Hits",
           x = "Developmental Timepoint",
           y = "Cell Type",
           fill = "Hits") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 12, face = "bold"),
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 9),
            legend.position = "right")
  } else {
    # Create empty placeholder if no data
    panel_C <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = "No isoTWAS-DE hits for this trait",
               size = 6) +
      theme_void()
  }

  # Combine panels using patchwork
  library(patchwork)
  combined_plot <- (panel_A | wrap_elements(panel_B_grob)) / panel_C +
    plot_layout(heights = c(1, 1.5)) +
    plot_annotation(
      title = paste0(trait_name, " - isoTWAS-DE Summary"),
      theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
    )

  # Print to current device
  print(combined_plot)

  cat("    ✓ Summary page created\n")
}

# Create heatmaps for each condition, celltype and timepoint----
generate_DE_heatmap <- function(seu.obj,
                                DE_df,
                                sampleID,
                                assay = "iso",
                                trait_name = NULL){
    if(assay == "RNA") {
      featureList = CollapseIsoDEdf2Gene(DE_df)
    } else {
      featureList = DE_df$gene
    }

    cat("Generating", assay, "heatmap for", sampleID, "with",
        length(featureList), "features\n")

    # Extract pseudobulked expression values for each consensus cell type
    # Use AggregateExpression (Seurat v5 recommended) with explicit group.by
    pseudobulk.Expression.list <- AggregateExpression(seu.obj,
                                                       features = featureList,
                                                       assays = assay,
                                                       group.by = "consensus_cell_type",
                                                       return.seurat = FALSE)

    # Extract the matrix for the specific assay
    pseudobulk.Expression.mat <- pseudobulk.Expression.list[[assay]] %>%
      as.data.frame()

    cat("  Pseudobulk matrix dimensions:", nrow(pseudobulk.Expression.mat), "features x",
        ncol(pseudobulk.Expression.mat), "cell types\n")

    # Check if features were dropped during aggregation
    n_dropped <- length(featureList) - nrow(pseudobulk.Expression.mat)
    if(n_dropped > 0) {
      cat("  ⚠ Warning:", n_dropped, "features dropped during aggregation (likely all zeros)\n")
    }

    # VALIDATION: Check minimum feature count
    if(nrow(pseudobulk.Expression.mat) < 3) {
      warning("Too few features (", nrow(pseudobulk.Expression.mat),
              ") for meaningful heatmap. Minimum 3 required. Skipping plot.")
      return(NULL)
    }

    # VALIDATION: Check for variance
    row_vars <- apply(pseudobulk.Expression.mat, 1, var, na.rm = TRUE)
    n_variable <- sum(row_vars > 0, na.rm = TRUE)
    cat("  Features with variance:", n_variable, "/", nrow(pseudobulk.Expression.mat), "\n")

    if(n_variable < 2) {
      warning("Insufficient variance in features for clustering (",
              n_variable, " variable features). Skipping plot.")
      return(NULL)
    }

    cat("  Cell types (x-axis labels):", paste(colnames(pseudobulk.Expression.mat), collapse = ", "), "\n")

  # Wrap pheatmap in tryCatch for robust error handling
  DEG_heatmap <- tryCatch({
    pheatmap::pheatmap(pseudobulk.Expression.mat,
                       cluster_rows = TRUE,
                       show_rownames = TRUE,
                       color = viridis::inferno(100),
                       breaks = seq(-3, 3, length.out = 101),
                       border_color = NA,
                       fontsize = 10,
                       scale = "row",
                       fontsize_row = 6,
                       fontsize_col = 9,
                       height = 20,
                       angle_col = 45,
                       main = paste0(sampleID, " - isoTWAS-DE ", trait_name, " GWAS sig Overlap\n(Row-scaled Z-score)"))
  }, error = function(e) {
    warning("pheatmap failed for ", sampleID, ": ", e$message)
    cat("  ✗ Heatmap generation failed - error:", e$message, "\n")
    return(NULL)
  }, warning = function(w) {
    # Suppress warnings but continue
    suppressWarnings(
      pheatmap::pheatmap(pseudobulk.Expression.mat,
                         cluster_rows = TRUE,
                         show_rownames = TRUE,
                         color = viridis::inferno(100),
                         breaks = seq(-3, 3, length.out = 101),
                         border_color = NA,
                         fontsize = 10,
                         scale = "row",
                         fontsize_row = 6,
                         fontsize_col = 9,
                         height = 20,
                         angle_col = 45,
                         main = paste0(sampleID, " - isoTWAS-DE ", trait_name, " GWAS sig Overlap\n(Row-scaled Z-score)"))
    )
  })

  if(is.null(DEG_heatmap)) {
    cat("  ✗ Heatmap generation returned NULL\n")
  } else {
    cat("  ✓ Heatmap generated successfully\n")
  }

  return(DEG_heatmap)
}

# Export transcript lists for each trait-timepoint combination----
export_transcript_lists <- function(split_list,
                                    timepoint_label,
                                    trait_name,
                                    output_dir = "./output_files/isoTWAS_transcript_lists") {

  # Check if trait exists in split list
  if (!(trait_name %in% names(split_list))) {
    cat("  No data for", trait_name, "at", timepoint_label, "\n")
    return(NULL)
  }

  # Get the dataframe for this trait
  de_df <- split_list[[trait_name]]

  # Check if dataframe has data
  if (is.null(de_df) || nrow(de_df) == 0) {
    cat("  Empty dataframe for", trait_name, "at", timepoint_label, "\n")
    return(NULL)
  }

  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }

  # Prepare gene symbol lookup from TWAS database
  gene_symbol_lookup <- isoform_TWAS_df %>%
    dplyr::select(Transcript, HGNC) %>%
    dplyr::distinct()

  # Create export dataframe with selected columns
  export_df <- de_df %>%
    # Join with TWAS database to get gene symbols
    dplyr::left_join(gene_symbol_lookup, by = c("trimmed_TxID" = "Transcript")) %>%
    # Select and rename columns for clarity
    dplyr::select(
      transcript_id = trimmed_TxID,
      gene_symbol = HGNC,
      cell_type_of_origin = cluster,
      trait = Trait,
      p_val_adj,
      avg_log2FC
    ) %>%
    # Add timepoint column
    dplyr::mutate(timepoint = timepoint_label) %>%
    # Reorder columns for better readability
    dplyr::select(transcript_id, gene_symbol, cell_type_of_origin, trait, timepoint, p_val_adj, avg_log2FC) %>%
    # Arrange by adjusted p-value for easier review
    dplyr::arrange(p_val_adj)

  # Create filename
  output_file <- file.path(output_dir, paste0(trait_name, "_", timepoint_label, "_transcript_list.csv"))

  # Write CSV
  write.csv(export_df, output_file, row.names = FALSE, quote = FALSE)

  cat("  ✓ Exported", nrow(export_df), "transcripts to:", output_file, "\n")

  return(export_df)
}

# Sample heatmap generation removed - heatmaps are now generated within createPDFSummary()
# Each trait gets its own PDF with timepoint-specific heatmaps on separate pages

# Create PDF summary containing heatmaps by TWAS trait (one PDF per trait, separate pages per timepoint)
createPDFSummary <- function(isoLvl = TRUE,
                             output_dir = "./output_files/isoTWAS_heatmaps"){

  # Create output directory if it doesn't exist
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }

  for(trait_name in traits) {
    cat("\n========================================\n")
    cat("Processing trait:", trait_name, "\n")
    cat("========================================\n")

    # Check if trait exists in any timepoint data
    has_6M <- trait_name %in% names(TWAS_DE_6M_Org.split)
    has_3M <- trait_name %in% names(TWAS_DE_3M_Org.split)
    has_1M <- trait_name %in% names(TWAS_DE_1M_Org.split)

    cat("  Data availability - 6M:", has_6M, "3M:", has_3M, "1M:", has_1M, "\n")

    if(!has_6M && !has_3M && !has_1M) {
      cat("  Skipping", trait_name, "- no data available\n")
      next
    }

    # Determine output file name (PDF instead of PNG)
    if(isoLvl) {
      assayType = "iso"
      fileName = file.path(output_dir, paste0(trait_name, "_GWAS_filtered_ISO.pdf"))
    } else {
      assayType = "RNA"
      fileName = file.path(output_dir, paste0(trait_name, "_GWAS_filtered_GENE.pdf"))
    }

    # Open PDF device (each heatmap will be on a separate page)
    pdf(fileName, width = 14, height = 10)

    # NEW: Generate summary page as Page 1
    total_gwas_for_trait <- sum(isoform_TWAS_df$Trait == trait_name)
    create_trait_summary_page(
      trait_name = trait_name,
      split_1M = TWAS_DE_1M_Org.split,
      split_3M = TWAS_DE_3M_Org.split,
      split_6M = TWAS_DE_6M_Org.split,
      total_gwas_records = total_gwas_for_trait
    )

    heatmap_count <- 1  # Summary page is page 1

    # Generate heatmaps for available timepoints (each on separate page)
    if(has_6M && nrow(TWAS_DE_6M_Org.split[[trait_name]]) > 0) {
      cat("  Plotting 6M_Org heatmap (page", heatmap_count + 1, ")...\n")
      result_6M <- generate_DE_heatmap(seu.obj = six_month_org,
                                       DE_df = TWAS_DE_6M_Org.split[[trait_name]],
                                       sampleID = "6M_Org",
                                       assay = assayType,
                                       trait_name = trait_name)
      if(!is.null(result_6M)) {
        heatmap_count <- heatmap_count + 1
        # Export transcript list for this heatmap
        export_transcript_lists(TWAS_DE_6M_Org.split, "6M", trait_name)
      } else {
        cat("  ⚠ Skipped 6M_Org - plotting failed\n")
      }
    }

    if(has_3M && nrow(TWAS_DE_3M_Org.split[[trait_name]]) > 0) {
      cat("  Plotting 3M_Org heatmap (page", heatmap_count + 1, ")...\n")
      result_3M <- generate_DE_heatmap(seu.obj = three_month_org,
                                       DE_df = TWAS_DE_3M_Org.split[[trait_name]],
                                       sampleID = "3M_Org",
                                       assay = assayType,
                                       trait_name = trait_name)
      if(!is.null(result_3M)) {
        heatmap_count <- heatmap_count + 1
        # Export transcript list for this heatmap
        export_transcript_lists(TWAS_DE_3M_Org.split, "3M", trait_name)
      } else {
        cat("  ⚠ Skipped 3M_Org - plotting failed\n")
      }
    }

    if(has_1M && nrow(TWAS_DE_1M_Org.split[[trait_name]]) > 0) {
      cat("  Plotting 1M_Org heatmap (page", heatmap_count + 1, ")...\n")
      result_1M <- generate_DE_heatmap(seu.obj = one_month_org,
                                       DE_df = TWAS_DE_1M_Org.split[[trait_name]],
                                       sampleID = "1M_Org",
                                       assay = assayType,
                                       trait_name = trait_name)
      if(!is.null(result_1M)) {
        heatmap_count <- heatmap_count + 1
        # Export transcript list for this heatmap
        export_transcript_lists(TWAS_DE_1M_Org.split, "1M", trait_name)
      } else {
        cat("  ⚠ Skipped 1M_Org - plotting failed\n")
      }
    }

    # Close PDF device
    dev.off()

    if(heatmap_count > 0) {
      cat("  ✓ Saved", heatmap_count, "heatmap page(s) to:", fileName, "\n")
    } else {
      # Remove empty PDF if no heatmaps were created
      if(file.exists(fileName)) {
        file.remove(fileName)
      }
      cat("  ✗ No valid heatmaps generated for", trait_name, "- PDF removed\n")
    }
  }

  cat("\n========================================\n")
  cat("PDF generation complete!\n")
  cat("========================================\n")
}

# Create trait-specific heatmaps using GWAS-filtered data
cat("\nCreating GWAS-filtered trait-specific heatmaps...\n")


cat("\nGenerating ISOFORM-level heatmaps...\n")
createPDFSummary(isoLvl = TRUE)

cat("\n========================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("========================================\n")
cat("GWAS-filtered heatmaps saved to: ./output_files/isoTWAS_heatmaps/\n")
cat("Each trait has its own PDF with timepoint-specific heatmaps on separate pages.\n")

# Print summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
cat("Input data summary:\n")
for(trait_name in traits) {
  trait_count <- sum(isoform_TWAS_df$Trait == trait_name)
  cat("  ", trait_name, ":", trait_count, 
      "GWAS-filtered records\n")
}

cat("\nDifferential expression overlaps:\n")
cat("6M organoids:", nrow(TWAS_DE_6M_Org),
    "TWAS-DE overlaps\n")
cat("3M organoids:", nrow(TWAS_DE_3M_Org),
    "TWAS-DE overlaps\n")
cat("1M organoids:", nrow(TWAS_DE_1M_Org),
    "TWAS-DE overlaps\n")

# ===== STANDALONE TRANSCRIPT LIST EXPORT =====
# This section manually exports transcript lists without re-running the full analysis.
# Uncomment the code below if you need to re-export CSVs:

# cat("\n========================================\n")
# cat("MANUAL TRANSCRIPT LIST EXPORT\n")
# cat("========================================\n")
#
# # Export all trait-timepoint combinations
# for (trait_name in traits) {
#   cat("\nExporting transcript lists for", trait_name, "...\n")
#
#   # Export 6M
#   export_transcript_lists(TWAS_DE_6M_Org.split, "6M", trait_name)
#
#   # Export 3M
#   export_transcript_lists(TWAS_DE_3M_Org.split, "3M", trait_name)
#
#   # Export 1M
#   export_transcript_lists(TWAS_DE_1M_Org.split, "1M", trait_name)
# }
#
# cat("\n========================================\n")
# cat("MANUAL EXPORT COMPLETE!\n")
# cat("========================================\n")

#################

# ===== Write "isovis" CSVs with columns: gene_id, transcript_id, <CellType...> =====
# Values are pseudobulk expression (AverageExpression on 'iso') of isoforms that are DE & TWAS-overlapping.

# Output directory
isovis_outdir <- "./output_files/isovis"
if (!dir.exists(isovis_outdir)) dir.create(isovis_outdir, recursive = TRUE)

# --- helpers to parse feature names like "ENST00000335137.4-GENE1" ---
extract_txid_full <- function(x) sub("-.*$", "", x)                  # keep version; drop -GENE
extract_parent_symbol <- function(x) ifelse(grepl("-", x),
                                            sub(".*-", "", x),
                                            NA_character_)
# try to discover a gene_id column in your TWAS table (optional)
detect_gene_id_col <- function(df) {
  candidates <- c("gene_id", "EnsemblGene", "ENSG", "Ensembl_GeneID", "GeneID", "Ensembl", "Gene") # last resort
  intersect(candidates, names(df))[1]
}
twas_gene_id_col <- detect_gene_id_col(isoform_TWAS_df)  # may be NULL

# Build one timepoint isovis data.frame for a trait (expression, not logFC)
build_isovis_expr_one_timepoint <- function(seu.obj, twas_de_split_list, trait_name, timepoint_label) {
  if (!(trait_name %in% names(twas_de_split_list))) return(NULL)
  de_df <- twas_de_split_list[[trait_name]]
  if (is.null(de_df) || nrow(de_df) == 0) return(NULL)
  
  # feature names present in the iso assay (e.g., "ENST...-GENE")
  feat_names <- unique(de_df$gene)

  # Average expression per cluster (pseudobulk) - using AggregateExpression for consistency
  # This ensures isovis CSVs have same cell type columns as heatmaps
  avg_list <- AggregateExpression(seu.obj,
                                   features = feat_names,
                                   assays = "iso",
                                   group.by = "consensus_cell_type",
                                   return.seurat = FALSE)
  avg <- as.data.frame(avg_list$iso)
  # Column names are already correct consensus cell types (no prefix to remove)
  
  # add transcript_id (full, with version) and gene_id (from TWAS if available, else parent symbol)
  df <- avg %>%
    tibble::rownames_to_column(var = "feature_name") %>%
    dplyr::mutate(
      transcript_id = extract_txid_full(feature_name),
      parent_symbol = extract_parent_symbol(feature_name)
    )
  
  # try to map to Ensembl gene_id using TWAS table (by transcript)
  # We normalize TWAS Transcript to transcript without version for a tolerant join;
  # if you have exact full IDs in TWAS, this still works because we coalesce.
  twas_map <- isoform_TWAS_df %>%
    dplyr::mutate(
      Transcript_full = Transcript,                              # as provided
      Transcript_nov  = sub("\\..*$", "", Transcript)            # trimmed version
    ) %>%
    {
      if (!is.null(twas_gene_id_col)) {
        dplyr::select(., !!twas_gene_id_col, Transcript_full, Transcript_nov, HGNC)
      } else {
        dplyr::select(., Transcript_full, Transcript_nov, HGNC)
      }
    } %>%
    dplyr::distinct()
  
  df <- df %>%
    dplyr::mutate(transcript_nov = sub("\\..*$", "", transcript_id)) %>%
    # left join twice to increase chances: full match first, then nov
    dplyr::left_join(twas_map, by = c("transcript_id" = "Transcript_full")) %>%
    dplyr::left_join(twas_map, by = c("transcript_nov" = "Transcript_nov"), suffix = c("", ".nov")) %>%
    dplyr::mutate(
      # pick gene_id: prefer Ensembl if available in TWAS; else fall back to parent symbol
      gene_id = dplyr::coalesce(!!rlang::sym(twas_gene_id_col %||% ""), !!rlang::sym(paste0(twas_gene_id_col, ".nov") %||% ""), parent_symbol)
    )
  
  # Final tidy table:
  # gene_id, transcript_id, then only cell-type columns (no timepoint prefix in headers)
  celltype_cols <- setdiff(names(df), c("feature_name","transcript_id","parent_symbol",
                                        "transcript_nov","HGNC","HGNC.nov","Transcript_full",
                                        "Transcript_nov","gene_id",
                                        if (!is.null(twas_gene_id_col)) c(twas_gene_id_col, paste0(twas_gene_id_col,".nov")) else character(0)))
  out <- df %>%
    dplyr::select(gene_id, transcript_id, dplyr::all_of(celltype_cols))
  
  # Optional: add a tiny epsilon to avoid exact zeros (if needed for log-scaling downstream)
  # out[celltype_cols] <- lapply(out[celltype_cols], function(v) ifelse(is.na(v), 0, v))
  
  attr(out, "timepoint_label") <- timepoint_label
  out
}

# Write per-trait, per-timepoint CSVs (columns: gene_id, transcript_id, <CellTypes...>)
write_isovis_for_trait <- function(trait_name) {
  cat("Creating expression isovis CSVs for", trait_name, "...\n")
  
  tab_6M <- build_isovis_expr_one_timepoint(six_month_org,   TWAS_DE_6M_Org.split, trait_name, "6M")
  tab_3M <- build_isovis_expr_one_timepoint(three_month_org, TWAS_DE_3M_Org.split, trait_name, "3M")
  tab_1M <- build_isovis_expr_one_timepoint(one_month_org,   TWAS_DE_1M_Org.split, trait_name, "1M")
  
  lst <- list(`6M`=tab_6M, `3M`=tab_3M, `1M`=tab_1M)
  for (nm in names(lst)) {
    tb <- lst[[nm]]
    if (is.null(tb) || nrow(tb) == 0) next
    # ensure numeric columns are numeric and NAs -> 0 (optional)
    ct_cols <- setdiff(names(tb), c("gene_id","transcript_id"))
    tb[ct_cols] <- lapply(tb[ct_cols], function(v) { v[is.na(v)] <- 0; as.numeric(v) })
    
    outfile <- file.path(isovis_outdir, paste0(trait_name, "_isovis_", nm, ".csv"))
    readr::write_csv(tb, outfile)
    cat("  Saved:", outfile, "with", nrow(tb), "isoforms and", length(ct_cols), "cell types.\n")
  }
  
  # Optional: also write an aggregated (across timepoints) file by taking the max per cell type
  if (!all(vapply(lst, is.null, logical(1)))) {
    # harmonize cell-type columns across timepoints by name, then take rowwise max
    # first, bind rows with a 'timepoint' tag, then aggregate
    bind_and_agg <- function(xlist) {
      xlist <- xlist[!vapply(xlist, is.null, logical(1))]
      if (length(xlist) == 0) return(NULL)
      # full join on gene_id + transcript_id so we don't lose isoforms unique to a TP
      out <- Reduce(function(a,b){
        dplyr::full_join(a, b, by = c("gene_id","transcript_id"), suffix = c(".a",".b"))
      }, xlist)
      
      # find cell-type columns (any non id columns)
      ct_cols <- setdiff(names(out), c("gene_id","transcript_id"))
      # If the same cell type appears multiple times (from different TPs), take rowwise max
      # Collapse by base name (strip any ".a"/".b" suffixes introduced by joins)
      base <- sub("\\.[ab]$", "", ct_cols)
      agg <- out[, c("gene_id","transcript_id"), drop = FALSE]
      for (nm in unique(base)) {
        cols <- ct_cols[base == nm]
        mat  <- as.matrix(out[, cols, drop = FALSE])
        agg[[nm]] <- apply(mat, 1, function(v) suppressWarnings(max(as.numeric(v), na.rm = TRUE)))
        agg[[nm]][!is.finite(agg[[nm]])] <- 0
      }
      agg
    }
    
    agg_tab <- bind_and_agg(lst)
    if (!is.null(agg_tab)) {
      outfile <- file.path(isovis_outdir, paste0(trait_name, "_isovis_aggregated.csv"))
      readr::write_csv(agg_tab, outfile)
      cat("  Saved aggregated:", outfile, "with", nrow(agg_tab), "isoforms and",
          ncol(agg_tab) - 2, "cell types.\n")
    }
  }
}

# Only run for traits present, and for focus on BP, MDD, and SCZ
traits_of_interest <- intersect(c("BP", "MDD", "SCZ"), traits)
for (tr in traits_of_interest) write_isovis_for_trait(tr)

cat("Done. Isovis CSVs are in:", isovis_outdir, "\n")

