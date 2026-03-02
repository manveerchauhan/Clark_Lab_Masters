# ============================================================================
# Canonical vs Non-Canonical Transcript Analysis for isoTWAS-DE Hits
# Author: Manveer Chauhan
# Date: 2025-01-28
# Description: Analyzes the distribution of canonical vs non-canonical
#              transcripts in TWAS-DE hits across developmental timepoints
# ============================================================================

# REQUIRED LIBRARIES --------------------------------------------------------
library(tidyverse)
library(AnnotationHub)
library(ensembldb)
library(ggplot2)
library(gridExtra)
library(grid)
library(patchwork)
library(scales)

# GLOBAL PARAMETERS ----------------------------------------------------------
BASE_DIR <- "/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space"
setwd(BASE_DIR)

# Input directory containing transcript list CSVs
INPUT_DIR <- "./output_files/isoTWAS_transcript_lists"

# Output directory for canonical analysis
OUTPUT_DIR <- "./output_files/isoTWAS_canonical_analysis"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Color scheme for visualizations
COLORS <- list(
  canonical = "#2E86AB",      # Blue
  non_canonical = "#F77F00",  # Orange
  unknown = "#6C757D"         # Grey
)

# SECTION 1: LOAD TRANSCRIPT METADATA FROM ANNOTATIONHUB --------------------
cat("\n========================================\n")
cat("LOADING TRANSCRIPT METADATA\n")
cat("========================================\n\n")

cat("Querying AnnotationHub for EnsDb.Hsapiens.v113 (GENCODE v47)...\n")

# Initialize AnnotationHub
ah <- AnnotationHub()

# Query for human EnsDb databases
human_ens <- query(ah, c("Homo sapiens", "EnsDb"))

# Find EnsDb v113 (corresponds to GENCODE v47)
# Note: The exact AH ID may vary - we'll search for v113
edb_candidates <- human_ens[grepl("v113", human_ens$title, ignore.case = TRUE)]

if (length(edb_candidates) > 0) {
  # Use the first v113 match
  edb <- edb_candidates[[1]]
  cat("✓ Loaded EnsDb:", edb_candidates$title[1], "\n")
  cat("  AH ID:", names(edb_candidates)[1], "\n")
} else {
  # Fallback: use latest available
  cat("⚠ EnsDb v113 not found. Using latest available EnsDb...\n")
  latest_idx <- length(human_ens)
  edb <- human_ens[[latest_idx]]
  cat("✓ Loaded EnsDb:", human_ens$title[latest_idx], "\n")
  cat("  AH ID:", names(human_ens)[latest_idx], "\n")
}

# Extract transcript-level metadata
cat("\nExtracting transcript metadata...\n")
tx_metadata_raw <- transcripts(edb, return.type = "data.frame")

cat("  Total transcripts in EnsDb:", nrow(tx_metadata_raw), "\n")

# Process metadata: select relevant columns and strip version numbers
tx_metadata <- tx_metadata_raw %>%
  dplyr::select(tx_id, tx_is_canonical, tx_biotype) %>%
  dplyr::mutate(tx_id_clean = str_remove(tx_id, "\\..*")) %>%  # Remove version
  dplyr::select(tx_id_clean, tx_is_canonical, tx_biotype)

cat("  Canonical transcripts:", sum(tx_metadata$tx_is_canonical == 1, na.rm = TRUE), "\n")
cat("  Non-canonical transcripts:", sum(tx_metadata$tx_is_canonical == 0, na.rm = TRUE), "\n")

# SECTION 2: LOAD TRANSCRIPT LIST CSVS --------------------------------------
cat("\n========================================\n")
cat("LOADING TRANSCRIPT LIST CSVS\n")
cat("========================================\n\n")

# Find only SCZ transcript list CSV files
csv_files <- list.files(INPUT_DIR, pattern = "SCZ.*_transcript_list\\.csv$", full.names = TRUE)

if (length(csv_files) == 0) {
  stop("❌ Error: No SCZ transcript list CSV files found in ", INPUT_DIR,
       "\n   Please run ch3_isoTWAS_DE_heatmaps.R first to generate transcript lists.")
}

cat("Found", length(csv_files), "SCZ transcript list files:\n")
for (f in basename(csv_files)) {
  cat("  -", f, "\n")
}

# Read all CSV files and combine
cat("\nReading CSV files...\n")
transcript_data_list <- lapply(csv_files, function(file) {
  df <- read.csv(file, stringsAsFactors = FALSE)
  df$source_file <- basename(file)
  return(df)
})

# Combine all dataframes
transcript_data <- bind_rows(transcript_data_list)

cat("✓ Loaded", nrow(transcript_data), "SCZ transcript records\n")
cat("  Unique transcripts:", length(unique(transcript_data$transcript_id)), "\n")
cat("  Trait: SCZ\n")
cat("  Timepoints:", paste(sort(unique(transcript_data$timepoint)), collapse = ", "), "\n")

# SECTION 3: ANNOTATE WITH CANONICAL STATUS ---------------------------------
cat("\n========================================\n")
cat("ANNOTATING WITH CANONICAL STATUS\n")
cat("========================================\n\n")

# Join with transcript metadata
cat("Matching transcript IDs with EnsDb metadata...\n")

annotated_data <- transcript_data %>%
  dplyr::left_join(tx_metadata, by = c("transcript_id" = "tx_id_clean"))

# Calculate matching statistics
n_total <- nrow(annotated_data)
n_matched <- sum(!is.na(annotated_data$tx_is_canonical))
n_unmatched <- sum(is.na(annotated_data$tx_is_canonical))
pct_matched <- round(100 * n_matched / n_total, 1)

cat("  ✓ Matched:", n_matched, "transcripts (", pct_matched, "%)\n", sep = "")
cat("  ⚠ Unmatched:", n_unmatched, "transcripts (", round(100 - pct_matched, 1), "%)\n", sep = "")

# Classify transcripts into canonical/non-canonical/unknown
annotated_data <- annotated_data %>%
  dplyr::mutate(
    canonical_status = case_when(
      is.na(tx_is_canonical) ~ "unknown",
      tx_is_canonical == 1 | tx_is_canonical == TRUE ~ "canonical",
      tx_is_canonical == 0 | tx_is_canonical == FALSE ~ "non-canonical",
      TRUE ~ "unknown"
    ),
    canonical_status = factor(canonical_status,
                              levels = c("canonical", "non-canonical", "unknown"))
  )

# Summary by canonical status
cat("\nCanonical status distribution:\n")
status_summary <- annotated_data %>%
  group_by(canonical_status) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percentage = round(100 * count / sum(count), 1))

print(status_summary, n = Inf)

# Save annotated data for potential downstream use
saveRDS(annotated_data, file.path(OUTPUT_DIR, "annotated_transcript_data.rds"))
cat("\n✓ Saved annotated data to:", file.path(OUTPUT_DIR, "annotated_transcript_data.rds"), "\n")

# SECTION 4: VISUALIZATION FUNCTIONS -----------------------------------------
cat("\n========================================\n")
cat("CREATING VISUALIZATIONS\n")
cat("========================================\n\n")

# Function to create summary statistics table
create_summary_table <- function(data, timepoint_label) {
  summary_df <- data %>%
    dplyr::filter(timepoint == timepoint_label) %>%
    group_by(trait, canonical_status) %>%
    summarise(count = n(), .groups = "drop") %>%
    pivot_wider(names_from = canonical_status, values_from = count, values_fill = 0)

  # Add missing canonical_status columns with 0 if they don't exist
  if(!"canonical" %in% names(summary_df)) summary_df$canonical <- 0
  if(!"non-canonical" %in% names(summary_df)) summary_df$`non-canonical` <- 0

  summary_df <- summary_df %>%
    mutate(
      total = canonical + `non-canonical`,
      pct_canonical = round(100 * canonical / total, 1),
      pct_noncanonical = round(100 * `non-canonical` / total, 1)
    )

  return(summary_df)
}

# Function to create overall pie chart for a timepoint
create_pie_chart <- function(data, timepoint_label) {
  plot_data <- data %>%
    dplyr::filter(timepoint == timepoint_label) %>%
    group_by(canonical_status) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(
      percentage = round(100 * count / sum(count), 1),
      label = paste0(canonical_status, "\n", count, " (", percentage, "%)")
    )

  p <- ggplot(plot_data, aes(x = "", y = count, fill = canonical_status)) +
    geom_bar(stat = "identity", width = 1, color = "white", size = 1) +
    coord_polar("y", start = 0) +
    geom_text(aes(label = label),
              position = position_stack(vjust = 0.5),
              size = 4, fontface = "bold") +
    scale_fill_manual(values = c(
      "canonical" = COLORS$canonical,
      "non-canonical" = COLORS$non_canonical,
      "unknown" = COLORS$unknown
    )) +
    labs(title = paste0(timepoint_label, " Overall Distribution"),
         fill = "Status") +
    theme_void() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11, face = "bold")
    )

  return(p)
}

# Function to create per-cell-type horizontal bar chart
create_celltype_barplot <- function(data, timepoint_label) {
  plot_data <- data %>%
    dplyr::filter(timepoint == timepoint_label) %>%
    group_by(cell_type_of_origin, canonical_status) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(cell_type_of_origin) %>%
    mutate(total = sum(count)) %>%
    ungroup() %>%
    arrange(desc(total)) %>%
    mutate(cell_type_of_origin = factor(cell_type_of_origin,
                                         levels = rev(unique(cell_type_of_origin))))

  p <- ggplot(plot_data, aes(x = cell_type_of_origin, y = count, fill = canonical_status)) +
    geom_bar(stat = "identity", position = "stack", color = "white", size = 0.3) +
    geom_text(aes(label = ifelse(count >= 3, count, "")),
              position = position_stack(vjust = 0.5),
              size = 2.8, color = "white", fontface = "bold") +
    coord_flip() +
    scale_fill_manual(values = c(
      "canonical" = COLORS$canonical,
      "non-canonical" = COLORS$non_canonical,
      "unknown" = COLORS$unknown
    )) +
    labs(
      title = paste0(timepoint_label, " SCZ Per-Cell-Type Breakdown (All Cell Types)"),
      x = "Cell Type of Origin",
      y = "Number of Transcripts",
      fill = "Status"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 9),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )

  return(p)
}

# Function to create per-biotype horizontal bar chart
create_biotype_barplot <- function(data, timepoint_label) {
  plot_data <- data %>%
    dplyr::filter(timepoint == timepoint_label) %>%
    # Handle missing biotypes
    mutate(tx_biotype = ifelse(is.na(tx_biotype), "unknown", tx_biotype)) %>%
    group_by(tx_biotype, canonical_status) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(tx_biotype) %>%
    mutate(total = sum(count)) %>%
    ungroup() %>%
    arrange(desc(total)) %>%
    mutate(tx_biotype = factor(tx_biotype,
                                levels = rev(unique(tx_biotype))))

  p <- ggplot(plot_data, aes(x = tx_biotype, y = count, fill = canonical_status)) +
    geom_bar(stat = "identity", position = "stack", color = "white", size = 0.3) +
    geom_text(aes(label = ifelse(count >= 3, count, "")),
              position = position_stack(vjust = 0.5),
              size = 2.8, color = "white", fontface = "bold") +
    coord_flip() +
    scale_fill_manual(values = c(
      "canonical" = COLORS$canonical,
      "non-canonical" = COLORS$non_canonical,
      "unknown" = COLORS$unknown
    )) +
    labs(
      title = paste0(timepoint_label, " SCZ Transcript Biotype Breakdown"),
      x = "Transcript Biotype",
      y = "Number of Transcripts",
      fill = "Status"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 9),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )

  return(p)
}

# Function to create summary statistics table plot
create_stats_table_plot <- function(data, timepoint_label) {
  summary_df <- create_summary_table(data, timepoint_label)

  # Format for display
  display_df <- summary_df %>%
    dplyr::select(
      Trait = trait,
      Canonical = canonical,
      `Non-Canonical` = `non-canonical`,
      Total = total,
      `% Canonical` = pct_canonical,
      `% Non-Canonical` = pct_noncanonical
    )

  # Create table grob
  table_grob <- gridExtra::tableGrob(
    display_df,
    rows = NULL,
    theme = gridExtra::ttheme_minimal(
      base_size = 10,
      padding = unit(c(4, 3), "mm"),
      core = list(
        fg_params = list(fontface = 1),
        bg_params = list(fill = c("grey95", "grey99"))
      ),
      colhead = list(
        fg_params = list(fontface = "bold"),
        bg_params = list(fill = "grey80")
      )
    )
  )

  return(table_grob)
}

# SECTION 5: GENERATE MULTI-PAGE PDF ----------------------------------------
cat("\nGenerating multi-page PDF...\n")

pdf_file <- file.path(OUTPUT_DIR, "SCZ_canonical_vs_noncanonical_analysis.pdf")
pdf(pdf_file, width = 11, height = 8.5)  # Landscape orientation

# PAGE 1: Overall Summary ----
cat("  Page 1: Overall Summary\n")

# All data is already SCZ-only (no filter needed)
scz_data <- annotated_data

overall_summary <- scz_data %>%
  group_by(canonical_status) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percentage = round(100 * count / sum(count), 1))

overall_bar <- ggplot(overall_summary, aes(x = canonical_status, y = count, fill = canonical_status)) +
  geom_bar(stat = "identity", color = "black", size = 0.5) +
  geom_text(aes(label = paste0(count, "\n(", percentage, "%)")),
            vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c(
    "canonical" = COLORS$canonical,
    "non-canonical" = COLORS$non_canonical,
    "unknown" = COLORS$unknown
  )) +
  labs(
    title = "SCZ Canonical vs Non-Canonical Distribution",
    subtitle = "SCZ TWAS-DE Hits Across All Timepoints",
    x = "Canonical Status",
    y = "Number of Transcripts"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "none"
  ) +
  ylim(0, max(overall_summary$count) * 1.15)

# Create summary statistics text
summary_text <- paste0(
  "SCZ Dataset Summary\n\n",
  "Total TWAS-DE Hits: ", nrow(scz_data), "\n",
  "Unique Transcripts: ", length(unique(scz_data$transcript_id)), "\n",
  "Timepoints: ", paste(sort(unique(scz_data$timepoint)), collapse = ", "), "\n\n",
  "Canonical Status:\n",
  "  Canonical: ", overall_summary$count[overall_summary$canonical_status == "canonical"],
  " (", overall_summary$percentage[overall_summary$canonical_status == "canonical"], "%)\n",
  "  Non-Canonical: ", overall_summary$count[overall_summary$canonical_status == "non-canonical"],
  " (", overall_summary$percentage[overall_summary$canonical_status == "non-canonical"], "%)\n",
  "  Unknown: ", overall_summary$count[overall_summary$canonical_status == "unknown"],
  " (", overall_summary$percentage[overall_summary$canonical_status == "unknown"], "%)"
)

summary_grob <- textGrob(summary_text,
                         x = 0.5, y = 0.5,
                         just = "center",
                         gp = gpar(fontsize = 12, fontfamily = "mono"))

# Combine plots
grid.arrange(
  overall_bar,
  summary_grob,
  ncol = 1,
  heights = c(2, 1)
)

# PAGES 2-4: Per-Timepoint Analysis ----
timepoints <- c("1M", "3M", "6M")

for (tp in timepoints) {
  page_num <- which(timepoints == tp) + 1
  cat("  Page", page_num, ":", tp, "Timepoint Analysis\n")

  # Check if data exists for this timepoint
  tp_data <- annotated_data %>% dplyr::filter(timepoint == tp)

  if (nrow(tp_data) == 0) {
    # Create placeholder page
    plot.new()
    text(0.5, 0.5, paste0("No data available for ", tp, " timepoint"),
         cex = 2, font = 2)
    next
  }

  # Create all panels
  pie_plot <- create_pie_chart(annotated_data, tp)
  celltype_plot <- create_celltype_barplot(annotated_data, tp)
  biotype_plot <- create_biotype_barplot(annotated_data, tp)
  stats_table <- create_stats_table_plot(annotated_data, tp)

  # Combine using patchwork - 4-panel layout
  combined_plot <- (pie_plot | celltype_plot) / biotype_plot / wrap_elements(stats_table) +
    plot_layout(heights = c(1.5, 0.8, 0.7)) +
    plot_annotation(
      title = paste0(tp, " Timepoint: SCZ Canonical vs Non-Canonical Analysis"),
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )

  print(combined_plot)
}

# Close PDF device
dev.off()

cat("\n✓ PDF saved to:", pdf_file, "\n")

# SECTION 6: OPTIONAL CSV EXPORT ---------------------------------------------
cat("\n========================================\n")
cat("EXPORTING ENHANCED SCZ CSV FILES\n")
cat("========================================\n\n")

# Export enhanced CSVs with canonical status (SCZ only)
for (tp in unique(annotated_data$timepoint)) {
  tr <- "SCZ"  # Fixed trait - all data is SCZ
  export_data <- annotated_data %>%
    dplyr::filter(timepoint == tp) %>%
    dplyr::select(transcript_id, gene_symbol, cell_type_of_origin, trait, timepoint,
           p_val_adj, avg_log2FC, canonical_status, tx_biotype) %>%
    arrange(p_val_adj)

  if (nrow(export_data) > 0) {
    output_file <- file.path(OUTPUT_DIR, paste0(tr, "_", tp, "_with_canonical.csv"))
    write.csv(export_data, output_file, row.names = FALSE, quote = FALSE)
    cat("  ✓ Exported:", basename(output_file), "\n")
  }
}

# COMPLETION MESSAGE ---------------------------------------------------------
cat("\n========================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("========================================\n\n")
cat("Output files:\n")
cat("  📊 PDF Report:", pdf_file, "\n")
cat("  📁 Enhanced CSVs:", OUTPUT_DIR, "\n")
cat("  💾 Annotated data:", file.path(OUTPUT_DIR, "annotated_transcript_data.rds"), "\n\n")
