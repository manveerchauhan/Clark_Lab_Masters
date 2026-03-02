# =============================================================================
# COMPREHENSIVE BIOTYPE DISCOVERY CURVES - All Protocol/Data Type Combinations
# =============================================================================

library(ggplot2)
library(dplyr)
library(scales)

# Load our utils system
source("utils/data_navigator.R")
source("utils/biotype_analyzer.R")

theme_set(theme_minimal())

cat("🧬 COMPREHENSIVE BIOTYPE DISCOVERY ANALYSIS 🧬\n")
cat("Biotype expression discovery across all 7 protocol/data type combinations\n\n")

# -----------------------------------------------------------------------------
# Helpers: palettes, ordering, and consistent scales/themes
# -----------------------------------------------------------------------------

# Colorblind-friendly palettes (Okabe–Ito inspired)
get_biotype_palette <- function(mode = "major") {
  if (mode == "major") {
    # 3-category palette
    colors <- c(
      "Protein Coding" = "#0072B2",          # blue
      "Long Non-Coding RNA" = "#D55E00",    # vermillion
      "Others" = "#009E73"                   # green
    )
  } else {
    # 9-category palette for "all" biotypes
    colors <- c(
      "Protein Coding" = "#0072B2",            # blue
      "Long Non-Coding RNA" = "#D55E00",      # vermillion
      "Mitochondrial RNA" = "#009E73",        # green
      "Nonsense Mediated Decay" = "#E69F00",  # orange
      "Processed Transcript" = "#CC79A7",     # purple
      "Pseudogenes" = "#56B4E9",              # sky blue
      "Retained Intron" = "#F0E442",          # yellow
      "Short Non-Coding RNA" = "#000000",     # black
      "Other" = "#7F7F7F"                      # gray
    )
  }
  return(colors)
}

# Consistent order for legends and factors
get_biotype_order <- function(mode = "major") {
  if (mode == "major") return(c("Protein Coding", "Long Non-Coding RNA", "Others"))
  return(c("Protein Coding", "Long Non-Coding RNA", "Mitochondrial RNA",
           "Nonsense Mediated Decay", "Processed Transcript", "Pseudogenes",
           "Retained Intron", "Short Non-Coding RNA", "Other"))
}

# Consistent axis formatting for log scales
log_y_scale <- function(name) {
  scale_y_continuous(
    name = name,
    breaks = scales::log_breaks(n = 7),
    labels = scales::comma,
    trans = "log10"
  )
}

log_x_scale <- function() {
  scale_x_continuous(
    name = "# Valid Aligned Reads\n(oarfish 'valid_best_aln' reads)",
    breaks = c(1e5, 2e5, 5e5, 1e6, 2e6, 5e6, 1e7, 2e7, 5e7),
    labels = c("100K", "200K", "500K", "1M", "2M", "5M", "10M", "20M", "50M"),
    trans = "log10"
  )
}

clean_theme <- function() {
  theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "#F0F0F0"),
      panel.grid.minor = element_line(color = "#F5F5F5")
    )
}

# =============================================================================
# CONFIGURATION
# =============================================================================

# All 7 technology/data type combinations to analyze
COMBINATIONS <- list(
  list(tech = "dRNA", data_type = "bulk", name = "dRNA Bulk"),
  list(tech = "ONT", data_type = "bulk", name = "ONT Bulk"),
  list(tech = "ONT", data_type = "single_cell", name = "ONT Single-Cell"),
  list(tech = "ONT", data_type = "single_nucleus", name = "ONT Single-Nucleus"),
  list(tech = "PacBio", data_type = "bulk", name = "PacBio Bulk"),
  list(tech = "PacBio", data_type = "single_cell", name = "PacBio Single-Cell"),
  list(tech = "PacBio", data_type = "single_nucleus", name = "PacBio Single-Nucleus")
)

# Include all available cell lines for comprehensive visualization
CELL_LINES <- c("H146", "H1975", "H69", "H211", "H2228", "H526", "HCC827", "SHP77")

# Biotype category definitions
MAJOR_BIOTYPES <- c("Protein Coding", "Long Non-Coding RNA", "Others")
ALL_BIOTYPES <- c("Protein Coding", "Long Non-Coding RNA", "Mitochondrial RNA", 
                  "Nonsense Mediated Decay", "Other", "Processed Transcript", 
                  "Pseudogenes", "Retained Intron", "Short Non-Coding RNA")

# Expression threshold (same as rarefaction analysis)
COUNT_THRESHOLD <- 1

# =============================================================================
# DATA EXTRACTION
# =============================================================================

# Helper function to extract proportional biotype data directly from sample sheets
extract_proportional_biotype_data_direct <- function(combinations, cell_lines, biotype_mode = "major") {
  cat("Extracting proportional biotype data directly from sample sheets...\n")
  
  # Load both sample sheets
  #bulk_data <- read.csv("bulk_cleaned_sample_sheet_with_counts_with_groupings_with_read_counts.csv", 
  #                      stringsAsFactors = FALSE)
  #sc_sn_data <- read.csv("sc_sn_cleaned_sample_sheet_with_counts_with_groupings_with_read_counts.csv", 
  #                      stringsAsFactors = FALSE)

  # Longbench spec sample sheets:
  bulk_data <- read.csv("bulk_cleaned_sample_sheet_with_counts_with_groupings_with_valid_best_aln.csv", 
                        stringsAsFactors = FALSE)
  sc_sn_data <- read.csv("pseudobulk_sample_sheet_with_valid_best_aln_clean.csv", 
                         stringsAsFactors = FALSE)
  
  # Fix column compatibility (use common columns only)
  common_cols <- intersect(names(bulk_data), names(sc_sn_data))
  bulk_data <- bulk_data[, common_cols]
  sc_sn_data <- sc_sn_data[, common_cols]
  
  # Filter out "All" cell lines from sc_sn data
  sc_sn_data <- sc_sn_data[sc_sn_data$cell_line != "All", ]
  
  # Combine datasets
  all_data <- rbind(bulk_data, sc_sn_data)
  
  # Filter for proportional data only
  proportional_data <- all_data[all_data$measurement_type == "proportional" & 
                               !is.na(all_data$valid_best_aln), ]
  
  # Initialize results list to collect all rows
  all_biotype_rows <- list()
  
  # Process each combination
  for (combo in combinations) {
    combo_name <- combo$name
    tech <- combo$tech
    data_type <- combo$data_type
    
    # Filter data for this combination
    combo_data <- proportional_data[proportional_data$technology == tech & 
                                   proportional_data$data_type == data_type &
                                   proportional_data$cell_line %in% cell_lines, ]
    
    if (nrow(combo_data) == 0) next
    
    cat(sprintf("Processing %s proportional data (%d entries)...\n", combo_name, nrow(combo_data)))
    
    for (i in 1:nrow(combo_data)) {
      row <- combo_data[i, ]
      
      # Get biotype data for this file
      if (file.exists(row$file_path)) {
        biotype_summary <- BiotypeAnalyzer$get_expressed_by_biotype(
          technology = tech,
          data_type = data_type,
          cell_line = row$cell_line,
          sampling_depth = row$sampling_depth,
          count_threshold = COUNT_THRESHOLD
        )
        
        if (!is.null(biotype_summary) && nrow(biotype_summary) > 0) {
          # Process biotype data based on mode
          if (biotype_mode == "major") {
            # Group into major categories
            protein_coding <- sum(biotype_summary$expressed_transcripts[
              biotype_summary$broad_tx_biotype == "Protein Coding"], na.rm = TRUE)
            lnc_rna <- sum(biotype_summary$expressed_transcripts[
              biotype_summary$broad_tx_biotype == "Long Non-Coding RNA"], na.rm = TRUE)
            others <- sum(biotype_summary$expressed_transcripts[
              !biotype_summary$broad_tx_biotype %in% c("Protein Coding", "Long Non-Coding RNA")], na.rm = TRUE)
            
            # Create data rows for each major biotype
            for (biotype_cat in c("Protein Coding", "Long Non-Coding RNA", "Others")) {
              expressed_count <- switch(biotype_cat,
                "Protein Coding" = protein_coding,
                "Long Non-Coding RNA" = lnc_rna,
                "Others" = others
              )
              
              new_row <- data.frame(
                combination = as.character(combo_name),
                technology = as.character(tech),
                data_type = as.character(data_type),
                cell_line = as.character(row$cell_line),
                sampling_depth = as.character(row$sampling_depth),
                measurement_type = as.character("proportional"),
                valid_best_aln = as.numeric(row$valid_best_aln),
                biotype_category = as.character(biotype_cat),
                expressed_transcripts = as.numeric(expressed_count),
                sample_id = as.character(row$sample_id),
                stringsAsFactors = FALSE
              )
              
              all_biotype_rows[[length(all_biotype_rows) + 1]] <- new_row
            }
          } else {
            # Use all individual biotype categories
            for (j in 1:nrow(biotype_summary)) {
              biotype_row <- biotype_summary[j, ]
              
              new_row <- data.frame(
                combination = as.character(combo_name),
                technology = as.character(tech),
                data_type = as.character(data_type),
                cell_line = as.character(row$cell_line),
                sampling_depth = as.character(row$sampling_depth),
                measurement_type = as.character("proportional"),
                valid_best_aln = as.numeric(row$valid_best_aln),
                biotype_category = as.character(biotype_row$broad_tx_biotype),
                expressed_transcripts = as.numeric(biotype_row$expressed_transcripts),
                sample_id = as.character(row$sample_id),
                stringsAsFactors = FALSE
              )
              
              all_biotype_rows[[length(all_biotype_rows) + 1]] <- new_row
            }
          }
        }
      }
    }
  }
  
  # Combine all rows into a single data frame
  if (length(all_biotype_rows) > 0) {
    all_biotype_data <- do.call(rbind, all_biotype_rows)
  } else {
    # Return empty data frame with proper structure
    all_biotype_data <- data.frame(
      combination = character(0),
      technology = character(0),
      data_type = character(0),
      cell_line = character(0),
      sampling_depth = character(0),
      measurement_type = character(0),
      valid_best_aln = numeric(0),
      biotype_category = character(0),
      expressed_transcripts = numeric(0),
      sample_id = character(0),
      stringsAsFactors = FALSE
    )
  }
  
  cat(sprintf("Proportional biotype data extracted: %d total data points\n", nrow(all_biotype_data)))
  return(all_biotype_data)
}

# Main function to extract comprehensive biotype discovery data
extract_biotype_discovery_data <- function(combinations = COMBINATIONS, cell_lines = CELL_LINES, biotype_mode = "major") {
  cat("=== EXTRACTING COMPREHENSIVE BIOTYPE DISCOVERY DATA ===\n")
  cat(sprintf("Mode: %s biotypes\n", biotype_mode))
  cat(sprintf("Expression threshold: > %d counts\n", COUNT_THRESHOLD))
  cat(sprintf("Combinations: %d\n", length(combinations)))
  cat(sprintf("Cell lines: %d\n", length(cell_lines)))
  cat("\n")
  
  # Initialize main data frame with proper column structure
  all_biotype_data <- data.frame(
    combination = character(0),
    technology = character(0),
    data_type = character(0),
    cell_line = character(0),
    sampling_depth = character(0),
    measurement_type = character(0),
    valid_best_aln = numeric(0),
    biotype_category = character(0),
    expressed_transcripts = numeric(0),
    sample_id = character(0),
    stringsAsFactors = FALSE
  )
  
  # Extract proportional data first (direct method)
  proportional_data <- extract_proportional_biotype_data_direct(combinations, cell_lines, biotype_mode)
  all_biotype_data <- rbind(all_biotype_data, proportional_data)
  
  # Extract absolute data using DataNav (fixed method using rarefaction pattern)
  cat("Extracting absolute biotype data using DataNav...\n")
  
  # Initialize list to collect all absolute data rows
  all_absolute_rows <- list()
  
  for (combo in combinations) {
    combo_name <- combo$name
    tech <- combo$tech
    data_type <- combo$data_type
    
    cat(sprintf("Processing %s absolute data...\n", combo_name))
    
    # Get all available depths for absolute measurements using the correct DataNav pattern
    available_meta <- DataNav$get(
      technology = tech,
      data_type = data_type,
      measurement_type = "absolute",
      return_meta = TRUE
    )
    
    if (nrow(available_meta) > 0) {
      available_depths <- unique(available_meta$sampling_depth)
      # Sort depths properly
      numeric_depths <- suppressWarnings(as.numeric(gsub("M|full", "", available_depths)))
      sorted_indices <- order(numeric_depths, na.last = TRUE)
      available_depths <- available_depths[sorted_indices]
      
      cat(sprintf("  • Available depths: %s\n", paste(available_depths, collapse=", ")))
      
      # Extract absolute data for each cell line and depth
      for (cell_line in cell_lines) {
        for (depth in available_depths) {
          meta <- DataNav$get(
            technology = tech,
            cell_line = cell_line,
            data_type = data_type,
            measurement_type = "absolute",
            sampling_depth = depth,
            return_meta = TRUE
          )
          
          if (nrow(meta) > 0) {
            # Check if valid_best_aln column exists and has a valid value
            if ("valid_best_aln" %in% names(meta) && 
                !is.null(meta$valid_best_aln[1]) && 
                !is.na(meta$valid_best_aln[1])) {
              
              # Get biotype data for this file
              biotype_summary <- BiotypeAnalyzer$get_expressed_by_biotype(
                technology = tech,
                data_type = data_type,
                cell_line = cell_line,
                sampling_depth = depth,
                count_threshold = COUNT_THRESHOLD
              )
              
              if (!is.null(biotype_summary) && nrow(biotype_summary) > 0) {
                # Process biotype data based on mode
                if (biotype_mode == "major") {
                  # Group into major categories
                  protein_coding <- sum(biotype_summary$expressed_transcripts[
                    biotype_summary$broad_tx_biotype == "Protein Coding"], na.rm = TRUE)
                  lnc_rna <- sum(biotype_summary$expressed_transcripts[
                    biotype_summary$broad_tx_biotype == "Long Non-Coding RNA"], na.rm = TRUE)
                  others <- sum(biotype_summary$expressed_transcripts[
                    !biotype_summary$broad_tx_biotype %in% c("Protein Coding", "Long Non-Coding RNA")], na.rm = TRUE)
                  
                  # Create data rows for each major biotype
                  for (biotype_cat in c("Protein Coding", "Long Non-Coding RNA", "Others")) {
                    expressed_count <- switch(biotype_cat,
                      "Protein Coding" = protein_coding,
                      "Long Non-Coding RNA" = lnc_rna,
                      "Others" = others
                    )
                    
                    new_row <- data.frame(
                      combination = as.character(combo_name),
                      technology = as.character(tech),
                      data_type = as.character(data_type),
                      cell_line = as.character(cell_line),
                      sampling_depth = as.character(depth),
                      measurement_type = as.character("absolute"),
                      valid_best_aln = as.numeric(meta$valid_best_aln[1]),
                      biotype_category = as.character(biotype_cat),
                      expressed_transcripts = as.numeric(expressed_count),
                      sample_id = as.character(meta$sample_id[1]),
                      stringsAsFactors = FALSE
                    )
                    
                    all_absolute_rows[[length(all_absolute_rows) + 1]] <- new_row
                  }
                } else {
                  # Use all individual biotype categories
                  for (j in 1:nrow(biotype_summary)) {
                    biotype_row <- biotype_summary[j, ]
                    
                    new_row <- data.frame(
                      combination = as.character(combo_name),
                      technology = as.character(tech),
                      data_type = as.character(data_type),
                      cell_line = as.character(cell_line),
                      sampling_depth = as.character(depth),
                      measurement_type = as.character("absolute"),
                      valid_best_aln = as.numeric(meta$valid_best_aln[1]),
                      biotype_category = as.character(biotype_row$broad_tx_biotype),
                      expressed_transcripts = as.numeric(biotype_row$expressed_transcripts),
                      sample_id = as.character(meta$sample_id[1]),
                      stringsAsFactors = FALSE
                    )
                    
                    all_absolute_rows[[length(all_absolute_rows) + 1]] <- new_row
                  }
                }
              }
            } else {
              cat(sprintf("    Warning: %s %s %s absolute - missing valid_best_aln data\n", 
                         tech, cell_line, depth))
            }
          }
        }
      }
      
      # Count absolute data points for this combination so far
      combo_absolute_points <- length(all_absolute_rows)
      cat(sprintf("  • Absolute biotype data points collected so far: %d\n", combo_absolute_points))
    } else {
      cat(sprintf("  No absolute data found for %s\n", combo_name))
    }
    cat("\n")
  }
  
  # Combine all absolute rows and add to main data frame
  if (length(all_absolute_rows) > 0) {
    absolute_data <- do.call(rbind, all_absolute_rows)
    all_biotype_data <- rbind(all_biotype_data, absolute_data)
  }
  
  cat("\n=== BIOTYPE DISCOVERY DATA EXTRACTION SUMMARY ===\n")
  total_points <- nrow(all_biotype_data)
  proportional_points <- nrow(all_biotype_data[all_biotype_data$measurement_type == "proportional", ])
  absolute_points <- nrow(all_biotype_data[all_biotype_data$measurement_type == "absolute", ])
  
  cat(sprintf("Total biotype data points: %d\n", total_points))
  cat(sprintf("Proportional: %d\n", proportional_points))
  cat(sprintf("Absolute: %d\n", absolute_points))
  
  # Show biotype category breakdown
  if (biotype_mode == "major") {
    cat("\nMajor biotype categories:\n")
    biotype_counts <- table(all_biotype_data$biotype_category)
    for (biotype in names(biotype_counts)) {
      cat(sprintf("  %s: %d data points\n", biotype, biotype_counts[biotype]))
    }
  } else {
    cat(sprintf("\nAll biotype categories: %d unique categories\n", 
                length(unique(all_biotype_data$biotype_category))))
    biotype_counts <- table(all_biotype_data$biotype_category)
    for (biotype in names(biotype_counts)) {
      cat(sprintf("  %s: %d data points\n", biotype, biotype_counts[biotype]))
    }
  }
  
  cat("\n")
  return(all_biotype_data)
}

# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

# Create comprehensive biotype discovery plot (single plot, all protocols)
plot_biotype_discovery_comprehensive <- function(data, biotype_mode = "major", 
                                                 title = "Comprehensive Biotype Discovery Curves", 
                                                 filename = NULL) {
  
  # Sort data to ensure proper curve connections
  data <- data %>%
    arrange(combination, biotype_category, cell_line, valid_best_aln, expressed_transcripts)
  
  # Optional: summarize across cell lines (median line with IQR ribbon)
  # Set summarize_cell_lines = TRUE to enable summary mode
  summarize_cell_lines <- FALSE
  
  # Biotype palette and ordering
  biotype_colors <- get_biotype_palette(biotype_mode)
  biotype_levels <- get_biotype_order(biotype_mode)
  data$biotype_category <- factor(data$biotype_category, levels = biotype_levels)
  
  # Create the plot
  p <- ggplot(data, aes(x = valid_best_aln, y = expressed_transcripts)) +
    geom_smooth(aes(color = biotype_category, group = interaction(combination, biotype_category, cell_line)), 
                method = "loess", se = FALSE, span = 0.3, size = 0.8, alpha = 0.7) +
    geom_point(aes(color = biotype_category), size = 1.5, alpha = 0.6) +
    
    # Color and scales
    scale_color_manual(values = biotype_colors, name = "Biotype Category") +
    log_x_scale() +
    log_y_scale("Expressed Transcripts") +
    
    # Theming
    labs(
      title = title,
      subtitle = sprintf("Biotype expression discovery across all protocols (%s categories)", biotype_mode)
    ) +
    clean_theme() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      strip.text = element_text(size = 10, face = "bold")
    )
  
  # Save plot if filename provided
  if (!is.null(filename)) {
    plot_dir <- dirname(filename)
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir, recursive = TRUE)
      cat(sprintf("Created directory: %s\n", plot_dir))
    }
    ggsave(filename, plot = p, width = 14, height = 10, dpi = 300)
    cat(sprintf("Comprehensive biotype discovery plot saved as: %s\n", filename))
  }
  
  return(p)
}

# Create grid biotype discovery plot (one panel per protocol)
plot_biotype_discovery_grid <- function(data, biotype_mode = "major", 
                                        title = "Biotype Discovery by Protocol", 
                                        filename = NULL) {
  
  # Sort data to ensure proper curve connections
  data <- data %>%
    arrange(combination, biotype_category, cell_line, valid_best_aln, expressed_transcripts)
  
  # Toggle for summary mode (median line with IQR ribbon across cell lines)
  summarize_cell_lines <- FALSE
  
  # Biotype palette and ordering
  biotype_colors <- get_biotype_palette(biotype_mode)
  biotype_levels <- get_biotype_order(biotype_mode)
  data$biotype_category <- factor(data$biotype_category, levels = biotype_levels)
  
  # Create the grid plot (optionally summarized)
  if (summarize_cell_lines) {
    # Summarize across cell lines by sampling_depth (or binned aligned reads if needed)
    suppressPackageStartupMessages(require(dplyr))
    summary_df <- data %>%
      group_by(combination, biotype_category, sampling_depth) %>%
      summarise(
        valid_best_aln = median(valid_best_aln, na.rm = TRUE),
        y_med = median(expressed_transcripts, na.rm = TRUE),
        y_q1 = quantile(expressed_transcripts, 0.25, na.rm = TRUE),
        y_q3 = quantile(expressed_transcripts, 0.75, na.rm = TRUE),
        .groups = "drop"
      )
    p <- ggplot(summary_df, aes(x = valid_best_aln, y = y_med, color = biotype_category)) +
      geom_ribbon(aes(ymin = y_q1, ymax = y_q3, fill = biotype_category), alpha = 0.15, color = NA) +
      geom_line(size = 1.3, alpha = 0.95) +
      scale_fill_manual(values = biotype_colors, guide = "none")
  } else {
    p <- ggplot(data, aes(x = valid_best_aln, y = expressed_transcripts)) +
      geom_smooth(aes(color = biotype_category, group = interaction(biotype_category, cell_line)), 
                  method = "loess", se = FALSE, span = 0.3, size = 1.2, alpha = 0.9) +
      geom_point(aes(color = biotype_category), size = 1.0, alpha = 0.5)
  }
    
    # Grid layout and scales/theme appended to p
    p <- p +
      facet_wrap(~ combination, scales = "free", ncol = 2) +
      
      # Color and scales
      scale_color_manual(values = biotype_colors, name = "Biotype Category") +
      log_x_scale() +
      log_y_scale("Expressed Transcripts") +
      
      # Theming
      labs(
        title = title,
        subtitle = sprintf("Biotype expression discovery across protocols (%s categories)", biotype_mode)
      ) +
      clean_theme() +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  
  # Calculate dynamic height based on number of combinations
  n_combos <- length(unique(data$combination))
  n_rows <- ceiling(n_combos / 2)
  plot_height <- max(8, n_rows * 3)
  
  # Save plot if filename provided
  if (!is.null(filename)) {
    plot_dir <- dirname(filename)
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir, recursive = TRUE)
      cat(sprintf("Created directory: %s\n", plot_dir))
    }
    ggsave(filename, plot = p, width = 14, height = plot_height, dpi = 300)
    cat(sprintf("Grid biotype discovery plot saved as: %s\n", filename))
  }
  
  return(p)
}

# Create cell line-specific biotype discovery plots (one page per cell line)
plot_biotype_discovery_by_cell_line <- function(data, biotype_mode = "major", 
                                                title_prefix = "Biotype Discovery", 
                                                filename = NULL) {
  
  # Sort data to ensure proper curve connections
  data <- data %>%
    arrange(combination, biotype_category, cell_line, valid_best_aln, expressed_transcripts)
  
  # Define color palettes based on biotype mode
  if (biotype_mode == "major") {
    biotype_colors <- c(
      "Protein Coding" = "#1F78B4",
      "Long Non-Coding RNA" = "#E31A1C", 
      "Others" = "#33A02C"
    )
  } else {
    # Use a larger color palette for all biotypes
    biotype_colors <- c(
      "Protein Coding" = "#1F78B4",
      "Long Non-Coding RNA" = "#E31A1C",
      "Mitochondrial RNA" = "#33A02C",
      "Nonsense Mediated Decay" = "#FF7F00",
      "Other" = "#6A3D9A",
      "Processed Transcript" = "#B15928",
      "Pseudogenes" = "#FB9A99",
      "Retained Intron" = "#FDBF6F",
      "Short Non-Coding RNA" = "#CAB2D6"
    )
  }
  
  # Protocol colors for distinguishing between protocols
  protocol_colors <- c(
    "dRNA Bulk" = "#E31A1C",
    "ONT Bulk" = "#1F78B4", 
    "ONT Single-Cell" = "#33A02C",
    "ONT Single-Nucleus" = "#FF7F00",
    "PacBio Bulk" = "#6A3D9A",
    "PacBio Single-Cell" = "#B15928",
    "PacBio Single-Nucleus" = "#FB9A99"
  )
  
  # Create list to store individual plots
  cell_line_plots <- list()
  
  # Get unique cell lines
  unique_cell_lines <- unique(data$cell_line)
  
  # Create a plot for each cell line
  for (cell_line in unique_cell_lines) {
    cell_data <- data[data$cell_line == cell_line, ]
    
    if (nrow(cell_data) == 0) next
    
    # Create plot for this cell line
    p <- ggplot(cell_data, aes(x = valid_best_aln, y = expressed_transcripts)) +
      geom_smooth(aes(color = biotype_category, linetype = combination,
                      group = interaction(biotype_category, combination)), 
                  method = "loess", se = FALSE, span = 0.3, size = 1.0, alpha = 0.8) +
      geom_point(aes(color = biotype_category, shape = combination), 
                 size = 2, alpha = 0.7) +
      
      # Color and scales
      scale_color_manual(values = biotype_colors, name = "Biotype Category") +
      scale_linetype_manual(
        values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "solid"),
        name = "Protocol"
      ) +
      scale_shape_manual(
        values = c(16, 17, 15, 18, 8, 11, 12),
        name = "Protocol"
      ) +
      scale_x_continuous(
        name = "# Valid Aligned Reads",
        breaks = c(100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000, 20000000, 50000000),
        labels = c("100K", "200K", "500K", "1M", "2M", "5M", "10M", "20M", "50M"),
        trans = "log10"
      ) +
      scale_y_continuous(
        name = "Expressed Transcripts",
        labels = comma_format(),
        trans = "log10"
      ) +
      
      # Theming
      labs(
        title = sprintf("%s - %s Cell Line", title_prefix, cell_line),
        subtitle = sprintf("Biotype discovery across all protocols (%s categories)", biotype_mode)
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal"
      ) +
      guides(
        color = guide_legend(override.aes = list(size = 3), ncol = 3),
        linetype = guide_legend(override.aes = list(size = 1), ncol = 4),
        shape = guide_legend(override.aes = list(size = 3), ncol = 4)
      )
    
    cell_line_plots[[cell_line]] <- p
  }
  
  # Save as multi-page PDF if filename provided
  if (!is.null(filename)) {
    # Calculate appropriate dimensions
    plot_width <- 12
    plot_height <- 8
    
    plot_dir <- dirname(filename)
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir, recursive = TRUE)
      cat(sprintf("Created directory: %s\n", plot_dir))
    }
    pdf(filename, width = plot_width, height = plot_height)
    for (cell_line in names(cell_line_plots)) {
      print(cell_line_plots[[cell_line]])
    }
    dev.off()
    
    cat(sprintf("Cell line-specific biotype discovery plots saved as: %s (%d pages)\n", 
                filename, length(cell_line_plots)))
  }
  
  return(cell_line_plots)
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

# New: biotype summary by technology with IQR ribbons
plot_biotype_summary_by_technology <- function(data, biotype_mode = "all",
                                               title = "Biotype Discovery by Technology (Summary)",
                                               filename = NULL) {
  # Ensure ordering and palette
  biotype_colors <- get_biotype_palette(biotype_mode)
  biotype_levels <- get_biotype_order(biotype_mode)
  data$biotype_category <- factor(data$biotype_category, levels = biotype_levels)
  
  # Pretty labels for data_type for the legend
  data <- data %>% mutate(
    data_type_pretty = dplyr::recode(
      data_type,
      bulk = "Bulk",
      single_cell = "Single-Cell",
      single_nucleus = "Single-Nucleus"
    )
  )
  
  # Summarize across cell lines for each technology × data_type × biotype × depth
  summary_df <- data %>%
    group_by(technology, data_type_pretty, biotype_category, sampling_depth) %>%
    summarise(
      valid_best_aln = median(valid_best_aln, na.rm = TRUE),
      y_med = median(expressed_transcripts, na.rm = TRUE),
      y_q1 = quantile(expressed_transcripts, 0.25, na.rm = TRUE),
      y_q3 = quantile(expressed_transcripts, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Order biotypes in legend by overall median abundance (descending)
  biotype_order_by_median <- summary_df %>%
    group_by(biotype_category) %>%
    summarise(med = median(y_med, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(med)) %>% pull(biotype_category)
  summary_df$biotype_category <- factor(summary_df$biotype_category, levels = biotype_order_by_median)
  
  # Linetype mapping
  lty_values <- c("Bulk" = "solid", "Single-Cell" = "dashed", "Single-Nucleus" = "dotted")
  
  n_cell_lines <- length(unique(data$cell_line))
  n_depths <- length(unique(data$sampling_depth))
  
  p <- ggplot(summary_df, aes(x = valid_best_aln, y = y_med,
                              color = biotype_category, linetype = data_type_pretty,
                              group = interaction(biotype_category, data_type_pretty))) +
    geom_ribbon(aes(ymin = y_q1, ymax = y_q3, fill = biotype_category), alpha = 0.12, color = NA) +
    geom_line(size = 1.35) +
    facet_wrap(~ technology, scales = "free_x", ncol = 3) +
    scale_color_manual(values = biotype_colors, name = "Biotype Category") +
    scale_fill_manual(values = biotype_colors, guide = "none") +
    scale_linetype_manual(values = lty_values, name = "Data Type") +
    log_x_scale() +
    log_y_scale("Expressed Transcripts") +
    coord_cartesian(ylim = c(max(min(summary_df$y_q1, summary_df$y_med, na.rm = TRUE), 0.1),
                             max(summary_df$y_q3, summary_df$y_med, na.rm = TRUE))) +
    labs(
      title = title,
      subtitle = sprintf("Median ± IQR across %d cell lines; depths: %d levels", n_cell_lines, n_depths)
    ) +
    clean_theme() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      strip.text = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.spacing.y = grid::unit(0.1, "lines"),
      legend.key.height = grid::unit(0.5, "lines"),
      legend.key.width = grid::unit(0.8, "lines"),
      legend.box.spacing = grid::unit(0.2, "lines")
    ) +
    guides(
      color = guide_legend(ncol = 2, override.aes = list(size = 2)),
      linetype = guide_legend(override.aes = list(size = 1.6))
    )
  
  if (!is.null(filename)) {
    plot_dir <- dirname(filename)
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir, recursive = TRUE)
      cat(sprintf("Created directory: %s\n", plot_dir))
    }
    # Export as PDF for publication-quality vector output
    if (grepl("\\.pdf$", filename, ignore.case = TRUE)) {
      grDevices::pdf(filename, width = 14, height = 6)
      print(p)
      grDevices::dev.off()
    } else {
      ggsave(filename, plot = p, width = 14, height = 6, dpi = 300)
    }
    cat(sprintf("Biotype summary by technology plot saved as: %s\n", filename))
  }
  
  return(p)
}

# Runner for the summary plot
run_biotype_summary_by_technology <- function(biotype_mode = "all", save_plots = TRUE) {
  cat("📊 Creating biotype summary-by-technology plot...\n")
  biotype_data <- extract_biotype_discovery_data(biotype_mode = biotype_mode)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  out_file <- if (save_plots) sprintf("analysis_scripts/plots/biotype_summary_by_technology_%s_%s.pdf", biotype_mode, timestamp) else NULL
  p <- plot_biotype_summary_by_technology(
    biotype_data,
    biotype_mode = biotype_mode,
    title = sprintf("Biotype Discovery by Technology (Summary, %s)", tools::toTitleCase(biotype_mode)),
    filename = out_file
  )
  invisible(p)
}

# Function to run comprehensive biotype discovery analysis
run_biotype_discovery_analysis <- function(biotype_mode = "major", save_plots = TRUE) {
  cat("🚀 STARTING COMPREHENSIVE BIOTYPE DISCOVERY ANALYSIS 🚀\n")
  cat(sprintf("Mode: %s biotypes\n", biotype_mode))
  cat("==================================================\n\n")
  
  # Extract data
  biotype_data <- extract_biotype_discovery_data(biotype_mode = biotype_mode)
  
  if (nrow(biotype_data) == 0) {
    cat("❌ No biotype data extracted. Aborting analysis.\n")
    return(NULL)
  }
  
  # Generate timestamp for unique filenames
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # Create comprehensive plot
  cat("📊 Creating comprehensive biotype discovery plot...\n")
  comp_filename <- if (save_plots) sprintf("analysis_scripts/plots/biotype_discovery_comprehensive_%s_%s.png", biotype_mode, timestamp) else NULL
  comp_plot <- plot_biotype_discovery_comprehensive(
    biotype_data, 
    biotype_mode = biotype_mode,
    title = sprintf("Comprehensive Biotype Discovery Curves (%s)", tools::toTitleCase(biotype_mode)),
    filename = comp_filename
  )
  
  # Create grid plot
  cat("📊 Creating grid biotype discovery plot...\n")
  grid_filename <- if (save_plots) sprintf("analysis_scripts/plots/biotype_discovery_grid_%s_%s.pdf", biotype_mode, timestamp) else NULL
  grid_plot <- plot_biotype_discovery_grid(
    biotype_data, 
    biotype_mode = biotype_mode,
    title = sprintf("Biotype Discovery by Protocol (%s)", tools::toTitleCase(biotype_mode)),
    filename = grid_filename
  )
  
  # Create cell line-specific plots (NEW!)
  cat("📊 Creating cell line-specific biotype discovery plots...\n")
  cell_line_filename <- if (save_plots) sprintf("analysis_scripts/plots/biotype_discovery_by_cell_line_%s_%s.pdf", biotype_mode, timestamp) else NULL
  cell_line_plots <- plot_biotype_discovery_by_cell_line(
    biotype_data, 
    biotype_mode = biotype_mode,
    title_prefix = sprintf("Biotype Discovery (%s)", tools::toTitleCase(biotype_mode)),
    filename = cell_line_filename
  )
  
  # Save the data
  saveRDS(list(
    biotype_data = biotype_data,
    biotype_mode = biotype_mode,
    combinations_analyzed = COMBINATIONS,
    comprehensive_plot = comp_plot,
    grid_plot = grid_plot,
    cell_line_plots = cell_line_plots
  ), sprintf("analysis_scripts/biotype_discovery_data_%s.rds", biotype_mode))
  
  # Summary statistics
  cat("\n📈 BIOTYPE DISCOVERY ANALYSIS COMPLETE 📈\n")
  cat("==========================================\n")
  cat(sprintf("Total data points analyzed: %s\n", format(nrow(biotype_data), big.mark = ",")))
  cat(sprintf("Biotype categories: %d\n", length(unique(biotype_data$biotype_category))))
  cat(sprintf("Protocols analyzed: %d\n", length(unique(biotype_data$combination))))
  cat(sprintf("Cell lines included: %d\n", length(unique(biotype_data$cell_line))))
  
  if (save_plots) {
    cat("\n📁 Generated Files:\n")
    if (!is.null(comp_filename)) cat(sprintf("• %s\n", comp_filename))
    if (!is.null(grid_filename)) cat(sprintf("• %s\n", grid_filename))
    if (!is.null(cell_line_filename)) cat(sprintf("• %s (one page per cell line)\n", cell_line_filename))
    cat(sprintf("• analysis_scripts/biotype_discovery_data_%s.rds\n", biotype_mode))
  }
  
  cat("\n✅ Biotype discovery analysis complete! 🎉\n")
  
  # Return the data and plots for further analysis
  return(list(
    data = biotype_data,
    comprehensive_plot = comp_plot,
    grid_plot = grid_plot,
    cell_line_plots = cell_line_plots
  ))
}

# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

# Quick run with major biotypes (recommended default)
run_major_biotype_analysis <- function(save_plots = TRUE) {
  return(run_biotype_discovery_analysis(biotype_mode = "major", save_plots = save_plots))
}

# Quick run with all biotype categories
run_all_biotype_analysis <- function(save_plots = TRUE) {
  return(run_biotype_discovery_analysis(biotype_mode = "all", save_plots = save_plots))
}

cat("🧬 Biotype Discovery Analysis Ready! 🧬\n")
cat("=======================================\n")
cat("Quick start commands:\n")
cat("• run_major_biotype_analysis()     - Major biotypes (Protein Coding, lncRNA, Others)\n")
cat("• run_all_biotype_analysis()       - All 9 biotype categories\n") 
cat("• run_biotype_discovery_analysis(biotype_mode='major', save_plots=TRUE)\n")
cat("• run_biotype_summary_by_technology('all') - Summary by technology (PDF)\n")
cat("\nConfigurable biotype_mode: 'major' or 'all'\n")
cat("Plot types generated:\n")
cat("  📊 Comprehensive plot (all protocols & cell lines)\n")
cat("  📋 Grid plot (protocols in panels)\n")
cat("  📄 Cell line plots (separate pages per cell line)\n\n")

# Remove the automatic execution
# result <- run_all_biotype_analysis()
