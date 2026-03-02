# =============================================================================
# COMPREHENSIVE RAREFACTION CURVES - All Protocol/Data Type Combinations
# =============================================================================

library(ggplot2)
library(dplyr)
library(scales)

# Load our utils system
source("utils/data_navigator.R")
source("utils/biotype_analyzer.R")

theme_set(theme_minimal())

cat("🧬 COMPREHENSIVE RAREFACTION ANALYSIS 🧬\n")
cat("All 7 protocol/data type combinations on a single curve\n\n")

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

# =============================================================================
# DATA EXTRACTION
# =============================================================================

# Helper function to extract proportional data directly from sample sheets
extract_proportional_data_direct <- function(combinations, cell_lines) {
  cat("Extracting proportional data directly from sample sheets...\n")
  
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
  
  # Filter out "All" cell lines from sc_sn data
  sc_sn_data <- sc_sn_data[sc_sn_data$cell_line != "All", ]
  
  all_proportional_data <- data.frame()
  
  for (combo in combinations) {
    tech <- combo$tech
    data_type <- combo$data_type
    combo_name <- combo$name
    
    cat(sprintf("  Processing proportional data for %s...\n", combo_name))
    
    # Select the appropriate dataset
    if (data_type == "bulk") {
      current_data <- bulk_data
    } else {
      current_data <- sc_sn_data
    }
    
    # Filter for this technology, data type, and proportional measurement
    filtered_data <- current_data[
      current_data$technology == tech & 
      current_data$data_type == data_type & 
      current_data$measurement_type == "proportional" &
      current_data$cell_line %in% cell_lines, 
    ]
    
    if (nrow(filtered_data) > 0) {
      cat(sprintf("    Found %d proportional entries\n", nrow(filtered_data)))
      
      # Create standardized data frame for each entry
      for (i in 1:nrow(filtered_data)) {
        row <- filtered_data[i, ]
        
        # Only include entries with valid valid_best_aln
        if (!is.na(row$valid_best_aln) && !is.null(row$valid_best_aln) && row$valid_best_aln != "") {
          new_row <- data.frame(
            combination = combo_name,
            technology = tech,
            data_type = data_type,
            cell_line = row$cell_line,
            sampling_depth = row$sampling_depth,
            measurement_type = "proportional",
            valid_best_aln = as.numeric(row$valid_best_aln),
            expressed_transcripts = as.numeric(row$expressed_transcripts),
            total_transcripts = as.numeric(row$total_transcripts),
            expression_rate = as.numeric(row$expression_rate_percent),
            sample_id = row$sample_id,
            stringsAsFactors = FALSE
          )
          all_proportional_data <- rbind(all_proportional_data, new_row)
        }
      }
    } else {
      cat(sprintf("    No proportional data found for %s\n", combo_name))
    }
  }
  
  cat(sprintf("Total proportional data points extracted: %d\n\n", nrow(all_proportional_data)))
  return(all_proportional_data)
}

# Extract rarefaction data for all combinations
extract_comprehensive_rarefaction <- function(combinations, cell_lines) {
  all_rarefaction_data <- data.frame()
  
  # First, extract proportional data directly from sample sheets
  proportional_data <- extract_proportional_data_direct(combinations, cell_lines)
  all_rarefaction_data <- rbind(all_rarefaction_data, proportional_data)
  
  # Then, extract absolute data using the existing DataNav method
  cat("Extracting absolute data using DataNav...\n")
  
  for (combo in combinations) {
    tech <- combo$tech
    data_type <- combo$data_type
    combo_name <- combo$name
    
    cat(sprintf("Processing absolute data for %s...\n", combo_name))
    
    # Get all available depths for absolute measurements only
    available_meta <- DataNav$get(
      technology = tech,
      data_type = data_type,
      measurement_type = "absolute",  # Only absolute for this method
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
            measurement_type = "absolute",  # Only absolute
            sampling_depth = depth,
            return_meta = TRUE
          )
          
          if (nrow(meta) > 0) {
            # Check if valid_best_aln column exists and has a valid value
            if ("valid_best_aln" %in% names(meta) && 
                !is.null(meta$valid_best_aln[1]) && 
                !is.na(meta$valid_best_aln[1])) {
              
              new_row <- data.frame(
                combination = combo_name,
                technology = tech,
                data_type = data_type,
                cell_line = cell_line,
                sampling_depth = depth,
                measurement_type = "absolute",
                valid_best_aln = meta$valid_best_aln[1],
                expressed_transcripts = meta$expressed_transcripts[1],
                total_transcripts = meta$total_transcripts[1],
                expression_rate = meta$expression_rate_percent[1],
                sample_id = meta$sample_id[1],
                stringsAsFactors = FALSE
              )
              all_rarefaction_data <- rbind(all_rarefaction_data, new_row)
            } else {
              cat(sprintf("    Warning: %s %s %s absolute - missing valid_best_aln data\n", 
                         tech, cell_line, depth))
            }
          }
        }
      }
      
      # Count absolute data points for this combination
      combo_absolute_points <- nrow(all_rarefaction_data[
        all_rarefaction_data$combination == combo_name & 
        all_rarefaction_data$measurement_type == "absolute", 
      ])
      cat(sprintf("  • Absolute data points collected: %d\n", combo_absolute_points))
    }
    cat("\n")
  }
  
  return(all_rarefaction_data)
}

# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

# Create comprehensive rarefaction plot
plot_comprehensive_rarefaction <- function(data, title = "Comprehensive Rarefaction Curves", filename = NULL) {
  
  # Sort data to ensure proper line connections
  data <- data %>%
    arrange(combination, cell_line, valid_best_aln, expressed_transcripts)
  
  # Ensure a consistent legend/draw order for combinations
  combination_levels <- c(
    "dRNA Bulk",
    "ONT Bulk",
    "ONT Single-Cell",
    "ONT Single-Nucleus",
    "PacBio Bulk",
    "PacBio Single-Cell",
    "PacBio Single-Nucleus"
  )
  data <- data %>% mutate(combination = factor(combination, levels = combination_levels))
  
  # Color mapping for technologies
  color_mapping <- c(
    "dRNA Bulk" = "#00A651",      # Green for ONT dRNA
    "ONT Bulk" = "#003F5C",       # Deep blue for ONT
    "ONT Single-Cell" = "#003F5C", # Deep blue for ONT
    "ONT Single-Nucleus" = "#003F5C", # Deep blue for ONT
    "PacBio Bulk" = "#DA1884",    # Magenta for PacBio
    "PacBio Single-Cell" = "#DA1884", # Magenta for PacBio
    "PacBio Single-Nucleus" = "#DA1884" # Magenta for PacBio
  )
  
  # Line type mapping for sample types
  linetype_mapping <- c(
    "dRNA Bulk" = "solid",
    "ONT Bulk" = "solid",
    "ONT Single-Cell" = "dashed",
    "ONT Single-Nucleus" = "dotted", 
    "PacBio Bulk" = "solid",
    "PacBio Single-Cell" = "dashed",
    "PacBio Single-Nucleus" = "dotted"
  )
  
  # Create the plot with better grouping
  p <- ggplot(data, aes(x = valid_best_aln, y = expressed_transcripts)) +
    geom_point(aes(color = combination, shape = cell_line), 
               size = 1.6, alpha = 0.55) +
    geom_smooth(aes(color = combination, linetype = combination, 
                    group = interaction(combination, cell_line)), 
                method = "loess", se = FALSE, span = 0.3, linewidth = 1.4, alpha = 0.9) +
    scale_x_continuous(
      name = "# Valid Aligned Reads\n(oarfish 'valid_best_aln' reads)",
      breaks = c(100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000, 20000000, 50000000),
      labels = c("100K", "200K", "500K", "1M", "2M", "5M", "10M", "20M", "50M"),
      trans = "log10"
    ) +
    scale_y_continuous(
      name = "Expressed Transcripts Discovered (count > 1)",
      breaks = scales::log_breaks(n = 7),
      labels = scales::comma,
      trans = "log10"
    ) +
    scale_color_manual(
      name = "Protocol/Data Type",
      values = color_mapping
    ) +
    scale_linetype_manual(
      name = "Protocol/Data Type", 
      values = linetype_mapping
    ) +
    scale_shape_manual(
      name = "Cell Line",
      values = c(16, 17, 15, 18, 8, 11, 12, 13)
    ) +
    labs(
      title = title,
      subtitle = paste("Using actual aligned read counts • All measurement types • Smoothed curves •", 
                      length(unique(data$combination)), "protocol combinations •",
                      length(unique(data$cell_line)), "cell lines")
    ) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_line(color = "#F0F0F0"),
      panel.grid.minor = element_line(color = "#F5F5F5")
    ) +
    guides(
      color = guide_legend(override.aes = list(shape = 16, size = 3)),
      linetype = guide_legend(override.aes = list(size = 1.5)),
      shape = guide_legend(override.aes = list(color = "black", size = 3))
    )
  
  # Save plot if filename provided
  if (!is.null(filename)) {
    # Create directory if it doesn't exist
    plot_dir <- dirname(filename)
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir, recursive = TRUE)
      cat(sprintf("Created directory: %s\n", plot_dir))
    }
    ggsave(filename, plot = p, width = 16, height = 10, dpi = 300)
  }
  
  return(p)
}

# Create grid plot by cell line
plot_grid_by_cell_line <- function(data, title = "Rarefaction Curves by Cell Line", filename = NULL) {
  
  # Sort data to ensure proper line connections
  data <- data %>%
    arrange(combination, cell_line, valid_best_aln, expressed_transcripts)
  
  # Ensure a consistent legend/draw order for combinations
  combination_levels <- c(
    "dRNA Bulk",
    "ONT Bulk",
    "ONT Single-Cell",
    "ONT Single-Nucleus",
    "PacBio Bulk",
    "PacBio Single-Cell",
    "PacBio Single-Nucleus"
  )
  data <- data %>% mutate(combination = factor(combination, levels = combination_levels))
  
  # Color mapping for technologies
  color_mapping <- c(
    "dRNA Bulk" = "#00A651",      # Green for ONT dRNA
    "ONT Bulk" = "#003F5C",       # Deep blue for ONT
    "ONT Single-Cell" = "#003F5C", # Deep blue for ONT
    "ONT Single-Nucleus" = "#003F5C", # Deep blue for ONT
    "PacBio Bulk" = "#DA1884",    # Magenta for PacBio
    "PacBio Single-Cell" = "#DA1884", # Magenta for PacBio
    "PacBio Single-Nucleus" = "#DA1884" # Magenta for PacBio
  )
  
  # Line type mapping for sample types
  linetype_mapping <- c(
    "dRNA Bulk" = "solid",
    "ONT Bulk" = "solid",
    "ONT Single-Cell" = "dashed",
    "ONT Single-Nucleus" = "dotted", 
    "PacBio Bulk" = "solid",
    "PacBio Single-Cell" = "dashed",
    "PacBio Single-Nucleus" = "dotted"
  )
  
  # Create the grid plot with smoothed curves
  p <- ggplot(data, aes(x = valid_best_aln, y = expressed_transcripts)) +
    geom_point(aes(color = combination), 
               size = 1.2, alpha = 0.55) +
    geom_smooth(aes(color = combination, linetype = combination,
                    group = combination), 
                method = "loess", se = FALSE, span = 0.4, linewidth = 1.2, alpha = 0.9) +
    scale_x_continuous(
      name = "# Valid Aligned Reads\n(oarfish 'valid_best_aln' reads)",
      breaks = c(100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000, 20000000, 50000000),
      labels = c("100K", "200K", "500K", "1M", "2M", "5M", "10M", "20M", "50M"),
      trans = "log10"
    ) +
    scale_y_continuous(
      name = "Expressed Transcripts Discovered",
      breaks = scales::log_breaks(n = 7),
      labels = scales::comma,
      trans = "log10"
    ) +
    scale_color_manual(
      name = "Protocol/Data Type",
      values = color_mapping
    ) +
    scale_linetype_manual(
      name = "Protocol/Data Type", 
      values = linetype_mapping
    ) +
    facet_wrap(~ cell_line, scales = "free", ncol = 2) +
    labs(
      title = title,
      subtitle = paste("Each panel shows one cell line • All measurement types • Smoothed curves •", 
                      length(unique(data$combination)), "protocol combinations •",
                      length(unique(data$cell_line)), "cell lines")
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      strip.text = element_text(size = 10, face = "bold"),
      panel.grid.major = element_line(color = "#F0F0F0"),
      panel.grid.minor = element_line(color = "#F5F5F5"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9)
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 2), ncol = 4),
      linetype = guide_legend(override.aes = list(size = 1), ncol = 4)
    )
  
  # Save plot if filename provided
  if (!is.null(filename)) {
    # Create directory if it doesn't exist
    plot_dir <- dirname(filename)
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir, recursive = TRUE)
      cat(sprintf("Created directory: %s\n", plot_dir))
    }
    
    # Calculate appropriate dimensions based on number of cell lines
    n_cell_lines <- length(unique(data$cell_line))
    n_cols <- 2
    n_rows <- ceiling(n_cell_lines / n_cols)
    
    plot_width <- 12
    plot_height <- max(8, 3 + (n_rows * 3))
    
    ggsave(filename, plot = p, width = plot_width, height = plot_height, dpi = 300)
    cat(sprintf("Grid plot saved: %s (%d cell lines, %.1f x %.1f inches)\n", 
                filename, n_cell_lines, plot_width, plot_height))
  }
  
  return(p)
}

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

cat("Extracting rarefaction data for all 7 protocol/data type combinations...\n\n")

# Extract comprehensive rarefaction data
comprehensive_data <- extract_comprehensive_rarefaction(COMBINATIONS, CELL_LINES)

cat("=== DATA EXTRACTION SUMMARY ===\n")
cat(sprintf("Total data points collected: %d\n", nrow(comprehensive_data)))
cat(sprintf("Combinations: %d\n", length(unique(comprehensive_data$combination))))
cat(sprintf("Cell lines: %d\n", length(unique(comprehensive_data$cell_line))))
cat(sprintf("Unique depths: %d\n", length(unique(comprehensive_data$sampling_depth))))

# Show extraction method breakdown
extraction_method_summary <- comprehensive_data %>%
  group_by(measurement_type) %>%
  summarise(
    data_points = n(),
    combinations = length(unique(combination)),
    cell_lines = length(unique(cell_line)),
    unique_depths = length(unique(sampling_depth)),
    read_count_range = paste(scales::number(range(valid_best_aln), scale = 1e-6, suffix = "M", accuracy = 0.1), collapse = " - "),
    .groups = 'drop'
  )

cat(sprintf("\nExtraction method results:\n"))
cat(sprintf("  • Proportional data (direct CSV): %d points\n", 
            sum(comprehensive_data$measurement_type == "proportional")))
cat(sprintf("  • Absolute data (DataNav): %d points\n", 
            sum(comprehensive_data$measurement_type == "absolute")))

# Show breakdown by combination
combination_summary <- comprehensive_data %>%
  group_by(combination, measurement_type) %>%
  summarise(
    data_points = n(),
    depths = length(unique(sampling_depth)),
    read_count_range = paste(scales::number(range(valid_best_aln), scale = 1e-6, suffix = "M", accuracy = 0.1), collapse = " - "),
    .groups = 'drop'
  )

cat("\nData points by combination and measurement type:\n")
print(combination_summary)

# Show overall measurement type breakdown
cat("\nOverall breakdown by measurement type:\n")
print(extraction_method_summary)

# Generate comprehensive plot
cat("\n=== GENERATING COMPREHENSIVE RAREFACTION PLOT ===\n")

# Create the comprehensive plot
comprehensive_plot <- plot_comprehensive_rarefaction(
  comprehensive_data,
  "Comprehensive Rarefaction Curves - All Protocols",
  "analysis_scripts/plots/comprehensive_rarefaction_all_protocols.png"
)

# Generate grid plot by cell line
cat("\n=== GENERATING GRID PLOT BY CELL LINE ===\n")

# Create the grid plot
grid_plot <- plot_grid_by_cell_line(
  comprehensive_data,
  "Rarefaction Curves by Cell Line - All Protocols",
  "analysis_scripts/plots/rarefaction_curves_by_cell_line.pdf"
)

# Save the data
saveRDS(list(
  comprehensive_data = comprehensive_data,
  combination_summary = combination_summary,
  combinations_analyzed = COMBINATIONS,
  comprehensive_plot = comprehensive_plot,
  grid_plot = grid_plot
), "analysis_scripts/comprehensive_rarefaction_data.rds")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Generated plots:\n")
cat("  📊 Comprehensive plot: analysis_scripts/plots/comprehensive_rarefaction_all_protocols.png\n")
cat("  📋 Grid plot by cell line: analysis_scripts/plots/rarefaction_curves_by_cell_line.pdf\n")
cat("  💾 Data: analysis_scripts/comprehensive_rarefaction_data.rds\n\n")

cat("🎯 COMPREHENSIVE RAREFACTION CURVES COMPLETE! 🎯\n")
cat("All 7 protocol/data type combinations plotted using actual aligned read counts!\n")
cat("Enhanced data extraction: direct CSV querying for proportional + DataNav for absolute data.\n")
cat("Includes both absolute and proportional measurement types for maximum data coverage.\n")
cat("Using smoothed curves to ensure proper rarefaction behavior (no backwards lines).\n")
cat("Grid plot shows each cell line in a separate panel for detailed comparison.\n")

# Display the plots
cat("\n=== DISPLAYING PLOTS ===\n")
print(comprehensive_plot)
cat("\nGrid plot by cell line:\n")
print(grid_plot)
