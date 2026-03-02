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

# Focus on key cell lines for cleaner visualization
CELL_LINES <- c("H146", "H1975", "H69", "H211")

# Convert sampling depth to thousands of reads for plotting
depth_to_thousands <- function(depth_string) {
  if (depth_string == "full") return(100000)  # Assign high value for "full"
  numeric_value <- as.numeric(gsub("M", "", depth_string))
  return(numeric_value * 1000)  # Convert millions to thousands (e.g., 0.1M = 100K)
}

# =============================================================================
# DATA EXTRACTION
# =============================================================================

# Extract rarefaction data for all combinations
extract_comprehensive_rarefaction <- function(combinations, cell_lines) {
  all_rarefaction_data <- data.frame()
  
  for (combo in combinations) {
    tech <- combo$tech
    data_type <- combo$data_type
    combo_name <- combo$name
    
    cat(sprintf("Extracting data for %s...\n", combo_name))
    
    # Get all available depths for this specific combination
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
      
      # Extract data for each cell line and depth
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
            new_row <- data.frame(
              combination = combo_name,
              technology = tech,
              data_type = data_type,
              cell_line = cell_line,
              sampling_depth = depth,
              depth_thousands = depth_to_thousands(depth),
              expressed_transcripts = meta$expressed_transcripts[1],
              total_transcripts = meta$total_transcripts[1],
              expression_rate = meta$expression_rate_percent[1],
              sample_id = meta$sample_id[1],
              stringsAsFactors = FALSE
            )
            all_rarefaction_data <- rbind(all_rarefaction_data, new_row)
          }
        }
      }
      
      # Count data points for this combination
      combo_points <- nrow(all_rarefaction_data[all_rarefaction_data$combination == combo_name, ])
      cat(sprintf("  • Data points collected: %d\n", combo_points))
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
  
  # Remove "full" depth for cleaner visualization
  plot_data <- data[data$sampling_depth != "full", ]
  
  # Create color and line type mappings
  color_mapping <- c(
    "dRNA Bulk" = "#E31A1C",
    "ONT Bulk" = "#1F78B4", 
    "ONT Single-Cell" = "#33A02C",
    "ONT Single-Nucleus" = "#FF7F00",
    "PacBio Bulk" = "#6A3D9A",
    "PacBio Single-Cell" = "#B15928",
    "PacBio Single-Nucleus" = "#FB9A99"
  )
  
  linetype_mapping <- c(
    "dRNA Bulk" = "solid",
    "ONT Bulk" = "solid",
    "ONT Single-Cell" = "dashed",
    "ONT Single-Nucleus" = "dotted", 
    "PacBio Bulk" = "solid",
    "PacBio Single-Cell" = "dashed",
    "PacBio Single-Nucleus" = "dotted"
  )
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = depth_thousands, y = expressed_transcripts)) +
    geom_line(aes(color = combination, linetype = combination, group = interaction(combination, cell_line)), 
              linewidth = 1.2, alpha = 0.8) +
    geom_point(aes(color = combination, shape = cell_line), size = 2, alpha = 0.7) +
    scale_x_continuous(
      name = "Total Reads (Thousands)",
      breaks = c(100, 200, 500, 1000, 5000, 10000, 20000, 30000, 60000),
      labels = c("100K", "200K", "500K", "1M", "5M", "10M", "20M", "30M", "60M"),
      trans = "log10"
    ) +
    scale_y_continuous(
      name = "Expressed Transcripts Discovered (count > 1)",
      labels = comma,
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
      values = c(16, 17, 15, 18)
    ) +
    labs(
      title = title,
      subtitle = paste("All available sampling depths •", 
                      length(unique(plot_data$combination)), "protocol combinations •",
                      length(unique(plot_data$cell_line)), "cell lines")
    ) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    ) +
    guides(
      color = guide_legend(override.aes = list(shape = 16, size = 3)),
      linetype = guide_legend(override.aes = list(linewidth = 1.5)),
      shape = guide_legend(override.aes = list(color = "black", size = 3))
    )
  
  # Save plot if filename provided
  if (!is.null(filename)) {
    ggsave(filename, plot = p, width = 16, height = 10, dpi = 300)
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

# Show breakdown by combination
combination_summary <- comprehensive_data %>%
  group_by(combination) %>%
  summarise(
    data_points = n(),
    depths = length(unique(sampling_depth)),
    depth_range = paste(range(depth_thousands), collapse = " - "),
    .groups = 'drop'
  )

cat("\nData points by combination:\n")
print(combination_summary)

# Generate comprehensive plots
cat("\n=== GENERATING COMPREHENSIVE RAREFACTION PLOTS ===\n")

# Log scale plot (better for wide range)
log_plot <- plot_comprehensive_rarefaction(
  comprehensive_data,
  "Comprehensive Rarefaction Curves - All Protocols (Log Scale)",
  "analysis_scripts/plots/comprehensive_rarefaction_log.png"
)

# Save the data
saveRDS(list(
  comprehensive_data = comprehensive_data,
  combination_summary = combination_summary,
  combinations_analyzed = COMBINATIONS
), "analysis_scripts/comprehensive_rarefaction_data.rds")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Generated plots:\n")
cat("  📊 Log scale: analysis_scripts/plots/comprehensive_rarefaction_log.png\n")
cat("  💾 Data: analysis_scripts/comprehensive_rarefaction_data.rds\n\n")

cat("🎯 COMPREHENSIVE RAREFACTION CURVES COMPLETE! 🎯\n")
cat("All 7 protocol/data type combinations plotted with ALL available sampling depths!\n")

# Display the plot
print(log_plot) 