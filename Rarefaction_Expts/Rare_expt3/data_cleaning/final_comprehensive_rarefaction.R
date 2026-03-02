# Final comprehensive rarefaction plot - ALL 7 combinations, ALL depths
library(ggplot2)
library(dplyr)
library(scales)
source('utils/data_navigator.R')

cat('🧬 FINAL COMPREHENSIVE RAREFACTION ANALYSIS 🧬\n')
cat('All 7 protocol combinations with ALL available sampling depths\n\n')

# All 7 combinations
combinations <- list(
  c('dRNA', 'bulk', 'dRNA Bulk'),
  c('ONT', 'bulk', 'ONT Bulk'),
  c('ONT', 'single_cell', 'ONT Single-Cell'),
  c('ONT', 'single_nucleus', 'ONT Single-Nucleus'),
  c('PacBio', 'bulk', 'PacBio Bulk'),
  c('PacBio', 'single_cell', 'PacBio Single-Cell'),
  c('PacBio', 'single_nucleus', 'PacBio Single-Nucleus')
)

# Use all 4 cell lines
cell_lines <- c('H146', 'H1975', 'H69', 'H211')

# Extract ALL data
comprehensive_data <- data.frame()

for (combo in combinations) {
  tech <- combo[1]
  dtype <- combo[2] 
  name <- combo[3]
  
  cat(sprintf('Extracting %s...\n', name))
  
  # Get ALL available depths for this combination
  meta <- DataNav$get(technology=tech, data_type=dtype, measurement_type='absolute', return_meta=TRUE)
  if (nrow(meta) > 0) {
    all_depths <- unique(meta$sampling_depth)
    all_depths <- all_depths[all_depths != 'full']  # Remove 'full' for cleaner plot
    
    # Sort depths numerically
    numeric_depths <- as.numeric(gsub('M', '', all_depths))
    all_depths <- all_depths[order(numeric_depths)]
    
    cat(sprintf('  • %d depths: %s\n', length(all_depths), paste(all_depths, collapse=', ')))
    
    # Extract for all cell lines and all depths
    for (cell_line in cell_lines) {
      for (depth in all_depths) {
        sample_meta <- DataNav$get(
          technology = tech,
          cell_line = cell_line,
          data_type = dtype,
          measurement_type = 'absolute',
          sampling_depth = depth,
          return_meta = TRUE
        )
        
        if (nrow(sample_meta) > 0) {
          comprehensive_data <- rbind(comprehensive_data, data.frame(
            combination = name,
            technology = tech,
            data_type = dtype,
            cell_line = cell_line,
            depth = depth,
            depth_numeric = as.numeric(gsub('M', '', depth)),
            expressed_transcripts = sample_meta$expressed_transcripts[1],
            total_transcripts = sample_meta$total_transcripts[1],
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
}

cat(sprintf('\nTotal data points: %d\n', nrow(comprehensive_data)))

# Create the comprehensive plot
color_mapping <- c(
  'dRNA Bulk' = '#E31A1C',
  'ONT Bulk' = '#1F78B4', 
  'ONT Single-Cell' = '#33A02C',
  'ONT Single-Nucleus' = '#FF7F00',
  'PacBio Bulk' = '#6A3D9A',
  'PacBio Single-Cell' = '#B15928',
  'PacBio Single-Nucleus' = '#FB9A99'
)

linetype_mapping <- c(
  'dRNA Bulk' = 'solid',
  'ONT Bulk' = 'solid',
  'ONT Single-Cell' = 'dashed',
  'ONT Single-Nucleus' = 'dotted', 
  'PacBio Bulk' = 'solid',
  'PacBio Single-Cell' = 'dashed',
  'PacBio Single-Nucleus' = 'dotted'
)

# Create comprehensive plot
p <- ggplot(comprehensive_data, aes(x = depth_numeric, y = expressed_transcripts)) +
  geom_line(aes(color = combination, linetype = combination, group = interaction(combination, cell_line)), 
            size = 1.2, alpha = 0.8) +
  geom_point(aes(color = combination, shape = cell_line), size = 2, alpha = 0.7) +
  scale_x_continuous(
    name = 'Total Reads (Million)',
    breaks = c(0.1, 0.5, 1, 5, 10, 20, 30, 60),
    labels = c('0.1M', '0.5M', '1M', '5M', '10M', '20M', '30M', '60M'),
    trans = 'log10'
  ) +
  scale_y_continuous(
    name = 'Expressed Transcripts Discovered (count > 1)',
    labels = comma,
    trans = 'log10'
  ) +
  scale_color_manual(
    name = 'Protocol/Data Type',
    values = color_mapping
  ) +
  scale_linetype_manual(
    name = 'Protocol/Data Type', 
    values = linetype_mapping
  ) +
  scale_shape_manual(
    name = 'Cell Line',
    values = c(16, 17, 15, 18)
  ) +
  labs(
    title = 'Comprehensive Rarefaction Curves - All Protocols & Data Types',
    subtitle = 'Transcript discovery across all available sampling depths (absolute subsampled data)'
  ) +
  theme_minimal() +
  theme(
    legend.position = 'bottom',
    legend.box = 'horizontal',
    plot.title = element_text(size = 16, face = 'bold'),
    plot.subtitle = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  guides(
    color = guide_legend(override.aes = list(shape = 16, size = 3)),
    linetype = guide_legend(override.aes = list(size = 1.5)),
    shape = guide_legend(override.aes = list(color = 'black', size = 3))
  )

# Save the plot
ggsave('analysis_scripts/plots/FINAL_comprehensive_rarefaction_all_protocols.png', 
       plot = p, width = 16, height = 10, dpi = 300)

# Save the data
saveRDS(comprehensive_data, 'analysis_scripts/FINAL_comprehensive_rarefaction_data.rds')

cat('\n=== ANALYSIS COMPLETE ===\n')
cat('📊 Plot: analysis_scripts/plots/FINAL_comprehensive_rarefaction_all_protocols.png\n')
cat('💾 Data: analysis_scripts/FINAL_comprehensive_rarefaction_data.rds\n')

# Show summary by combination
summary_stats <- comprehensive_data %>%
  group_by(combination) %>%
  summarise(
    data_points = n(),
    unique_depths = length(unique(depth_numeric)),
    min_depth = min(depth_numeric),
    max_depth = max(depth_numeric),
    .groups = 'drop'
  )

cat('\nSummary by combination:\n')
print(summary_stats)

cat('\n🎯 COMPREHENSIVE RAREFACTION COMPLETE! 🎯\n')
cat('All 7 protocol/data type combinations with ALL available sampling depths!\n')
