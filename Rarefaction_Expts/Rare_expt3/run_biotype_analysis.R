#!/usr/bin/env Rscript

# =============================================================================
# RUN BIOTYPE DISCOVERY ANALYSIS - ALL BIOTYPES
# =============================================================================

# Source the comprehensive biotype discovery script
source("analysis_scripts/comprehensive_biotype_discovery.R")

cat("🚀 Running ALL Biotypes Discovery Analysis...\n")
cat("This will generate 3 types of plots showing all 9 individual biotype categories:\n")
cat("  📊 Comprehensive plot (all protocols & cell lines)\n")
cat("  📋 Grid plot (protocols in panels)\n") 
cat("  📄 Cell line plots (separate pages per cell line)\n\n")

cat("🧬 ANALYZING ALL INDIVIDUAL BIOTYPE CATEGORIES:\n")
cat("  • Protein Coding\n")
cat("  • Long Non-Coding RNA\n")
cat("  • Mitochondrial RNA\n")
cat("  • Nonsense Mediated Decay\n")
cat("  • Other\n")
cat("  • Processed Transcript\n")
cat("  • Pseudogenes\n")
cat("  • Retained Intron\n")
cat("  • Short Non-Coding RNA\n\n")

# Run ALL biotype analysis (9 individual categories)
cat("=== RUNNING ALL BIOTYPE ANALYSIS ===\n")
result <- run_all_biotype_analysis(save_plots = TRUE)

# Summary
if (!is.null(result)) {
  cat("\n✅ ALL BIOTYPES ANALYSIS COMPLETED SUCCESSFULLY!\n")
  cat("📁 Check the analysis_scripts/plots/ directory for output files.\n")
  
  cat(sprintf("\n📊 Data Summary:\n"))
  cat(sprintf("  • Total data points: %s\n", format(nrow(result$data), big.mark = ",")))
  cat(sprintf("  • Protocols: %d\n", length(unique(result$data$combination))))
  cat(sprintf("  • Cell lines: %d\n", length(unique(result$data$cell_line))))
  cat(sprintf("  • Biotype categories: %d individual categories\n", 
              length(unique(result$data$biotype_category))))
  
  cat("\n📊 Individual biotype categories found:\n")
  for (biotype in sort(unique(result$data$biotype_category))) {
    cat(sprintf("  • %s\n", biotype))
  }
  
  cat("\n📄 Cell line plots generated for:\n")
  for (cell_line in names(result$cell_line_plots)) {
    cat(sprintf("  • %s\n", cell_line))
  }
  
  cat("\n📁 FILES GENERATED:\n")
  cat("  • biotype_discovery_comprehensive_all_[timestamp].png\n")
  cat("  • biotype_discovery_grid_all_[timestamp].pdf\n")
  cat("  • biotype_discovery_by_cell_line_all_[timestamp].pdf ← Cell line pages\n")
  
} else {
  cat("\n❌ Analysis failed - check for errors above.\n")
} 