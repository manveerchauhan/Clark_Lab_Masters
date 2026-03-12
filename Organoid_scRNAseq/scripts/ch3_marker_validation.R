## Marker Gene Validation Script - Refined Version
# Author: Manveer Chauhan
# This script validates developmental stage-specific cell type marker genes across
# integrated organoid timepoints using module scores, individual gene expression UMAPs,
# and consensus cell type-based dotplots

# Load essential packages
library(Seurat)
library(tidyverse)
library(patchwork)
library(gridExtra)
library(grid)

# Set working directory
setwd("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space")

# Create output directory
dir.create("./output_files/marker_validation", recursive = TRUE, showWarnings = FALSE)

cat("\n========================================\n")
cat("MARKER GENE VALIDATION ANALYSIS\n")
cat("Refined Version - Developmental Stage Specific\n")
cat("========================================\n\n")

## DEFINE TIMEPOINT-SPECIFIC MARKER GENE SETS------
# Organized hierarchically: Timepoint → Category → Marker sets
# Reflects developmental progression and stage-appropriate biology

markers_by_timepoint <- list(
  "1M_Org" = list(
    "Progenitors" = list(
      "Dorsal_Forebrain_Progenitors" = c("EMX1", "PAX6", "FOXG1", "SOX2")
    ),
    "Neurons" = list(
      "Pan_Neuronal" = c("MAP2")
    )
  ),

  "3M_Org" = list(
    "Progenitors" = list(
      "Radial_Glia" = c("PAX6", "SOX2"),
      "Outer_Radial_Glia" = c("HOPX", "TNC"),
      "Intermediate_Progenitors" = c("EOMES", "PPP1R17", "TMEM158"),
      "Cycling_Progenitors" = c("MKI67", "TOP2A", "BIRC5", "PCNA", "MCM6")
    ),
    "Projection_Neurons" = list(
      "Immature_Projection_Neurons" = c("TBR1"),
      "Callosal_Projection_Neurons" = c("SATB2", "INHBA", "FRMD4B"),
      "Corticothalamic_Projection_Neurons" = c("ZFPM2", "TBR1"),
      "Corticofugal_Projection_Neurons" = c("BCL11B", "CHRM1", "CHRM2", "CHRM4", "CHRM5", "TLE4")
    ),
    "Interneurons" = list(
      "Interneurons" = c("DLX1", "DLX2", "GAD2", "LGALS3"),
      "Cajal_Retzius" = c("RELN")
    )
  ),

  "6M_Org" = list(
    "Progenitors" = list(
      "Outer_Radial_Glia" = c("HOPX", "TNC")
    ),
    "Glia" = list(
      "Astroglia" = c("GFAP", "S100B")
    ),
    "Synaptic_Maturation" = list(
      "Synaptic_Markers" = c("SLC17A7", "DLG4")  # VGLUT1, PSD95
    ),
    "Oligodendrocyte_Lineage" = list(
      "OPC_Core" = c("PDGFRA", "OLIG2", "SOX10", "CSPG4"),
      "Committed_Precursors" = c("GPR17", "BCAS1"),
      "Mature_Oligodendrocytes" = c("MBP", "PLP1", "MOG")
    )
  )
)

# Define gene aliases for common synonyms
gene_aliases <- list(
  "VGLUT1" = "SLC17A7",
  "PSD95" = "DLG4",
  "CTIP2" = "BCL11B",
  "FOG2" = "ZFPM2"
)

cat("Defined timepoint-specific marker gene sets:\n")
for(tp in names(markers_by_timepoint)) {
  n_categories <- length(markers_by_timepoint[[tp]])
  n_marker_sets <- sum(sapply(markers_by_timepoint[[tp]], length))
  n_genes <- sum(sapply(markers_by_timepoint[[tp]], function(cat) {
    sum(sapply(cat, length))
  }))
  cat("  ", tp, ": ", n_categories, " categories, ", n_marker_sets,
      " marker sets, ", n_genes, " genes\n", sep = "")
}

## FUNCTION 1: Resolve gene symbol to full feature name with aliases------
resolve_gene_name_with_aliases <- function(gene_symbol, seurat_obj, verbose = FALSE) {

  # First try direct match
  pattern <- paste0('-', gene_symbol, '$')
  matches <- grep(pattern, rownames(seurat_obj), value = TRUE, ignore.case = TRUE)

  if(length(matches) > 0) {
    if(length(matches) == 1) {
      if(verbose) cat("    ", gene_symbol, " → ", matches[1], "\n", sep = "")
      return(matches[1])
    } else {
      warning("Multiple matches for ", gene_symbol, ": ", paste(matches, collapse = ", "))
      if(verbose) cat("    ", gene_symbol, " → ", matches[1], " (multiple matches)\n", sep = "")
      return(matches[1])
    }
  }

  # Try alias if no direct match
  if(gene_symbol %in% names(gene_aliases)) {
    alias <- gene_aliases[[gene_symbol]]
    if(verbose) cat("    Trying alias: ", gene_symbol, " → ", alias, "\n", sep = "")
    pattern_alias <- paste0('-', alias, '$')
    matches_alias <- grep(pattern_alias, rownames(seurat_obj), value = TRUE, ignore.case = TRUE)

    if(length(matches_alias) > 0) {
      if(verbose) cat("    ", gene_symbol, " (", alias, ") → ", matches_alias[1], "\n", sep = "")
      return(matches_alias[1])
    }
  }

  # No matches found
  if(verbose) cat("    ", gene_symbol, " → NOT FOUND\n", sep = "")
  return(NA)
}

## FUNCTION 2: Validate gene set against Seurat object------
validate_gene_set <- function(gene_symbols, seurat_obj, set_name, verbose = TRUE) {

  if(verbose) cat("\n  Validating", set_name, ":", length(gene_symbols), "genes\n")

  validation_results <- data.frame(
    gene_symbol = gene_symbols,
    full_name = NA_character_,
    found = FALSE,
    stringsAsFactors = FALSE
  )

  for(i in seq_along(gene_symbols)) {
    full_name <- resolve_gene_name_with_aliases(gene_symbols[i], seurat_obj, verbose = FALSE)
    validation_results$full_name[i] <- full_name
    validation_results$found[i] <- !is.na(full_name)
  }

  n_found <- sum(validation_results$found)
  n_missing <- sum(!validation_results$found)
  pct_found <- round(100 * n_found / length(gene_symbols), 1)

  if(verbose) cat("    Found:", n_found, "/", length(gene_symbols), "(", pct_found, "%)\n")

  if(n_missing > 0 && verbose) {
    missing_genes <- validation_results$gene_symbol[!validation_results$found]
    cat("    Missing:", paste(missing_genes, collapse = ", "), "\n")
  }

  return(validation_results)
}

## FUNCTION 3: Calculate top expressed genes from a set------
get_top_genes <- function(gene_full_names, seurat_obj, n = 5) {

  # Remove NAs
  valid_genes <- gene_full_names[!is.na(gene_full_names)]

  if(length(valid_genes) == 0) {
    return(character(0))
  }

  # Calculate mean expression for each gene (Seurat v5 compatible)
  mean_expr <- sapply(valid_genes, function(gene) {
    mean(expm1(GetAssayData(seurat_obj, assay = "RNA", layer = "data")[gene, ]))
  })

  # Sort and take top n
  top_genes <- names(sort(mean_expr, decreasing = TRUE))[1:min(n, length(mean_expr))]

  return(top_genes)
}

## FUNCTION 4: Create marker set validation plot------
create_marker_plot <- function(seurat_obj, marker_set_name, validation_results,
                               module_score_name, timepoint_name) {

  cat("\n  Creating plots for", marker_set_name, "\n")

  # Get valid genes
  valid_genes <- validation_results$full_name[validation_results$found]
  missing_genes <- validation_results$gene_symbol[!validation_results$found]

  # Check if we have enough genes
  if(length(valid_genes) == 0) {
    cat("    WARNING: No valid genes found - creating empty plot\n")
    empty_plot <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste0(marker_set_name, "\n\nNo valid genes found"),
               size = 6, color = "gray50") +
      theme_void()

    gene_list_text <- paste(validation_results$gene_symbol, collapse = ", ")
    return(list(
      title_info = list(
        marker_set_name = marker_set_name,
        gene_list = gene_list_text,
        missing_genes = missing_genes
      ),
      module_plot = empty_plot,
      cluster_plot = empty_plot,
      gene_plots = list()
    ))
  }

  if(length(valid_genes) < length(validation_results$gene_symbol) * 0.5) {
    cat("    WARNING: <50% genes found - proceeding with caution\n")
  }

  # Get top 5 genes by expression
  top_genes <- get_top_genes(valid_genes, seurat_obj, n = 5)
  cat("    Top genes:", paste(gsub(".*-", "", top_genes), collapse = ", "), "\n")

  # Determine max expression for consistent color scaling
  max_expr <- max(GetAssayData(seurat_obj, assay = "RNA", layer = "data")[valid_genes, ], na.rm = TRUE)

  # Plot 1: Module score
  module_plot <- FeaturePlot(seurat_obj,
                              features = module_score_name,
                              reduction = "umap.harmony") +
    scale_colour_gradient(low = "lightgrey", high = "#e31837",
                          name = "Module\nScore") +
    ggtitle(paste0(marker_set_name, " - Module Score")) +
    theme(plot.title = element_text(size = 10, face = "bold"))

  # Plot 2: DimPlot with consensus cell types
  cluster_plot <- DimPlot(seurat_obj,
                          reduction = "umap.harmony",
                          group.by = "consensus_cell_type",
                          label = TRUE,
                          label.size = 2.5,
                          repel = TRUE) +
    ggtitle("Consensus Cell Types") +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 7)
    )

  # Plots 3-7: Top individual genes
  gene_plots <- list()
  for(i in seq_along(top_genes)) {
    gene_symbol <- gsub(".*-", "", top_genes[i])

    gene_plots[[i]] <- FeaturePlot(seurat_obj,
                                    features = top_genes[i],
                                    reduction = "umap.harmony",
                                    min.cutoff = 'q10') +
      scale_colour_gradient(low = "lightgrey", high = "#e31837",
                            limits = c(0, max_expr),
                            na.value = "#e31837",
                            name = "Expression") +
      ggtitle(gene_symbol) +
      theme(plot.title = element_text(size = 10, face = "bold"))
  }

  # Prepare title information
  gene_list_text <- paste(validation_results$gene_symbol, collapse = ", ")

  return(list(
    title_info = list(
      marker_set_name = marker_set_name,
      gene_list = gene_list_text,
      missing_genes = missing_genes
    ),
    module_plot = module_plot,
    cluster_plot = cluster_plot,
    gene_plots = gene_plots
  ))
}

## FUNCTION 5: Create DotPlot for a marker set------
create_marker_dotplot <- function(seurat_obj, marker_set_name, category_name,
                                  validation_results, timepoint_name) {

  cat("    Creating DotPlot for", marker_set_name, "\n")

  # Get valid genes
  valid_genes <- validation_results$full_name[validation_results$found]

  if(length(valid_genes) == 0) {
    cat("      WARNING: No valid genes - creating empty plot\n")
    empty_plot <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste0(marker_set_name, "\n\nNo valid genes found"),
               size = 6, color = "gray50") +
      theme_void()
    return(empty_plot)
  }

  # Create DotPlot
  dotplot <- DotPlot(seurat_obj,
                     features = valid_genes,
                     group.by = "consensus_cell_type",
                     dot.scale = 8) +
    coord_flip() +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 11, face = "bold"),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5),
      legend.position = "right"
    ) +
    scale_color_gradient(low = "lightgrey", high = "#e31837", name = "Avg\nExpression") +
    labs(
      title = paste0(marker_set_name, " (", category_name, ")"),
      subtitle = paste0(timepoint_name, " - Genes: ",
                       paste(gsub(".*-", "", valid_genes), collapse = ", ")),
      x = "Genes",
      y = "Consensus Cell Type"
    )

  return(dotplot)
}

## FUNCTION 6: Create summary page------
create_summary_page <- function(timepoint_name, seurat_obj, all_validation_results,
                               marker_structure) {

  cat("\n=== Creating summary page for", timepoint_name, "===\n")

  # Collect all missing genes across all sets
  missing_summary <- data.frame(
    category = character(),
    marker_set = character(),
    gene_symbol = character(),
    stringsAsFactors = FALSE
  )

  for(category_name in names(marker_structure)) {
    for(set_name in names(marker_structure[[category_name]])) {
      validation_key <- paste0(category_name, "::", set_name)
      validation <- all_validation_results[[validation_key]]
      missing_genes <- validation$gene_symbol[!validation$found]

      if(length(missing_genes) > 0) {
        missing_summary <- rbind(
          missing_summary,
          data.frame(
            category = rep(category_name, length(missing_genes)),
            marker_set = rep(set_name, length(missing_genes)),
            gene_symbol = missing_genes,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }

  # Create summary statistics
  total_categories <- length(marker_structure)
  total_marker_sets <- sum(sapply(marker_structure, length))
  total_genes_tested <- sum(sapply(all_validation_results, nrow))
  total_genes_found <- sum(sapply(all_validation_results, function(x) sum(x$found)))
  total_genes_missing <- total_genes_tested - total_genes_found
  pct_found <- round(100 * total_genes_found / total_genes_tested, 1)

  # Create summary text plot
  summary_text <- paste0(
    "MARKER GENE VALIDATION SUMMARY\n",
    "Timepoint: ", timepoint_name, "\n\n",
    "Dataset Information:\n",
    "  Total Cells: ", ncol(seurat_obj), "\n",
    "  Number of Clusters: ", length(unique(seurat_obj$seurat_clusters)), "\n",
    "  Consensus Cell Types: ", length(unique(seurat_obj$consensus_cell_type)), "\n\n",
    "Marker Structure:\n",
    "  Categories: ", total_categories, "\n",
    "  Marker Sets: ", total_marker_sets, "\n",
    "  Total Genes Tested: ", total_genes_tested, "\n",
    "  Genes Found: ", total_genes_found, " (", pct_found, "%)\n",
    "  Genes Missing: ", total_genes_missing, "\n"
  )

  summary_plot <- ggplot() +
    annotate("text", x = 0, y = 1,
             label = summary_text,
             hjust = 0, vjust = 1,
             size = 5, family = "mono") +
    theme_void() +
    xlim(0, 1) + ylim(0, 1)

  # Create missing genes table plot if any
  if(nrow(missing_summary) > 0) {
    missing_table <- gridExtra::tableGrob(missing_summary, rows = NULL)

    combined <- grid.arrange(
      summary_plot,
      missing_table,
      nrow = 2,
      heights = c(1, 1.5),
      top = textGrob(paste(timepoint_name, "- Marker Validation Summary"),
                     gp = gpar(fontsize = 16, fontface = "bold"))
    )
  } else {
    combined <- grid.arrange(
      summary_plot,
      nrow = 1,
      top = textGrob(paste(timepoint_name, "- Marker Validation Summary"),
                     gp = gpar(fontsize = 16, fontface = "bold"))
    )
  }

  return(combined)
}

## FUNCTION 7: Process one timepoint------
process_timepoint <- function(timepoint_name, file_path, marker_structure) {

  cat("\n\n####################################################\n")
  cat("PROCESSING TIMEPOINT:", timepoint_name, "\n")
  cat("####################################################\n")

  # Load consensus-annotated Seurat object
  cat("\nLoading consensus-annotated object from:", file_path, "\n")
  seurat_obj <- readRDS(file_path)
  cat("  Loaded:", ncol(seurat_obj), "cells,", nrow(seurat_obj), "features\n")
  cat("  Clusters:", length(unique(seurat_obj$seurat_clusters)), "\n")
  cat("  Consensus cell types:", length(unique(seurat_obj$consensus_cell_type)), "\n")

  # Join layers for Seurat v5 compatibility
  cat("  Joining layers for v5 compatibility...\n")
  seurat_obj <- JoinLayers(seurat_obj)
  cat("  Layers joined successfully\n")

  # Validate all marker sets
  cat("\n=== Validating marker gene sets ===\n")
  all_validation_results <- list()

  for(category_name in names(marker_structure)) {
    cat("\nCategory:", category_name, "\n")
    for(set_name in names(marker_structure[[category_name]])) {
      genes <- marker_structure[[category_name]][[set_name]]
      validation_key <- paste0(category_name, "::", set_name)
      all_validation_results[[validation_key]] <- validate_gene_set(
        genes,
        seurat_obj,
        set_name
      )
    }
  }

  # Add module scores for each marker set
  cat("\n=== Adding module scores ===\n")
  module_score_names <- character()
  module_idx <- 1

  for(category_name in names(marker_structure)) {
    for(set_name in names(marker_structure[[category_name]])) {
      validation_key <- paste0(category_name, "::", set_name)
      validation <- all_validation_results[[validation_key]]
      valid_genes <- validation$full_name[validation$found]

      if(length(valid_genes) > 0) {
        module_name <- paste0("MarkerSet_", module_idx)
        cat("  Adding module score for", set_name, "\n")
        seurat_obj <- AddModuleScore(
          seurat_obj,
          features = list(valid_genes),
          name = module_name
        )
        module_score_names[validation_key] <- paste0(module_name, "1")
        module_idx <- module_idx + 1
      } else {
        cat("  Skipping module score for", set_name, "(no valid genes)\n")
        module_score_names[validation_key] <- NA
      }
    }
  }

  # Create plots for each marker set
  cat("\n=== Creating validation plots ===\n")
  all_plots <- list()

  for(category_name in names(marker_structure)) {
    for(set_name in names(marker_structure[[category_name]])) {
      validation_key <- paste0(category_name, "::", set_name)
      all_plots[[validation_key]] <- create_marker_plot(
        seurat_obj,
        set_name,
        all_validation_results[[validation_key]],
        module_score_names[validation_key],
        timepoint_name
      )
    }
  }

  # Generate main validation PDF
  cat("\n=== Generating main validation PDF ===\n")
  pdf_file <- paste0("./output_files/marker_validation/",
                     timepoint_name, "_marker_validation.pdf")
  pdf(pdf_file, width = 14, height = 10)

  # Page 1: Summary
  create_summary_page(timepoint_name, seurat_obj, all_validation_results, marker_structure)

  # Pages 2+: Marker set plots organized by category
  for(category_name in names(marker_structure)) {
    cat("\n  Processing category:", category_name, "\n")
    for(set_name in names(marker_structure[[category_name]])) {
      validation_key <- paste0(category_name, "::", set_name)
      cat("    Plotting", set_name, "\n")
      plot_data <- all_plots[[validation_key]]

      # Create title text
      title_info <- plot_data$title_info
      if(length(title_info$missing_genes) > 0) {
        title_text <- paste0(
          title_info$marker_set_name, " (", category_name, ")\n",
          "Genes: ", title_info$gene_list, "\n",
          "[Missing: ", paste(title_info$missing_genes, collapse = ", "), "]"
        )
      } else {
        title_text <- paste0(
          title_info$marker_set_name, " (", category_name, ")\n",
          "Genes: ", title_info$gene_list
        )
      }

      # Create title grob
      title_grob <- textGrob(
        title_text,
        gp = gpar(fontsize = 10, fontface = "bold")
      )

      # Page 1 for this marker set: Module score + consensus cell types
      grid.arrange(
        plot_data$module_plot,
        plot_data$cluster_plot,
        ncol = 2,
        nrow = 1,
        top = title_grob
      )

      # Page 2 for this marker set: Individual genes
      if(length(plot_data$gene_plots) > 0) {
        grid.arrange(
          grobs = plot_data$gene_plots,
          ncol = 3,
          top = textGrob(
            paste0(title_info$marker_set_name, " - Individual Gene Expression"),
            gp = gpar(fontsize = 10, fontface = "bold")
          )
        )
      }
    }
  }

  dev.off()
  cat("\nMain validation PDF saved:", pdf_file, "\n")

  # Generate separate DotPlot PDF
  cat("\n=== Generating DotPlot PDF ===\n")
  dotplot_pdf_file <- paste0("./output_files/marker_validation/",
                              timepoint_name, "_marker_dotplots.pdf")
  pdf(dotplot_pdf_file, width = 12, height = 8)

  for(category_name in names(marker_structure)) {
    cat("\n  Creating DotPlots for category:", category_name, "\n")
    for(set_name in names(marker_structure[[category_name]])) {
      validation_key <- paste0(category_name, "::", set_name)

      dotplot <- create_marker_dotplot(
        seurat_obj,
        set_name,
        category_name,
        all_validation_results[[validation_key]],
        timepoint_name
      )

      print(dotplot)
    }
  }

  dev.off()
  cat("\nDotPlot PDF saved:", dotplot_pdf_file, "\n")

  cat("\n====================================\n")
  cat("VALIDATION COMPLETE:", timepoint_name, "\n")
  cat("====================================\n")

  return(list(
    validation_results = all_validation_results,
    plots = all_plots,
    pdf_file = pdf_file,
    dotplot_pdf_file = dotplot_pdf_file
  ))
}

## MAIN EXECUTION: Process all timepoints------
cat("\n\n####################################################\n")
cat("STARTING MARKER VALIDATION FOR ALL TIMEPOINTS\n")
cat("####################################################\n\n")

# Define timepoint file paths (consensus-annotated objects)
timepoint_files <- list(
  "1M_Org" = "./output_files/integrated_objects/1M_Org_integrated_harmony_consensus.rds",
  "3M_Org" = "./output_files/integrated_objects/3M_Org_integrated_harmony_consensus.rds",
  "6M_Org" = "./output_files/integrated_objects/6M_Org_integrated_harmony_consensus.rds"
)

# Process each timepoint
results <- list()

for(timepoint_name in names(timepoint_files)) {
  # Get timepoint-specific marker structure
  marker_structure <- markers_by_timepoint[[timepoint_name]]

  results[[timepoint_name]] <- process_timepoint(
    timepoint_name,
    timepoint_files[[timepoint_name]],
    marker_structure
  )

  # Force garbage collection
  gc()
}

cat("\n\n####################################################\n")
cat("ALL MARKER VALIDATIONS COMPLETE\n")
cat("####################################################\n")
cat("\nOutputs saved in: ./output_files/marker_validation/\n")
cat("\nGenerated validation reports:\n")
for(tp in names(results)) {
  cat("  ", tp, ":\n", sep = "")
  cat("    Main PDF: ", basename(results[[tp]]$pdf_file), "\n", sep = "")
  cat("    DotPlot PDF: ", basename(results[[tp]]$dotplot_pdf_file), "\n", sep = "")
}
cat("\nScript completed successfully!\n")
