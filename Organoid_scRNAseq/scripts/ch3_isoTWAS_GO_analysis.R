#!/usr/bin/env Rscript
# ============================================================================
# GO Enrichment Analysis for SCZ isoTWAS-DE Hits: oRG to Astrocytes (6M)
# Author: Manveer Chauhan
# Date: 2025-10-28
# Description: Performs Gene Ontology enrichment analysis on SCZ isoTWAS-DE
#              gene-level hits for oRG to Astrocytes (Transitioning) at
#              6 months, using genes identified from isoform-level TWAS
#              differential expression analysis.
# ============================================================================

# REQUIRED LIBRARIES ---------------------------------------------------------
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(tidyverse)
  library(ggplot2)
})

# GLOBAL PARAMETERS ----------------------------------------------------------
BASE_DIR <- "/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space"
setwd(BASE_DIR)

# Input CSV file
INPUT_CSV <- "./output_files/isoTWAS_canonical_analysis/SCZ_6M_with_canonical.csv"

# Target cell type
TARGET_CELL_TYPE <- "oRG to Astrocytes (Transitioning)"

# Output directory
OUTPUT_DIR <- "./output_files/isoTWAS_GO_analysis"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Output prefix for all files
OUTPUT_PREFIX <- "SCZ_6M_oRG_to_Astrocytes"

# GO enrichment parameters
PVALUE_CUTOFF <- 0.05
QVALUE_CUTOFF <- 0.05
P_ADJUST_METHOD <- "BH"

# HELPER FUNCTIONS -----------------------------------------------------------

#' Convert gene symbols to ENTREZ IDs with error handling
#'
#' @param gene_symbols Character vector of gene symbols
#' @return data.frame with columns: gene_symbol, entrez_id, conversion_success
convert_symbols_to_entrez <- function(gene_symbols) {
  cat("Converting gene symbols to ENTREZ IDs...\n")
  cat("  Input genes:", length(gene_symbols), "\n")
  cat("  Unique genes:", length(unique(gene_symbols)), "\n\n")

  # Attempt conversion using bitr
  conversion_result <- tryCatch({
    bitr(
      geneID = unique(gene_symbols),
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
  }, error = function(e) {
    warning("bitr conversion failed: ", e$message)
    return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
  })

  # Rename columns for clarity
  colnames(conversion_result) <- c("gene_symbol", "entrez_id")

  # Keep only first ENTREZ ID per gene symbol to avoid duplicates
  conversion_result <- conversion_result %>%
    dplyr::group_by(gene_symbol) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  # Create full mapping dataframe
  input_df <- data.frame(
    gene_symbol = unique(gene_symbols),
    stringsAsFactors = FALSE
  )

  output_df <- left_join(input_df, conversion_result, by = "gene_symbol")
  output_df$conversion_success <- !is.na(output_df$entrez_id)

  # Calculate success rate
  success_rate <- sum(output_df$conversion_success) / nrow(output_df) * 100
  cat("  ✓ Conversion success rate:", round(success_rate, 2), "%\n")
  cat("  ✓ Successfully converted:", sum(output_df$conversion_success), "genes\n")
  cat("  ✗ Failed conversions:", sum(!output_df$conversion_success), "genes\n")

  if (sum(!output_df$conversion_success) > 0) {
    failed_genes <- output_df$gene_symbol[!output_df$conversion_success]
    cat("\n  Failed gene symbols:\n")
    cat("   ", paste(head(failed_genes, 10), collapse = ", "))
    if (length(failed_genes) > 10) {
      cat(" ... (", length(failed_genes) - 10, " more)", sep = "")
    }
    cat("\n")
  }
  cat("\n")

  return(output_df)
}

#' Perform GO enrichment analysis for a single ontology
#'
#' @param entrez_ids Character vector of ENTREZ gene IDs
#' @param ontology One of "BP", "MF", "CC"
#' @param pvalue_cutoff Numeric p-value cutoff (default: 0.05)
#' @param qvalue_cutoff Numeric q-value cutoff (default: 0.05)
#' @return enrichResult object or NULL if enrichment fails
run_go_enrichment <- function(entrez_ids, ontology,
                              pvalue_cutoff = PVALUE_CUTOFF,
                              qvalue_cutoff = QVALUE_CUTOFF) {

  ontology_names <- c(
    "BP" = "Biological Process",
    "MF" = "Molecular Function",
    "CC" = "Cellular Component"
  )

  cat("Running GO enrichment:", ontology_names[ontology], "\n")
  cat("  Input genes:", length(entrez_ids), "\n")
  cat("  p-value cutoff:", pvalue_cutoff, "\n")
  cat("  q-value cutoff:", qvalue_cutoff, "\n")

  go_result <- tryCatch({
    enrichGO(
      gene = entrez_ids,
      OrgDb = org.Hs.eg.db,
      ont = ontology,
      pAdjustMethod = P_ADJUST_METHOD,
      pvalueCutoff = pvalue_cutoff,
      qvalueCutoff = qvalue_cutoff,
      readable = TRUE  # Convert ENTREZ IDs back to gene symbols in results
    )
  }, error = function(e) {
    warning(ontology, " enrichment failed: ", e$message)
    return(NULL)
  })

  if (!is.null(go_result) && nrow(go_result@result) > 0) {
    cat("  ✓ Enriched terms found:", nrow(go_result@result), "\n")
  } else {
    cat("  ✗ No enriched terms found\n")
  }
  cat("\n")

  return(go_result)
}

#' Create dotplot for GO enrichment results (manual ggplot version)
#'
#' @param go_result enrichResult object
#' @param title Plot title
#' @param showCategory Number of top categories to show (default: 15)
#' @return ggplot object
create_go_dotplot <- function(go_result, title, showCategory = 15) {
  if (is.null(go_result) || nrow(go_result@result) == 0) {
    return(NULL)
  }

  tryCatch({
    # Extract top results
    top_results <- go_result@result %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = showCategory) %>%
      dplyr::mutate(
        Description = factor(Description, levels = rev(Description)),
        GeneRatio_numeric = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
      )

    # Create manual dotplot
    p <- ggplot(top_results, aes(x = GeneRatio_numeric, y = Description)) +
      geom_point(aes(size = Count, color = p.adjust)) +
      scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +
      scale_size_continuous(name = "Gene Count", range = c(3, 8)) +
      labs(
        title = title,
        x = "Gene Ratio",
        y = NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10, face = "bold"),
        panel.grid.major.y = element_line(color = "grey90"),
        panel.grid.minor = element_blank()
      )

    return(p)
  }, error = function(e) {
    warning("Manual dotplot creation failed: ", e$message)
    return(NULL)
  })
}

#' Create barplot for GO enrichment results (manual ggplot version)
#'
#' @param go_result enrichResult object
#' @param title Plot title
#' @param showCategory Number of top categories to show (default: 20)
#' @return ggplot object
create_go_barplot <- function(go_result, title, showCategory = 20) {
  if (is.null(go_result) || nrow(go_result@result) == 0) {
    return(NULL)
  }

  tryCatch({
    # Extract top results
    top_results <- go_result@result %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = showCategory) %>%
      dplyr::mutate(
        Description = factor(Description, levels = rev(Description))
      )

    # Create manual barplot
    p <- ggplot(top_results, aes(x = Count, y = Description, fill = p.adjust)) +
      geom_bar(stat = "identity") +
      scale_fill_gradient(low = "red", high = "blue", name = "p.adjust") +
      labs(
        title = title,
        x = "Gene Count",
        y = NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10, face = "bold"),
        panel.grid.major.y = element_line(color = "grey90"),
        panel.grid.minor = element_blank()
      )

    return(p)
  }, error = function(e) {
    warning("Manual barplot creation failed: ", e$message)
    return(NULL)
  })
}

#' Save GO enrichment results
#'
#' @param go_result enrichResult object
#' @param ontology One of "BP", "MF", "CC"
#' @param prefix Output file prefix
#' @param output_dir Output directory path
save_go_results <- function(go_result, ontology, prefix, output_dir) {
  if (is.null(go_result) || nrow(go_result@result) == 0) {
    cat("  Skipping", ontology, "- no results to save\n\n")
    return(invisible(NULL))
  }

  ontology_names <- c(
    "BP" = "Biological Process",
    "MF" = "Molecular Function",
    "CC" = "Cellular Component"
  )

  cat("Saving", ontology, "results...\n")

  # 1. Save dotplot
  dotplot_path <- file.path(output_dir, paste0(prefix, "_", ontology, "_dotplot.pdf"))
  p_dot <- create_go_dotplot(go_result,
                              title = paste0(prefix, " - ", ontology_names[ontology]))

  if (!is.null(p_dot)) {
    tryCatch({
      ggsave(
        filename = dotplot_path,
        plot = p_dot,
        width = 12,
        height = 10,
        device = "pdf"
      )
      cat("  ✓ Saved dotplot:", basename(dotplot_path), "\n")
    }, error = function(e) {
      warning("Failed to save dotplot: ", e$message)
      cat("  ✗ Dotplot save failed:", e$message, "\n")
    })
  } else {
    cat("  ✗ Dotplot creation failed (skipping)\n")
  }

  # 2. Save barplot
  barplot_path <- file.path(output_dir, paste0(prefix, "_", ontology, "_barplot.pdf"))
  p_bar <- create_go_barplot(go_result,
                              title = paste0(prefix, " - ", ontology_names[ontology]))

  if (!is.null(p_bar)) {
    tryCatch({
      ggsave(
        filename = barplot_path,
        plot = p_bar,
        width = 12,
        height = 10,
        device = "pdf"
      )
      cat("  ✓ Saved barplot:", basename(barplot_path), "\n")
    }, error = function(e) {
      warning("Failed to save barplot: ", e$message)
      cat("  ✗ Barplot save failed:", e$message, "\n")
    })
  } else {
    cat("  ✗ Barplot creation failed (skipping)\n")
  }

  # 3. Save enrichment table as CSV
  csv_path <- file.path(output_dir, paste0(prefix, "_", ontology, "_results.csv"))
  write.csv(go_result@result, csv_path, row.names = FALSE, quote = FALSE)
  cat("  ✓ Saved results table:", basename(csv_path), "\n")

  # 4. Save enrichResult object as RDS
  rds_path <- file.path(output_dir, paste0(prefix, "_", ontology, "_results.rds"))
  saveRDS(go_result, rds_path)
  cat("  ✓ Saved RDS object:", basename(rds_path), "\n")

  cat("\n")

  return(invisible(NULL))
}

# MAIN ANALYSIS --------------------------------------------------------------

cat("\n")
cat("================================================================================\n")
cat("  GO Enrichment Analysis: SCZ isoTWAS-DE Hits - oRG to Astrocytes (6M)\n")
cat("  Author: Manveer Chauhan\n")
cat("  Date:", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("================================================================================\n\n")

# SECTION 1: Load and Filter Data --------------------------------------------
cat("========================================\n")
cat("SECTION 1: Loading Data\n")
cat("========================================\n\n")

cat("Loading CSV file:", INPUT_CSV, "\n")

if (!file.exists(INPUT_CSV)) {
  stop("Error: Input CSV file not found at ", INPUT_CSV)
}

data <- read.csv(INPUT_CSV, stringsAsFactors = FALSE)

cat("  ✓ Loaded", nrow(data), "transcript records\n")
cat("  Cell types present:", length(unique(data$cell_type_of_origin)), "\n")
cat("\n")

# Filter for target cell type
cat("Filtering for target cell type:", TARGET_CELL_TYPE, "\n")
data_filtered <- data %>%
  dplyr::filter(cell_type_of_origin == TARGET_CELL_TYPE)

if (nrow(data_filtered) == 0) {
  stop("Error: No data found for cell type '", TARGET_CELL_TYPE, "'")
}

cat("  ✓ Filtered to", nrow(data_filtered), "transcript records\n")
cat("  Unique genes:", length(unique(data_filtered$gene_symbol)), "\n")
cat("  Canonical status breakdown:\n")
print(table(data_filtered$canonical_status))
cat("\n")

# Extract unique gene symbols (collapse multiple transcripts per gene)
gene_list <- unique(data_filtered$gene_symbol)
cat("Gene list for GO analysis:\n")
cat("  Total unique genes:", length(gene_list), "\n\n")

# SECTION 2: Gene ID Conversion ----------------------------------------------
cat("========================================\n")
cat("SECTION 2: Gene ID Conversion\n")
cat("========================================\n\n")

conversion_df <- convert_symbols_to_entrez(gene_list)

# Save gene list with conversion mapping
gene_list_path <- file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_gene_list.csv"))
write.csv(conversion_df, gene_list_path, row.names = FALSE, quote = FALSE)
cat("Saved gene list with conversion mapping to:", gene_list_path, "\n\n")

# Extract successfully converted ENTREZ IDs
entrez_ids <- conversion_df %>%
  dplyr::filter(conversion_success) %>%
  pull(entrez_id)

if (length(entrez_ids) == 0) {
  stop("Error: No genes could be converted to ENTREZ IDs. Cannot proceed with GO analysis.")
}

cat("Proceeding with", length(entrez_ids), "genes for GO enrichment\n\n")

# SECTION 3: GO Enrichment Analysis ------------------------------------------
cat("========================================\n")
cat("SECTION 3: GO Enrichment Analysis\n")
cat("========================================\n\n")

go_results <- list()

# Biological Process
go_results[["BP"]] <- run_go_enrichment(entrez_ids, "BP")

# Molecular Function
go_results[["MF"]] <- run_go_enrichment(entrez_ids, "MF")

# Cellular Component
go_results[["CC"]] <- run_go_enrichment(entrez_ids, "CC")

# SECTION 4: Save Results ----------------------------------------------------
cat("========================================\n")
cat("SECTION 4: Saving Results\n")
cat("========================================\n\n")

for (ont in c("BP", "MF", "CC")) {
  save_go_results(go_results[[ont]], ont, OUTPUT_PREFIX, OUTPUT_DIR)
}

# SECTION 5: Generate Summary Report -----------------------------------------
cat("========================================\n")
cat("SECTION 5: Summary Report\n")
cat("========================================\n\n")

# Collect summary statistics
summary_text <- c(
  "================================================================================",
  "  GO Enrichment Analysis Summary",
  "  SCZ isoTWAS-DE Hits: oRG to Astrocytes (Transitioning) at 6M",
  "================================================================================",
  "",
  "Analysis Date:", format(Sys.Date(), "%Y-%m-%d"),
  "",
  "INPUT DATA:",
  paste0("  CSV File: ", basename(INPUT_CSV)),
  paste0("  Cell Type: ", TARGET_CELL_TYPE),
  paste0("  Total Transcripts: ", nrow(data_filtered)),
  paste0("  Unique Genes: ", length(gene_list)),
  "",
  "GENE ID CONVERSION:",
  paste0("  Input Genes: ", nrow(conversion_df)),
  paste0("  Successfully Converted: ", sum(conversion_df$conversion_success),
         " (", round(sum(conversion_df$conversion_success) / nrow(conversion_df) * 100, 1), "%)"),
  paste0("  Failed Conversions: ", sum(!conversion_df$conversion_success)),
  "",
  "GO ENRICHMENT PARAMETERS:",
  paste0("  p-value cutoff: ", PVALUE_CUTOFF),
  paste0("  q-value cutoff: ", QVALUE_CUTOFF),
  paste0("  p-adjust method: ", P_ADJUST_METHOD),
  "",
  "ENRICHMENT RESULTS:",
  ""
)

ontology_names <- c(
  "BP" = "Biological Process",
  "MF" = "Molecular Function",
  "CC" = "Cellular Component"
)

for (ont in c("BP", "MF", "CC")) {
  go_result <- go_results[[ont]]

  summary_text <- c(
    summary_text,
    paste0(ontology_names[ont], " (", ont, "):"),
    ""
  )

  if (!is.null(go_result) && nrow(go_result@result) > 0) {
    n_terms <- nrow(go_result@result)
    summary_text <- c(
      summary_text,
      paste0("  ✓ Enriched Terms: ", n_terms),
      ""
    )

    # Add top 5 terms
    if (n_terms > 0) {
      top_terms <- head(go_result@result[order(go_result@result$p.adjust), ], 5)
      summary_text <- c(
        summary_text,
        "  Top 5 Enriched Terms:",
        ""
      )

      for (i in 1:min(5, nrow(top_terms))) {
        term_line <- sprintf("    %d. %s", i, top_terms$Description[i])
        pval_line <- sprintf("       p.adjust = %.2e, Count = %s",
                             top_terms$p.adjust[i],
                             top_terms$Count[i])
        summary_text <- c(summary_text, term_line, pval_line)
      }
      summary_text <- c(summary_text, "")
    }
  } else {
    summary_text <- c(
      summary_text,
      "  ✗ No enriched terms found",
      ""
    )
  }
}

summary_text <- c(
  summary_text,
  "OUTPUT FILES:",
  paste0("  Output Directory: ", OUTPUT_DIR),
  "",
  "  Generated files:",
  paste0("    - ", OUTPUT_PREFIX, "_gene_list.csv"),
  paste0("    - ", OUTPUT_PREFIX, "_[BP/MF/CC]_dotplot.pdf"),
  paste0("    - ", OUTPUT_PREFIX, "_[BP/MF/CC]_barplot.pdf"),
  paste0("    - ", OUTPUT_PREFIX, "_[BP/MF/CC]_results.csv"),
  paste0("    - ", OUTPUT_PREFIX, "_[BP/MF/CC]_results.rds"),
  paste0("    - ", OUTPUT_PREFIX, "_GO_summary.txt"),
  "",
  "================================================================================",
  ""
)

# Print summary to console
cat(paste(summary_text, collapse = "\n"))

# Save summary to file
summary_path <- file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_GO_summary.txt"))
writeLines(summary_text, summary_path)
cat("\n✓ Summary saved to:", summary_path, "\n")

# COMPLETION MESSAGE ---------------------------------------------------------
cat("\n")
cat("================================================================================\n")
cat("  Analysis Complete!\n")
cat("  Output directory:", OUTPUT_DIR, "\n")
cat("================================================================================\n\n")
