library(tidyverse)
library(knitr)
library(kableExtra)

#######################################################
# Script to create sequencing depth table across modalities and ranks
#######################################################

setwd("/data/gpfs/projects/punim2251/Aim1_LongBench/longbench-analysis-pipeline")

# Define input paths
SC_SN_SAMPLE_SHEET <- "sc_processed_matrices/sc_sn_sample_sheet_used.csv"
BULK_SAMPLE_SHEET <- "bulk_processed_matrices/bulk_sample_sheet_used.csv"

# Output directory
OUTPUT_DIR <- "sequencing_depth_tables"

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR)
}

#######################################################
# Load and process sample sheets
#######################################################

message("Loading sample sheets...")

# Load SC/SN sample sheet
sc_sn_data <- read.csv(SC_SN_SAMPLE_SHEET, stringsAsFactors = FALSE)

# Load bulk sample sheet  
bulk_data <- read.csv(BULK_SAMPLE_SHEET, stringsAsFactors = FALSE)

# Filter for transcript data only (since that's what we're analyzing)
sc_sn_transcript <- sc_sn_data %>%
  filter(data_type == "transcript") %>%
  select(technology, sample_type, rank_order, sequencing_depth)

bulk_transcript <- bulk_data %>%
  filter(data_type == "transcript") %>%
  select(technology, sample_type, rank_order, sequencing_depth)

message(paste("Found", nrow(sc_sn_transcript), "SC/SN transcript entries"))
message(paste("Found", nrow(bulk_transcript), "bulk transcript entries"))

#######################################################
# Create comprehensive sequencing depth table
#######################################################

# Combine all data
all_data <- bind_rows(sc_sn_transcript, bulk_transcript)

# Create a more descriptive modality column
all_data <- all_data %>%
  mutate(
    modality = case_when(
      technology == "ONT" & sample_type == "sc" ~ "ONT_SC",
      technology == "ONT" & sample_type == "sn" ~ "ONT_SN", 
      technology == "ONT" & sample_type == "bulk" ~ "ONT_Bulk",
      technology == "PacBio" & sample_type == "sc" ~ "PacBio_SC",
      technology == "PacBio" & sample_type == "sn" ~ "PacBio_SN",
      technology == "PacBio" & sample_type == "bulk" ~ "PacBio_Bulk",
      TRUE ~ paste(technology, sample_type, sep = "_")
    )
  )

# Pivot to wide format for table
depth_table <- all_data %>%
  select(rank_order, modality, sequencing_depth) %>%
  pivot_wider(
    names_from = modality,
    values_from = sequencing_depth,
    names_sort = TRUE
  ) %>%
  arrange(rank_order)

# Reorder columns logically
column_order <- c("rank_order", "ONT_SC", "ONT_SN", "ONT_Bulk", "PacBio_SC", "PacBio_SN", "PacBio_Bulk")
available_columns <- intersect(column_order, names(depth_table))
depth_table <- depth_table %>% select(all_of(available_columns))

# Rename columns for publication
final_table <- depth_table %>%
  rename(
    "Rank" = "rank_order",
    "ONT SC" = "ONT_SC",
    "ONT SN" = "ONT_SN", 
    "ONT Bulk" = "ONT_Bulk",
    "PacBio SC" = "PacBio_SC",
    "PacBio SN" = "PacBio_SN",
    "PacBio Bulk" = "PacBio_Bulk"
  )

# Display the table
message("\nSequencing Depth Table (Full Numbers):")
message("SC/SN columns show median reads per cell")
message("Bulk columns show total reads")
message(paste(rep("=", 80), collapse = ""))
print(final_table)

# Save the table
write.csv(final_table, file.path(OUTPUT_DIR, "sequencing_depths_table.csv"), row.names = FALSE)

#######################################################
# Create publication-ready PDF table
#######################################################

message("Creating publication-ready PDF table...")

if (requireNamespace("kableExtra", quietly = TRUE)) {
  
  # Format numbers with commas but keep full precision
  formatted_for_pdf <- final_table %>%
    mutate(
      across(-Rank, ~ format(.x, big.mark = ",", scientific = FALSE))
    )
  
  # Create LaTeX table for PDF
  latex_table <- formatted_for_pdf %>%
    kable(
      format = "latex",
      caption = "Sequencing Depth by Rank Order Across Long-Read Sequencing Modalities",
      col.names = c("Rank", "ONT SC", "ONT SN", "ONT Bulk", "PacBio SC", "PacBio SN", "PacBio Bulk"),
      align = c('c', 'r', 'r', 'r', 'r', 'r', 'r'),
      booktabs = TRUE,
      linesep = ""
    ) %>%
    kable_styling(
      latex_options = c("striped", "hold_position", "scale_down"),
      full_width = FALSE,
      font_size = 10
    ) %>%
    add_header_above(c(" " = 1, "Single Cell/Nucleus\\n(median reads per cell)" = 4, "Bulk\\n(total reads)" = 2), 
                     escape = FALSE) %>%
    add_header_above(c(" " = 1, "Oxford Nanopore" = 2, " " = 1, "PacBio" = 2, " " = 1)) %>%
    footnote(
      general = c(
        "SC = Single Cell, SN = Single Nucleus",
        "Sequencing depths represent median reads per cell for SC/SN data and total reads for bulk data"
      ),
      threeparttable = TRUE,
      general_title = "Notes:"
    )
  
  # Save LaTeX table to file
  latex_file <- file.path(OUTPUT_DIR, "sequencing_depths_publication_table.tex")
  writeLines(latex_table, latex_file)
  
  #######################################################
  # Create publication-ready HTML table
  #######################################################
  
  message("Creating publication-ready HTML table...")
  
  # Create HTML table with proper formatting
  html_table <- formatted_for_pdf %>%
    kable(
      format = "html",
      caption = "Table 1. Sequencing Depth by Rank Order Across Long-Read Sequencing Modalities",
      col.names = c("Rank", "ONT SC", "ONT SN", "ONT Bulk", "PacBio SC", "PacBio SN", "PacBio Bulk"),
      align = c('c', 'r', 'r', 'r', 'r', 'r', 'r'),
      table.attr = 'style="font-family: Arial, sans-serif; border-collapse: collapse; margin: 20px auto; width: 90%; min-width: 800px;"'
    ) %>%
    kable_styling(
      bootstrap_options = c("striped", "hover", "condensed"),
      full_width = TRUE,
      position = "center",
      font_size = 14
    ) %>%
    add_header_above(c(" " = 1, "Single Cell/Nucleus (median reads per cell)" = 4, "Bulk (total reads)" = 2)) %>%
    add_header_above(c(" " = 1, "Oxford Nanopore" = 2, " " = 1, "PacBio" = 2, " " = 1)) %>%
    footnote(
      general = c(
        "SC = Single Cell, SN = Single Nucleus",
        "Sequencing depths represent median reads per cell for SC/SN data and total reads for bulk data"
      ),
      general_title = "Notes:",
      footnote_as_chunk = TRUE
    )
  
  # Create a complete HTML document
  html_document <- paste0(
    '<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Sequencing Depth Summary Table</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 40px auto;
            max-width: 1200px;
            line-height: 1.6;
        }
        h1 {
            text-align: center;
            color: #333;
            margin-bottom: 30px;
        }
        .table-container {
            margin: 20px 0;
            text-align: center;
            overflow-x: auto;
        }
        table {
            margin: 0 auto;
            border-collapse: collapse;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            width: 90%;
            min-width: 800px;
        }
        th, td {
            white-space: nowrap;
            padding: 10px 15px !important;
            border: 1px solid #ddd;
        }
        th {
            background-color: #f8f9fa;
            font-weight: bold;
        }
        td {
            text-align: right;
        }
        .rank-col {
            text-align: center !important;
        }
        @media print {
            body { 
                margin: 10px; 
                max-width: none;
            }
            table { 
                font-size: 9pt;
                width: 100%;
            }
            th, td {
                padding: 6px 8px !important;
            }
        }
        @media screen and (max-width: 768px) {
            body {
                margin: 20px;
            }
            .table-container {
                overflow-x: scroll;
            }
        }
    </style>
</head>
<body>
    <h1>Sequencing Depth Summary Table</h1>
    <div class="table-container">
        ', html_table, '
    </div>
</body>
</html>'
  )
  
  # Save HTML document
  html_file <- file.path(OUTPUT_DIR, "sequencing_depths_publication_table.html")
  writeLines(html_document, html_file)
  
  message(paste("HTML table created:", html_file))
  message("To create PDF: Open the HTML file in your browser and print to PDF")
  
  # Create a standalone LaTeX document
  latex_document <- paste0(
    "\\documentclass[11pt]{article}\n",
    "\\usepackage[margin=1in]{geometry}\n",
    "\\usepackage{booktabs}\n",
    "\\usepackage{threeparttable}\n",
    "\\usepackage{array}\n",
    "\\begin{document}\n",
    "\\title{Sequencing Depth Summary Table}\n",
    "\\date{}\n",
    "\\maketitle\n\n",
    latex_table, "\n\n",
    "\\end{document}"
  )
  
  # Save the complete LaTeX document
  latex_doc_file <- file.path(OUTPUT_DIR, "sequencing_depths_document.tex")
  writeLines(latex_document, latex_doc_file)
  
  message("LaTeX files created successfully!")
  message(paste("LaTeX table:", latex_file))
  message(paste("LaTeX document:", latex_doc_file))
  
  # Try to compile to PDF if pdflatex is available
  tryCatch({
    system_result <- system(paste("pdflatex -output-directory", OUTPUT_DIR, latex_doc_file), 
                           ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    if (system_result == 0) {
      # Clean up auxiliary files
      aux_files <- list.files(OUTPUT_DIR, pattern = "\\.(aux|log)$", full.names = TRUE)
      file.remove(aux_files)
      
      pdf_file <- file.path(OUTPUT_DIR, "sequencing_depths_document.pdf")
      message(paste("PDF created successfully:", pdf_file))
    } else {
      message("pdflatex not available or compilation failed.")
      message("You can compile manually with: pdflatex sequencing_depths_document.tex")
    }
  }, error = function(e) {
    message("Could not compile PDF automatically.")
    message("LaTeX files created - compile manually with: pdflatex sequencing_depths_document.tex")
  })
  
} else {
  message("kableExtra package not available. Install with: install.packages('kableExtra')")
}

message(paste(rep("=", 80), collapse = ""))
message(paste("Files saved to:", OUTPUT_DIR))
message("• sequencing_depths_table.csv - Raw data table")
message("• sequencing_depths_publication_table.html - Publication-ready HTML table")
message("• sequencing_depths_publication_table.tex - LaTeX table only")
message("• sequencing_depths_document.tex - Complete LaTeX document")
message("• sequencing_depths_document.pdf - Publication-ready PDF (if compilation succeeded)")
message("")
message("=== TO CREATE PDF FROM HTML ===")
message("1. Open sequencing_depths_publication_table.html in your web browser")
message("2. Press Ctrl+P (or Cmd+P on Mac) to print")
message("3. Select 'Save as PDF' as the destination")
message("4. Adjust margins if needed and save") 