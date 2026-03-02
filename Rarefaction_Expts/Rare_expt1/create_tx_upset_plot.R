library(ggplot2)
library(ComplexUpset)
library(tidyverse)

#######################################################
# Define directories and setup
#######################################################

setwd("/data/gpfs/projects/punim2251/Aim1_LongBench/longbench-analysis-pipeline")

# Input directories with processed matrices
SC_INPUT_DIR <- "sc_processed_matrices_with_metadata"
BULK_INPUT_DIR <- "bulk_processed_matrices_with_metadata"

# Output directory for upset plots
OUTPUT_DIR <- "upset_plots"

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR)
}

# Control how many intersection sets to save per rank (set to 25 by default)
MAX_INTERSECTIONS_TO_SAVE <- 25

# Minimum size threshold for intersection sets (should match upset plot min_size)
MIN_INTERSECTION_SIZE <- 100

#######################################################
# Load transcript data from different modalities
#######################################################

message("Loading transcript data from all modalities...")

# Load single-cell transcript data
sc_tx_ont_sc <- readRDS(file.path(SC_INPUT_DIR, "filtered_isoPseudobulkDfs_ont_sc_tx_with_metadata.rds"))
sc_tx_pb_sc <- readRDS(file.path(SC_INPUT_DIR, "filtered_isoPseudobulkDfs_pacbio_sc_tx_with_metadata.rds"))

# Load single-nucleus transcript data
sc_tx_ont_sn <- readRDS(file.path(SC_INPUT_DIR, "filtered_isoPseudobulkDfs_ont_sn_tx_with_metadata.rds"))
sc_tx_pb_sn <- readRDS(file.path(SC_INPUT_DIR, "filtered_isoPseudobulkDfs_pacbio_sn_tx_with_metadata.rds"))

# Load bulk transcript data
bulk_tx_ont <- readRDS(file.path(BULK_INPUT_DIR, "bulk_isoDfs_ont_filtered_with_metadata.rds"))
bulk_tx_pb <- readRDS(file.path(BULK_INPUT_DIR, "bulk_isoDfs_pb_filtered_with_metadata.rds"))

#######################################################
# Determine available rank orders
#######################################################

# Extract rank orders from ALL datasets to get complete set
all_datasets <- list(sc_tx_ont_sc, sc_tx_pb_sc, sc_tx_ont_sn, sc_tx_pb_sn, bulk_tx_ont, bulk_tx_pb)
all_ranks <- c()

for(dataset in all_datasets) {
  if(!is.null(dataset)) {
    ranks_in_dataset <- sapply(names(dataset), function(x) {
      rank_match <- stringr::str_extract(x, "rank(\\d+)")
      if (!is.na(rank_match)) {
        as.numeric(stringr::str_extract(rank_match, "\\d+"))
      } else {
        NA
      }
    })
    all_ranks <- c(all_ranks, ranks_in_dataset)
  }
}

# Get unique ranks and remove NA values
available_ranks <- sort(unique(all_ranks[!is.na(all_ranks)]))

message(paste("Available rank orders:", paste(available_ranks, collapse = ", ")))

# Debug: Show rank counts per dataset
message("Rank counts per dataset:")
dataset_names <- c("sc_tx_ont_sc", "sc_tx_pb_sc", "sc_tx_ont_sn", "sc_tx_pb_sn", "bulk_tx_ont", "bulk_tx_pb")
for(i in 1:length(all_datasets)) {
  if(!is.null(all_datasets[[i]])) {
    ranks_in_this <- sapply(names(all_datasets[[i]]), function(x) {
      rank_match <- stringr::str_extract(x, "rank(\\d+)")
      if (!is.na(rank_match)) {
        as.numeric(stringr::str_extract(rank_match, "\\d+"))
      } else {
        NA
      }
    })
    unique_ranks <- sort(unique(ranks_in_this[!is.na(ranks_in_this)]))
    message(paste(dataset_names[i], "- ranks:", paste(unique_ranks, collapse = ", ")))
  }
}

######################################################
# Process each rank order
######################################################

for(rank_order in available_ranks){
  message(paste0("Processing rank order: ", rank_order, "..."))
  
  # Find the corresponding matrix from each dataset using rank order
  ont_sc_matrix <- NULL
  pb_sc_matrix <- NULL
  ont_sn_matrix <- NULL
  pb_sn_matrix <- NULL
  ont_bulk_matrix <- NULL
  pb_bulk_matrix <- NULL
  
  # Find matrices by rank order
  for(name in names(sc_tx_ont_sc)) {
    if(grepl(paste0("rank", rank_order), name)) {
      ont_sc_matrix <- sc_tx_ont_sc[[name]]
      break
    }
  }
  
  for(name in names(sc_tx_pb_sc)) {
    if(grepl(paste0("rank", rank_order), name)) {
      pb_sc_matrix <- sc_tx_pb_sc[[name]]
      break
    }
  }
  
  for(name in names(sc_tx_ont_sn)) {
    if(grepl(paste0("rank", rank_order), name)) {
      ont_sn_matrix <- sc_tx_ont_sn[[name]]
      break
    }
  }
  
  for(name in names(sc_tx_pb_sn)) {
    if(grepl(paste0("rank", rank_order), name)) {
      pb_sn_matrix <- sc_tx_pb_sn[[name]]
      break
    }
  }
  
  for(name in names(bulk_tx_ont)) {
    if(grepl(paste0("rank", rank_order), name)) {
      ont_bulk_matrix <- bulk_tx_ont[[name]]
      break
    }
  }
  
  for(name in names(bulk_tx_pb)) {
    if(grepl(paste0("rank", rank_order), name)) {
      pb_bulk_matrix <- bulk_tx_pb[[name]]
      break
    }
  }
  
  # Check that we found matrices for this rank
  matrices_found <- sum(!sapply(list(ont_sc_matrix, pb_sc_matrix, ont_sn_matrix, 
                                    pb_sn_matrix, ont_bulk_matrix, pb_bulk_matrix), is.null))
  
  if(matrices_found == 0) {
    message(paste("No matrices found for rank", rank_order, "- skipping"))
    next
  }
  
  message(paste("Found", matrices_found, "matrices for rank", rank_order))
  
  # Extract transcript IDs from each matrix (handle potential NULLs)
  ont_sc_tx_ids <- if(!is.null(ont_sc_matrix)) ont_sc_matrix$tx_id else character(0)
  pb_sc_tx_ids <- if(!is.null(pb_sc_matrix)) pb_sc_matrix$tx_id else character(0)
  ont_sn_tx_ids <- if(!is.null(ont_sn_matrix)) ont_sn_matrix$tx_id else character(0)
  pb_sn_tx_ids <- if(!is.null(pb_sn_matrix)) pb_sn_matrix$tx_id else character(0)
  ont_bulk_tx_ids <- if(!is.null(ont_bulk_matrix)) ont_bulk_matrix$tx_id else character(0)
  pb_bulk_tx_ids <- if(!is.null(pb_bulk_matrix)) pb_bulk_matrix$tx_id else character(0)
  
  # Create a data frame with all unique transcript IDs
  all_tx_ids <- unique(c(ont_sc_tx_ids, pb_sc_tx_ids, ont_sn_tx_ids, 
                        pb_sn_tx_ids, ont_bulk_tx_ids, pb_bulk_tx_ids))
  
  if(length(all_tx_ids) == 0) {
    message(paste("No transcript IDs found for rank", rank_order, "- skipping"))
    next
  }
  
  # Create a comprehensive presence/absence matrix that preserves ALL original metadata
  # Initialize the core data frame
  all_data <- data.frame(
    tx_id = all_tx_ids,
    ONT_SC = all_tx_ids %in% ont_sc_tx_ids,
    PB_SC = all_tx_ids %in% pb_sc_tx_ids,
    ONT_SN = all_tx_ids %in% ont_sn_tx_ids,
    PB_SN = all_tx_ids %in% pb_sn_tx_ids,    
    ONT_Bulk = all_tx_ids %in% ont_bulk_tx_ids,
    PB_Bulk = all_tx_ids %in% pb_bulk_tx_ids,
    stringsAsFactors = FALSE
  )
  
  # Get all possible metadata columns from the first available matrix
  sample_matrix <- NULL
  if(!is.null(ont_sc_matrix)) sample_matrix <- ont_sc_matrix
  else if(!is.null(pb_sc_matrix)) sample_matrix <- pb_sc_matrix
  else if(!is.null(ont_sn_matrix)) sample_matrix <- ont_sn_matrix
  else if(!is.null(pb_sn_matrix)) sample_matrix <- pb_sn_matrix
  else if(!is.null(ont_bulk_matrix)) sample_matrix <- ont_bulk_matrix
  else if(!is.null(pb_bulk_matrix)) sample_matrix <- pb_bulk_matrix
  
  # Add all metadata columns from the sample matrix
  if(!is.null(sample_matrix)) {
    metadata_cols <- setdiff(colnames(sample_matrix), c("tx_id", "counts", "sampleID", "targetMedian", "rarefactionDesc"))
    for(col in metadata_cols) {
      all_data[[col]] <- NA
    }
  }
  
  # Add summary columns for counts and expression
  all_data$counts <- NA
  all_data$log2_cpm_plus1 <- NA
  all_data$sequencing_depth <- NA
  
  # Modified approach to calculate average counts across modalities
  for(i in 1:nrow(all_data)) {
    tx <- all_data$tx_id[i]
    counts_values <- c()
    log2_cpm_values <- c()
    
    # Collect counts and metadata from all available modalities
    if(all_data$ONT_SC[i] && !is.null(ont_sc_matrix)) {
      idx <- which(ont_sc_matrix$tx_id == tx)
      if(length(idx) > 0) {
        # Fill all metadata columns if not already filled
        for(col in metadata_cols) {
          if(col %in% colnames(ont_sc_matrix) && is.na(all_data[[col]][i])) {
            all_data[[col]][i] <- ont_sc_matrix[[col]][idx[1]]
          }
        }
        # Set sequencing depth if not set
        if(is.na(all_data$sequencing_depth[i]) && "targetMedian" %in% colnames(ont_sc_matrix)) {
          all_data$sequencing_depth[i] <- ont_sc_matrix$targetMedian[idx[1]]
        }
        counts_values <- c(counts_values, ont_sc_matrix$counts[idx[1]])
        log2_cpm_values <- c(log2_cpm_values, ont_sc_matrix$log2_cpm_plus1[idx[1]])
      }
    }
    
    if(all_data$PB_SC[i] && !is.null(pb_sc_matrix)) {
      idx <- which(pb_sc_matrix$tx_id == tx)
      if(length(idx) > 0) {
        # Fill all metadata columns if not already filled
        for(col in metadata_cols) {
          if(col %in% colnames(pb_sc_matrix) && is.na(all_data[[col]][i])) {
            all_data[[col]][i] <- pb_sc_matrix[[col]][idx[1]]
          }
        }
        # Set sequencing depth if not set
        if(is.na(all_data$sequencing_depth[i]) && "targetMedian" %in% colnames(pb_sc_matrix)) {
          all_data$sequencing_depth[i] <- pb_sc_matrix$targetMedian[idx[1]]
        }
        counts_values <- c(counts_values, pb_sc_matrix$counts[idx[1]])
        log2_cpm_values <- c(log2_cpm_values, pb_sc_matrix$log2_cpm_plus1[idx[1]])
      }
    }
    
    if(all_data$ONT_SN[i] && !is.null(ont_sn_matrix)) {
      idx <- which(ont_sn_matrix$tx_id == tx)
      if(length(idx) > 0) {
        # Fill all metadata columns if not already filled
        for(col in metadata_cols) {
          if(col %in% colnames(ont_sn_matrix) && is.na(all_data[[col]][i])) {
            all_data[[col]][i] <- ont_sn_matrix[[col]][idx[1]]
          }
        }
        # Set sequencing depth if not set
        if(is.na(all_data$sequencing_depth[i]) && "targetMedian" %in% colnames(ont_sn_matrix)) {
          all_data$sequencing_depth[i] <- ont_sn_matrix$targetMedian[idx[1]]
        }
        counts_values <- c(counts_values, ont_sn_matrix$counts[idx[1]])
        log2_cpm_values <- c(log2_cpm_values, ont_sn_matrix$log2_cpm_plus1[idx[1]])
      }
    }
    
    if(all_data$PB_SN[i] && !is.null(pb_sn_matrix)) {
      idx <- which(pb_sn_matrix$tx_id == tx)
      if(length(idx) > 0) {
        # Fill all metadata columns if not already filled
        for(col in metadata_cols) {
          if(col %in% colnames(pb_sn_matrix) && is.na(all_data[[col]][i])) {
            all_data[[col]][i] <- pb_sn_matrix[[col]][idx[1]]
          }
        }
        # Set sequencing depth if not set
        if(is.na(all_data$sequencing_depth[i]) && "targetMedian" %in% colnames(pb_sn_matrix)) {
          all_data$sequencing_depth[i] <- pb_sn_matrix$targetMedian[idx[1]]
        }
        counts_values <- c(counts_values, pb_sn_matrix$counts[idx[1]])
        log2_cpm_values <- c(log2_cpm_values, pb_sn_matrix$log2_cpm_plus1[idx[1]])
      }
    }
    
    if(all_data$ONT_Bulk[i] && !is.null(ont_bulk_matrix)) {
      idx <- which(ont_bulk_matrix$tx_id == tx)
      if(length(idx) > 0) {
        # Fill all metadata columns if not already filled
        for(col in metadata_cols) {
          if(col %in% colnames(ont_bulk_matrix) && is.na(all_data[[col]][i])) {
            all_data[[col]][i] <- ont_bulk_matrix[[col]][idx[1]]
          }
        }
        # Set sequencing depth if not set
        if(is.na(all_data$sequencing_depth[i]) && "targetMedian" %in% colnames(ont_bulk_matrix)) {
          all_data$sequencing_depth[i] <- ont_bulk_matrix$targetMedian[idx[1]]
        }
        counts_values <- c(counts_values, ont_bulk_matrix$counts[idx[1]])
        log2_cpm_values <- c(log2_cpm_values, ont_bulk_matrix$log2_cpm_plus1[idx[1]])
      }
    }
    
    if(all_data$PB_Bulk[i] && !is.null(pb_bulk_matrix)) {
      idx <- which(pb_bulk_matrix$tx_id == tx)
      if(length(idx) > 0) {
        # Fill all metadata columns if not already filled
        for(col in metadata_cols) {
          if(col %in% colnames(pb_bulk_matrix) && is.na(all_data[[col]][i])) {
            all_data[[col]][i] <- pb_bulk_matrix[[col]][idx[1]]
          }
        }
        # Set sequencing depth if not set
        if(is.na(all_data$sequencing_depth[i]) && "targetMedian" %in% colnames(pb_bulk_matrix)) {
          all_data$sequencing_depth[i] <- pb_bulk_matrix$targetMedian[idx[1]]
        }
        counts_values <- c(counts_values, pb_bulk_matrix$counts[idx[1]])
        log2_cpm_values <- c(log2_cpm_values, pb_bulk_matrix$log2_cpm_plus1[idx[1]])
      }
    }
    
    # Calculate average counts and log2_cpm if values were found
    if(length(counts_values) > 0) {
      all_data$counts[i] <- mean(counts_values)
      all_data$log2_cpm_plus1[i] <- mean(log2_cpm_values)
    }
  }
  
    # Set row names for upset function
  rownames(all_data) <- all_data$tx_id
  
  # Create upset plot with annotations for tx_len, normalized counts, and exon count
  upset_plot <- upset(
    all_data,
    c("ONT_SC", "PB_SC", "ONT_SN", "PB_SN", "ONT_Bulk", "PB_Bulk"),
    name = 'sequencing_modality',
    annotations = list(
      # Transcript length distribution\n(Log10 scale in bp)
      'Transcript Length (log10 bp)' = list(
        aes = aes(x = intersection, y = log10(get(ifelse("tx_len" %in% colnames(all_data), "tx_len", "transcript_length")) + 1)),
        geom = list(
          geom_boxplot(na.rm = TRUE, fill = "lightblue", alpha = 0.7, width = 0.3),
          geom_violin(na.rm = TRUE, fill = "lightblue", alpha = 0.4)
        )
      ),
      # Normalized counts distribution (capped between 5th and 95th percentile)
      'Log10(CPM + 1)\n(5th-95th percentile)' = list(
        aes = aes(x = intersection, y = pmax(pmin(log10((2^log2_cpm_plus1 - 1) + 1), quantile(log10((2^log2_cpm_plus1 - 1) + 1), 0.95, na.rm = TRUE)), quantile(log10((2^log2_cpm_plus1 - 1) + 1), 0.05, na.rm = TRUE))),
        geom = list(
          geom_boxplot(na.rm = TRUE, fill = "lightyellow", alpha = 0.7, width = 0.3),
          geom_violin(na.rm = TRUE, fill = "lightyellow", alpha = 0.4)
        )
      ),
      # Number of exons distribution (capped at 95th percentile for better visualization)
      'Number of Exons\n(capped at 95th percentile)' = list(
        aes = aes(x = intersection, y = pmin(get(ifelse("nexon" %in% colnames(all_data), "nexon", ifelse("num_exons" %in% colnames(all_data), "num_exons", "nexon"))), quantile(get(ifelse("nexon" %in% colnames(all_data), "nexon", ifelse("num_exons" %in% colnames(all_data), "num_exons", "nexon"))), 0.95, na.rm = TRUE))),
        geom = list(
          geom_boxplot(na.rm = TRUE, fill = "lightcoral", alpha = 0.7, outlier.alpha = 0.3)
        )
      )
    ),
    min_size = MIN_INTERSECTION_SIZE,  # Only show intersections with at least this many transcripts
    width_ratio = 0.15,  # Increase width ratio for more space between bars
    sort_sets = FALSE,  # Keep the original order of sets
    stripes = "white",
    base_annotations = list(
      'Intersection size' = intersection_size(
        text = list(size = 0),  # Hide text by setting size to 0
        bar_number_threshold = Inf  # Never show numbers
      ) + 
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10))  # More y-axis breaks
    )
  )
  
  # Get sequencing depth for file naming
  seq_depth <- if(!is.null(ont_sc_matrix)) {
    ont_sc_matrix$sequencing_depth[1]
  } else if(!is.null(pb_sc_matrix)) {
    pb_sc_matrix$sequencing_depth[1]
  } else if(!is.null(ont_sn_matrix)) {
    ont_sn_matrix$sequencing_depth[1]
  } else if(!is.null(pb_sn_matrix)) {
    pb_sn_matrix$sequencing_depth[1]
  } else if(!is.null(ont_bulk_matrix)) {
    ont_bulk_matrix$sequencing_depth[1]
  } else if(!is.null(pb_bulk_matrix)) {
    pb_bulk_matrix$sequencing_depth[1]
  } else {
    "unknown"
  }
  
  # Save the plot with higher resolution and dimensions for the annotations
  # Save as PNG
  output_file_png <- file.path(OUTPUT_DIR, paste0("transcript_overlap_rank", rank_order, "_depth", seq_depth, ".png"))
  png(output_file_png, width = 2800, height = 1600, res = 150)
  print(upset_plot)
  dev.off()
  
  # Save as PDF
  output_file_pdf <- file.path(OUTPUT_DIR, paste0("transcript_overlap_rank", rank_order, "_depth", seq_depth, ".pdf"))
  pdf(output_file_pdf, width = 2800/150, height = 1600/150)  # Convert pixels to inches for PDF
  print(upset_plot)
  dev.off()
  
  # Save the data for further analysis
  saveRDS(all_data, file.path(OUTPUT_DIR, paste0("transcript_intersection_data_rank", rank_order, "_depth", seq_depth, ".rds")))
  
  message(paste0("Completed upset plot for rank ", rank_order, " (depth: ", seq_depth, ")"))
}

# Create a summary table of transcript counts per modality and shared across modalities
summary_data <- data.frame(
  Rank = integer(),
  Sequencing_Depth = character(),
  ONT_SC_Total = integer(),
  PB_SC_Total = integer(),
  ONT_SN_Total = integer(),
  PB_SN_Total = integer(),
  ONT_Bulk_Total = integer(),
  PB_Bulk_Total = integer(),
  All_Modalities_Shared = integer(),
  stringsAsFactors = FALSE
)

# Add additional summary data about transcript properties
transcript_properties <- data.frame(
  Rank = integer(),
  Sequencing_Depth = character(),
  Intersection = character(),
  Count = integer(),
  Median_Length = numeric(),
  Median_Counts = numeric(),
  Median_Log2CPM = numeric(),
  stringsAsFactors = FALSE
)

# Create a nested list to store transcript IDs for each intersection set
tx_intersection_ids <- list()

for(rank_order in available_ranks) {
  # Read the saved intersection data
  data_files <- list.files(OUTPUT_DIR, pattern = paste0("transcript_intersection_data_rank", rank_order, "_"), full.names = TRUE)
  
  if(length(data_files) == 0) {
    message(paste("No data files found for rank", rank_order))
    next
  }
  
  intersection_data <- readRDS(data_files[1])
  
  # Extract sequencing depth from filename
  seq_depth <- stringr::str_extract(data_files[1], "depth([^.]+)")
  seq_depth <- stringr::str_replace(seq_depth, "depth", "")
  
  # Create a list for this rank
  tx_intersection_ids[[paste0("rank", rank_order)]] <- list()
  
  # Print summary statistics for transcript lengths
  message(paste0("\nSummary statistics for transcript lengths at rank ", rank_order, " (depth: ", seq_depth, "):"))
  print(summary(intersection_data$tx_len))
  message(paste0("Number of transcripts: ", nrow(intersection_data)))
  
  # Count transcripts in each modality
  ont_sc_count <- sum(intersection_data$ONT_SC, na.rm = TRUE)
  pb_sc_count <- sum(intersection_data$PB_SC, na.rm = TRUE)
  ont_sn_count <- sum(intersection_data$ONT_SN, na.rm = TRUE)
  pb_sn_count <- sum(intersection_data$PB_SN, na.rm = TRUE)
  ont_bulk_count <- sum(intersection_data$ONT_Bulk, na.rm = TRUE)
  pb_bulk_count <- sum(intersection_data$PB_Bulk, na.rm = TRUE)
  
  # Count transcripts shared across all modalities
  shared_count <- sum(intersection_data$ONT_SC & intersection_data$PB_SC & 
                     intersection_data$ONT_SN & intersection_data$PB_SN &
                     intersection_data$ONT_Bulk & intersection_data$PB_Bulk, na.rm = TRUE)
  
  # Add to summary data
  summary_data <- rbind(summary_data, data.frame(
    Rank = rank_order,
    Sequencing_Depth = seq_depth,
    ONT_SC_Total = ont_sc_count,
    PB_SC_Total = pb_sc_count,
    ONT_SN_Total = ont_sn_count,
    PB_SN_Total = pb_sn_count,
    ONT_Bulk_Total = ont_bulk_count,
    PB_Bulk_Total = pb_bulk_count,
    All_Modalities_Shared = shared_count,
    stringsAsFactors = FALSE
  ))
  
     # Store transcript IDs for each modality (with consistent data structure)
   modality_names <- c("ONT_SC_Total", "PB_SC_Total", "ONT_SN_Total", "PB_SN_Total", "ONT_Bulk_Total", "PB_Bulk_Total")
   modality_filters <- list(
     intersection_data$ONT_SC,
     intersection_data$PB_SC,
     intersection_data$ONT_SN,
     intersection_data$PB_SN,
     intersection_data$ONT_Bulk,
     intersection_data$PB_Bulk
   )
   
   for(i in 1:length(modality_names)) {
     subset_data <- intersection_data[modality_filters[[i]], ]
     if(nrow(subset_data) > 0) {
       tx_intersection_ids[[paste0("rank", rank_order)]][[modality_names[i]]] <- list(
         tx_ids = subset_data$tx_id,
         metadata = subset_data[, !colnames(subset_data) %in% c("ONT_SC", "PB_SC", "ONT_SN", "PB_SN", "ONT_Bulk", "PB_Bulk")]
       )
     }
   }
   
   # All modalities shared
   shared_filter <- intersection_data$ONT_SC & intersection_data$PB_SC & 
                   intersection_data$ONT_SN & intersection_data$PB_SN &
                   intersection_data$ONT_Bulk & intersection_data$PB_Bulk
   shared_subset_data <- intersection_data[shared_filter, ]
   if(nrow(shared_subset_data) > 0) {
     tx_intersection_ids[[paste0("rank", rank_order)]][["All_Modalities_Shared"]] <- list(
       tx_ids = shared_subset_data$tx_id,
       metadata = shared_subset_data[, !colnames(shared_subset_data) %in% c("ONT_SC", "PB_SC", "ONT_SN", "PB_SN", "ONT_Bulk", "PB_Bulk")]
     )
   }
  
  # Calculate ALL possible intersection sets and find the top 25 largest
  message(paste0("Calculating all intersection sets for rank ", rank_order, "..."))
  
  modality_cols <- c("ONT_SC", "PB_SC", "ONT_SN", "PB_SN", "ONT_Bulk", "PB_Bulk")
  
  # Get all possible intersection combinations (2^6 - 1 = 63 combinations)
  all_combinations <- list()
  for(i in 1:length(modality_cols)) {
    combinations_i <- combn(modality_cols, i, simplify = FALSE)
    all_combinations <- c(all_combinations, combinations_i)
  }
  
  # Calculate intersection sizes and create proper names
  intersection_info <- data.frame(
    combination = character(),
    size = integer(),
    intersection_name = character(),
    stringsAsFactors = FALSE
  )
  
  for(combo in all_combinations) {
    # Create boolean vector for this exact intersection
    # Must be present in ALL modalities in the combination
    # AND absent from all other modalities
    in_combo <- rep(TRUE, nrow(intersection_data))
    out_combo <- rep(TRUE, nrow(intersection_data))
    
    # Must be present in all modalities in the combination
    for(mod in combo) {
      in_combo <- in_combo & intersection_data[[mod]]
    }
    
    # Must be absent from all modalities NOT in the combination
    for(mod in setdiff(modality_cols, combo)) {
      out_combo <- out_combo & !intersection_data[[mod]]
    }
    
    # Final intersection: present in combo modalities AND absent from others
    final_intersection <- in_combo & out_combo
    intersection_size <- sum(final_intersection, na.rm = TRUE)
    
    # Create a readable name for this intersection
    intersection_name <- if(length(combo) == 1) {
      paste0(combo[1], "_only")
    } else if(length(combo) == length(modality_cols)) {
      "All_Modalities_Shared"
    } else if(length(combo) == 2) {
      paste(combo, collapse = "_and_")
    } else if(length(combo) == 3) {
      paste(combo, collapse = "_and_")
    } else if(length(combo) == 4) {
      paste(combo, collapse = "_and_")
    } else if(length(combo) == 5) {
      # For 5-way intersections, name by what's excluded
      excluded <- setdiff(modality_cols, combo)
      paste0("All_except_", excluded[1])
    } else {
      paste(combo, collapse = "_and_")
    }
    
    if(intersection_size >= MIN_INTERSECTION_SIZE) {
      intersection_info <- rbind(intersection_info, data.frame(
        combination = paste(combo, collapse = "&"),
        size = intersection_size,
        intersection_name = intersection_name,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Sort by size and get top N (controlled by MAX_INTERSECTIONS_TO_SAVE)
  intersection_info <- intersection_info[order(intersection_info$size, decreasing = TRUE), ]
  top_intersections <- head(intersection_info, min(MAX_INTERSECTIONS_TO_SAVE, nrow(intersection_info)))
  
  message(paste("Found", nrow(intersection_info), "intersection sets >=", MIN_INTERSECTION_SIZE, "transcripts for rank", rank_order))
  message(paste("Saving top", nrow(top_intersections), "largest intersection sets:"))
  print(top_intersections[, c("intersection_name", "size")])
  
  # Store transcript IDs for all top intersection sets
  for(i in 1:nrow(top_intersections)) {
    combo <- strsplit(top_intersections$combination[i], "&")[[1]]
    intersection_name <- top_intersections$intersection_name[i]
    
    # Recreate the boolean vector for this intersection
    in_combo <- rep(TRUE, nrow(intersection_data))
    out_combo <- rep(TRUE, nrow(intersection_data))
    
    for(mod in combo) {
      in_combo <- in_combo & intersection_data[[mod]]
    }
    
    for(mod in setdiff(modality_cols, combo)) {
      out_combo <- out_combo & !intersection_data[[mod]]
    }
    
    final_intersection <- in_combo & out_combo
    subset_data <- intersection_data[final_intersection, ]
    
    if(nrow(subset_data) > 0) {
      # Add to transcript properties summary
      transcript_properties <- rbind(transcript_properties, data.frame(
        Rank = rank_order,
        Sequencing_Depth = seq_depth,
        Intersection = intersection_name,
        Count = nrow(subset_data),
        Median_Length = median(subset_data$tx_len, na.rm = TRUE),
        Median_Counts = median(subset_data$counts, na.rm = TRUE),
        Median_Log2CPM = median(subset_data$log2_cpm_plus1, na.rm = TRUE),
        stringsAsFactors = FALSE
      ))
      
      # Store transcript IDs AND full metadata for this intersection
      tx_intersection_ids[[paste0("rank", rank_order)]][[intersection_name]] <- list(
        tx_ids = subset_data$tx_id,
        metadata = subset_data[, !colnames(subset_data) %in% c("ONT_SC", "PB_SC", "ONT_SN", "PB_SN", "ONT_Bulk", "PB_Bulk")]
      )
    }
  }
  
  # Also save the intersection summary for this rank
  write.csv(top_intersections, 
           file.path(OUTPUT_DIR, paste0("intersection_summary_rank", rank_order, "_depth", seq_depth, ".csv")), 
           row.names = FALSE)
}

# Save summary tables
write.csv(summary_data, file.path(OUTPUT_DIR, "transcript_counts_summary.csv"), row.names = FALSE)
write.csv(transcript_properties, file.path(OUTPUT_DIR, "transcript_properties_summary.csv"), row.names = FALSE)

# Save the transcript IDs for each intersection set
saveRDS(tx_intersection_ids, file.path(OUTPUT_DIR, "transcript_intersection_ids.rds"))
message(paste0("Saved transcript IDs for top ", MAX_INTERSECTIONS_TO_SAVE, " intersection sets per rank to transcript_intersection_ids.rds"))
message("Individual intersection summaries saved as CSV files for each rank")

# Create summary plots
# 1. Bar plot of transcript counts by modality and rank
summary_long <- tidyr::pivot_longer(
  summary_data,
  cols = c("ONT_SC_Total", "PB_SC_Total", "ONT_SN_Total", "PB_SN_Total", 
          "ONT_Bulk_Total", "PB_Bulk_Total", "All_Modalities_Shared"),
  names_to = "Modality",
  values_to = "Count"
)

# Order the Rank factor
summary_long$Rank <- factor(summary_long$Rank, levels = sort(unique(summary_long$Rank)))

counts_plot <- ggplot(summary_long, aes(x = Rank, y = Count, fill = Modality)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(
    title = "Transcript Counts by Sequencing Modality and Rank Order",
    x = "Sequencing Depth Rank",
    y = "Number of Transcripts"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the summary plot
ggsave(file.path(OUTPUT_DIR, "transcript_counts_summary.png"), 
       counts_plot, width = 12, height = 8, dpi = 150)

# Create a comprehensive intersection summary across all ranks
create_comprehensive_intersection_summary <- function() {
  message("Creating comprehensive intersection summary across all ranks...")
  
  all_intersection_files <- list.files(OUTPUT_DIR, pattern = "intersection_summary_rank.*\\.csv$", full.names = TRUE)
  
  if(length(all_intersection_files) == 0) {
    message("No intersection summary files found")
    return()
  }
  
  # Combine all intersection summaries
  comprehensive_summary <- data.frame()
  
  for(file in all_intersection_files) {
    # Extract rank and depth from filename
    filename <- basename(file)
    rank_num <- stringr::str_extract(filename, "rank(\\d+)")
    rank_num <- stringr::str_replace(rank_num, "rank", "")
    depth_str <- stringr::str_extract(filename, "depth([^.]+)")
    depth_str <- stringr::str_replace(depth_str, "depth", "")
    
    # Read the file
    rank_summary <- read.csv(file, stringsAsFactors = FALSE)
    rank_summary$Rank <- as.numeric(rank_num)
    rank_summary$Depth <- depth_str
    
    comprehensive_summary <- rbind(comprehensive_summary, rank_summary)
  }
  
  # Save comprehensive summary
  write.csv(comprehensive_summary, 
           file.path(OUTPUT_DIR, "comprehensive_intersection_summary_all_ranks.csv"), 
           row.names = FALSE)
  
  message(paste("Comprehensive intersection summary saved with", nrow(comprehensive_summary), "total intersection sets across all ranks"))
  
  # Create a summary of intersection types across ranks
  intersection_type_summary <- comprehensive_summary %>%
    group_by(intersection_name) %>%
    summarise(
      appears_in_n_ranks = n(),
      min_size = min(size),
      max_size = max(size),
      mean_size = round(mean(size), 1),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_size))
  
  write.csv(intersection_type_summary, 
           file.path(OUTPUT_DIR, "intersection_types_summary_across_ranks.csv"), 
           row.names = FALSE)
  
  message("Intersection types summary across ranks saved")
}

# Generate the comprehensive summary
create_comprehensive_intersection_summary()

message("All upset plots, summary data, and summary plots have been created successfully!")
message(paste0("✓ Top ", MAX_INTERSECTIONS_TO_SAVE, " intersection sets (≥", MIN_INTERSECTION_SIZE, " transcripts) saved per rank"))
message("✓ Saved intersection sets match exactly what's visualized in upset plots")
message("✓ FULL metadata preserved for all transcript IDs (including all original annotations)")
message("✓ Comprehensive intersection summaries created")
message("✓ All intersection data properly organized by rank order")

# Display what metadata is available
if(length(available_ranks) > 0) {
  sample_rank <- names(tx_intersection_ids)[1]
  sample_intersection <- names(tx_intersection_ids[[sample_rank]])[1] 
  if(is.list(tx_intersection_ids[[sample_rank]][[sample_intersection]]) && 
     "metadata" %in% names(tx_intersection_ids[[sample_rank]][[sample_intersection]])) {
    metadata_cols <- colnames(tx_intersection_ids[[sample_rank]][[sample_intersection]]$metadata)
    message("\n📊 Available metadata columns for each transcript:")
    message(paste("   ", paste(metadata_cols, collapse = ", ")))
    message("\n💡 Access example: tx_intersection_ids[['rank1']][['ONT_SC_Total']]$metadata")
    message("💡 Get transcript IDs: tx_intersection_ids[['rank1']][['ONT_SC_Total']]$tx_ids")
    message("💡 For intersection sets: tx_intersection_ids[['rank1']][['All_Modalities_Shared']]$metadata")
  }
}

#######################################################
# Custom function to create upset plot for specific rank and top X intersections
#######################################################

create_custom_upset_plot <- function(rank_order, top_x_intersections = 15, 
                                     min_size = 10, save_plot = TRUE, 
                                     show_numbers = FALSE) {
  
  # Validate inputs
  if (!rank_order %in% available_ranks) {
    stop(paste("Rank", rank_order, "not available. Available ranks:", paste(available_ranks, collapse = ", ")))
  }
  
  if (top_x_intersections < 1) {
    stop("top_x_intersections must be at least 1")
  }
  
  message(paste0("Creating custom upset plot for rank ", rank_order, " with top ", top_x_intersections, " intersection sets..."))
  if (show_numbers) {
    message("Numbers will be displayed on intersection bars")
  }
  
  # Load the intersection data for the specified rank (same as original)
  data_files <- list.files(OUTPUT_DIR, pattern = paste0("transcript_intersection_data_rank", rank_order, "_"), full.names = TRUE)
  
  if(length(data_files) == 0) {
    stop(paste("No data files found for rank", rank_order))
  }
  
  # Use the same data structure as the original code
  all_data <- readRDS(data_files[1])
  
  # Extract sequencing depth from filename
  seq_depth <- stringr::str_extract(data_files[1], "depth([^.]+)")
  seq_depth <- stringr::str_replace(seq_depth, "depth", "")
  
  # Set row names for upset function (same as original)
  rownames(all_data) <- all_data$tx_id
  
  # First, create a temporary upset plot to identify intersection sizes
  # Use the exact same parameters as the original, but with min_size = 1 to capture all intersections
  temp_upset <- upset(
    all_data,
    c("ONT_SC", "PB_SC", "ONT_SN", "PB_SN", "ONT_Bulk", "PB_Bulk"),
    min_size = 1,  # Capture all intersections
    sort_intersections_by = "cardinality"  # Sort by size
  )
  
  # Extract intersection information from the temp plot
  # This is a bit tricky - we need to get the intersection data from ComplexUpset
  
  # Alternative approach: Calculate intersections manually but using ComplexUpset logic
  modality_cols <- c("ONT_SC", "PB_SC", "ONT_SN", "PB_SN", "ONT_Bulk", "PB_Bulk")
  
  # Get all possible intersection combinations
  all_combinations <- list()
  for(i in 1:length(modality_cols)) {
    combinations_i <- combn(modality_cols, i, simplify = FALSE)
    all_combinations <- c(all_combinations, combinations_i)
  }
  
  # Calculate intersection sizes using ComplexUpset's logic (presence in specified sets only)
  intersection_info <- data.frame(
    combination = character(),
    size = integer(),
    transcript_indices = I(list()),
    stringsAsFactors = FALSE
  )
  
  for(combo in all_combinations) {
    # Create boolean vector for this exact intersection
    # Must be present in ALL modalities in the combination
    # AND absent from all other modalities
    in_combo <- rep(TRUE, nrow(all_data))
    out_combo <- rep(TRUE, nrow(all_data))
    
    # Must be present in all modalities in the combination
    for(mod in combo) {
      in_combo <- in_combo & all_data[[mod]]
    }
    
    # Must be absent from all modalities NOT in the combination
    for(mod in setdiff(modality_cols, combo)) {
      out_combo <- out_combo & !all_data[[mod]]
    }
    
    # Final intersection: present in combo modalities AND absent from others
    final_intersection <- in_combo & out_combo
    intersection_size <- sum(final_intersection, na.rm = TRUE)
    
    if(intersection_size >= min_size) {
      intersection_info <- rbind(intersection_info, data.frame(
        combination = paste(combo, collapse = "&"),
        size = intersection_size,
        transcript_indices = I(list(which(final_intersection))),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Sort by size and get top X intersections
  intersection_info <- intersection_info[order(intersection_info$size, decreasing = TRUE), ]
  top_intersections <- head(intersection_info, top_x_intersections)
  
  message(paste("Found", nrow(intersection_info), "intersection sets >=", min_size, "transcripts"))
  message(paste("Showing top", nrow(top_intersections), "intersection sets:"))
  print(top_intersections[, c("combination", "size")])
  
  # Get all transcript indices that belong to the top intersections
  top_transcript_indices <- unique(unlist(top_intersections$transcript_indices))
  
  # Filter the data to only include these transcripts
  filtered_all_data <- all_data[top_transcript_indices, ]
  
  # Set row names for the filtered data
  rownames(filtered_all_data) <- filtered_all_data$tx_id
  
  message(paste("Filtered to", nrow(filtered_all_data), "transcripts for plotting"))
  
  # Create the upset plot using the EXACT same structure as the original
  upset_plot <- upset(
    filtered_all_data,
    c("ONT_SC", "PB_SC", "ONT_SN", "PB_SN", "ONT_Bulk", "PB_Bulk"),
    name = 'sequencing_modality',
    annotations = list(
      # Transcript length distribution (same as original)
      'Transcript Length (log10 bp)' = list(
        aes = aes(x = intersection, y = log10(tx_len + 1)),
        geom = list(
          geom_boxplot(na.rm = TRUE, fill = "lightblue", alpha = 0.7, width = 0.3),
          geom_violin(na.rm = TRUE, fill = "lightblue", alpha = 0.4)
        )
      ),
      # Normalized counts distribution (same as original)
      'Log10(CPM + 1)\n(5th-95th percentile)' = list(
        aes = aes(x = intersection, y = pmax(pmin(log10((2^log2_cpm_plus1 - 1) + 1), quantile(log10((2^log2_cpm_plus1 - 1) + 1), 0.95, na.rm = TRUE)), quantile(log10((2^log2_cpm_plus1 - 1) + 1), 0.05, na.rm = TRUE))),
        geom = list(
          geom_boxplot(na.rm = TRUE, fill = "lightyellow", alpha = 0.7, width = 0.3),
          geom_violin(na.rm = TRUE, fill = "lightyellow", alpha = 0.4)
        )
      ),
      # Number of exons distribution (same as original)
      'Number of Exons\n(capped at 95th percentile)' = list(
        aes = aes(x = intersection, y = pmin(nexon, quantile(nexon, 0.95, na.rm = TRUE))),
        geom = list(
          geom_boxplot(na.rm = TRUE, fill = "lightcoral", alpha = 0.7, outlier.alpha = 0.3)
        )
      )
    ),
    min_size = min_size,  # Use the provided min_size parameter
    width_ratio = 0.15,  # Same as original
    sort_sets = FALSE,  # Same as original
    stripes = "white",  # Same as original
    base_annotations = list(
      'Intersection size' = intersection_size(
        text = if(show_numbers) list(size = 2, color = "black") else list(size = 0),  # Toggleable text display
        bar_number_threshold = if(show_numbers) 1 else Inf  # Toggleable number threshold
      ) + 
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10))  # Same as original
    )
  )
  
  # Save the plot if requested (same structure as original)
  if(save_plot) {
    # Save as PNG (same dimensions as original)
    output_file_png <- file.path(OUTPUT_DIR, paste0("custom_transcript_overlap_rank", rank_order, 
                                                   "_depth", seq_depth, "_top", top_x_intersections, ".png"))
    png(output_file_png, width = 2800, height = 1600, res = 150)
    print(upset_plot)
    dev.off()
    
    # Save as PDF (same dimensions as original)
    output_file_pdf <- file.path(OUTPUT_DIR, paste0("custom_transcript_overlap_rank", rank_order, 
                                                   "_depth", seq_depth, "_top", top_x_intersections, ".pdf"))
    pdf(output_file_pdf, width = 2800/150, height = 1600/150)
    print(upset_plot)
    dev.off()
    
    # Save the filtered data for further analysis (same as original)
    saveRDS(filtered_all_data, file.path(OUTPUT_DIR, paste0("custom_transcript_intersection_data_rank", rank_order, "_depth", seq_depth, "_top", top_x_intersections, ".rds")))
    
    message(paste0("Custom upset plot saved to: ", output_file_png, " and ", output_file_pdf))
  }
  
  # Return the plot and summary data
  return(list(
    plot = upset_plot,
    intersection_summary = top_intersections[, c("combination", "size")],
    filtered_data = filtered_all_data,
    total_transcripts = nrow(filtered_all_data),
    all_intersection_info = intersection_info[, c("combination", "size")]
  ))
}

#######################################################
# Example usage of the custom function
#######################################################

# Example: Create upset plot for rank 6 with top 15 intersection sets
# result <- create_custom_upset_plot(rank_order = 6, top_x_intersections = 15, min_size = 10)
# print(result$plot)

# Example: Create upset plot with numbers displayed on intersection bars

result <- create_custom_upset_plot(rank_order = 6, 
                                   top_x_intersections = 19, 
                                   show_numbers = TRUE)
