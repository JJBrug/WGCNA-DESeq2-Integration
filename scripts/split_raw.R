# ==============================================================================
#   Helper: Split Raw Counts & Metadata by Category
# ==============================================================================
#   Usage: 
#   Run this once to generate the input files needed to run General WGCNA
#   script on individual datasets (e.g., split by year, tissue, treatment).

library(tidyverse)

# --- CONFIGURATION ---
input_counts  <- "rice_raw_counts.csv"
input_traits  <- "rice_metadata.csv"
output_folder <- "split_data"
split_col     <- "Treatment"      # The column in metadata to split by
file_prefix   <- "rice_"      # Optional prefix for output files

# Create output directory
if (!dir.exists(output_folder)) dir.create(output_folder)

# 1. Load Data
# Auto-detect file type for counts
if (endsWith(input_counts, ".csv")) {
  counts <- read.csv(input_counts, header = TRUE, row.names = 1, check.names = FALSE)
} else if (endsWith(input_counts, ".tsv") || endsWith(input_counts, ".txt")) {
  counts <- read.table(input_counts, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
} else {
  stop("Unsupported file format for counts. Please provide a .csv, .tsv, or .txt file.")
}

# Auto-detect file type for traits
if (endsWith(input_traits, ".csv")) {
  traits <- read.csv(input_traits, row.names = 1, check.names = FALSE)
} else if (endsWith(input_traits, ".tsv") || endsWith(input_traits, ".txt")) {
  traits <- read.table(input_traits, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
} else {
  stop("Unsupported file format for traits. Please provide a .csv, .tsv, or .txt file.")
}

# 2. Identify Unique Categories
if (!split_col %in% colnames(traits)) stop(paste("Column", split_col, "not found in metadata!"))
categories <- unique(traits[[split_col]])

cat(paste("Found", length(categories), "groups in column", split_col, ":", paste(categories, collapse=", "), "\n"))

# 3. Loop and Split
for (cat_val in categories) {
  
  # Identify samples for this category
  samples_in_cat <- rownames(traits)[traits[[split_col]] == cat_val]
  
  # Intersect with counts (ensure we only grab samples that exist in both)
  common_samples <- intersect(samples_in_cat, colnames(counts))
  
  if (length(common_samples) == 0) {
    warning(paste("No matching samples found for group", cat_val, "- Skipping."))
    next
  }
  
  # Subset
  counts_sub <- counts[, common_samples]
  traits_sub <- traits[common_samples, ]
  
  # Define Filenames 
  out_name <- paste0(file_prefix, cat_val) 
  
  file_counts <- file.path(output_folder, paste0("raw_counts_", out_name, ".tsv"))
  file_traits <- file.path(output_folder, paste0("metadata_", out_name, ".csv"))
  
  # Write to disk
  write.table(counts_sub, file_counts, sep="\t", quote=FALSE, col.names=NA)
  write.csv(traits_sub, file_traits, quote=FALSE)
  
  cat(paste("  -> Saved:", out_name, "| Samples:", length(common_samples), "\n"))
}
