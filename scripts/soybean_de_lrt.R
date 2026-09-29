# ==============================================================================
#   A Given Example of Potential Differential Expression Analysis Methods
#   Replicates the HSP Paper Methodology (LRT) for WGCNA Integration
# ==============================================================================

library(DESeq2)
library(tidyverse)
options(stringsAsFactors = FALSE)

# --- CONFIGURATION ---
file_raw_counts <- "raw_counts.tsv" # DESeq2 requires RAW, un-normalized counts
file_metadata   <- "metadata.csv"

col_region <- "region"
col_year   <- "year"
group_east <- "east" # Control/Denominator
group_west <- "west" # Treatment/Numerator

# --- LOAD & FORMAT DATA ---
counts_raw <- read.table(file_raw_counts, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
meta       <- read.csv(file_metadata, row.names = 1, check.names = FALSE)

# Ensure samples match
common_samples <- intersect(colnames(counts_raw), rownames(meta))
counts_raw     <- counts_raw[, common_samples]
meta           <- meta[common_samples, ]

# Ensure metadata columns are factors
meta[[col_region]] <- factor(meta[[col_region]], levels = c(group_east, group_west))
meta[[col_year]]   <- factor(meta[[col_year]])

# --- PRE-FILTERING ---
# Filter out low-abundance genes (>= 10 counts in at least 3 samples) as per manuscript
keep <- rowSums(counts_raw >= 10) >= 3
counts_filt <- counts_raw[keep, ]

# --- DESEQ2 LIKELIHOOD RATIO TEST (LRT) ---
# Full Model: ~ Year + Region (Tests the effect of Region while controlling for Year)
dds <- DESeqDataSetFromMatrix(countData = counts_filt, 
                              colData = meta, 
                              design = as.formula(paste0("~ ", col_year, " + ", col_region)))

# Run LRT comparing the full model to the reduced model (~ Year)
dds <- DESeq(dds, test = "LRT", reduced = as.formula(paste0("~ ", col_year)))

# Extract Results specifically for the Region contrast (West vs East)
res <- results(dds, contrast = c(col_region, group_west, group_east))

# --- FORMAT OUTPUT FOR ENRICHMENT SCRIPT ---
df_res <- as.data.frame(res)
df_res$GeneID <- rownames(df_res)

# Select and reorder only the necessary columns (matching Wilcoxon format)
df_res <- df_res %>%
  select(GeneID, log2FoldChange, pvalue, padj)

# Save to CSV
write.csv(df_res, "de_results/DE_results_LRT_West_vs_East.csv", row.names = FALSE)

cat("\nDESeq2 LRT analysis complete. Found",
    sum(df_res$padj < 0.01 & abs(df_res$log2FoldChange) > 1,
        na.rm = TRUE), "significant genes based on strict thresholds.\n")
