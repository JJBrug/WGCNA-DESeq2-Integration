# ==============================================================================
#   Interaction Differential Expression Analysis (LRT)
#   Isolating Stage-Specific Rewiring for WGCNA Integration
# ==============================================================================

library(DESeq2)
library(tidyverse)
options(stringsAsFactors = FALSE)

# --- CONFIGURATION ---
file_raw_counts <- "rice_raw_counts.csv" 
file_metadata   <- "rice_metadata.csv"

col_treatment <- "Treatment"
col_time      <- "Time"
group_control <- "Control" 
group_severe  <- "Severe_Drought" 

# --- LOAD & FORMAT DATA ---
counts_raw <- read.csv(file_raw_counts, row.names = 1, check.names = FALSE)
meta       <- read.csv(file_metadata, row.names = 1, check.names = FALSE)

# Ensure samples match
common_samples <- intersect(colnames(counts_raw), rownames(meta))
counts_raw     <- counts_raw[, common_samples]
meta           <- meta[common_samples, ]

# Ensure metadata columns are factors
meta[[col_treatment]] <- factor(meta[[col_treatment]], levels = c(group_control, group_severe))
meta[[col_time]]      <- factor(meta[[col_time]])

# --- PRE-FILTERING ---
# Filter out low-abundance genes (>= 10 counts in at least 3 samples)
keep <- rowSums(counts_raw >= 10) >= 3
counts_filt <- counts_raw[keep, ]

# --- DESEQ2 INTERACTION LIKELIHOOD RATIO TEST (LRT) ---
dds <- DESeqDataSetFromMatrix(countData = counts_filt, 
                              colData = meta, 
                              design = as.formula(paste0("~ ", col_time, " + ", col_treatment, " + ", col_time, ":", col_treatment)))

# Run LRT comparing the full model to the reduced model (main effects only)
dds <- DESeq(dds, test = "LRT", reduced = as.formula(paste0("~ ", col_time, " + ", col_treatment)))

# Extract Results 
res <- results(dds)

# --- FORMAT OUTPUT FOR ENRICHMENT SCRIPT ---
df_res <- as.data.frame(res)
df_res$GeneID <- rownames(df_res)

# Select and reorder (Note: LFC in an interaction model represents the effect at a specific reference level, but padj governs significance)
df_res <- df_res %>%
  select(GeneID, log2FoldChange, pvalue, padj)

# Save to CSV
if (!dir.exists("de_results")) dir.create("de_results")
write.csv(df_res, "de_results/DE_results_LRT_Interaction.csv", row.names = FALSE)

cat("\nDESeq2 Interaction LRT analysis complete. Found",
    sum(df_res$padj < 0.05, na.rm = TRUE), "genes with significantly altered temporal trajectories.\n")