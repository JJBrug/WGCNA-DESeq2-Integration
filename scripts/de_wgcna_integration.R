# ==============================================================================
#   WGCNA & Differential Expression (DE) Integration 
#   (Module Overrepresentation Analysis via Hypergeometric Test)
# ==============================================================================
#   
#   Purpose: 
#   Calculates whether specific WGCNA modules contain a significantly higher 
#   proportion of Differentially Expressed (DE) genes than expected by chance.
#
#   Usage:
#   1. Run your DE analysis (DESeq2, edgeR, etc.) entirely separately.
#   2. Export the DE results as a CSV.
#   3. Point this script to that CSV and your previously generated WGCNA outputs.
#
# ==============================================================================

library(tidyverse)
library(ggplot2)
options(stringsAsFactors = FALSE)

# ==============================================================================
#   1. CONFIGURATION
# ==============================================================================

output_dir    <- "integration_results"
output_prefix <- "WGCNA_DE_Enrichment"
if (!dir.exists(output_dir)) dir.create(output_dir)

# --- WGCNA Inputs ---
is_consensus <- FALSE  # Set to FALSE if analyzing a single-environment network

# Point to your expression and color .rds files
file_wgcna_expr   <- "general_results/rice_Severe_Drought_datExpr.rds"
file_wgcna_colors <- "general_results/rice_Severe_Drought_moduleColors.rds"

# --- DE Inputs ---
# Point to your DESeq2/edgeR/limma results table.
file_de_results   <- "de_results/DE_results_LRT_Interaction.csv"

# Tell the script what your DE columns are named
col_gene_id <- "GeneID"          
col_padj    <- "padj"            
col_lfc     <- "log2FoldChange"  

# Define what constitutes a "Significant DE Gene"
padj_threshold <- 0.05
lfc_threshold  <- 1.0  

# ==============================================================================
#   2. DATA INGESTION & FORMATTING
# ==============================================================================

# 1. Load WGCNA Data
if(!file.exists(file_wgcna_expr)) stop("WGCNA expression file missing")
if(!file.exists(file_wgcna_colors)) stop("WGCNA colors file missing")

expr_obj     <- readRDS(file_wgcna_expr)
moduleColors <- readRDS(file_wgcna_colors)

# Dynamically extract gene names based on network type
if (is_consensus) {
  # Consensus structure: list of lists
  wgcna_genes <- colnames(expr_obj[[1]]$data)
} else {
  # Single network structure: flat matrix
  wgcna_genes <- colnames(expr_obj)
}

# Combine into a dataframe
df_wgcna <- data.frame(
  GeneID = wgcna_genes,
  Module = moduleColors
)

# 2. Load DE Results
if(!file.exists(file_de_results)) stop("DE results file missing")
df_de <- read.csv(file_de_results)

# Ensure necessary DE columns exist
required_cols <- c(col_gene_id, col_padj, col_lfc)
missing_cols <- setdiff(required_cols, colnames(df_de))
if(length(missing_cols) > 0) {
  stop(paste("Missing DE columns:", paste(missing_cols, collapse = ", ")))
}

# ==============================================================================
#   3. DEFINE THE STATISTICAL BACKGROUND
# ==============================================================================
# CRITICAL: The hypergeometric test background (the "urn") must ONLY contain 
# genes that were tested in BOTH WGCNA and DESeq2. 

# Merge datasets
df_merged <- merge(df_wgcna, df_de, by.x = "GeneID", by.y = col_gene_id)

if(nrow(df_merged) < 1000) {
  warning("Very few common genes found between WGCNA and DE results. Check Gene IDs.")
}

# Clean out NAs from DE testing (genes dropped by DESeq2 independent filtering)
df_merged <- df_merged %>% filter(!is.na(!!sym(col_padj)))

# Define binary DE status
df_merged <- df_merged %>%
  mutate(
    is_DE = ifelse(!!sym(col_padj) < padj_threshold & abs(!!sym(col_lfc)) > lfc_threshold, 
                   TRUE, FALSE)
  )

total_background <- nrow(df_merged)
total_DE         <- sum(df_merged$is_DE)
total_non_DE     <- total_background - total_DE

cat(paste("\nBackground Universe:", total_background, "genes"))
cat(paste("\nTotal DE Genes in Universe:", total_DE, "\n"))

# ==============================================================================
#   4. HYPERGEOMETRIC ENRICHMENT TESTING
# ==============================================================================

modules <- unique(df_merged$Module)
enrichment_results <- list()

for (mod in modules) {
  # Skip the grey (unassigned) module
  if (mod %in% c("grey", "gold")) next
  
  # Module specific stats
  mod_genes <- df_merged %>% filter(Module == mod)
  mod_size  <- nrow(mod_genes)
  mod_DE    <- sum(mod_genes$is_DE)
  
  # Hypergeometric Test Math:
  # q = successes in sample (minus 1 for lower.tail = FALSE)
  # m = total successes in population
  # n = total failures in population
  # k = sample size
  
  if (mod_DE > 0) {
    p_val <- phyper(q = mod_DE - 1, 
                    m = total_DE, 
                    n = total_non_DE, 
                    k = mod_size, 
                    lower.tail = FALSE)
  } else {
    p_val <- 1
  }
  
  enrichment_results[[mod]] <- data.frame(
    Module         = mod,
    Module_Size    = mod_size,
    DE_in_Module   = mod_DE,
    Expected_DE    = (total_DE / total_background) * mod_size,
    Enrichment_P   = p_val
  )
}

# Compile and FDR correct
df_enrich <- do.call(rbind, enrichment_results)
df_enrich$Enrichment_FDR <- p.adjust(df_enrich$Enrichment_P, method = "BH")

# Sort by significance
df_enrich <- df_enrich %>% arrange(Enrichment_FDR)

write.csv(df_enrich, file.path(output_dir, paste0(output_prefix, "_Results.csv")), row.names = FALSE)

# ==============================================================================
#   5. VISUALIZATION
# ==============================================================================

# Filter for plotting (only show modules with at least 1 DE gene)
df_plot <- df_enrich %>% filter(DE_in_Module > 0)

# Calculate -log10 FDR for scaling
df_plot$MinusLog10FDR <- -log10(df_plot$Enrichment_FDR)

# Force actual colors onto the plot
mod_colors_map <- setNames(df_plot$Module, df_plot$Module)

p <- ggplot(df_plot, aes(x = reorder(Module, MinusLog10FDR), y = MinusLog10FDR, fill = Module)) +
  geom_bar(stat = "identity", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  coord_flip() +
  scale_fill_manual(values = mod_colors_map) +
  labs(
    title = "WGCNA Module Enrichment for DE Genes",
    subtitle = "Red dashed line indicates FDR = 0.05",
    x = "Module",
    y = "-log10(FDR Adjusted P-value)"
  ) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(file.path(output_dir, paste0(output_prefix, "_Barplot.pdf")), p, width = 8, height = 10)
