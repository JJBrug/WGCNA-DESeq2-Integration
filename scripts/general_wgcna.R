# ==============================================================================
#   Generic WGCNA Pipeline for RNA-seq Data (Raw Counts -> Network -> Hubs)
# ==============================================================================
#
#   Usage:
#   1. Place your raw counts file (genes as rows, samples as cols) in the working dir.
#   2. Place your metadata/traits file (samples as rows, traits as cols) in the
#      working dir.
#   3. Edit the "CONFIGURATION" section below to match your file names and
#      preferences.
#   4. Run the script.
#
# ==============================================================================

# ==============================================================================
#   1. SETUP & CONFIGURATION
# ==============================================================================

# --- User Parameters (EDIT THESE) ---
input_counts_file <- "split_data/raw_counts_rice_Severe_Drought.tsv"    # Path to raw counts (CSV/TSV)
input_traits_file <- "split_data/metadata_rice_Severe_Drought.csv"      # Path to traits/metadata (CSV)
output_prefix     <- "rice_Severe_Drought"    # Prefix for all output files
output_dir        <- "general_results"

target_trait      <- "Time"  # column name in metadata to correlate with hubs

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Filtering Parameters
min_counts        <- 10       # Minimum counts per gene
min_samples       <- 5        # Minimum samples required to have 'min_counts'
variance_filter   <- 0.25     # Remove bottom % of low-variance genes (0.25 = remove bottom 25%)
gene_id_pattern   <- ""       # Optional regex to keep specific genes (e.g., "^GLYMA_"). Leave "" for all.

# Network Construction Parameters
net_type          <- "signed"   # 'signed' is recommended for biological data
cor_type          <- "bicor"    # 'bicor' is more robust to outliers than 'pearson'
block_size        <- 4000       # Max genes per block (lower this if you have low RAM)
soft_power_cut    <- 0.75       # Target R^2 for scale-free topology
min_module_size   <- 30         # Minimum number of genes to form a module
merge_cut_height  <- 0.25       # Merge modules with > 0.75 correlation (1 - 0.25)

# Parallel Processing
enable_threads    <- FALSE
n_threads         <- 4          # Adjust based on your CPU

# ------------------------------------

# Load required packages
required_packages <- c("WGCNA", "DESeq2", "tidyverse", "corrr", "gridExtra")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg == "DESeq2") {
      BiocManager::install("DESeq2")
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

# WGCNA Options
options(stringsAsFactors = FALSE)
if (enable_threads) enableWGCNAThreads(nThreads = n_threads)

# ==============================================================================
#   2. DATA INGESTION & NORMALIZATION
# ==============================================================================
# Purpose: WGCNA requires normalized data (that follows a normal distribution).
# Raw counts are negative binomial; we use DESeq2's Variance Stabilizing Transformation (VST)
# to make them homoscedastic (variance is independent of the mean).

# Load Counts
# Detect separator based on file extension
sep_counts <- if (grepl("\\.tsv$", input_counts_file)) "\t" else ","
counts_raw <- read.table(input_counts_file, header = TRUE, row.names = 1,
                         sep = sep_counts, check.names = FALSE)

# Load Traits
sep_traits <- if (grepl("\\.tsv$", input_traits_file)) "\t" else ","
traits_all <- read.table(input_traits_file, header = TRUE, row.names = 1,
                         sep = sep_traits, check.names = FALSE)

# Match Samples: Ensure samples in counts exist in traits
common_samples <- intersect(colnames(counts_raw), rownames(traits_all))

if (length(common_samples) < ncol(counts_raw)) {
  warning(paste("Dropped", ncol(counts_raw) - length(common_samples),
                "samples missing from metadata."))
}

counts_raw <- counts_raw[, common_samples]
traits_all <- traits_all[common_samples, ]

# Pre-filtering: Remove low count genes
keep_genes <- rowSums(counts_raw >= min_counts) >= min_samples
counts_filt <- counts_raw[keep_genes, ]

# Optional: Filter for specific gene ID patterns
if (gene_id_pattern != "") {
  keep_rows <- grepl(gene_id_pattern, rownames(counts_filt))
  counts_filt <- counts_filt[keep_rows, ]
}

# VST Normalization
dds <- DESeqDataSetFromMatrix(countData = counts_filt,
                              colData = traits_all,
                              design = ~ 1) # Blind design, we just want normalization
vst_obj <- vst(dds, blind = TRUE)
datExpr0 <- assay(vst_obj) # This is the normalized matrix (Genes x Samples)

# Transpose for WGCNA (Needs Samples x Genes)
datExpr0 <- t(datExpr0)

# Variance Filtering
# Purpose: Genes that don't vary across samples cannot be co-expressed with anything.
# We remove the genes with the lowest variance to reduce computational load and noise.
vars <- apply(datExpr0, 2, var)
cutoff_var <- quantile(vars, probs = variance_filter)
datExpr <- datExpr0[, vars > cutoff_var]


# ==============================================================================
#   3. SAMPLE OUTLIER DETECTION
# ==============================================================================
# Purpose: Outlier samples (e.g., failed library prep) can skew the entire network.
# We use clustering to visualize and potentially remove them.

# Outlier Removal Paramters
remove_outliers   <- TRUE   # Set TRUE to auto-remove bad samples
sample_cut_height <- 70      # Cut height for the sample dendrogram (check plot first!)
min_cluster_size  <- 10       # Minimum size to keep a cluster

# Check for missing values/zero variance before clustering
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  # Optionally auto-remove bad genes/samples here if any slipped through
  if (sum(!gsg$goodGenes) > 0) printFlush(paste("Removing genes:",
                                                paste(names(datExpr)[!gsg$goodGenes],
                                                      collapse = ", ")))
  if (sum(!gsg$goodSamples) > 0) printFlush(paste("Removing samples:",
                                                  paste(rownames(datExpr)[!gsg$goodSamples],
                                                        collapse = ", ")))
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  traits_all <- traits_all[gsg$goodSamples, ]
}

# Build the sample tree
sampleTree <- hclust(dist(datExpr), method = "average")

# Plot the tree
pdf(file.path(output_dir, paste0(output_prefix, "_SampleClustering.pdf")),
    width = 18, height = 9)
par(cex = 0.6); par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# If filtering is enabled, draw the cut line
if (remove_outliers) {
  abline(h = sample_cut_height, col = "red")
}
dev.off()

# Apply the Cut !!! Check PDF plot first !!!
if (remove_outliers) {
  # cutreeStatic determines which cluster each sample belongs to
  # Cluster 0 is unassigned (outliers), Cluster 1 is the main group
  clust <- cutreeStatic(sampleTree, cutHeight = sample_cut_height,
                        minSize = min_cluster_size)
  
  # Keep samples in Cluster 1 (the largest cluster)
  keepSamples <- (clust == 1)
  
  # Filter Expression Data
  n_genes_before <- nrow(datExpr)
  n_samples_before <- nrow(datExpr)
  datExpr <- datExpr[keepSamples, ]
  
  # Filter Traits Data to match
  traits_all <- traits_all[keepSamples, ]
  
  # Report results
  n_removed <- sum(!keepSamples)
  
  # Re-check genes after sample removal (sometimes removing samples makes genes constant)
  # This is a safe double-check
  gsg_post <- goodSamplesGenes(datExpr, verbose = 0)
  if (!gsg_post$allOK) {
    datExpr <- datExpr[gsg_post$goodSamples, gsg_post$goodGenes]
    traits_all <- traits_all[gsg_post$goodSamples, ]
  }
}


# ==============================================================================
#   4. SOFT THRESHOLDING (Scale-Free Topology)
# ==============================================================================
# Purpose: Biological networks are "scale-free" (few hubs, many non-hubs).
# WGCNA raises the correlation matrix to a power (beta) to enforce this topology.
# We check which power gives us a high Scale-Free Fit Index (R^2 > 0.8 or 0.9).

powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(
  datExpr, 
  powerVector = powers, 
  verbose = 5,
  networkType = net_type, 
  blockSize = 2000,
  corFnc = cor_type,                      # "bicor"
  corOptions = list(maxPOutliers = 0.05)  # Matches your network construction settings
)

# Plot on PDF for visualizing Scale Free Topology
pdf(file.path(output_dir, paste0(output_prefix, "_SoftThreshold.pdf")),
    width = 9, height = 5)
par(mfrow = c(1, 2))
# Scale-free topology fit index
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = 0.9, col = "red")
abline(h = soft_power_cut, col = "red")

# Mean connectivity
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers,
     cex = 0.9, col = "red")
dev.off()

# Automatically select the lowest power that passes the cutoff
softPower <- sft$fitIndices %>% 
  filter(SFT.R.sq >= soft_power_cut) %>% 
  pull(Power) %>% 
  min()

if(is.infinite(softPower)) {
  warning("No power reached the R^2 cutoff. Defaulting to power 6 (check your data!).")
  softPower <- 6
}
cat(paste("Selected Soft Threshold Power:", softPower, "\n"))


# ==============================================================================
#   5. NETWORK CONSTRUCTION (Blockwise Modules)
# ==============================================================================
# Purpose: This steps calculates the Adjacency Matrix -> Topological Overlap Matrix (TOM).
# It then clusters genes based on TOM dissimilarity to find "Modules" (co-expressed gene clusters).
# 'blockwiseModules' is used to handle large datasets by splitting them into manageable blocks.

net <- blockwiseModules(datExpr,
                        power = softPower,
                        maxBlockSize = block_size,
                        TOMType = net_type,
                        networkType = net_type,
                        minModuleSize = min_module_size,
                        reassignThreshold = 0,
                        mergeCutHeight = merge_cut_height,
                        numericLabels = TRUE,
                        deepSplit = 4,   # Int from 0 to 4 controlling sensitivity
                        pamRespectsDendro = FALSE,
                        saveTOMs = FALSE,
                        saveTOMFileBase = file.path(output_dir, output_prefix),
                        corType = cor_type,
                        maxPOutliers = 0.05, # Robustness against outliers
                                             # if bicor is used
                        verbose = 3)


# Convert labels to colors for plotting
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]

# Plot Dendrogram
n_blocks <- length(net$dendrograms)

pdf(file.path(output_dir, paste0(output_prefix, "_Dendrograms.pdf")), width = 12, height = 9)

for (b in 1:n_blocks) {
  plotDendroAndColors(net$dendrograms[[b]], 
                      moduleColors[net$blockGenes[[b]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = paste("Gene dendrogram and module colors - Block", b))
}

dev.off()


# ==============================================================================
#   6. MODULE-TRAIT CORRELATION
# ==============================================================================
# Purpose: Identify which modules are biologically relevant by correlating
# Module Eigengenes (MEs) with external traits (Yield, Disease resistance, etc.).
# The 'Eigengene' is the 1st Principal Component of the module (representative profile).

# Ensure traits are numeric
traits_numeric <- traits_all %>% select_if(is.numeric)

# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

# Correlate MEs with Traits
moduleTraitCor <- cor(MEs, traits_numeric, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Visualize Heatmap
pdf(file.path(output_dir, paste0(output_prefix, "_ModuleTraitHeatmap.pdf")),
    width = 10, height = 10)

# Format text for the heatmap (Correlation + p-value)
textMatrix <- paste(signif(moduleTraitCor, 2),
                    "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

par(mar = c(6, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traits_numeric),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = paste("Module-trait relationships"))
dev.off()


# ==============================================================================
#   7. GENE SIGNIFICANCE (GS) & MODULE MEMBERSHIP (kME)
# ==============================================================================
# Purpose:
# 1. Gene Significance (GS): Correlation of Gene Expression <-> Trait.
#    (Does this gene go up when Yield goes up?)
# 2. Module Membership (kME): Correlation of Gene Expression <-> Module Eigengene.
#    (Is this gene a "hub" or central player in the module?)
#
# Finding Hubs: We look for genes with HIGH GS (relevant to trait) and HIGH kME (central to module).

# Define thresholds. Criteria: |GS| > 0.5 and |kME| > 0.8 (Standard strict thresholds)
kME_cut <- 0.8
GS_cut  <- 0.5

# Calculate kME (Signed Module Membership)
datKME <- signedKME(datExpr, MEs, outputColumnName = "kME_")

# Select a Trait of Interest for Hub Analysis
selected_trait_name <- target_trait 

# Safety check
if (!selected_trait_name %in% colnames(traits_numeric)) {
  stop(paste("Trait", selected_trait_name, "not found in numeric metadata! Check spelling."))
}

trait_vect <- traits_numeric[[selected_trait_name]]

# Calculate GS for this trait
geneTraitSignificance <- as.data.frame(cor(datExpr, trait_vect, use = "p"))
names(geneTraitSignificance) <- paste0("GS.", selected_trait_name)

# Create a Summary Table of all Genes
geneInfo <- data.frame(GeneID = colnames(datExpr),
                       ModuleColor = moduleColors,
                       geneTraitSignificance,
                       datKME)

# Filter for Top Hub Candidates
hub_genes <- geneInfo %>%
  rowwise() %>%
  filter(
    abs(get(paste0("GS.", selected_trait_name))) > GS_cut,
    
    abs(get(paste0("kME_", ModuleColor))) > kME_cut
  ) %>%
  ungroup() %>% 
  arrange(desc(abs(get(paste0("GS.", selected_trait_name)))))

write.csv(hub_genes, 
          file = file.path(output_dir, paste0(output_prefix,
                                              "_HubGenes_Candidates.csv")), 
          row.names = FALSE)


# ==============================================================================
#   8. PLOT GS vs kME (The "Hub Plot")
# ==============================================================================
# Purpose: Visual confirmation. A good module for a trait will show a strong correlation
# between GS and kME. (Genes that define the module are also the ones affecting the trait).

modules_to_plot <- unique(hub_genes$ModuleColor) # Only plot modules with hubs

pdf(file.path(output_dir, paste0(output_prefix, "_GS_vs_kME.pdf")))

for (module in modules_to_plot) {
  # Match the module name to the column index in MEs
  # Assumes MEs names are like "MEturquoise" -> substring removes "ME"
  column <- match(module, substring(names(MEs), 3))
  moduleGenes <- (moduleColors == module)
  
  # Create the base scatterplot
  verboseScatterplot(abs(datKME[moduleGenes, paste0("kME_", module)]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership (kME) in", module, "module"),
                     ylab = paste("Gene Significance (GS) for", selected_trait_name),
                     main = paste("Membership vs. Significance\n", module, "module"),
                     cex.main = 1.2, cex.lab = 1.2, pch = 21, col = "black",
                     bg = module)
  
  # --- Add Cutoff Lines ---
  # Vertical line for Module Membership (kME)
  abline(v = kME_cut, col = "grey60", lty = 2, lwd = 1.5)
  
  # Horizontal line for Gene Significance (GS)
  abline(h = GS_cut, col = "grey60", lty = 2, lwd = 1.5)
  
  # Optional: Add a legend or text to identify the quadrant
  # legend("bottomright", legend = c("Hub Thresholds"), col = "red", lty = 2, bty = "n")
}
dev.off()

# ==============================================================================
#   9. SAVE DATA AS RDS
# ==============================================================================

# Save the critical objects individually
saveRDS(datExpr,      file = file.path(output_dir,
                                       paste0(output_prefix, "_datExpr.rds")))
saveRDS(traits_all,   file = file.path(output_dir,
                                       paste0(output_prefix, "_traits.rds")))
saveRDS(moduleColors, file = file.path(output_dir,
                                       paste0(output_prefix, "_moduleColors.rds")))
saveRDS(net,          file = file.path(output_dir,
                                       paste0(output_prefix, "_net.rds")))
saveRDS(softPower,    file = file.path(output_dir,
                                       paste0(output_prefix, "_softPower.rds")))
