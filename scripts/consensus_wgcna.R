# ==============================================================================
#   Generic Consensus WGCNA Pipeline
#   (Multi-Set Analysis for Consistent Module Detection)
# ==============================================================================
#
#   Purpose: 
#   Construct a consensus gene co-expression network across multiple datasets 
#   (e.g., different years, tissues, or treatments) to find conserved biological modules.
#
#   Usage:
#   1. Edit the "CONFIGURATION" section.
#   2. Choose your input mode:
#      - LOAD_FROM_RAW = TRUE: Load one big count matrix and split it by a metadata column.
#      - LOAD_FROM_RAW = FALSE: Load separate .rds files (e.g. from the previous single-network script).
#
# ==============================================================================

# ==============================================================================
#   1. SETUP & CONFIGURATION
# ==============================================================================

# --- INPUT MODE SWITCH ---
# TRUE  = Load one large raw count matrix and split it (Data Ingestion -> VST -> Split).
# FALSE = Load separate, pre-normalized datExpr files (Data Ingestion -> Merge).
LOAD_FROM_RAW <- FALSE

# --- CONFIGURATION: SCENARIO A (LOAD_FROM_RAW = TRUE) ---
input_counts_file <- "raw_counts.tsv"      # Genes (rows) x Samples (cols)
input_traits_file <- "metadata.csv"        # Samples (rows) x Traits (cols)
split_column      <- "year"                # Column in metadata to split datasets by
gene_id_pattern   <- ""                    # Optional regex to keep specific genes (e.g., "^GLYMA_"). Leave "" for all.
sample_cut_height <- 50000                 # Cut height for outlier removal (adjust based on data)

# --- CONFIGURATION: SCENARIO B (LOAD_FROM_RAW = FALSE) ---
# List paths to your normalized expression matrices (datExpr) and traits.
# Names (L18, L19...) will be used as dataset label examples.
preprocessed_files <- list(
  L18 = list(expr = "general_results/L18_datExpr.rds", traits = "general_results/L18_traits.rds"),
  L19 = list(expr = "general_results/L19_datExpr.rds", traits = "general_results/L19_traits.rds"),
  L20 = list(expr = "general_results/L20_datExpr.rds", traits = "general_results/L20_traits.rds"),
  L21 = list(expr = "general_results/L21_datExpr.rds", traits = "general_results/L21_traits.rds")
)

# --- GENERAL PARAMETERS ---
output_dir        <- "consensus_results_"
output_prefix     <- "Soybean_Consensus_"
net_type          <- "signed"
cor_type          <- "bicor"
max_block_size    <- 4000
soft_power_cut    <- 0.7
consensus_split   <- 2
min_module_size   <- 30
merge_cut_height  <- 0.25
variance_filter   <- 0.25 # Remove bottom 25% of low-variance genes

# Parallel Processing
enable_threads    <- FALSE
n_threads         <- 4

# ------------------------------------
# Load Packages
packages <- c("WGCNA", "DESeq2", "tidyverse", "gridExtra")
for (p in packages) {
  if (!require(p, character.only = T)) {
    if (p == "DESeq2") BiocManager::install(p) else install.packages(p)
    library(p, character.only = T)
  }
}

options(stringsAsFactors = FALSE)

if (enable_threads) enableWGCNAThreads(nThreads = n_threads)

if (!dir.exists(output_dir)) dir.create(output_dir)


# ==============================================================================
#   2. DATA INGESTION
# ==============================================================================

multiExpr   <- list() # Container for expression data
multiTraits <- list() # Container for trait data

if (LOAD_FROM_RAW) {
  
  # 1. Load Data
  counts_raw <- read.table(input_counts_file, header = T, row.names = 1, sep = "\t", check.names = F)
  traits_all <- read.csv(input_traits_file, row.names = 1, check.names = F)
  
  # 2. Filter Gene IDs (Optional)
  if (gene_id_pattern != "") {
    counts_raw <- counts_raw[grepl(gene_id_pattern, rownames(counts_raw)), ]
  }
  
  # 3. Match Samples
  common <- intersect(colnames(counts_raw), rownames(traits_all))
  counts_raw <- counts_raw[, common]
  traits_all <- traits_all[common, ]
  
  # 4. Filter Low Counts
  keep <- rowSums(counts_raw >= 10) >= 5
  counts_filt <- counts_raw[keep, ]
  
  # 5. VST Normalization
  dds <- DESeqDataSetFromMatrix(countData = counts_filt, colData = traits_all, design = ~1)
  vst_mat <- assay(vst(dds, blind = TRUE)) # Genes x Samples
  
  # 6. Calculate variance across ALL samples before splitting
  vars <- apply(vst_mat, 1, var)
  cutoff_var <- quantile(vars, probs = variance_filter)
  vst_mat <- vst_mat[vars > cutoff_var, ]
  
  # 7. Split Data by Metadata Column
  if (!split_column %in% colnames(traits_all)) stop("Split column not found in metadata!")
  
  sets <- unique(traits_all[[split_column]])
  
  for (set in sets) {
    # Identify samples for this set
    samples_in_set <- rownames(traits_all)[traits_all[[split_column]] == set]
    
    # Subset Data
    expr_subset <- t(vst_mat[, samples_in_set]) # Transpose to Samples x Genes
    traits_subset <- traits_all[samples_in_set, ]
    
    # Basic Outlier Removal (Clustering)
    tree <- hclust(dist(expr_subset), method = "average")
    clust <- cutreeStatic(tree, cutHeight = sample_cut_height, minSize = 10) 
    keep_samp <- (clust == 1)
    
    # Store in multi-list structure required by WGCNA
    multiExpr[[as.character(set)]] <- list(data = expr_subset[keep_samp, ])
    multiTraits[[as.character(set)]] <- list(data = traits_subset[keep_samp, ])
  }
  
} else {
  
  sets <- names(preprocessed_files)
  
  for (set in sets) {
    if (!file.exists(preprocessed_files[[set]]$expr))
      stop(paste("File missing:", preprocessed_files[[set]]$expr))
    
    # Load and force to matrix
    expr_data <- readRDS(preprocessed_files[[set]]$expr)
    trait_data <- readRDS(preprocessed_files[[set]]$traits)
    
    multiExpr[[set]] <- list(data = as.matrix(expr_data))
    multiTraits[[set]] <- list(data = as.data.frame(trait_data))
  }
  
  # CRITICAL: Ensure all sets have the exact same genes
  genes_list <- lapply(multiExpr, function(x) colnames(x$data))
  common_genes <- Reduce(intersect, genes_list)
  
  # Subset all matrices to common genes
  for (set in names(multiExpr)) {
    multiExpr[[set]]$data <- multiExpr[[set]]$data[, common_genes]
  }
}

# Final Structure Check
checkSets(multiExpr)


# ==============================================================================
#   3. SOFT THRESHOLDING (Consensus)
# ==============================================================================

# We pick a power that works well for ALL sets.
# Strategy: Calculate scale-free fit for each set, then pick the lowest power
# where ALL sets meet the threshold (e.g. R^2 > 0.85).

powers <- c(seq(4,10,by=1), seq(12,20,by=2))
powerTables <- vector(mode = "list", length = length(multiExpr))

for (i in 1:length(multiExpr)) {
  sft <- pickSoftThreshold(
    multiExpr[[i]]$data, 
    powerVector = powers, 
    networkType = net_type,
    corFnc = cor_type,
    corOptions = list(maxPOutliers = 0.05),
    blockSize = 2000,
    verbose = 5
  )
  powerTables[[i]] <- sft$fitIndices
}

# Visualization
pdf(file.path(output_dir, paste0(output_prefix, "_Consensus_SoftPower.pdf")),
    height=6, width=10)
par(mfrow = c(1, 2))

# Plot R^2
colors <- 1:length(multiExpr)
plot(powers, powerTables[[1]]$SFT.R.sq, type="n", ylim=c(0,1),
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit (R^2)",
     main="Scale independence")
abline(h=soft_power_cut, col="red")

for (i in 1:length(multiExpr)) {
  points(powers, powerTables[[i]]$SFT.R.sq, col=colors[i], type="b", pch=19)
}
legend("bottomright", legend=names(multiExpr), col=colors, pch=19, cex=0.8)

# Plot Connectivity
plot(powers, powerTables[[1]]$mean.k., type="n",
     ylim=c(0,max(powerTables[[1]]$mean.k.)),
     xlab="Soft Threshold (power)", ylab="Mean Connectivity",
     main="Mean Connectivity")
for (i in 1:length(multiExpr)) {
  points(powers, powerTables[[i]]$mean.k., col=colors[i], type="b", pch=19)
}
dev.off()

# Automatic Selection (Lowest power where ALL sets > cutoff)
valid_powers <- powers
for (i in 1:length(multiExpr)) {
  set_valid <- powerTables[[i]]$Power[powerTables[[i]]$SFT.R.sq > soft_power_cut]
  valid_powers <- intersect(valid_powers, set_valid)
}

if (length(valid_powers) > 0) {
  softPower <- min(valid_powers)
} else {
  warning("No power reached consensus cutoff. Defaulting to 12. Check plots!")
  softPower <- 12
}

# ==============================================================================
#   4. CONSENSUS NETWORK CONSTRUCTION
# ==============================================================================
# blockwiseConsensusModules constructs the network for all sets simultaneously.

netConsensus <- blockwiseConsensusModules(
  multiExpr, 
  power = softPower, 
  minModuleSize = min_module_size, 
  deepSplit = consensus_split, 
  pamRespectsDendro = FALSE, 
  mergeCutHeight = merge_cut_height, 
  numericLabels = TRUE,
  saveConsensusTOMs = FALSE,
  saveIndividualTOMs = FALSE,
  verbose = 5,
  networkType = net_type,
  TOMType = net_type,
  corType = cor_type,
  maxPOutliers = 0.05,
  maxBlockSize = max_block_size
)

# Extract Results
consMEs <- netConsensus$multiMEs
moduleLabels <- netConsensus$colors
moduleColors <- labels2colors(moduleLabels)
nSets <- length(multiExpr)

# Save Network Object
saveRDS(netConsensus, file.path(output_dir, paste0(output_prefix, "_Network.rds")))

# Plot Consensus Dendrogram
pdf(file.path(output_dir, paste0(output_prefix, "_Dendrogram.pdf")),
    width=12, height=9)
plotDendroAndColors(netConsensus$dendrograms[[1]], 
                    moduleColors[netConsensus$blockGenes[[1]]],
                    "Module Colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus Gene Dendrogram")
dev.off()


# ==============================================================================
#   5. EIGENGENE NETWORK ANALYSIS
# ==============================================================================
# This compares the relationship between modules ACROSS the different sets.
# Are the modules correlated with each other in the same way in L18 vs L19?

pdf(file.path(output_dir, paste0(output_prefix, "_EigengeneNetworks.pdf")),
    width=30, height=30)

# 1. Preservation Plot (Multi-set)
# We set plotHeatmaps = TRUE because the preservation grid IS a heatmap.
plotEigengeneNetworks(
  consMEs, 
  setLabels = names(multiExpr), 
  plotDendrograms = FALSE, 
  plotHeatmaps = TRUE, 
  plotPreservation = "standard", 
  marHeatmap = c(5, 5, 5, 5)
)

# 2. Individual Networks (Loop)
# We use $data to pass a data frame, not a list.
for (i in 1:nSets) {
  plotEigengeneNetworks(
    consMEs[[i]]$data, 
    setLabels = names(multiExpr)[i],
    plotDendrograms = TRUE, 
    plotHeatmaps = TRUE, 
    signed = TRUE,
    marHeatmap = c(10, 12, 5, 5),
    marDendro = c(0, 4, 4, 0)
  )
}

dev.off()


# ==============================================================================
#   6. MODULE-TRAIT RELATIONSHIPS (Consensus)
# ==============================================================================
# We calculate the correlation of Consensus MEs with Traits for EACH set.
# Then we display a "Consensus" heatmap showing relationships that are consistent.

# Initialize storage
moduleTraitCor <- list()
moduleTraitPvalue <- list()

for (i in 1:nSets) {
  # Align traits to MEs
  curr_traits <- multiTraits[[i]]$data
  curr_traits_num <- curr_traits %>% select_if(is.numeric)
  
  # Correlate
  res <- corAndPvalue(consMEs[[i]]$data, curr_traits_num, use = "p")
  moduleTraitCor[[i]] <- res$cor
  moduleTraitPvalue[[i]] <- res$p
}

# --- CONSENSUS CALCULATION ---
# A relationship is "Consensus" if it has the SAME SIGN in all sets.
# If signs match, we take the minimum correlation magnitude.
# If signs differ, the consensus correlation is 0.

nMods <- ncol(consMEs[[1]]$data)
nTraits <- ncol(moduleTraitCor[[1]])
consensusCor <- matrix(0, nrow = nMods, ncol = nTraits)
consensusPvalue <- matrix(1, nrow = nMods, ncol = nTraits)

for (m in 1:nMods) {
  for (t in 1:nTraits) {
    # Extract corrs for this module-trait pair across all sets
    corrs <- sapply(moduleTraitCor, function(x) x[m, t])
    pvals <- sapply(moduleTraitPvalue, function(x) x[m, t])
    
    # Check for NA first. If any set is NA, we cannot claim consensus.
    if (any(is.na(corrs))) {
      consensusCor[m, t] <- 0
      consensusPvalue[m, t] <- 1  # Not significant
    } 
    # Check sign consistency (only if no NAs)
    else if (all(corrs > 0) | all(corrs < 0)) {
      # Keep the weakest correlation (conservative estimate)
      consensusCor[m, t] <- min(abs(corrs)) * sign(corrs[1])
      # Take the least significant p-value (conservative)
      consensusPvalue[m, t] <- max(pvals) 
    } else {
      # Signs disagree
      consensusCor[m, t] <- 0
    }
  }
}

rownames(consensusCor) <- substring(names(consMEs[[1]]$data), 3) # Remove "ME"
colnames(consensusCor) <- colnames(moduleTraitCor[[1]])

# Plot Heatmap
pdf(file.path(output_dir, paste0(output_prefix, "_ModuleTraitHeatmap.pdf")),
    width=10, height=20)

textMatrix <- paste(signif(consensusCor, 2),
                    "\n(", signif(consensusPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(consensusCor)

labeledHeatmap(Matrix = consensusCor,
               xLabels = colnames(consensusCor),
               yLabels = paste("ME", rownames(consensusCor), sep=""),
               ySymbols = paste("ME", rownames(consensusCor), sep=""),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Consensus Module-Trait Relationships\n(Consistent correlations across",
                            nSets, "sets)"))
dev.off()

# ==============================================================================
#   7. CONSENSUS HUB GENE EXTRACTION
# ==============================================================================
# Purpose: Extract genes that act as stable hubs ACROSS ALL environments.
# A true consensus hub must have high Module Membership (kME) AND high 
# Gene Significance (GS) consistently in every dataset, with no sign-flipping.

# Define Hub Thresholds
# These will depend on preferred strictness and the data used
target_trait <- "ayield" # Replace with exact metadata column name
kME_cut      <- 0.7
GS_cut       <- 0.3

# Verify trait exists in the first set (assuming uniform metadata)
if (!target_trait %in% colnames(multiTraits[[1]]$data)) {
  stop(paste("Trait", target_trait, "not found in metadata! Check spelling."))
}

nGenes <- ncol(multiExpr[[1]]$data)
gene_names <- colnames(multiExpr[[1]]$data)

# Initialize storage for kME and GS across all sets
kME_list <- list()
GS_list  <- list()

for (i in 1:nSets) {
  # Calculate kME for set i
  kME_list[[i]] <- signedKME(multiExpr[[i]]$data, consMEs[[i]]$data, outputColumnName = "kME_")
  
  # Calculate GS for set i
  trait_vect <- multiTraits[[i]]$data[[target_trait]]
  gs_df <- as.data.frame(cor(multiExpr[[i]]$data, trait_vect, use = "p"))
  colnames(gs_df) <- "GS"
  GS_list[[i]]  <- gs_df
}

# --- CONSENSUS CORRELATION FUNCTION ---
# Logic: If signs match across all sets, keep the minimum magnitude (conservative).
# If signs flip (e.g., positive in L18, negative in L19), the consensus is 0.
get_consensus_cor <- function(cor_list) {
  # Bind matrices into a 3D array
  arr <- array(unlist(cor_list), 
               dim = c(nrow(cor_list[[1]]), ncol(cor_list[[1]]), length(cor_list)))
  
  # Find max and min across the sets (3rd dimension)
  max_val <- apply(arr, c(1,2), max, na.rm = TRUE)
  min_val <- apply(arr, c(1,2), min, na.rm = TRUE)
  
  # Evaluate sign consistency
  consensus <- ifelse(max_val < 0, max_val,            # All negative: keep the one closest to 0 (max)
                      ifelse(min_val > 0, min_val, 0)) # All positive: keep the one closest to 0 (min)
  return(consensus)
}

# Apply consensus logic
consensus_kME <- get_consensus_cor(kME_list)
colnames(consensus_kME) <- colnames(kME_list[[1]])
rownames(consensus_kME) <- gene_names

consensus_GS <- get_consensus_cor(GS_list)
colnames(consensus_GS) <- paste0("Consensus_GS.", target_trait)
rownames(consensus_GS) <- gene_names

# Compile Summary Table
geneInfoConsensus <- data.frame(
  GeneID = gene_names,
  ModuleColor = moduleColors,
  ModuleLabel = moduleLabels, # Use the numeric labels from WGCNA
  Consensus_GS = consensus_GS[, 1]
)
geneInfoConsensus <- cbind(geneInfoConsensus, consensus_kME)

# Safely extract the exact kME value matching each gene's numeric module label
kme_col_names <- paste0("kME_", geneInfoConsensus$ModuleLabel)
col_idx <- match(kme_col_names, colnames(geneInfoConsensus))

# Use matrix indexing to pull the diagonal values and force them to numeric
geneInfoConsensus$Module_kME <-
  as.numeric(geneInfoConsensus[cbind(1:nrow(geneInfoConsensus), col_idx)])

# Filter for Consensus Hub Candidates
consensus_hubs <- geneInfoConsensus %>%
  filter(
    abs(Consensus_GS) > GS_cut,
    abs(Module_kME) > kME_cut
  ) %>%
  arrange(desc(abs(Consensus_GS)))

# Save to CSV
write.csv(consensus_hubs, 
          file = file.path(output_dir,
                           paste0(output_prefix, "_Consensus_HubGenes.csv")), 
          row.names = FALSE)

cat(paste("\nFound", nrow(consensus_hubs),
          "strict consensus hub genes for trait:", target_trait, "\n"))
# ==============================================================================
#   8. EXPORT DATA
# ==============================================================================

# Save useful variables for downstream plotting or inspection
saveRDS(multiExpr, file.path(output_dir,
                             paste0(output_prefix, "_multiExpr.rds")))
saveRDS(multiTraits, file.path(output_dir,
                               paste0(output_prefix, "_multiTraits.rds")))
saveRDS(consMEs, file.path(output_dir,
                           paste0(output_prefix, "_ConsensusMEs.rds")))
saveRDS(moduleColors, file.path(output_dir,
                                paste0(output_prefix, "_ModuleColors.rds")))
