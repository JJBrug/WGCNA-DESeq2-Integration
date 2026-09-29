# ==============================================================================
#   WGCNA Module Preservation
# ==============================================================================
#   
#   Purpose: 
#   Validates modules from a "Reference" network against multiple "Test" datasets.
#
#   Usage:
#   1. Define your Reference files (Network + Expression).
#   2. Define your Test files (Expression only).
#   3. Run.
#
# ==============================================================================

library(WGCNA); library(tidyverse); library(ggplot2)
options(stringsAsFactors = FALSE); #enableWGCNAThreads(nThreads = 4)

# --- CONFIGURATION ------------------------------------------------------------
output_dir     <- "preservation_results"
output_prefix  <- "rice_Module_Preservation"
if (!dir.exists(output_dir)) dir.create(output_dir)

# 1. THE REFERENCE (Source of Truth / Hypothesis)
label_ref       <- "Reference_Set"
file_ref_net    <- "general_results/rice_Severe_Drought_net.rds"
file_ref_expr   <- "general_results/rice_Severe_Drought_datExpr.rds"

# 2. THE TEST DATASETS (Validation Sets)
#    Define as many or as few as you need in this list. 
#    Format: "Label" = "Path/to/datExpr.rds"
test_datasets <- list(
  Test_Set_1 = "general_results/rice_Control_datExpr.rds"#,
  #Test_Set_2 = "general_results/L20_datExpr.rds",
  #Test_Set_3 = "general_results/L21_datExpr.rds"
)

n_permutations  <- 200

# ==============================================================================
#   1. LOAD DATA & DEFINE COMMON GENES
# ==============================================================================

# Load Reference
if(!file.exists(file_ref_net)) stop("Reference Network file missing")
net_ref <- readRDS(file_ref_net)

if(!file.exists(file_ref_expr)) stop("Reference Expr file missing")
d_ref <- readRDS(file_ref_expr)

# Load Tests into a list
d_tests <- list()
for (t_label in names(test_datasets)) {
  if(!file.exists(test_datasets[[t_label]])) stop(paste("Test file missing:", test_datasets[[t_label]]))
  d_tests[[t_label]] <- readRDS(test_datasets[[t_label]])
}

# Intersect to find common genes across ALL sets
common <- unlist(colnames(d_ref))
for (t_label in names(d_tests)) {
  common <- intersect(common, colnames(d_tests[[t_label]]))
}

# Safety Check
if (length(common) < 500) {
  stop("Error: Too few common genes found. Check column name formats.")
}

# ==============================================================================
#   2. BUILD WGCNA INPUTS
# ==============================================================================

multiExpr  <- list()
multiColor <- list()

# Extract Colors from Reference Network
refColors_raw <- labels2colors(net_ref$colors)
names(refColors_raw) <- unlist(colnames(d_ref))
refColors_common <- refColors_raw[common]

# Helper to subset data and FORCE column names
format_data <- function(mat, genes) {
  sub_mat <- mat[, genes]
  colnames(sub_mat) <- genes 
  return(sub_mat)
}

# 1. Reference Assignment
multiExpr[[label_ref]]  <- list(data = format_data(d_ref, common))
multiColor[[label_ref]] <- refColors_common

# 2. Dynamic Test Set Assignments
for (t_label in names(d_tests)) {
  multiExpr[[t_label]] <- list(data = format_data(d_tests[[t_label]], common))
  multiColor[[t_label]] <- refColors_common
}

# ==============================================================================
#   3. RUN PRESERVATION STATS
# ==============================================================================

# referenceNetworks = 1 means "Use the first dataset in multiExpr as the Reference"
# (which corresponds to label_ref above)

mp <- modulePreservation(
  multiExpr, 
  multiColor,
  referenceNetworks = 1,
  nPermutations = n_permutations,
  randomSeed = 1,
  verbose = 3,
  networkType = "signed", 
  corFnc = "bicor",
  corOptions = "maxPOutliers = 0.05",
  savePermutedStatistics = FALSE
)

saveRDS(mp, file.path(output_dir, paste0(output_prefix, "_mpObject.rds")))


# ==============================================================================
#   4. EXTRACT & PLOT RESULTS
# ==============================================================================

stats_list <- list()
test_sets <- names(test_datasets)

# 1. Locate the Reference Slot
# WGCNA names the slot "ref.{NameOfReferenceSet}"
ref_slot_name <- paste0("ref.", label_ref)

if (!ref_slot_name %in% names(mp$preservation$Z)) {
  stop(paste("Error: Could not find reference slot:", ref_slot_name, 
             "\nAvailable slots:", paste(names(mp$preservation$Z), collapse=", ")))
}

ref_stats <- mp$preservation$Z[[ref_slot_name]]

# 2. Loop through Test Sets
for (test in test_sets) {
  
  # WGCNA names the comparison slot "inColumnsAlsoPresentIn.{NameOfTestSet}"
  target_slot <- paste0("inColumnsAlsoPresentIn.", test)
  
  # Extract Z-stats
  z_stats <- ref_stats[[target_slot]]
  
  if (!is.null(z_stats)) {
    df <- data.frame(
      Module     = rownames(z_stats),
      ModuleSize = z_stats$moduleSize,
      Zsummary   = z_stats$Zsummary.pres,
      TestSet    = test
    )
    # Remove noise modules
    stats_list[[test]] <- df %>% filter(!Module %in% c("grey", "gold"))
  } else {
    warning(paste("Stats missing for:", target_slot))
  }
}

# 3. Bind & Save
full_stats <- do.call(rbind, stats_list)

if (is.null(full_stats) || nrow(full_stats) == 0) {
  stop("Error: full_stats is empty. Check if test set names match the slots in 'mp'.")
}

write.csv(full_stats, file.path(output_dir, paste0(output_prefix, "_Zsummary.csv")), row.names=F)

# 4. Plot
unique_mods <- unique(full_stats$Module)
mod_colors_map <- setNames(unique_mods, unique_mods)

p1 <- ggplot(full_stats, aes(x = ModuleSize, y = Zsummary, color = Module, label = Module)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text(vjust = -0.5, size = 3, check_overlap = TRUE, colour = "black") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "darkgreen") +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  scale_color_manual(values = mod_colors_map) +
  scale_x_log10() +
  facet_wrap(~TestSet) +
  labs(
    title = paste("Preservation of", label_ref, "Modules"),
    subtitle = "Z > 10: Preserved | Z < 2: Non-Preserved",
    x = "Module Size", y = "Z-Summary"
  ) +
  theme_bw() + theme(legend.position = "none")

ggsave(file.path(output_dir, paste0(output_prefix, "_Zsummary.pdf")), p1, width=12, height=8)

# ==============================================================================
#   5. MEDIAN RANK PLOT (Robust to Module Size)
# ==============================================================================

rank_list <- list()
# Re-define standard slot name variables
ref_slot_name <- paste0("ref.", label_ref)

for (test in test_sets) {
  target_slot <- paste0("inColumnsAlsoPresentIn.", test)
  
  # Note: Median Rank is stored in the 'observed' slot, not 'Z'
  obs_stats <- mp$preservation$observed[[ref_slot_name]][[target_slot]]
  
  if (!is.null(obs_stats)) {
    df <- data.frame(
      Module     = rownames(obs_stats),
      MedianRank = obs_stats$medianRank.pres,
      TestSet    = test
    )
    # Remove noise
    rank_list[[test]] <- df %>% filter(!Module %in% c("grey", "gold"))
  }
}

full_rank <- do.call(rbind, rank_list)
write.csv(full_rank, file.path(output_dir, paste0(output_prefix, "_MedianRank.csv")), row.names=F)

# Plotting
# We reverse the Y-axis because Rank 1 is "Best" (Top of chart)
p2 <- ggplot(full_rank, aes(x = MedianRank, y = reorder(Module, -MedianRank), color = Module)) +
  geom_point(size = 4) +
  scale_color_manual(values = mod_colors_map) +
  facet_wrap(~TestSet) +
  labs(
    title = paste("Median Rank of", label_ref, "Modules"),
    subtitle = "Lower Rank (Left) = Better Preservation",
    x = "Median Rank",
    y = "Module"
  ) +
  theme_bw() + 
  theme(legend.position = "none")

ggsave(file.path(output_dir, paste0(output_prefix, "_MedianRank.pdf")), p2, width=10, height=8)
