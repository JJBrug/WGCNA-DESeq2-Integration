# WGCNA-DESeq2-Integration
An R-based framework unifying DESeq2 generalized linear modeling with WGCNA network topology.

## Overview
Standard bioinformatics pipelines often treat secondary biological variables as nuisance batch effects or execute differential expression and network topology as entirely disjointed operations. This repository provides customized accessible R scripts that bridge this interpretive gap. It allows researchers to explicitly model complex variance using likelihood ratio tests (LRT) and systematically map those dynamic targets to co-expression modules via hypergeometric integration.

## Getting Started
* **To run this integration pipeline:** Navigate to the `/scripts` folder. Ensure your input files (`raw_counts.tsv` and `metadata.csv`) are placed in your working directory and configured in the top 20 lines of the scripts.
* **To view results:** Check the generated `general_results`, `consensus_results`, or `preservation_results` subfolders for PDF heatmaps, dendrograms, and the final `HubGenes_Candidates.csv` files.

## Pipeline Execution
WGCNA is highly memory-intensive. These scripts use the `max_block_size` parameter to manage RAM constraints automatically. For standard laptops (8GB-16GB RAM), leave `max_block_size <- 4000`. For 32GB+ workstations, set it as high as possible to force WGCNA to process all genes in a single unified block.

### 1. `split_raw.R` (Data Partitioning)
Splits a master counts matrix into independent datasets based on a specified metadata column (e.g., Year, Tissue, Treatment). Generates the formatted inputs required for subsequent steps.

### 2. `[species]_de_lrt.R` (Differential Expression Modeling)
Leverages DESeq2 to execute generalized linear modeling (e.g., additive models or Interaction Likelihood Ratio Tests) to isolate genes whose trajectories fundamentally deviate due to applied stress or environmental variance.

### 3. `general_wgcna.R` (Single Network Construction)
Constructs an independent co-expression network for a single dataset. 
* Performs automated Variance Stabilizing Transformation (VST) via DESeq2. 
* Removes outlier samples via hierarchical clustering. 
* Calculates module eigengenes, correlates them with a user-defined continuous trait, and extracts module-specific hub genes.

### 4. `consensus_wgcna.R` (Consensus Network Construction)
Builds a consensus network across multiple independent datasets simultaneously to isolate immutable genetic signatures across noisy environments. 
* Applies a strict sign-matching penalty to trait correlations. 
* Extracts "Consensus Hub Genes" that maintain strict Module Membership (kME) and Gene Significance (GS) thresholds across every tested dataset.

### 5. `preservation_wgcna.R` (Network Preservation Analysis)
Tests the robustness of a "Reference" network against dynamic "Test" datasets using permutation testing to map how systems collapse and rebuild. 
* Outputs a Z-summary plot to visualize which biological networks are conserved (Z > 10) and which are *de novo* or collapsed (Z < 2). 

### 6. `de_wgcna_integration.R` (Hypergeometric Integration)
The capstone integration script. Executes a mathematically strict hypergeometric test to determine if dynamically regulated genes isolated by the LRT are overrepresented in specific topological modules, adjusted via FDR.

## Data Availability
The full transcriptomic datasets used in the original publication case studies (*Glycine max* multi-environment field trials (https://doi.org/10.1186/s12870-026-08603-w) and *Oryza sativa* drought time-courses (https://doi.org/10.1186/s12870-025-07175-5)) are available via their respective DOIs.
