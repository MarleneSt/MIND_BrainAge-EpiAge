################################################################################
# 0: BASIC QC + NORMALIZATION 
# Authors: Isabel Schuurmans, Matthew Suderman 
# Analyst: YOUR NAME
# Cohort: YOUR COHORT
################################################################################

# IMPORTANT: you should run this script for each batch separately
# --- Indicate Batch and Timepoint ---
batch <- "BATCH_NAME"        # e.g. "450K"
timepoint <- "TIMEPOINT_NAME" # e.g. "T1" or "age10"

# Create suffix for file names
suffix <- paste0("_", batch, "_", timepoint)

# --- Load packages ---
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("remotes", quietly = TRUE)) install.packages("remotes")
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!require("tidyr", quietly = TRUE))   install.packages("tidyr")
library(ggplot2); library(tidyr)

suppressPackageStartupMessages({
  if (!require("meffil", quietly = TRUE)) {
    remotes::install_github("perishky/meffil")
    library(meffil) ## https://github.com/perishky/meffil
  } else library(meffil)
  
  if (!require("EpiDISH", quietly = TRUE)) install.packages("EpiDISH")
  library(EpiDISH)
  
  if (!require("matrixStats", quietly = TRUE)) install.packages("matrixStats")
  library(matrixStats)
  
  if (!require("dplyr", quietly = TRUE)) install.packages("dplyr")
  library(dplyr)
})

#-------------------------------------------------------------------------------

# now, let's get the folder structure right

# Define the output folder where you want to store all folders for this study
################################################################################
# Create full directory structure for brain_epi project
# Author: Isabel Schuurmans
# Date: 2025-10-09
################################################################################

# Define the base project folder
base_dir <- "/home/i.schuurmans/EpiAge_BrainAge/brain_epi/"

#  test plate - DELETE!
rawdnam.dir  <- file.path(base_dir, "epi/preprocessed") # where the normalized data is
share.dir   <- file.path(base_dir, "epi/descriptives_and_plots") # where we want to save output that we are going to share
paths.script <- file.path(base_dir, "epi/scripts") # where saved the scripts
cellcount.dir <- file.path(base_dir, "epi/preprocessed") # where we want to save the cell types

# --- Functions ---
source(file.path(paths.script, "/functions/1.functions_2025-10-09.R"))

#-------------------------------------------------------------------------------
# IMPORTANT NOTE: you can skip step 1-2 if you have the following cell counts already available:
# "combined cord blood" (Salas/Gervin); "blood gse35069" (Reinius/Houseman);
# "saliva gse48472" (saliva/buccal)
#-------------------------------------------------------------------------------

# ---- Step 1: Read in the normalized object ---
# Read object from Script 0 with the SAME suffix

#####################################
# IMPORTANT PART BELOW! READ IF YOU DONT HAVE IDAT DATA
####################################

# The file 'norm.objects.rds' is created by meffil during preprocessing (Script 0). 
# It contains the quality-controlled and functionally normalized methylation data 
# for all samples, and is required for downstream cell count estimation. 
# If you don't have this file yet, please run Script 0 first. 
# If you don't want to run script 0 (not preferred) / don't have IDAT files available 
# than please skip step 1 and 2 but go directly to step 3, and unhash the cell count calculation with note 
# 'UNHASH IF STEP 1 COULD NOT BE DONE' 
# And finally, if the cell-count of interest is already available in your dataset 
# you could also just read in this available set. 

norm.objects <- readRDS(file.path(rawdnam.dir, paste0("norm.objects", suffix, ".rds")))

# ---- Step 2: Calculate the counts (based on meffil) ---
counts <- t(meffil.cell.count.estimates(norm.objects))  # samples x cell types
cell_counts_meffil <- as.data.frame(counts * 100)
cell_counts_meffil$ID <- rownames(cell_counts_meffil)

saveRDS(cell_counts_meffil, file.path(cellcount.dir, paste0("cell_counts_from_norm", suffix, ".rds")))
write.csv(cell_counts_meffil, file.path(cellcount.dir, paste0("cell_counts_from_norm", suffix, ".csv")), row.names = FALSE)

# ---- Step 3: Read in the EWAS object (beta matrix or list with $betas) ---
meth  <- readRDS(file.path(rawdnam.dir, paste0("meth_new", suffix, ".rds"))) 
# you can change this if you are reading in your own DNAm data directly. 
# Please use a functional normalized set

betas <- safe_as_matrix(meth)    # CpGs x Samples
# ensure samples are columns; EpiDISH expects CpGs x Samples
if (nrow(betas) < ncol(betas) && grepl("^cg", rownames(betas)[1])) {
  # likely already CpG x sample; keep
} else if (ncol(betas) < nrow(betas) && grepl("^cg", colnames(betas)[1])) {
  betas <- t(betas)
}

# UNHASH IF STEP 1 COULD NOT BE DONE
#counts <- meffil.estimate.cell.counts.from.betas(betas, cell.type.reference = "blood gse35069", verbose = TRUE) # "combined cord blood" (cord); "blood gse35069" (peripheral blood); "saliva gse48472" (saliva/buccal)
#cell_counts_meffil <- as.data.frame(counts * 100)
#cell_counts_meffil$ID <- rownames(cell_counts_meffil)

# ---- Step 4: EpiDISH counts ----
# 4a. UniLIFE for everyone (life-course reference)
data(centUniLIFE.m, package = "EpiDISH")
epid_unilife <- EpiDISH::epidish(beta.m = betas, ref.m = centUniLIFE.m, method = "RPC")
cell_counts_unilife <- as.data.frame(epid_unilife$estF * 100)
cell_counts_unilife$ID <- rownames(cell_counts_unilife)
saveRDS(cell_counts_unilife, file.path(cellcount.dir, paste0("cell_counts_unilife", suffix, ".rds")))
write.csv(cell_counts_unilife, file.path(cellcount.dir, paste0("cell_counts_unilife", suffix, ".csv")), row.names = FALSE)

# 4b–d. Tissue-specific deconvolution (only ONE will run)

# ==== ONE-TISSUE DATASET SWITCH ====
# Set this for the whole dataset:
dataset_tissue <- "peripheral_blood"  # or "cord_blood" or "saliva"
array <- 'EPIC' # or "450K"

cell_counts_epidish_cord    <- NULL
cell_counts_epidish_blood   <- NULL
cell_counts_hepidish_saliva <- NULL

.array <- toupper(trimws(array))
if (!array %in% c("EPIC", "450K")) {
  stop("[ERROR] 'array' must be 'EPIC' or '450K'. You provided: ", array)
}

if (dataset_tissue == "cord_blood") {
  # --- Cord blood ---
  message("[INFO] dataset_tissue='cord_blood' but no ref_cord provided; skipping cord epiDISH.")
  
} else if (dataset_tissue == "peripheral_blood") {
  # --- Peripheral blood ---
  ref_name <- if (.array == "EPIC") "cent12CT.m" else "cent12CT450k.m"
  
  # Load from env if present; otherwise try to load from EpiDISH package
  if (!exists(ref_name, inherits = TRUE)) {
    suppressWarnings({
      if (.array == "EPIC") {
        data(cent12CT.m, package = "EpiDISH", envir = environment())
      } else {
        data(cent12CT450k.m, package = "EpiDISH", envir = environment())
      }
    })
  }
  
  if (!exists(ref_name, inherits = TRUE)) {
    stop("[ERROR] Could not find reference '", ref_name,
         "'. Ensure EpiDISH is installed and the reference is available for ", .array, ".")
  }
  
  ref_blood <- get(ref_name, inherits = TRUE)
  message("[INFO] Using blood reference: ", ref_name, " (array=", .array, ").")
  
  epid_blood <- EpiDISH::epidish(beta.m = betas, ref.m = ref_blood, method = "RPC")
  cell_counts_epidish_blood <- as.data.frame(epid_blood$estF * 100)
  cell_counts_epidish_blood$ID <- rownames(cell_counts_epidish_blood)
  
  saveRDS(cell_counts_epidish_blood, file.path(cellcount.dir, paste0("cell_counts_epidish_blood", suffix, ".rds")))
  write.csv(cell_counts_epidish_blood, file.path(cellcount.dir, paste0("cell_counts_epidish_blood", suffix, ".csv")), row.names = FALSE)
  
  
} else if (dataset_tissue == "saliva") {
  
  # --- Saliva/Buccal (hepiDISH two-step) ---
  data(centEpiFibIC.m, package = "EpiDISH")
  data(centBloodSub.m,  package = "EpiDISH")
  ref2_name <- if (.array == "EPIC") "centBloodSub.m" else "centBloodSub450k.m"
  frac_saliva <- EpiDISH::hepidish(
    beta.m   = betas,
    ref1.m   = centEpiFibIC.m,  # step 1: epithelial, fibroblast, immune
    ref2.m   = centBloodSub.m,  # step 2: decompose immune
    h.CT.idx = 3,
    method   = "RPC"
  )
  cell_counts_hepidish_saliva <- as.data.frame(frac_saliva * 100)
  cell_counts_hepidish_saliva$ID <- rownames(cell_counts_hepidish_saliva)
  saveRDS(cell_counts_hepidish_saliva, file.path(cellcount.dir, paste0("cell_counts_hepidish_saliva", suffix, ".rds")))
  write.csv(cell_counts_hepidish_saliva, file.path(cellcount.dir, paste0("cell_counts_hepidish_saliva", suffix, ".csv")), row.names = FALSE)
}

# ---- Step 5: Descriptives of cell counts ---
table_with_descriptives_cellcounts_meffil   <- desc_counts(cell_counts_meffil)
table_with_descriptives_cellcounts_unilife  <- desc_counts(cell_counts_unilife)

table_with_descriptives_cellcounts_epidish_cord    <- if (!is.null(cell_counts_epidish_cord))    desc_counts(cell_counts_epidish_cord)    else NULL
table_with_descriptives_cellcounts_epidish_blood   <- if (!is.null(cell_counts_epidish_blood))   desc_counts(cell_counts_epidish_blood)   else NULL
table_with_descriptives_cellcounts_hepidish_saliva <- if (!is.null(cell_counts_hepidish_saliva)) desc_counts(cell_counts_hepidish_saliva) else NULL

# ---- Step 6: Correlations of cell counts ----
correlations_cellcounts_meffil_unilife <- compute_pairwise_corr(cell_counts_meffil, cell_counts_unilife, "meffil", "unilife")
correlations_cellcounts_unilife_epidish_blood  <- if (!is.null(cell_counts_epidish_blood))   compute_pairwise_corr(cell_counts_unilife, cell_counts_epidish_blood, "unilife", "epidish_blood") else NULL
correlations_cellcounts_meffil_epidish_blood   <- if (!is.null(cell_counts_epidish_blood))   compute_pairwise_corr(cell_counts_meffil, cell_counts_epidish_blood,  "meffil",  "epidish_blood") else NULL
correlations_cellcounts_unilife_hepidish_saliva <- if (!is.null(cell_counts_hepidish_saliva)) compute_pairwise_corr(cell_counts_unilife, cell_counts_hepidish_saliva, "unilife", "hepidish_saliva") else NULL
correlations_cellcounts_meffil_hepidish_saliva  <- if (!is.null(cell_counts_hepidish_saliva)) compute_pairwise_corr(cell_counts_meffil,  cell_counts_hepidish_saliva,  "meffil",  "hepidish_saliva") else NULL

# ---- Step 7: plots cell counts ----
plot_cellcount_boxplot(cell_counts_meffil,          share.dir, paste0("meffil", suffix))
plot_cellcount_boxplot(cell_counts_unilife,         share.dir, paste0("unilife", suffix))
plot_cellcount_boxplot(cell_counts_epidish_blood,   share.dir, paste0("epidish_blood", suffix))
plot_cellcount_boxplot(cell_counts_epidish_cord,    share.dir, paste0("epidish_cord", suffix))
plot_cellcount_boxplot(cell_counts_hepidish_saliva, share.dir, paste0("hepidish_saliva", suffix))

# ---- Step 8: Descriptives of beta values (CpGs) ----
table_with_descriptives_cpgs <- descriptives(betas)

# ---- Step 9: Bundle + save (to SHARE folder) ----
to_share <- Filter(Negate(is.null), mget(c(
  "table_with_descriptives_cellcounts_meffil",
  "table_with_descriptives_cellcounts_unilife",
  "table_with_descriptives_cellcounts_epidish_cord",
  "table_with_descriptives_cellcounts_epidish_blood",
  "table_with_descriptives_cellcounts_hepidish_saliva",
  "correlations_cellcounts_meffil_unilife",
  "correlations_cellcounts_unilife_epidish_blood",
  "correlations_cellcounts_meffil_epidish_blood",
  "correlations_cellcounts_unilife_hepidish_saliva",
  "correlations_cellcounts_meffil_hepidish_saliva",
  "table_with_descriptives_cpgs"
), ifnotfound = list(NULL)))

# Save bundle
save(to_share, file = file.path(share.dir, paste0("descriptives", suffix, ".RData")))
save(to_share, file = file.path(cellcount.dir, paste0("descriptives", suffix, ".RData")))
write.csv(table_with_descriptives_cpgs, file.path(share.dir, paste0("descriptives_cpgs", suffix, ".csv")), row.names = FALSE)
