################################################################################
# 0: BASIC QC + NORMALIZATION 
# Authors: Isabel Schuurmans, Matthew Suderman 
# Analyst: YOUR NAME
# Cohort: YOUR COHORT
################################################################################

# YOU DONT HAVE TO RUN THIS SCRIPT IF YOUR DATA IS ALREADY PREPROCESSED
# BUT WE DO RECOMMEND READING THROUGH THE SCRIPT TO CHECK FOR ANY MAJOR INCONSISTENCIES

# IMPOERANT: you should run this script for each batch separately
# --- Indicate Batch and Timepoint ---
batch <- "BATCH_NAME"        # e.g. "450K"
timepoint <- "TIMEPOINT_NAME" # e.g. "T1" or "age10"

# Create suffix for file names
suffix <- paste0("_", batch, "_", timepoint)

# --- Load packages ---
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("remotes"))
  install.packages("remotes") 
if (!require("meffil"))
  remotes::install_github("perishky/meffil")

library(meffil) ## https://github.com/perishky/meffil

# --- Main paths ---

# Define the base project folder
base_dir <- "/home/i.schuurmans/EpiAge_BrainAge/brain_epi/"

idats.dir <- "PATH WHERE IDATS ARE SAVED" 
rawdnam.dir <- file.path(base_dir, "epi/preprocessed")

# test plate - DELETE!
idats.dir <- '/home/i.schuurmans/DNAmChecks/April2025/Plate1_2025-08-14/EPIC2023-400_plate1_NEW/idats/'

# --- Step 1: Reload sample sheet ---
samplesheet = meffil.create.samplesheet(idats.dir, recursive=TRUE)
head(samplesheet)
# IF SEX IS NA, PLEASE ADD SEX TO THIS DATASET 

# SELECT THE RIGHT REFERENCE HERE BASED ON YOUR TISSUE!
# --- Step 2: Run QC ---
qc.objects = meffil.qc(samplesheet, verbose=T,
                       cell.type.reference=
                         #  "combined cord blood" # cord blood
                         #  "blood gse35069" # peripheral blood
                         #  "saliva gse48472" # saliva/buccal 
) 

qc.parameters <- meffil.qc.parameters(
  beadnum.samples.threshold = 0.1,
  detectionp.samples.threshold = 0.1,
  detectionp.cpgs.threshold = 0.1,
  beadnum.cpgs.threshold = 0.1,
  sex.outlier.sd = 5,
  snp.concordance.threshold = 0.95,
  sample.genotype.concordance.threshold = 0.8
)

qc.summary <- meffil.qc.summary(qc.objects, parameters = qc.parameters)

# --- Step 3: Save original QC objects ---
saveRDS(qc.objects, file.path(rawdnam.dir, paste0("qc.objects_new", suffix, ".rds")))
saveRDS(qc.summary, file.path(rawdnam.dir, paste0("qc.summary_new", suffix, ".rds")))
meffil.qc.report(qc.summary, file.path(rawdnam.dir, paste0("qc-report_new", suffix, ".html")))

# --- Step 4: Filter bad samples ---
qc.objects <- meffil.remove.samples(qc.objects, qc.summary$bad.samples$sample.name) 

# --- Step 5: Quantile normalization ---
## select number of control probe PCs for normalization 
## based on the following plot
pc.fit = meffil.plot.pc.fit(qc.objects, n.cross=3)
ggsave(pc.fit$plot, filename = file.path(rawdnam.dir, paste0("pc-fit", suffix, ".pdf")), height=6, width=6)
number.pcs = 10 # THIS IS THE DEFAULT FOR MIND

## normalize the data
norm.objects <- meffil.normalize.quantiles(
  qc.objects,
  number.pcs = number.pcs,
  random.effects = "Slide"
)

# if you have only one slide, use this code instead.
# unsure? if your sample sheet has 96 rows, likely you only have one slide
#norm.objects <- meffil.normalize.quantiles(
#  qc.objects,
#  number.pcs = number.pcs
#)

saveRDS(norm.objects, file.path(rawdnam.dir, paste0("norm.objects", suffix, ".rds")))

# ---- Step 6: Create beta matrix ----
meth <- meffil.normalize.samples(
  norm.objects,
  just.beta = TRUE,
  cpglist.remove = qc.summary$bad.cpgs$name,
  verbose = TRUE
)
saveRDS(meth, file.path(rawdnam.dir, paste0("meth_new", suffix, ".rds")))

## generate normalization report
pcs <- meffil.methylation.pcs(meth)
saveRDS(pcs, file.path(rawdnam.dir, paste0("pcs_all_new", suffix, ".rds")))

batch.variables <- c("Slide", "sentrix_row", "sentrix_col")
parameters <- meffil.normalization.parameters(norm.objects, variables = batch.variables)
summary <- meffil.normalization.summary(norm.objects, pcs, parameters)

# ---- Step 7: Normalization report  ---
meffil.normalization.report(
  summary,
  output.file = file.path(rawdnam.dir, paste0("norm-report_new", suffix, ".html")),
  author = "YOUR NAME",
  study = "YOUR COHORT"
)

