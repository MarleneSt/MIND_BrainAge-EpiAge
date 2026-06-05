# Date: 2025/03/12
# Author: VB

# This script obtains epi clock predictions

# _______________________________________________#
#                       SET UP
# _______________________________________________#

# load packages 
library(foreign)
library(tidyverse)
# install.packages("writexl") 
library(writexl)
library(openxlsx)

# ! set working directory to the "brain_epi/epi" folder we provided 
wd = "/campaign/VB-FM5HPC-001/Vilte/Projects/Marlene-epi-clocks/brain_epi/epi/" # (! update with your path to brain_epi/epi/ folder)
setwd(wd)

# set paths 
postprocess_path = paste0(wd, "postprocessed/imputed/") # path to 'postprocessed/imputed' folder (where you stored imputed and winsorized DNAm .rds files)
# ! make sure the "postprocessed/imputed/" folder also contains the corresponding samplesheet for your methylation data
output = paste0(wd, "models/") # specify path to your output folder (to store the epi predictions) 
functions_path = paste0(wd, "scripts/functions/") # specify path to "functions" folder provided by us

# if 'output' path does not exist, create it
if (!dir.exists(output)) { dir.create(output, recursive = TRUE) }

# _______________________________________________#
#             Install methylclock
# _______________________________________________#
# trying on 4.3
# module load CMake/3.26.3-GCCcore-12.3.0

# install.packages("nloptr", type = "source")
# install.packages("quantreg", dependencies = TRUE, configure.args = "--disable-gsl")
# install.packages("https://cran.r-project.org/src/contrib/Archive/MatrixModels/MatrixModels_0.4-1.tar.gz", repos = NULL, type = "source")
# 
# devtools::install_github("isglobal-brge/methylclock@main", lib = "/campaign/VB-FM5HPC-001/Vilte/R/x86_64-pc-linux-gnu-library/4.3")

library(methylclock)

# add the name of your imputed DNAm rds file and specify your timepoint 
# (for ALSPAC we have 3 time points, one of which is analysed twice - once for 450k and once for epic array)
# if you only have one time point, remove or comment out the irrelevant lines in this script 
# IMPORTANT: if you are using the birth time point, please set timepoint="birth" below. Other time points can be sample-specific (i.e., anything you like) 
source(paste0(functions_path, "processDNAmAge.R"))
results_450k_cord <- processDNAmAge(DNAm_imputed_rds="DNAm.imputed_wins_cord_450k.rds", timepoint="birth", output_dir=postprocess_path) 
results_450k_F7 <- processDNAmAge(DNAm_imputed_rds="DNAm.imputed_wins_F7_450k.rds", timepoint="F7", output_dir=postprocess_path)
results_450k_15up <- processDNAmAge(DNAm_imputed_rds="DNAm.imputed_wins_up15_450k.rds", timepoint="15up", output_dir=postprocess_path)
results_epic_15up <- processDNAmAge(DNAm_imputed_rds="DNAm.imputed_wins_up15_EPIC.rds", timepoint="15up", output_dir=postprocess_path)

# save each file
# !! update with the appropriate file names for your timepoint and array 
source(paste0(functions_path, "save_results_to_excel.R"))
save_results_to_excel(results_450k_cord, paste0(output, "epi_ages_450k_cord_pred.xlsx"))
save_results_to_excel(results_450k_F7, paste0(output, "epi_ages_450k_F7_pred.xlsx"))
save_results_to_excel(results_450k_15up, paste0(output, "epi_ages_450k_15up_pred.xlsx"))
save_results_to_excel(results_epic_15up, paste0(output, "epi_ages_epic_15up_pred.xlsx"))


# _______________________________________________#
#             Install dnaMethyAge 
# _______________________________________________#

## Make sure 'devtools' is installed in your R
# install.packages("devtools")
#devtools::install_github("yiluyucheng/dnaMethyAge")
library(dnaMethyAge)


# !! IMPORTANT: dnaMethyAge needs age information in addition to methylation data,
# please replace samplesheet[].csv with your covariate file that contains a 'Sample_Name' (i.e. sample ID), 'age', and 'sex' columns
# the file must be a .csv file 
# 'age' will be in years, except for birth, where age will be in gestational weeks
# 'sex' will be coded as "Female" (or anything that starts with an "F") and "Male" (or anything that starts with an "M")
source(paste0(functions_path, "compute_epigenetic_ages.R"))

# cord 450k 
epi_ages_450k_cord_dnaMethyAge <- compute_epigenetic_ages(
  samplesheet_path = paste0(postprocess_path, "samplesheet_cord_450k.csv"),
  methylation_rds_path = paste0(postprocess_path, "DNAm.imputed_wins_cord_450k.rds"),
  output_path = output,
  clocks = c("PCGrimAge", "ShirebyG2020", "BernabeuE2023c"),
  edited_methyAge_path = paste0(functions_path, "methyAge_edited.R"),
  array = "450k"
)

# F7 450k
epi_ages_450k_F7_dnaMethyAge <- compute_epigenetic_ages(
  samplesheet_path = paste0(postprocess_path, "samplesheet_F7_450k.csv"),
  methylation_rds_path = paste0(postprocess_path, "DNAm.imputed_wins_F7_450k.rds"),
  output_path = output,
  clocks = c("PCGrimAge", "ShirebyG2020", "BernabeuE2023c"),
  edited_methyAge_path = paste0(functions_path, "methyAge_edited.R"),
  array = "450k"
)

# 15up 450k
epi_ages_450k_15up_dnaMethyAge <- compute_epigenetic_ages(
  samplesheet_path = paste0(postprocess_path, "samplesheet_up15_450k.csv"),
  methylation_rds_path = paste0(postprocess_path, "DNAm.imputed_wins_up15_450k.rds"),
  output_path = output,
  clocks = c("PCGrimAge", "ShirebyG2020", "BernabeuE2023c"),
  edited_methyAge_path = paste0(functions_path, "methyAge_edited.R"),
  array = "450k"
)

# 15up epic
epi_ages_epic_15up_dnaMethyAge <- compute_epigenetic_ages(
  samplesheet_path = paste0(postprocess_path, "samplesheet_up15_epic.csv"),
  methylation_rds_path = paste0(postprocess_path, "DNAm.imputed_wins_up15_EPIC.rds"),
  output_path = output,
  clocks = c("PCGrimAge", "ShirebyG2020", "BernabeuE2023c"),
  edited_methyAge_path = paste0(functions_path, "methyAge_edited.R"),
  array = "epic"
)

# save each file
writexl::write_xlsx(epi_ages_450k_cord_dnaMethyAge, paste0(output, "epi_ages_450k_cord_pred_dnaMethyAge.xlsx"))
writexl::write_xlsx(epi_ages_450k_F7_dnaMethyAge, paste0(output, "epi_ages_450k_F7_pred_dnaMethyAge.xlsx"))
writexl::write_xlsx(epi_ages_450k_15up_dnaMethyAge, paste0(output, "epi_ages_450k_15up_pred_dnaMethyAge.xlsx"))
writexl::write_xlsx(epi_ages_epic_15up_dnaMethyAge, paste0(output, "epi_ages_epic_15up_pred_dnaMethyAge.xlsx"))


# _______________________________________________#
#             Install calcPCBrainAge 
# _______________________________________________#
#devtools::install_github("MorganLevineLab/calcPCBrainAge")
library(calcPCBrainAge)
## VB: github version does not work, needed some changes, so sourcing the code instead
source(paste0(functions_path, "calcPCBrainAge_edited.R"))

# load methylation data
DNAm_450k_cord <- readRDS(paste0(postprocess_path, "DNAm.imputed_wins_cord_450k.rds"))
DNAm_450k_F7 <- readRDS(paste0(postprocess_path, "DNAm.imputed_wins_F7_450k.rds"))
DNAm_450k_15up <- readRDS(paste0(postprocess_path, "DNAm.imputed_wins_up15_450k.rds"))
DNAm_epic_15up <- readRDS(paste0(postprocess_path, "DNAm.imputed_wins_up15_EPIC.rds"))

# transpose meth data as required by calcPCBrainAge
DNAm_450k_cord_T <- as.data.frame(t(DNAm_450k_cord))
DNAm_450k_F7_T <- as.data.frame(t(DNAm_450k_F7))
DNAm_450k_15up_T <- as.data.frame(t(DNAm_450k_15up))
DNAm_epic_15up_T <- as.data.frame(t(DNAm_epic_15up))

# calculate PCBrainAge
result_cord <- calcPCBrainAge(DNAm_450k_cord_T)
result_F7 <- calcPCBrainAge(DNAm_450k_F7_T)
result_15up <- calcPCBrainAge(DNAm_450k_15up_T)
result_15up_epic <- calcPCBrainAge(DNAm_epic_15up_T)

# save to excel
dataset_name = "450k_cord"
write.xlsx(result_cord, file = paste0(output, "epi_ages_", dataset_name, "_pred_calcPCBrainAge.xlsx"), rowNames = FALSE)
dataset_name = "450k_F7"
write.xlsx(result_F7, file = paste0(output, "epi_ages_", dataset_name, "_pred_calcPCBrainAge.xlsx"), rowNames = FALSE)
dataset_name = "450k_15up"
write.xlsx(result_15up, file = paste0(output, "epi_ages_", dataset_name, "_pred_calcPCBrainAge.xlsx"), rowNames = FALSE)
dataset_name = "epic_15up"
write.xlsx(result_15up_epic, file = paste0(output, "epi_ages_", dataset_name, "_pred_calcPCBrainAge.xlsx"), rowNames = FALSE)

