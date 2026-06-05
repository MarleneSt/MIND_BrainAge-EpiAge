## This script preprocesses methylation data separately for each time point and array

# _______________________________________________#
#                       SET UP
# _______________________________________________#
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("impute")
library(impute) # for impute.knn

## load packages needed for postprocessing
#devtools::install_github("markgene/maxprobes")
#library(maxprobes)  # for identifying cross-reactive probes
#install.packages("matrixStats")
library(matrixStats)  # for rowIQR
# source("http://bioconductor.org/biocLite.R")
# install.packages("devtools") # if the devtools package is not installed
# library(devtools)
# install_github("perishky/meffil")
library(meffil)

# install matrixStats package for cpg descriptives
#install.packages("matrixStats")
library(matrixStats)

# set working directory to the "brain_epi/epi" folder we provided 
wd = "/campaign/VB-FM5HPC-001/Vilte/Projects/Marlene-epi-clocks/brain_epi/epi/" # (! update with your path to brain_epi/epi/ folder)
setwd(wd)

# specify path to the folder that contains normalised methylation data (e.g., functionally or quantile normalised) 
# with CpGs in rows and samples in columns
meth_path = paste0(wd, "preprocessed/") 

# _______________________________________________#
#       PREPROCESS DATA: 450K and epic
# _______________________________________________#

# read in meth data for each time point and array
meth_cord_450k <- readRDS(paste0(meth_path,  "meth_cord_450k.rds"))
meth_F7_450k <- readRDS(paste0(meth_path, "meth_F7_450k.rds"))
meth_up15_450k <- readRDS(paste0(meth_path,"meth_up15_450k.rds"))
meth_up15_epic <- readRDS(paste0(meth_path,  "meth_up15_epic.rds"))

# preprocess_methylation() will preprocess meth dataset and 
# save the imputed and winsorized methylation matrix as "DNAm.imputed_wins_array_type_time_point.rds"
source("scripts/functions/preprocess_methylation.R") # function to impute and winsorize meth data
source("scripts/functions/cpg_descriptives.R") # function to obtain per-CpG descriptives (rows) before and after preprocessing

# ! note: 
# define output folder (**relative to working directory, which should be the 'brain_epi/epi' folder we provided**)
output_folder <- file.path(wd, "postprocessed")

# create output folder if it doesn't exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}


## cord 450k 
meth_cord_450k_processed <- preprocess_methylation(meth_matrix = meth_cord_450k, # should be matrix (not data frame), with rows as CpGs and columns as individuals
                                                   remove_sexchr = FALSE, # option available but should be set to FALSE for epi clocks 
                                                   remove_crossreactive = FALSE, # option available but should be set to FALSE for epi clocks 
                                                   array_type = "450k", # either "450k" or "EPIC"
                                                   time_point = "cord", # name your time point 
                                                   output_prefix = file.path(output_folder, "DNAm.imputed_wins"), # will save as paste0(output_prefix, "_", array_type, "_", time_point, ".rds")
                                                   epic_chunk_size = NA, # for EPIC array specify chunk size, e.g. "100000" or "200000" instead of "NA" (not relevant for 405k, as 450k array is imputed all at once)
                                                   log_file = file.path(output_folder, "450k_cord_postprocess_log.txt")) # update with the relevant name for your array and time point 
                                                   #log_file = file.path(output_folder, "logs", "450k_cord_postprocess_log.txt"))
                                                   

## F7 450k 
meth_F7_450k_processed <- preprocess_methylation(meth_matrix = meth_F7_450k, 
                                                   remove_sexchr = FALSE, 
                                                   remove_crossreactive = FALSE, 
                                                   array_type = "450k", 
                                                   time_point = "F7",
                                                   output_prefix = file.path(output_folder, "DNAm.imputed_wins"),
                                                   epic_chunk_size = NA,
                                                   log_file = file.path(output_folder, "450k_F7_postprocess_log.txt"))

## 15up 450k 
meth_up15_450k_processed <- preprocess_methylation(meth_matrix = meth_up15_450k, 
                                                 remove_sexchr = FALSE, 
                                                 remove_crossreactive = FALSE,
                                                 array_type = "450k", 
                                                 time_point = "up15", 
                                                 output_prefix = file.path(output_folder, "DNAm.imputed_wins"), 
                                                 epic_chunk_size = NA,
                                                 log_file = file.path(output_folder, "450k_up15_postprocess_log.txt"))

## 15up epic  
meth_up15_epic_processed <- preprocess_methylation(meth_matrix = meth_up15_epic, 
                                                   remove_sexchr = FALSE,  
                                                   remove_crossreactive = FALSE,
                                                   array_type = "EPIC", 
                                                   time_point = "up15", 
                                                   output_prefix = file.path(output_folder, "DNAm.imputed_wins"),
                                                   epic_chunk_size = 200000, # will apply knn in chunks of 200,000 CpGs for epic array (helps if experiencing memory issues, can adapt chunksize as needed)
                                                   log_file = file.path(output_folder, "EPIC_up15_postprocess_log.txt"))

