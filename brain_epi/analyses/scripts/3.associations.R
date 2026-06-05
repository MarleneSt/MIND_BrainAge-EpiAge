# _______________________________________________#
#         3.associations.R
# _______________________________________________#
# Date: 2025/08/07; last updated 2025/12/10 
#
# Author: MS
#
# Purpose: 
#   - Run association models between brain and epigenetic age outcomes
#   - Loop over 3 outcome types (_predage, _PAD, _PAR) and 3 models 
#       - Model 1: brain age ~ epi age 
#       - Model 2: Model 1 + sex + age_brain + age_epi + batch 
#                               only age_brain included if corr ≥ 0.8
#       - Model 3: Model 2 + cell type 
#   - Produce model estimates using both OLS (lm) and robust regression (rlm)
#   - Save model availability logs and association results (option to create 
#     diagnostic plots is currently disable)
#
# Note: In the pre-registration, we specified using PAR, the other outcomes are only
#   considered so that we have them in case reviewers ask 
#
# Input: 
#   - overlap_brain_epi_<timepoint>_wins_PAR_PAD.rds (from performance script)
#   - configuration settings in site_config (define later in the script)
#
# Output:
#   - Model availability log under /results/associations/
#   - Model outputs as CSVs under /results/associations/
#
# General instructions:
#   - Follow script section by section and make necessary adjustments
#   - the below code is for one timepoint, if you have multiple timepoints copy
#     the whole chunk from '1st TP with brain + epi' onwards and adapt for next 
#     timepoint
#
# _______________________________________________#
# ---- Set UP  ----
# _______________________________________________#

# ---- Set WD to main folder ----
setwd("/brain_epi/analyses/") # main analysis folder (brain_epi/analyses)

# ---- Load packages ----
if (!require(sandwich)) install.packages("sandwich"); library(sandwich)
if (!require(lmtest)) install.packages("lmtest"); library(lmtest)
if (!require(MASS)) install.packages("MASS"); library(MASS) 
if (!require(car)) install.packages("car"); library(car)
if (!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)

# ---- Source scripts ----
source("scripts/functions/association_helpers.R")

# ---- General configurations ----

var_types <- c("_predage", "_PAR", "_PAD") #defining the 3 different outcome variables
model_nums <- 1:3 # 3 models will be run, model_nums used for output files 

dir.create("results/associations/", recursive = TRUE, showWarnings = FALSE) #create associations subfolder if it doesn't exist already

# ---- load data  ----
# please load data for all timepoints with brain & epi data 
# these data should have been saved during the 2.performance script 
# naming: overlap_brain_epi_<timepoint>_wins_PAR_PAD.rds
# saved under /data 

# the following example is for ALSPAC which only has one timepoint for associations
# but for cohorts with multiple timepoints, please load all data here 
# DF name should always be dat_<timepoint>

dat_15up <- readRDS("data/overlap_brain_epi_15up_wins_PAR_PAD.rds")

#dat_<timepoint> <- readRDS("data/overlap_brain_epi_<timepoint>_wins_PAR_PAD.rds")

# ---- check var type accuracy  ----

# while this should have been handled by previous steps, please check if variable
# types are correct via str(dat_<timepoint)
# i.e., numeric vars = numeric (e.g., all the age variables)
# categorical vars = factors (e.g., SEX)
# this will ensure that the functions to generate standardised betas work 
str(dat_15up)

# _______________________________________________#
# ---- 1st TP with brain + epi ----
# _______________________________________________#

# we will now commence with the association analyses for the first timepoint 

# ---- assign data frame  ----
data <- dat_15up

# ---- site configuration ----

# site_config helps to configure all variables for the 3 models 
# that are run in the association analyses as the naming of batch 
# and cell-type variables will differ between sites 

# please use the below to define your sample & timepoint specific 
#   - batch variables (batch_vars)
#   - cell types (cell_type_vars), to circumvent issues with
#     multicollinearity, please remove one of the cell types from the
#     chosen panel
#   - extra covariates that need to be adjusted for; these should only be any of
#     the following and should only be included if 1) relevant to your sample 
#     and 2) have not already been adjusted for 
#         - study - if imaging or DNAm or imaging data comes from different 
#           substudies that need to be accounted for 
#         - scanner or coil or similar - if imaging data comes from different 
#           scanners/scanner settings etc. 
#         !!!! Please do not include any other additional covariates without 
#         !!!! Discussion with Marlene (ms2290@bath.ac.uk)
#   - timepoint_var - string variable to indicate TP name (used for output naming)
#   - sample_name - string variable to indicate sample name (used for output naming)

site_config <- list(
  batch_vars = c("V1", "V4", "V5", "V7", "V8", "V9"), #include your batch variables here (if needed)
  cell_type_vars = c("blood.idoloptimised.CD4T", "blood.idoloptimised.CD8T" , "blood.idoloptimised.Mono", "blood.idoloptimised.Neu","blood.idoloptimised.NK"), #cell type minus 1 category 
  extra_covariates = NULL, #additional study-specific technical covariates can be included here e.g., extra_covariates = c("scanner"), these wull be included in model 2+3
  timepoint_var = "15up",  
  sample_name = "ALSPAC" 
)

# ---- log model availability ----

# the function check_model_availability_agediff produces a log under 
# /results/associations which indicates which models are available & differences
# between age and MRI and DNAm assessment

# no changes required
# but please check the logs for issues 
# (e.g., a model indicatesd as NA when you think it should be there or implausible
# age differences)
# the function's outputs can then be used to define the available brain & epi 
# variables for analyses 

# output: log file under results/associations

availability <- check_model_availability_agediff(data, site_config)

# ---- define brain & epi vars for analysis  ----

# the outputs of check_model_availability can be used to define the available 
# brain & epi variables for analyses 

brain_vars <- availability$brain_vars
epi_vars   <- availability$epi_vars

# ---- run models ----

# the following loop will perform the association analyses across all available 
# brain age ~ epi age combinations 

# IMPORTANT - no changes are required to the loop it can just be run as is 

# 3 models will be run (see above) using 1) lm and 2) rlm
# outputs:
#   - CSV files with model outputs per model & outcome var (_predage, _PAR, _PAD)
#           (under results/associations)
#   - diagnostic plots for the LM models under results/associations/diagnostics_lm

for (var_type in var_types) {
  message("📦 Starting association analysis for variable type: ", var_type)
  for (model_num in model_nums) {
    message("🔍 Running Model ", model_num, " for outcome type: ", var_type)
    
    all_results <- list()
    
    lm_models <- list()  # Collect lm models for diagnostics
    
    #subset to available brain/epi vars + Dunedin models that have different suffixes
    brain_subset <- add_unsuffixed(grep(var_type, brain_vars, value=TRUE), var_type, data, is_brain = TRUE)
    epi_subset   <- add_unsuffixed(grep(var_type, epi_vars, value=TRUE), var_type, data, is_brain = FALSE)
    
    for (brain_var in brain_subset) {
      for (epi_var in epi_subset) {
        model_formula <- create_model_formula(
          outcome = brain_var,
          predictor = epi_var,
          model_num = model_num,
          config = site_config,
          data = data,
          var_type = var_type
        )
        
        #extract model formula for tracking & later adding to CSV
        model_string <- gsub("\\s+", " ", paste(deparse(model_formula), collapse = " "))
        
        # Run both lm and rlm
        lm_res <- run_lm_model_with_full_stats(model_formula, data, model_string)
        rlm_res <- run_rlm_model(model_formula, data)
        
        lm_models[[length(lm_models) + 1]] <- list(
          fit = lm_res$fit,
          brain = brain_var,
          epi = epi_var
        )
        
        # Merge results by Predictor
        full_res <- full_join(lm_res$results, rlm_res, by = "Predictor")
        
        # tag model info for CSV
        full_res$Outcome      <- brain_var
        full_res$EpiPredictor <- epi_var
        full_res$Model        <- paste0("Model_", model_num)
        full_res$VariableType <- var_type
        full_res$Timepoint <- site_config$timepoint_var
        full_res$ModelFormula <- model_string
        full_res <- full_res[, c("Outcome", "EpiPredictor", setdiff(names(full_res), c("Outcome", "EpiPredictor")))]
        
        all_results[[length(all_results)+1]] <- full_res
      }
    }

    combined <- bind_rows(all_results)
    
    #save results 
    filename <- paste0("results/associations/", paste0(site_config$sample_name,"_",site_config$timepoint_var, "_model", model_num, "_", gsub("_", "", var_type)), ".csv")
    write.csv(combined, filename, row.names = FALSE)
    message("💾 Done. Results saved to: ", filename)
    
    # Save diagnostic plots to PDF - this is commented out as not needed for all sites
   #  diagnostic_file <- paste0("results/associations/diagnostics_lm/", site_config$sample_name,"_",site_config$timepoint_var,"_" ,"diagnostics_model", model_num, "_", gsub("_", "", var_type), ".pdf")
   # save_combined_lm_diagnostics(lm_models, diagnostic_file)
   # message("📈 Diagnostic plots saved to: ", diagnostic_file)
   
   message("✅ Completed Model ", model_num, " for type: ", var_type, "\n")
  }
  message("🎉 All models completed for variable type: ", var_type, "\n")
}

# _______________________________________________#
# ---- 2nd TP with brain + epi ----
# _______________________________________________#

# please copy over and repeat the above with any additional timepoint

#### DONE WITH 3.associations.R ####

