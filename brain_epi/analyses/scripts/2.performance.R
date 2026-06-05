# _______________________________________________#
#         2.performance.R
# _______________________________________________#
# Date: 2025/08/07; last updated 2025/12/10 
#
# Author: MS
#
# Purpose: 
#   1. Add PAD and PAR metrics to each dataset
#   2. Compute performance metrics (e.g., MAE, RMSE, R²) for all prediction models
#   3. Generate correlation and covariance matrices
#   4. Produce summary plots (boxplots, scatterplots, corrplots)
#   5. Repeat steps 2-4  on merged data for timepoints with both brain + epi data
#
# Input: 
#   - RDS files from 1.epi_brain_clean_and_describe.R (per timepoint)
#
# Output:
#   - Performance metrics (.csv under /results/performance)
#   - Correlation and covariance matrices (.csv under /results/performance/corr_cov)
#   - Plots (.png under /results/performance/plots)
#   - same for overlap under results/performance/overlap
#   - Updated datasets with PAD & PAR (.rds under /data)
#
# General instructions:
#   - Follow script section by section and make necessary adjustments
#   - Run each module per timepoint and domain unless otherwise instructed
#   - Adjust covariates & stratification vars as needed
#
# _______________________________________________#
# ---- Set UP  ----
# _______________________________________________#

# ---- Set WD to main folder ----

setwd("/brain_epi/analyses/") #this should be the main analysis folder (brain_epi/analyses), which has subfolders like /data, /results etc.

# ---- Load packages ----

#for performance 
if (!require(caret)) install.packages("caret"); library(caret)
if (!require(boot)) install.packages("boot"); library(boot)

#for corr / cov matrices
if (!require(psych)) install.packages("psych"); library(psych)
if (!require(Hmisc)) install.packages("Hmisc"); library(Hmisc)
if (!require(dplyr)) install.packages("dplyr"); library(dplyr)
if (!require(tidyr)) install.packages("tidyr"); library(tidyr)

# for plotting 
if (!require(ggplot2)) install.packages("ggplot2"); library(ggplot2)
if (!require(GGally)) install.packages("GGally"); library(GGally)
if (!require(corrplot)) install.packages("corrplot"); library(corrplot)
if (!require(viridis)) install.packages("viridis"); library(viridis)

# ---- Source scripts ----

#adapt paths if needed
source("scripts/functions/add_PAD_PAR.R")
source("scripts/functions/compute_performance_stats.R")
source("scripts/functions/generate_corr_cov_matrices.R")
source("scripts/functions/generate_model_plots.R")
source("scripts/functions/generate_corrplots_overlap.R") #this will only be applied to the overlap data 

# ---- Load data  ----

# load the RDS files saved by 1.epi_brain_clean_and_describe
# please name the DF <domain>_<timepoint> as done below

# brain 
brain_15up <- readRDS("data/brain_15up_wins.rds")

# epi 
epi_cord <- readRDS("data/epi_cord_wins.rds")

epi_F7 <- readRDS("data/epi_F7_wins.rds")

epi_15up <- readRDS("data/epi_15up_wins.rds")

# ---- some data manipulation  ----

# exclude variables that we don't want to include

# for epi cord we will remove the non-converted gestational clocks
excluded_models <- c("Bohlin_predage", "EPIC_predage", "Knight_predage")
epi_cord <- epi_cord[, !(names(epi_cord) %in% excluded_models)]

# _______________________________________________________#
# ---- Assessing performance per domain & timepoint  ----
# _______________________________________________________#

# ---- Module 1: create PAD/PAR  ----

# apply function add_PAD_PAR separately for each time point,
# age_var = AGE_brain for brain ages, AGE_epi for epi ages 
# output: appends PAD and PAR measures to the DFs except for DunedinPACE & PACNI

brain_15up <- add_PAD_PAR(brain_15up, age_var = "AGE_brain")
epi_cord <- add_PAD_PAR(epi_cord, age_var = "AGE_epi") 
epi_F7 <- add_PAD_PAR(epi_F7, age_var = "AGE_epi")
epi_15up <- add_PAD_PAR(epi_15up, age_var = "AGE_epi")

# ---- Module 2: performance metrics  ----

# please run function compute_performance_stats once per domain and timepoint
  # & adjust the arguments accordingly

# outputs CSV files with performance metrics under /results/performance 
# for each var defined in stratify_vars and additional CSV file will be created 

# stratify_vars: please retain SEX (_brain/_epi) and add the following variables
  # if they are relevant for your study ( & have not been adjusted for so far)
  # brain age: study or scanner/coil if you data is from different substudies,
    # scanners or different scanning parameters 
  # epi age: array (if multiple were used at a TP), batch (if it's a categorical var)
# you don't necessarily need to include these variables here if they have already
  # been adjusted for but having the stratified output would still be useful 

# brain timepoints

# brain 15up
compute_performance_stats(
  df = brain_15up, #specify the df usually <domain>_<timepoint>
  age_var = "AGE_brain", #AGE_brain for brain ages, AGE_epi for epi ages 
  domain = "brain", #brain or epi 
  timepoint = "15up", #specify timepoint (cohort specific naming)
  output_dir = "results/performance", #output dir, assumes you’re in brain_epi/analyses
  stratify_vars = c("SEX_brain", "study") #provide categorical variables that you want to produce additional stratified results for, please retain SEX 
)

# epi timepoints

# epi cord
compute_performance_stats(
  df = epi_cord,
  age_var = "AGE_epi",
  domain = "epi",
  timepoint = "cord",
  output_dir = "results/performance",
  stratify_vars = c("SEX_epi")
)

# epi F7
compute_performance_stats(
  df = epi_F7,
  age_var = "AGE_epi",
  domain = "epi",
  timepoint = "F7",
  output_dir = "results/performance",
  stratify_vars = c("SEX_epi")
)

# epi 15up
compute_performance_stats(
  df = epi_15up,
  age_var = "AGE_epi",
  domain = "epi",
  timepoint = "15up",
  output_dir = "results/performance",
  stratify_vars = c("SEX_epi")
)

# ---- Module 3: Corr / Cov Matrices  ----

# please run function generate_corr_cov_matrices once per domain and timepoint &
  # adjust the arguments accordingly 

# outputs CSV files with correlation & covariance matrices + long format version 
  # with p-values, all under /results/performance/corr_cov

# brain timepoints

# brain 15up
generate_corr_cov_matrices(
  df = brain_15up, #specify the df usually <domain>_<timepoint> 
  age_var = "AGE_brain", #age variable (AGE_brain for brain, AGE_epi for epi)
  domain= "brain", #brain or epi 
  timepoint = "15up", #timepoint (cohort specific)
  output_dir = "results/performance/corr_cov" #output dir
)

# epi timepoints 

# epi cord
generate_corr_cov_matrices(
  df = epi_cord,
  age_var = "AGE_epi",
  domain= "epi",
  timepoint = "cord",
  output_dir = "results/performance/corr_cov"
)

# epi F7
generate_corr_cov_matrices(
  df = epi_F7,
  age_var = "AGE_epi",
  domain= "epi",
  timepoint = "F7",
  output_dir = "results/performance/corr_cov"
)

# epi 15up
generate_corr_cov_matrices(
  df = epi_15up,
  age_var = "AGE_epi",
  domain= "epi",
  timepoint = "15up",
  output_dir = "results/performance/corr_cov"
)

# ---- Module 4: Plots  ----

# please run function generate_model_plots once per domain and timepoint &
  # adjust the argument accordingly 

# outputs boxplots, corrplots & scatterplots under /results/performance/plots

# brain timepoints

# brain 15up
generate_model_plots(
  df = brain_15up, #specify the df usually <domain>_<timepoint>
  age_var = "AGE_brain", #age variable (AGE_brain for brain, AGE_epi for epi)
  domain = "brain", #brain or epi
  timepoint = "15up", #timepoint (cohort specific)
  output_dir = "results/performance/plots", #output dir
  include_categorical = "SEX_brain", #categorical variables to include in additional corrplots 
  sample_name = "ALSPAC" #sample name for plot titles 
)

# epi cord
generate_model_plots(
  df = epi_cord,
  age_var = "AGE_epi",
  domain = "epi",
  timepoint = "cord",
  output_dir = "results/performance/plots",
  include_categorical = "SEX_epi",
  sample_name = "ALSPAC"
)

# epi F7
generate_model_plots(
  df = epi_F7,
  age_var = "AGE_epi",
  domain = "epi",
  timepoint = "F7",
  output_dir = "results/performance/plots",
  include_categorical = "SEX_epi",
  sample_name = "ALSPAC"
)

# epi 15up
generate_model_plots(
  df = epi_15up,
  age_var = "AGE_epi",
  domain = "epi",
  timepoint = "15up",
  output_dir = "results/performance/plots",
  include_categorical = "SEX_epi",
  sample_name = "ALSPAC"
)

# __________________________________________________________________#
# ---- Assessing performance for subsample with brain & epi data ----
# __________________________________________________________________#

# In the next step we will apply the performance, corr/cov, and an adapted
  # corr plot scripts to the data that will go into the association analyses
  # i.e., timepoints with both brain and epi data 
# example: one timepoint (15up) in ALSPAC

# ---- Data prep -----

# merge data initially for subsetting 
dat_15up <- merge(brain_15up, epi_15up, by = "SUBJID")

# subset brain and epi DFs to participants across both for performance metrics
brain_15up_overlap <- brain_15up[brain_15up$SUBJID %in% dat_15up$SUBJID, ]
epi_15up_overlap <- epi_15up[epi_15up$SUBJID %in% dat_15up$SUBJID, ]

# ---- Module 1: Add PAD/PAR (overlap)  -----

# need to regenerated PAR for the subsample as PAR is a measure relative to the whole sample 

brain_15up_overlap <- add_PAD_PAR(brain_15up_overlap, age_var = "AGE_brain")

epi_15up_overlap <- add_PAD_PAR(epi_15up_overlap, age_var = "AGE_epi")


# ---- Inbetween step - merge files again  -----

# epi and brain data per timepoint now need to be combined again after new PAR values have been added 

# merge epi and brain DFs for timepoints with both 
# please name them dat_<timepoint>
# for ALSPAC this is only 1 timepoint but please repeat if you have multiple 
dat_15up <- merge(brain_15up_overlap, epi_15up_overlap, by = "SUBJID")

#inspect 
colnames(dat_15up)

#AGE_gest var might exist in both if so run the below code
dat_15up$AGE_gest <- dat_15up$AGE_gest.x
dat_15up$AGE_gest.x <- NULL
dat_15up$AGE_gest.y <- NULL

# ---- Module 2: Performance metrics (overlap)  -----

# performance metrics are now regenerated separately per domain and timepoint
# but for the subsample that will go into the association analyses 
# saved under /results/performance/overlap
# please perform this step PER domain on <domain>_<timepoint>_overlap, not on the 
# combined dataset dat_<timepoint>

# for brain at 15up
compute_performance_stats(
  df = brain_15up_overlap, #specify the df usually <domain>_<timepoint>_overlap, NOT dat_<timepoint>_overlap
  age_var = "AGE_brain", #AGE_brain for brain ages, AGE_epi for epi ages 
  domain = "brain", #brain or epi 
  timepoint = "15up", #specify timepoint (cohort specific naming)
  output_dir = "results/performance/overlap", #output dir, assumes you’re in brain_epi/analyses
  stratify_vars = c("SEX_brain", "study") #provide categorical variables that you want to produce additional stratified results for 
)

# for epi at 15up
compute_performance_stats(
  df = epi_15up_overlap, 
  age_var = "AGE_epi", 
  domain = "epi", 
  timepoint = "15up", 
  output_dir = "results/performance/overlap", 
  stratify_vars = c("SEX_epi")
)

# ---- Module 3: Corr / Cov matrices (overlap)  -----

# correlation & covariance matrices are now regenerated for the subsample that 
# will go into the association analyses based on the merged dataset dat_<timepoint>
# saved under /results/performance/overlap/corr_cov
# please perform this step on the merged dataframe dat_<timepoint>

generate_corr_cov_matrices(
  df = dat_15up, #specify the merged df, usually dat_<timepoint> 
  age_var = c("AGE_brain","AGE_epi"), #both age variables 
  domain= "brain_epi_overlap", #brain or epi 
  timepoint = "15up", #timepoint (cohort specific)
  output_dir = "results/performance/overlap/corr_cov" #output dir
)

# ---- Module 4: Plots adapted (overlap)  -----

# corrplots are now regenerated for the subsample that will go into the 
# association analyses based on the merged dataset dat_<timepoint>
# saved under /results/performance/overlap/plots
# please perform this step on the merged dataframe dat_<timepoint>

generate_corrplots_overlap(
  df = dat_15up, #specify the merged df, usually dat_<timepoint> 
  age_var = c("AGE_brain","AGE_epi"), #both age variables 
  timepoint = "15up", #timepoint (cohort specific)
  output_dir = "results/performance/overlap/plots", #output dir
  include_categorical = "SEX_brain", #categorical variables to include in additional corrplots 
  sample_name = "ALSPAC" #sample name for plot titles 
)

# __________________________________________________________________#
# ---- Save files with PAD & PAR and combined timepoints ----
# __________________________________________________________________#

# save domain specific 
saveRDS(brain_15up, file = "data/brain_15up_wins_PAR_PAD.rds")
saveRDS(epi_cord, file = "data/epi_cord_wins_PAR_PAD.rds")
saveRDS(epi_F7, file = "data/epi_F7_wins_PAR_PAD.rds")
saveRDS(epi_15up, file = "data/epi_15up_wins_PAR_PAD.rds")

# save data from timepoints with both 
saveRDS(dat_15up, file = "data/overlap_brain_epi_15up_wins_PAR_PAD.rds")

#### DONE WITH 2.performance.R ####
