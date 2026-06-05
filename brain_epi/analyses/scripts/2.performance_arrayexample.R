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
# if multiple tissue types or arrays were use this must be added to: 
# <domain>_<timepoint>_<array> or <domain>_<timepoint>_<tissue>
# to faciliate merging across modalities, epi specifiers (array or
# must also be included in the brain files)

# brain 
brain_15up_450k <- readRDS("data/brain_15up_450k_wins.rds") #loading 15up timepoint twice (once as 450k, once as epic) for later merging with epi data, the content of the two DFs would be exactly the same 

brain_15up_epic <- readRDS("data/brain_15up_epic_wins.rds")

# epi 
epi_cord <- readRDS("data/epi_cord_wins.rds")

epi_F7 <- readRDS("data/epi_F7_wins.rds")

epi_15up_450k <- readRDS("data/epi_15up_450k_wins.rds") #two 15up DFs due to 2 different arrays used

epi_15up_epic <- readRDS("data/epi_15up_epic_wins.rds")

# ---- some data manipulation  ----

# exclude variables that we don't want to include

# for epi cord we will remove the non-converted gestational clocks
excluded_models <- c("Bohlin_predage", "EPIC_predage", "Knight_predage")
epi_cord <- epi_cord[, !(names(epi_cord) %in% excluded_models)]

# note that some performance metrics are not meaningful for some age 
# vars, e.g., MAE is not interpretable for Dunedin clocks & DNAmTL

# _______________________________________________________#
# ---- Assessing performance per domain & timepoint  ----
# _______________________________________________________#

# ---- Module 1: create PAD/PAR  ----

# apply function add_PAD_PAR separately for each time point,
# age_var = AGE_brain for brain ages, AGE_epi for epi ages 
# output: appends PAD and PAR measures to the DFs except for DunedinPACE & PACNI

brain_15up_450k <- add_PAD_PAR(brain_15up_450k, age_var = "AGE_brain")
brain_15up_epic <- add_PAD_PAR(brain_15up_epic, age_var = "AGE_brain")

epi_cord <- add_PAD_PAR(epi_cord, age_var = "AGE_epi") 
epi_F7 <- add_PAD_PAR(epi_F7, age_var = "AGE_epi")
epi_15up_450k <- add_PAD_PAR(epi_15up_450k, age_var = "AGE_epi")
epi_15up_epic <- add_PAD_PAR(epi_15up_epic, age_var = "AGE_epi")

# ---- Module 2: performance metrics  ----

# please run function compute_performance_stats once per domain and timepoint
  # & adjust the arguments accordingly

# outputs CSV files with performance metrics under /results/performance 
# for each var defined in stratify_vars and additional CSV file will be created 

# stratify_vars: please retain SEX (_brain/_epi) and add the following variables
  # if they are relevant for your study ( & have not been adjusted for so far)
  # brain age: study or scanner/coil if you data is from different substudies,
  # scanners or different scanning parameters

# if your epi data had different tissue types or arrays, you need to run this 
# function once per subset but for brain data you can run it once only, as we've 
# created e.g., array specific datasets for later merging with the epi data 

# brain timepoints

# brain 15up - as brain data is the same across the _450k and _epic files we only need to run the script on one of them
compute_performance_stats(
  df = brain_15up_450k, #specify the df usually <domain>_<timepoint>
  age_var = "AGE_brain", #AGE_brain for brain ages, AGE_epi for epi ages 
  domain = "brain", #brain or epi 
  timepoint = "15up_450k", #specify timepoint (cohort specific naming)
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

# epi 15up - 450k
# two timepoint 15up epi chunks because 2 arrays were used 
compute_performance_stats(
  df = epi_15up_450k,
  age_var = "AGE_epi",
  domain = "epi",
  timepoint = "15up_450k",
  output_dir = "results/performance",
  stratify_vars = c("SEX_epi")
)

# epi 15up - epic
compute_performance_stats(
  df = epi_15up_epic,
  age_var = "AGE_epi",
  domain = "epi",
  timepoint = "15up_epic",
  output_dir = "results/performance",
  stratify_vars = c("SEX_epi")
)

# ---- Module 3: Corr / Cov Matrices  ----

# please run function generate_corr_cov_matrices once per domain and timepoint &
  # adjust the arguments accordingly 

# outputs CSV files with correlation & covariance matrices + long format version 
  # with p-values, all under /results/performance/corr_cov

# if your epi data had different tissue types or arrays, you need to run this 
# function once per subset but for brain data you can run it once only, as we've 
# created e.g., array specific datasets for later merging with the epi data 

# brain timepoints

# brain 15up - as brain data is the same across the _450k and _epic files we only need to run the script on one of them
generate_corr_cov_matrices(
  df = brain_15up_450k, #specify the df usually <domain>_<timepoint> 
  age_var = "AGE_brain", #age variable (AGE_brain for brain, AGE_epi for epi)
  domain= "brain", #brain or epi 
  timepoint = "15up_450k", #timepoint (cohort specific)
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

# epi 15up - 450k
# two timepoint 15up epi chunks because 2 arrays were used 
generate_corr_cov_matrices(
  df = epi_15up_450k,
  age_var = "AGE_epi",
  domain= "epi",
  timepoint = "15up_450k",
  output_dir = "results/performance/corr_cov"
)

# epi 15up - epic
generate_corr_cov_matrices(
  df = epi_15up_epic,
  age_var = "AGE_epi",
  domain= "epi",
  timepoint = "15up_epic",
  output_dir = "results/performance/corr_cov"
)

# ---- Module 4: Plots  ----

# please run function generate_model_plots once per domain and timepoint &
  # adjust the argument accordingly 

# outputs boxplots, corrplots & scatterplots under /results/performance/plots

# see notes on tissue & array for the previous modules

# brain timepoints

# brain 15up - as brain data is the same across the _450k and _epic files we only need to run the script on one of them
generate_model_plots(
  df = brain_15up_450k, #specify the df usually <domain>_<timepoint>
  age_var = "AGE_brain", #age variable (AGE_brain for brain, AGE_epi for epi)
  domain = "brain", #brain or epi
  timepoint = "15up_450k", #timepoint (cohort specific)
  output_dir = "results/performance/plots", #output dir
  include_categorical = "SEX_brain", #categorical variables to include in additional corrplots 
  sample_name = "ALSPAC" #sample name for plot titles 
)

# epi timepoints 

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

# epi 15up - 450k
# two timepoint 15up epi chunks because 2 arrays were used 
generate_model_plots(
  df = epi_15up_450k,
  age_var = "AGE_epi",
  domain = "epi",
  timepoint = "15up_450k",
  output_dir = "results/performance/plots",
  include_categorical = "SEX_epi",
  sample_name = "ALSPAC"
)

# epi 15up - epic
generate_model_plots(
  df = epi_15up_epic,
  age_var = "AGE_epi",
  domain = "epi",
  timepoint = "15up_epic",
  output_dir = "results/performance/plots",
  include_categorical = "SEX_epi",
  sample_name = "ALSPAC"
)

# __________________________________________________________________#
# ---- Assessing performance for subsample with brain & epi data ----
# __________________________________________________________________#

# In the next step we will apply the performance, corr/cov, and an adapted
  # corr plot script to the data that will go into the association analyses
  # i.e., timepoints with both brain and epi data 

# ---- Data prep -----

# 1) if you have a timepoint WITHOUT further subsetting by array (or tissue)

# merge data initially for subsetting 
#dat_15up <- merge(brain_15up, epi_15up, by = "SUBJID")

# subset brain and epi DFs to participants across both for performance metrics
#brain_15up_overlap <- brain_15up[brain_15up$SUBJID %in% dat_15up$SUBJID, ]
#epi_15up_overlap <- epi_15up[epi_15up$SUBJID %in% dat_15up$SUBJID, ]

# 2) if you have a timepoint WITH further subsetting by array (or tissue)

# merge data initially for subsetting 
dat_15up_450k <- merge(brain_15up_450k, epi_15up_450k, by = "SUBJID")
dat_15up_epic <- merge(brain_15up_epic, epi_15up_epic, by = "SUBJID")

# subset brain and epi DFs to participants across both for performance metrics
brain_15up_450k_overlap <- brain_15up_450k[brain_15up_450k$SUBJID %in% dat_15up_450k$SUBJID, ]
brain_15up_epic_overlap <- brain_15up_epic[brain_15up_epic$SUBJID %in% dat_15up_epic$SUBJID, ]

epi_15up_450k_overlap <- epi_15up_450k[epi_15up_450k$SUBJID %in% dat_15up_450k$SUBJID, ]
epi_15up_epic_overlap <- epi_15up_epic[epi_15up_epic$SUBJID %in% dat_15up_epic$SUBJID, ]

# ---- Module 1: Add PAD/PAR (overlap)  -----

# need to regenerated PAR for the subsample as PAR is a measure relative to the whole sample 

# 1) if you have a timepoint WITHOUT further subsetting by array (or tissue)

# brain_15up_overlap <- add_PAD_PAR(brain_15up_overlap, age_var = "AGE_brain")

# epi_15up_overlap <- add_PAD_PAR(epi_15up_overlap, age_var = "AGE_epi")

# 2) if you have a timepoint WITH further subsetting by array (or tissue)

brain_15up_450k_overlap <- add_PAD_PAR(brain_15up_450k_overlap, age_var = "AGE_brain")
brain_15up_epic_overlap <- add_PAD_PAR(brain_15up_epic_overlap, age_var = "AGE_brain")

epi_15up_450k_overlap <- add_PAD_PAR(epi_15up_450k_overlap, age_var = "AGE_epi")
epi_15up_epic_overlap <- add_PAD_PAR(epi_15up_epic_overlap, age_var = "AGE_epi")


# ---- Inbetween step - merge files again  -----

# epi and brain data per timepoint now need to be combined again after new PAR values have been added 

# merge epi and brain DFs for timepoints with both 
# please name them dat_<timepoint>
# and add additional suffixes for array and/or tissue if needed 
# for ALSPAC this is only 1 timepoint but please repeat if you have multiple 

# 1) if you have a timepoint WITHOUT further subsetting by array (or tissue)

# dat_15up <- merge(brain_15up_overlap, epi_15up_overlap, by = "SUBJID")
# 
# #inspect 
# colnames(dat_15up)
# 
# #AGE_gest var might exist in both if so run the below code
# dat_15up$AGE_gest <- dat_15up$AGE_gest.x
# dat_15up$AGE_gest.x <- NULL
# dat_15up$AGE_gest.y <- NULL

# 2) if you have a timepoint WITH further subsetting by array (or tissue)

# 450k 
dat_15up_450k <- merge(brain_15up_450k_overlap, epi_15up_450k_overlap, by = "SUBJID")

#inspect 
colnames(dat_15up_450k)

#AGE_gest var might exist in both if so run the below code
dat_15up_450k$AGE_gest <- dat_15up_450k$AGE_gest.x
dat_15up_450k$AGE_gest.x <- NULL
dat_15up_450k$AGE_gest.y <- NULL

# epic
dat_15up_epic <- merge(brain_15up_epic_overlap, epi_15up_epic_overlap, by = "SUBJID")

#inspect 
colnames(dat_15up_epic)

#AGE_gest var might exist in both if so run the below code
dat_15up_epic$AGE_gest <- dat_15up_epic$AGE_gest.x
dat_15up_epic$AGE_gest.x <- NULL
dat_15up_epic$AGE_gest.y <- NULL

# ---- Module 2: Performance metrics (overlap)  -----

# performance metrics are now regenerated separately per domain and timepoint
# but for the subsample that will go into the association analyses 
# saved under /results/performance/overlap
# please perform this step PER domain on <domain>_<timepoint>_overlap, not on the 
# combined dataset dat_<timepoint>

# 1) if you have a timepoint WITHOUT further subsetting by array (or tissue)

# # for brain at 15up
# compute_performance_stats(
#   df = brain_15up_overlap, #specify the df usually <domain>_<timepoint>_overlap, NOT dat_<timepoint>_overlap
#   age_var = "AGE_brain", #AGE_brain for brain ages, AGE_epi for epi ages 
#   domain = "brain", #brain or epi 
#   timepoint = "15up", #specify timepoint (cohort specific naming)
#   output_dir = "results/performance/overlap", #output dir, assumes you’re in brain_epi/analyses
#   stratify_vars = c("SEX_brain", "study") #provide categorical variables that you want to produce additional stratified results for 
# )
# 
# # for epi at 15up
# compute_performance_stats(
#   df = epi_15up_overlap, 
#   age_var = "AGE_epi", 
#   domain = "epi", 
#   timepoint = "15up", 
#   output_dir = "results/performance/overlap", 
#   stratify_vars = c("SEX_epi")
# )

# 2) if you have a timepoint WITH further subsetting by array (or tissue)

# for brain at 15up combined with epi 450k
compute_performance_stats(
  df = brain_15up_450k_overlap, #specify the df usually <domain>_<timepoint>_overlap, NOT dat_<timepoint>_overlap
  age_var = "AGE_brain", #AGE_brain for brain ages, AGE_epi for epi ages 
  domain = "brain", #brain or epi 
  timepoint = "15up_450k", #specify timepoint (cohort specific naming)
  output_dir = "results/performance/overlap", #output dir, assumes you’re in brain_epi/analyses
  stratify_vars = c("SEX_brain", "study") #provide categorical variables that you want to produce additional stratified results for 
)

# for brain at 15up combined with epi epic 
compute_performance_stats(
  df = brain_15up_epic_overlap, #specify the df usually <domain>_<timepoint>_overlap, NOT dat_<timepoint>_overlap
  age_var = "AGE_brain", #AGE_brain for brain ages, AGE_epi for epi ages 
  domain = "brain", #brain or epi 
  timepoint = "15up_epic", #specify timepoint (cohort specific naming)
  output_dir = "results/performance/overlap", #output dir, assumes you’re in brain_epi/analyses
  stratify_vars = c("SEX_brain", "study") #provide categorical variables that you want to produce additional stratified results for 
)

# for epi at 15up
compute_performance_stats(
  df = epi_15up_450k_overlap, 
  age_var = "AGE_epi", 
  domain = "epi", 
  timepoint = "15up_450k", 
  output_dir = "results/performance/overlap", 
  stratify_vars = c("SEX_epi")
)

compute_performance_stats(
  df = epi_15up_epic_overlap, 
  age_var = "AGE_epi", 
  domain = "epi", 
  timepoint = "15up_epic", 
  output_dir = "results/performance/overlap", 
  stratify_vars = c("SEX_epi")
)

# ---- Module 3: Corr / Cov matrices (overlap)  -----

# correlation & covariance matrices are now regenerated for the subsample that 
# will go into the association analyses based on the merged dataset dat_<timepoint>
# saved under /results/performance/overlap/corr_cov
# please perform this step on the merged dataframe dat_<timepoint>

# if array was no consideration you would include only one of the chunks
#with dat_15up as the DF in this specific example 

generate_corr_cov_matrices(
  df = dat_15up_450k, #specify the merged df, usually dat_<timepoint> 
  age_var = c("AGE_brain","AGE_epi"), #both age variables 
  domain= "brain_epi_overlap", #brain or epi 
  timepoint = "15up_450k", #timepoint (cohort specific)
  output_dir = "results/performance/overlap/corr_cov" #output dir
)

generate_corr_cov_matrices(
  df = dat_15up_epic, #specify the merged df, usually dat_<timepoint> 
  age_var = c("AGE_brain","AGE_epi"), #both age variables 
  domain= "brain_epi_overlap", #brain or epi 
  timepoint = "15up_epic", #timepoint (cohort specific)
  output_dir = "results/performance/overlap/corr_cov" #output dir
)

# ---- Module 4: Plots adapted (overlap)  -----

# corrplots are now regenerated for the subsample that will go into the 
# association analyses based on the merged dataset dat_<timepoint>
# saved under /results/performance/overlap/plots
# please perform this step on the merged dataframe dat_<timepoint>

# if array was no consideration you would include only one of the chunks
#with dat_15up as the DF in this specific example 

generate_corrplots_overlap(
  df = dat_15up_450k, #specify the merged df, usually dat_<timepoint> 
  age_var = c("AGE_brain","AGE_epi"), #both age variables 
  timepoint = "15up_450k", #timepoint (cohort specific)
  output_dir = "results/performance/overlap/plots", #output dir
  include_categorical = "SEX_brain", #categorical variables to include in additional corrplots 
  sample_name = "ALSPAC" #sample name for plot titles 
)

generate_corrplots_overlap(
  df = dat_15up_epic, #specify the merged df, usually dat_<timepoint> 
  age_var = c("AGE_brain","AGE_epi"), #both age variables 
  timepoint = "15up_epic", #timepoint (cohort specific)
  output_dir = "results/performance/overlap/plots", #output dir
  include_categorical = "SEX_brain", #categorical variables to include in additional corrplots 
  sample_name = "ALSPAC" #sample name for plot titles 
)

# __________________________________________________________________#
# ---- Save files with PAD & PAR and combined timepoints ----
# __________________________________________________________________#

# save domain specific 
saveRDS(brain_15up_450k, file = "data/brain_15up_450k_wins_PAR_PAD.rds")
saveRDS(brain_15up_epic, file = "data/brain_15up_epic_wins_PAR_PAD.rds") 

saveRDS(epi_cord, file = "data/epi_cord_wins_PAR_PAD.rds")
saveRDS(epi_F7, file = "data/epi_F7_wins_PAR_PAD.rds")
saveRDS(epi_15up_450k, file = "data/epi_15up_450k_wins_PAR_PAD.rds")
saveRDS(epi_15up_epic, file = "data/epi_15up_epic_wins_PAR_PAD.rds")

# save data from timepoints with both 
saveRDS(dat_15up_450k, file = "data/overlap_brain_epi_15up_450k_wins_PAR_PAD.rds")
saveRDS(dat_15up_epic, file = "data/overlap_brain_epi_15up_epic_wins_PAR_PAD.rds")


#### DONE WITH 2.performance.R ####
