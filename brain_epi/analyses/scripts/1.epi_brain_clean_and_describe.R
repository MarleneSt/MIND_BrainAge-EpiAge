# _______________________________________________#
#         1.epi_brain_clean_and_describe.R
# _______________________________________________#
# Date: 2025/08/07; last updated 2025/12/10 
#
# Author: MS
#
# Purpose: 
#   - Perform some simple cleaning tasks on epi and brain age data 
#   - extract descriptives pre and post winsorisation 
#   - perform winsorisation
#
# Input: CSV files per timepoint for epi & brain ages
#           incl. the covs outlined in the instructions
#
# Output:
#   1) RDS files pre+post winsorisation with correct var types (/data)
#   2) descriptives & log per timepoint (/results/descriptives)
#
# General instructions: 
#   - Follow script section by section and make necessary adjustments
#   - Run each module per timepoint and domain unless otherwise instructed
#   - Adjust covariates & stratification vars as needed

# _______________________________________________#
# ---- Set UP  ----
# _______________________________________________#

# ---- Set WD to main folder ----

#set WD to main analysis folder, which has the epi and brain age data in subfolders, e.g., /brain_epi/analyses

setwd("/brain_epi/analyses/")

# ---- libraries ----

packages <- c("dplyr", "purrr", "summarytools", "gt")

for (pkg in packages) {
  if (!suppressWarnings(require(pkg, character.only = TRUE))) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ---- source functions  ----

# align path if needed 
source("scripts/functions/prep_and_descriptives_epi_function.R")
source("scripts/functions/prep_and_descriptives_brain_function.R")
source("scripts/functions/descriptives_brain_epi_overlap_function.R")

# _______________________________________________#
# ---- EPI DATA  ----
# _______________________________________________#

# ----- defining epi variables -----

#these are the included epi ages, including the gestational ones
#names align with the names in the instructions
# please make sure your names align and do not change the names here 

common_epi_vars <- c(
  'PCBrainAge_predage', 'PCGrimAge_predage', 'CorticalClock_predage', 'cAge_predage', 'AltumAge_predage', 'Horvath2013_predage', 'DamAge_predage', 'AdaptAge_predage', 'DunedinPACE', 'Hannum_predage', 'PhenoAge_predage', 'skinHorvath_predage', 'PedBE_predage', 'Wu_predage', 'DNAmTL_predage', 'ZhangBLUP_predage', 'ZhangEN_predage', 'Bohlin_predage', 'EPIC_predage', 'Knight_predage', 'Bohlin_conv_predage', 'EPIC_conv_predage', 'Knight_conv_predage'
)

# ---- Configuration per timepoint ----

# here we define paths to the data from different time points 
# and provide the time point names (e.g., ALSPAC = cord, F7 and 15up) 
# we also define which epi variables to exclude from the winsorisation
# for cord this should always be the gestational clocks prior to conversion 
# "Bohlin_predage", "EPIC_predage", "Knight_predage"
# for the other time points, all gestational clocks should be excluded 

# the below example is for the ALSPAC epi data for 3 timepoints (cord/birth, F7, 15up)
# if you have more timepoints, add more lists into timepoint_configs

timepoint_configs <- list(
  cord = list( #cord = name of time point
    file = "data/epi_ages/epi_ages_covs_cord.csv", #path to CSV files
    exclude_vars = c("Bohlin_predage", "EPIC_predage", "Knight_predage")
  ), # name variables to exclude 
  F7 = list(
    file = "data/epi_ages/epi_ages_covs_F7.csv",
    exclude_vars = c("Bohlin_conv_predage", "EPIC_conv_predage", "Knight_conv_predage", "Bohlin_predage", "EPIC_predage", "Knight_predage") 
  ),
  `15up` = list(
    file = "data/epi_ages/epi_ages_covs_15up.csv",
    exclude_vars = c("Bohlin_conv_predage", "EPIC_conv_predage", "Knight_conv_predage", "Bohlin_predage", "EPIC_predage", "Knight_predage")
  )
) 

# ---- Run cleaning & descriptives function  ----

# the function prep_and_descr_timepoint_epi has to be run per timepoint 
# it will do the following:
  # rename variables for prettier plotting 
  # set variable types 
  # provide descriptives pre+post winsorisation + log file (saved under /results/descriptives) 
  #winsorise epi ages 

# !!!! IMPORTANT REQUIRED CHANGES !!!!!
  # the script assumes the covariates and epi clocks names from the instructions
  # document. Expected categorical variables (i.e., SEX_epi, batch) are set to factor
  # while all other vars are set to numeric
  # if you inlcuded a numeric batch variable (e.g., SVs), please open the
  # function prep_and_descr_timepoint_epi and remove 'batch' from the line
  # cat_candidates <- c("batch") - i.e., cat_candidates <- c("batch")
  # IF you have included additional technical covariates (after discussion with
  # Marlene), these ca be added in the cat_candidates vector

# Warning messages may appear if some columns are completely NA 
  # Warning messages:
  # 1: In min(x, na.rm = TRUE) :
  # no non-missing arguments to min; returning Inf
  # this happens when using summarytools if some of your variables are
  # completely NA (e.g., a model that could not be calculated) --> can usually be ignored

prep_and_descr_timepoint_epi("cord", common_epi_vars, timepoint_configs)

prep_and_descr_timepoint_epi("F7", common_epi_vars, timepoint_configs)

prep_and_descr_timepoint_epi("15up", common_epi_vars, timepoint_configs)

# _______________________________________________#
# ---- BRAIN DATA  ----
# _______________________________________________#

# The workflow for the brain data closely follows that for the epi clocks 

# ---- defining brain variables  ----

#these are the included brain ages
#names align with the names in the instructions
# please make sure your names align and do not change the names here 

common_brain_vars <- c('Centile2_predage', 'DevBrainAge_predage', 'DBN_predage', 'DunedinPACNI', 'ENIGMA_predage', 'Kaufmann_predage', 'PyBrainAge_predage', 'Pyment_predage'
)

# ---- Configuration per timepoint ----

# here we define paths to the data from different time points 
# and provide the time point names (e.g., ALSPAC has one brain timepoint = 15up) 
# we also define which brain variables to exclude from the winsorisation
# for brain ages, typically none should be excluded 
# the below example is for the ALSPAC brain data for 1 timepoint
# if you have more timepoints, add more lists into timepoint_configs

brain_timepoint_configs <- list(
  `15up` = list(
    file = "data/brain_ages/brain_ages_covs_15up.csv",
    exclude_vars = c()  # Add exclusions here if needed later
  )
)

# ---- Run cleaning & descriptives function ----

# the function prep_and_descr_timepoint_brain has to be run per timepoint 
  # it will do the following:
  # rename variables for prettier plotting 
  # set variable types 
  # provide descriptives pre+post winsorisation + log file (saved under /results/descriptives) 
  #winsorise epi ages 

# !!!! IMPORTANT REQUIRED CHANGES !!!!!
  # the script assumes the covariates and brain age model names from the instructions
  # document. Expected categorical variables (i.e., SEX_brain) and potentially included
  # variables ('study', 'scanner' and 'coil') are set to factor variables BUT any 
  # additional covariates are set to numeric
  # IF you have included additional covariates (after discussion with 
  # Marlene) that are not numeric, please open the function
  # prep_and_descr_timepoint_brain and make changes to this line:
  # cat_candidates <- c("study", "scanner", "coil") by adding additional variable 
  # names to the vector 

# Warning messages may appear if some columns are completely NA 
# Warning messages:
# 1: In min(x, na.rm = TRUE) :
# no non-missing arguments to min; returning Inf
# this happens when using summarytools if some of your variables are
# completely NA --> can usually be ignored

prep_and_descr_timepoint_brain("15up", common_brain_vars, brain_timepoint_configs)

# _______________________________________________#
# ---- TIMEPOINTS WITH EPI / BRAIN OVERLAP ----
# _______________________________________________#

# we will now run the descriptives again on merged data for timepoint with 
# both epi and brain data using the function descr_overlap 

# 1) Define  timepoints where both epi and brain ages are available
      # For ALSPAC, this is only 15up 
joint_timepoints <- c("15up") 

# 2) Run joint descriptives for each timepoint with overlap
for (tp in joint_timepoints) {
  descr_overlap(tp)
}

#### DONE WITH 1.epi_brain_clean_and_describe.R ####