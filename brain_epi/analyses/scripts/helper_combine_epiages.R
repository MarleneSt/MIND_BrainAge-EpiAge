#-------------------------------------------------------------------------#
#### load, merge & manipulate different epi ages  -----
#
# this is a helper script we used in the ALSPAC data to compile
# epi ages and covariates per timepoint and might require editing
# 
# has to be applied per timepoint
#
# !the examples paths provided should reflect the location of files 
# !based on the provided instructions 
#-------------------------------------------------------------------------#

#### setup ----
# set working directory to base folder with epi ages e.g., /brain_epi/epi/models/, this is where the results from the epi age scripts are saved, most importantly the all_predictions files per timepoint
setwd("path/to/models/folder")

# library
library(dplyr)
library(readxl) #to read in excel spreadsheet

#save a renaming map for later remapping (left = old, right = new) 
rename_map <- c(
  "ShirebyG2020"    = "CorticalClock",
  "BernabeuE2023c"  = "cAge",
  "altumage"        = "AltumAge",
  "horvath2013"     = "Horvath2013",
  "yingdamage"      = "DamAge",
  "yingadaptage"    = "AdaptAge",
  "dunedinpace"     = "DunedinPACE",
  "Levine"          = "PhenoAge",
  "TL"              = "DNAmTL",
  "BLUP"            = "ZhangBLUP",
  "EN"              = "ZhangEN",
  "Bohlin"          = "Bohlin"
)

#### birth timepoint ----

# example is 450k

##### load epi data for birth timepoint -----

dat_epi_cord_nongest <- read.csv("all_predictions_450k_cord.csv")

# gestational clocks are saved in a different file, 'cord' in the filename might have to be subsituted with the name you gave to this timepoints
dat_epi_cord_gest <- read_excel("../450k_cord/epi_ages_450k_cord_pred.xlsx", sheet = 2) # add in the gestational clocls 

names(dat_epi_cord_gest)[names(dat_epi_cord_gest) == "id"] <- "SUBJID" #or Sample_Name or whatever ID name is in your all_predictions file

#merge in the gestational clocks we're including 
dat_epi_cord <- merge(dat_epi_cord_nongest, dat_epi_cord_gest[, c("Bohlin","EPIC", "Knight", "SUBJID")], by="SUBJID" )

##### data cleaning and renaming ----

dat_epi_cord$Horvath <- NULL #take out Horvath2013 clock without BMIQ (we're keeping the BMIQ version)

colnames(dat_epi_cord) <- gsub("_mAge", "", colnames(dat_epi_cord)) #remove _mAge to shorten names

# rename variables for prettier plotting (only applies it to variables that exist)
existing_vars <- intersect(names(rename_map), names(dat_epi_cord))

if (length(existing_vars) > 0) {
  rename_map_subset <- rename_map[existing_vars]
  dat_epi_cord <- dat_epi_cord %>% dplyr::rename(!!!setNames(existing_vars, rename_map_subset))
}

cols_to_rename <- colnames(dat_epi_cord) # Get column names
cols_to_rename <- cols_to_rename[!(cols_to_rename %in% c("SUBJID", "DunedinPACE"))] # Identify columns to rename: exclude SUBJID and dunedinepace 
colnames(dat_epi_cord)[colnames(dat_epi_cord) %in% cols_to_rename] <- 
  paste0(cols_to_rename, "_predage") # Rename selected columns by appending "_predage"

#convert gestational weeks to years
dat_epi_cord$Bohlin_conv_predage <- 7 * (dat_epi_cord$Bohlin_predage  - 40) / 365
dat_epi_cord$EPIC_conv_predage <- 7 * (dat_epi_cord$EPIC_predage  - 40) / 365
dat_epi_cord$Knight_conv_predage <- 7 * (dat_epi_cord$Knight_predage  - 40) / 365

#### other non-birth timepoints ----

##### load epi data for other timepoints ----
dat_epi_F7 <- read.csv("all_predictions_450k_F7.csv")

# renaming and cleaning 

dat_epi_F7$Horvath <- NULL #take out Horvath2013 clock without BMIQ (we're keeping the BMIQ version)

colnames(dat_epi_F7) <- gsub("_mAge", "", colnames(dat_epi_F7)) #remove _mAge to shorten names

# rename variables for prettier plotting (only applies it to variables that exist)
existing_vars <- intersect(names(rename_map), names(dat_epi_F7))

if (length(existing_vars) > 0) {
  rename_map_subset <- rename_map[existing_vars]
  dat_epi_F7 <- dat_epi_F7 %>% dplyr::rename(!!!setNames(existing_vars, rename_map_subset))
}

cols_to_rename <- colnames(dat_epi_F7) # Get column names
cols_to_rename <- cols_to_rename[!(cols_to_rename %in% c("SUBJID", "DunedinPace"))] # Identify columns to rename: exclude SUBJID and dunedinepace 
colnames(dat_epi_F7)[colnames(dat_epi_F7) %in% cols_to_rename] <- 
  paste0(cols_to_rename, "_predage") # Rename selected columns by appending "_predage"


#### load covariate spreadsheets ----

covs_cord <- read.csv("/path/to/covariates/covariates_cord_epi.csv")

covs_F7 <- read.csv("/path/to/covariates/covariates_F7_epi.csv")

# merge 

dat_epi_cord_covs <- merge(dat_epi_cord, covs_cord, by = "SUBJID")

dat_epi_F7_covs <- merge(dat_epi_F7, covs_F7, by = "SUBJID")

#### perform some checks ----

# duplication of SUBJIDs?
sum(duplicated(dat_epi_cord_covs$SubjID)) #0
sum(duplicated(dat_epi_F7_covs$SubjID)) #0

# correct columns included? 
colnames(dat_epi_cord_covs)
colnames(dat_epi_F7_covs)

# for non-birth timepoints, columns include 
# "SUBJID"
# "SEX_epi"
# "AGE_epi"
# "AGE_gest"
# "PCBrainAge_predage"
# "PCGrimAge_predage"
# "CorticalClock_predage"
# "cAge_predage"
# "AltumAge_predage"
# "Horvath2013_predage"
# "DamAge_predage"
# "AdaptAge_predage"
# "DunedinPACE"
# "Hannum_predage"
# "PhenoAge_predage"
# "skinHorvath_predage"
# "PedBE_predage"
# "Wu_predage"
# "DNAmTL_predage"
# "ZhangBLUP_predage"
# "ZhangEN_predage"
# "Bcell"
# "CD4T"
# "CD8T"
# "Gran"
# "Mono"
# "NK"
# "V1"
# "V4"
# "V5"
# "V7"
# "V8"
# "V9"

# optional step: align order in line with the template (first main covs, then epi ages, then cell type, then batch) but not essential

# saving files 
write.csv(dat_epi_cord_covs, "/brain_epi/analyses/data/epi_ages/epi_ages_covs_cord.csv", row.names = FALSE)

write.csv(dat_epi_F7_covs, "/brain_epi/analyses/data/epi_ages/epi_ages_covs_F7.csv", row.names = FALSE)
