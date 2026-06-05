#-------------------------------------------------------------------------#
#### load, merge & manipulate different brain ages  -----
#
# this is a helper script we used in the ALSPAC data to compile
# brain ages and covariates per timepoint and might require editing
# 
# has to be applied per timepoint
#
# !the examples paths provided should reflect the location of files 
# !based on the provided instructions 
#-------------------------------------------------------------------------#

#### set paths ----

# set working directory to base folder with brain age folders e.g., /brain_epi/mri/models/
setwd("/brain_epi/mri/models/")

# read in covariates for first timepoint
covs <- read.csv("path/to/covariates") #should include SUBJID, SEX_brain, AGE_brain, AGE_gest and additional technical variables if needed

#### Centile ----
centile_male <- read.csv("centile/centile_male.csv") 
centile_female <- read.csv("centile/centile_female.csv")

centile <- rbind(centile_male, centile_female)

str(centile) #confirm 1st column = 'SUBJID' & 2nd column = 'Centile2_predage' 

#### DevBrainAge ----
devbrainage <- read.csv("devbrainage/DevelopmentalBrainAge_output.csv")

str(devbrainage) #confirm 1st column = 'SUBJID' & 2nd column = 'DevBrainAge_predage'

#### DunedinPACNI ----
dpacni <- read.csv("dunedinpacni/dunedinpacni_brainage.csv") #or dunedinpacni/imputed/dunedinpacni_brainage.csv if imputation was applied 

str(dpacni) #confirm 1st column = 'SUBJID' & 2nd column = 'DunedinPACNI'

#### DBN ----
dbn <- read.csv("DBN/dbn_predictions.csv")

str(dbn) #confirm 1st column = 'SUBJID' & 2nd column = 'DBN_predage'

#### ENIGMA ----
enigma <- read.csv("ENIGMA/all_ENIGMA_output.csv")

str(enigma) #confirm 1st column = 'SUBJID' & 2nd column = 'ENIGMA_predage'

#### Kaufmann ----
kaufmann <- read.csv("kaufmann/kaufmann_brainage.csv") #or kaufmann/imputed/kaufmann_brainage.csv if imputation was applied 

str(kaufmann) #confirm 1st column = 'SUBJID' & 2nd column = 'Kaufmann_predage'

#### PyBrainAge ----
pybrainage <- read.csv("pybrainage/PyBrainAge_Output.csv") #or pybrainage/imputed/PyBrainAge_Output.csv if imputation was applied 

# rename vars
names(pybrainage)[names(pybrainage) == "ID"] <- "SUBJID"
names(pybrainage)[names(pybrainage) == "BrainAge"] <- "PyBrainAge_predage"

# remove all vars except the renamed ones 
pybrainage <- pybrainage[, c("SUBJID", "PyBrainAge_predage")]

str(pybrainage) #confirm 1st column = 'SUBJID' & 2nd column = 'PyBrainAge_predage'

#### Pyment ----
pyment <- read.csv("pyment/pyment_predictions.csv")

str(pyment) #confirm 1st column = 'SUBJID' & 2nd column = 'Pyment_predage'

#### merge all brain ages ----

# list of data frames to merge
merge_list <- list(
  centile[, c("SUBJID", "Centile2_predage")],
  devbrainage[, c("SUBJID", "DevBrainAge_predage")],
  dbn[, c("SUBJID", "DBN_predage")],
  dpacni[, c("SUBJID", "DunedinPACNI")],
  enigma[, c("SUBJID", "ENIGMA_predage")],
  kaufmann[, c("SUBJID", "Kaufmann_predage")],
  pybrainage[, c("SUBJID", "PyBrainAge_predage")],
  pyment[, c("SUBJID", "Pyment_predage")]
)

dat_merged <- Reduce(function(x, y) merge(x, y, by = "SUBJID", all.x = TRUE), merge_list, init = covs)

#### some recommended checks ----

str(dat_merged)
colnames(dat_merged)

# the following columns should be included, please also include a brain age model 
# if it was not included in your sample but set all values to NA

# "SUBJID"
# "SEX_brain"
# "AGE_brain"
# "AGE_gest"
# "Centile2_predage"
# "DevBrainAge_predage"
# "DBN_predage"
# "DunedinPACNI"
# "ENIGMA_predage"
# "Kaufmann_predage"
# "PyBrainAge_predage"
# "Pyment_predage"
# "scanner" (only if applicable i.e., multiple scanners used)     

# how many Subject IDs are duplicated
sum(duplicated(dat_merged$SUBJID))

# who is duplicated 
dat_merged[duplicated(dat_merged$SUBJID),][, "SUBJID"]

# remove duplicates if needed

#### save combined brain age & covariate file

# example for one timepoint
write.csv(dat_merged, "brain_epi/analyses/data/brain_ages/brain_ages_covs_<timepoint>.csv", row.names = FALSE)

#### repeat for other timepoints ----
#...