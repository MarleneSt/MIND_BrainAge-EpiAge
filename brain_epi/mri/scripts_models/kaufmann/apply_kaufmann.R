###############################################################################
########### Script to apply the Kaufmann brain age model to new data ##########
###############################################################################

###### created by Marlene Staginnus based on information 
###### from https://github.com/tobias-kaufmann/brainage
###### Please first follow the instructions to generate the glasser parcellation

#### setting paths  ----
setwd("/brain_mri/mri/models/kaufmann/") # path where input features of the brain age model are stored (or subfolders if working with multiple TP)

script_path <- "/brain_epi/mri/scipts_models/kaufmann/" # path to scripts

#### load library ----
load.lib <- function(x) {
  for (i in x) {
    if (!require(i, character.only = TRUE)) {
      install.packages(i, dependencies = TRUE)
      library(i, character.only = TRUE)
    }
  }
}

load.lib(c("xgboost")) 

#### load data (covariates & imaging features) ----
covariates <- read.csv("../../Covariates.csv")

imaging <- read.csv("imputed/kaufmann_features_imputed.csv") 
# imaging <- read.csv("kaufmann_features.csv")  # use this if no data was imputed 

#### combine imaging & covs  ----
covs <- covariates[, c("SUBJID", "SEX_brain", "AGE_brain")]

dat <- merge(imaging, covs, by="SUBJID")

#### separate females and males ----
female <- dat[ which(dat$SEX_brain=='1'), ]
male <- dat[ which(dat$SEX_brain=='0'), ]

features_females <- female
features_females[ , c("SUBJID","SEX_brain","AGE_brain")] <- list(NULL)

features_males <- male
features_males[ , c("SUBJID","SEX_brain","AGE_brain")] <- list(NULL)

#### run model ----

##### load functions -----
load(file.path(script_path,"github-brainage-master/brainageModels.RData"))

##### apply kaufmann model -----
brainAge_females <- predict(mdl_agepred_female, as.matrix(features_females))
brainAge_males <- predict(mdl_agepred_male, as.matrix(features_males))

female$Kaufmann_predage <- brainAge_females
male$Kaufmann_predage <- brainAge_males

# retain only SUBJID & brainage values for saving 
female_save <- female[,c("SUBJID","Kaufmann_predage")]
male_save <- male[,c("SUBJID","Kaufmann_predage")]

# make one DF with both sexes 
combined <- rbind(female_save, male_save)

# save files 
write.csv(female_save, "kaufmann_female_brainage.csv", row.names = FALSE)
write.csv(female_save, "kaufmann_male_brainage.csv", row.names = FALSE)
write.csv(combined, "kaufmann_brainage.csv", row.names = FALSE)