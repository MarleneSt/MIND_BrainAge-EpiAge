#_____________________________________________________________________________#
######## Script to apply imputation to different input feature sets   ##########
#_____________________________________________________________________________#
## created by Marlene Staginnus 

#### set WD ----

setwd("/brain_epi/mri/models/kaufmann/") # path where  input features are saved (or subfolders if working with multiple TP)

#### load packages ----
load.lib <- function(x) {
  for (i in x) {
    if (!require(i, character.only = TRUE)) {
      install.packages(i, dependencies = TRUE)
      library(i, character.only = TRUE)
    }
  }
}
load.lib(c("skimr", "psych", "missForest")) 

#### load function ----
source("/brain_epi/mri/scripts_models/imputation_helpers.R")

#### load data (covariates & imaging features) ----
covariates <- read.csv("../../Covariates.csv")
cerebellum_sub <- read.csv("glasser_cerebellum_subcortical_enigma_aligned.csv", check.names = FALSE)
CT <- read.csv("glasser_ThickAvg.csv", check.names = FALSE)
SA <- read.csv("glasser_SurfAvg.csv", check.names = FALSE)
Vol <- read.csv("glasser_VolAvg.csv", check.names = FALSE)

# !!!!!IMPORTANT, if you observed issues with read.csv whereby the SUBJID column is mistakenly read in as the row header, you can try using readr::read_csv or resaving the input csv files

#### combine dataframes (order is CT, Vol, SA, rest) ----

# in this step subjects would be dropped that are entirely missing from one of the dataframes
# this is because we are expecting region-based missingness, not missingness of 
# a full outcome (e.g., thickness) - if this is the case for you, please get in touch 
# so that we can discuss the appropriateness of imputation

dat <- merge(CT, Vol, by ="SUBJID")
dat <- merge(dat, SA, by ="SUBJID")
dat <- merge(dat, cerebellum_sub, by ="SUBJID")

covs <- covariates[, c("SUBJID", "SEX_brain", "AGE_brain")]
names(covs) <- toupper(names(covs))  # Capitalize all for standardization

dat <- merge(dat, covs, by="SUBJID")

df_model <- dat 

#### generate descriptives & impute ----

# kaufmann features have to be read in due to the many features/number of characters
kaufmann_features <- scan("/brain_mri/mri/scripts_models/kaufmann/kaufmann_features.txt", what = character(), quiet = TRUE)

imputation_obj <- impute_brain_df(
  df         = df_model,
  id_col     = "SUBJID",
  sex_col    = "SEX_BRAIN",
  age_col    = "AGE_BRAIN",
  feature_cols = kaufmann_features,  #possible to feed in vector directly but kaufmann has too many features 
  seed       = 21021995,
  out_dir    = "imputed", #will create subdirectory for imputed values 
  out_prefix = "kaufmann", # files saved under imputed/kaufmann_...
  log_path   = "imputed/kaufmann_imputation.log",
  save_imputed = TRUE
)
