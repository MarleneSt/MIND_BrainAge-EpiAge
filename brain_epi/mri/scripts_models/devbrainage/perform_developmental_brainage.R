###############################################################################
######## Script to apply the Developmental Brain Age model to new data ########
###############################################################################

###### created by Marlene Staginnus based on https://github.com/GitDro/DevelopmentalBrainAge 

cat("Loading tidyverse, tidymodels, remotes and xgboost v1.0.0.2\nIf there are issues with this script, try installing these before running\n\n")

#### library ----

# Set CRAN mirror explicitly
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Function to install and load a package
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

install_if_missing("tidyverse")
install_if_missing("tidymodels")
install_if_missing("remotes")

# Install specific xgboost version if not already installed
if (!("xgboost" %in% installed.packages()[, "Package"]) || packageVersion("xgboost") != "1.0.0.2") {
  remotes::install_version("xgboost", version = "1.0.0.2", upgrade = "never")
}
library(xgboost)

cat("Loading data from current folder:\n
    - DevelopmentalBrainAge_input.csv\n
    - DevelopmentalBrainAge_input_withIDs.csv\n\n")

#### load data ----
dat <- read.csv("DevelopmentalBrainAge_input.csv")
dat_ID <- read.csv("DevelopmentalBrainAge_input_withIDs.csv")
your_df <- dat

cat("Loading model xgboost_9to19_brain_age_mod.rds \n\n")

#### load model ----
xgb_mod <- readRDS(
  file = file.path(
    "/brain_epi/mri/scripts_models/devbrainage/xgboost_9to19_brain_age_mod.rds" #this path has to be adjusted
  )
)

cat("Run model \n\n")

#### run & save model ----
test_out <- predict(xgb_mod$fit, newdata = as.matrix(your_df))

#### save variable ----
dat$DevBrainAge_predage <- test_out

dat$SUBJID <- dat_ID$SUBJID
dat <- dat[, c(190, 1:189)]
dat_nofeatures <- dat[, c(1,190)]

cat("Save model results with and without features\n
    -DevelopmentalBrainAge_output_withfeatures.csv\n
    -DevelopmentalBrainAge_output.csv\n\n")

# save file 
write.csv(dat, file="DevelopmentalBrainAge_output_withfeatures.csv",row.names = FALSE)
write.csv(dat_nofeatures, file="DevelopmentalBrainAge_output.csv",row.names = FALSE)

