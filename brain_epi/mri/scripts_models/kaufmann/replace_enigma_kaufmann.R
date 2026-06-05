#_____________________________________________________________________________#
######## Script to replace Glasser ROIs with QC'ed ENIGMA ROIS   ##########
#_____________________________________________________________________________#
## created by Marlene Staginnus 

#### setwd ----
setwd("/brain_epi/mri/models/kaufmann") # path where kaufmann input features are saved

#### load packages ----

packages <- c("dplyr")

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

#### load data ----
glasser_sub <- read.csv("glasser_cerebellum_subcortical.csv", check.names = FALSE) #check.names = FALSE as for some models columns include hyphenated names
enigma_sub <- read.csv("../../extracted_features/LandRvolumes.csv", check.names = FALSE)

#### source function ----
source("/brain_mri/mri/scripts_models/replace_enigma_helper.R") #adjust to correct script location

#### map glasser_sub - enigma_sub ----
map <- tibble::tribble(
  ~target,                       ~enigma,
  "Left.Lateral.Ventricle",       "LLatVent",
  "Left.Thalamus.Proper",         "Lthal",
  "Left.Caudate",                 "Lcaud",
  "Left.Putamen",                 "Lput",
  "Left.Pallidum",                "Lpal",
  "Left.Hippocampus",             "Lhippo",
  "Left.Amygdala",                "Lamyg",
  "Left.Accumbens.area",          "Laccumb",
  "Right.Lateral.Ventricle",      "RLatVent",
  "Right.Thalamus.Proper",        "Rthal",
  "Right.Caudate",                "Rcaud",
  "Right.Putamen",                "Rput",
  "Right.Pallidum",               "Rpal",
  "Right.Hippocampus",            "Rhippo",
  "Right.Amygdala",               "Ramyg",
  "Right.Accumbens.area",         "Raccumb",
  "EstimatedTotalIntraCranialVol","ICV"
)

#### replace NAs ----

out <- replace_from_enigma(
  target_df     = glasser_sub,   
  enigma_df     = enigma_sub, 
  id_col_target = "SUBJID",
  id_col_enigma = "SUBJID",
  mapping       = map,
  mode          = "overwrite_if_enigma_na",
  output_path   = "glasser_cerebellum_subcortical_enigma_aligned.csv"
)
