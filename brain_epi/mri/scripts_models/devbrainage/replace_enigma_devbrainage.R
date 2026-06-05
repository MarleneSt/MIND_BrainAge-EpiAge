#_____________________________________________________________________________#
######## Script to replace PyBrainAge ROIs with QC'ed ENIGMA ROIS   ##########
#_____________________________________________________________________________#
## created by Marlene Staginnus 

#### setwd ----
setwd("/brain_epi/mri/models/devbrainage") # path where pybrainage input features are saved

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
devbrainage <- read.csv("DevelopmentalBrainAge_input_withIDs.csv", check.names = FALSE) #check.names = FALSE as some column names are hyphenated
enigma_SA <- read.csv("../../extracted_features/CorticalMeasuresENIGMA_SurfAvg.csv", check.names = FALSE)
enigma_sub <- read.csv("../../extracted_features/LandRvolumes.csv", check.names = FALSE)

# merge enigma files into one 
enigma_all <- merge(enigma_SA,enigma_sub, by = "SUBJID" )

#### source function ----
source("/brain_epi/mri/scripts_models/replace_enigma_helper.R")

#### map devbrainage age - enigma_sub ----

# due to the man rows we are reading this one in 
map <- read.csv("/brain_epi/mri/scripts_models/devbrainage/devbrainage_feature_map_enigma.csv") %>%
  tibble::as_tibble()
# note: cortical volume & grey/white ratio regions get mapped to the SA NAs

#### replace NAs ----

out <- replace_from_enigma(
  target_df     = devbrainage,   
  enigma_df     = enigma_all, 
  id_col_target = "SUBJID",
  id_col_enigma = "SUBJID",
  mapping       = map,
  mode          = "overwrite_if_enigma_na",
  output_path   = "devbrainage_enigma_aligned_withIDs.csv"
)

devbrainage.no.id <- out[-1]
write.csv(devbrainage.no.id, "devbrainage_enigma_aligned.csv", row.names = FALSE)

