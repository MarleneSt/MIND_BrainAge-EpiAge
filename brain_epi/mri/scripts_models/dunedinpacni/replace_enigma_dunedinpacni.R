#_____________________________________________________________________________#
######## Script to replace PyBrainAge ROIs with QC'ed ENIGMA ROIS   ##########
#_____________________________________________________________________________#
## created by Marlene Staginnus 

#### setwd ----
setwd("/brain_epi/mri/models/dunedinpacni/") # path where dunedinpacni input features are saved

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
dunedinpacni <- read.csv("dunedinpacni_features.csv", check.names = FALSE) #check.names = FALSE as some column names for some models are hyphenated
enigma_CT <- read.csv("../../extracted_features/CorticalMeasuresENIGMA_ThickAvg.csv", check.names = FALSE)
enigma_SA <- read.csv("../../extracted_features/CorticalMeasuresENIGMA_SurfAvg.csv", check.names = FALSE)
enigma_sub <- read.csv("../../extracted_features/LandRvolumes.csv", check.names = FALSE)

# merge enigma files into one 

enigma_all <- merge(enigma_CT,enigma_SA, by = "SUBJID" )
enigma_all <- merge(enigma_all,enigma_sub, by = "SUBJID" )

#### source function ----
source("/brain_epi/mri/scripts_models/replace_enigma_helper.R") #adjust to path to helper script if needed

#### map dunedinpacni age - enigma_sub ----

# due to the man rows we are reading this one in 
map <- read.csv("/brain_epi/mri/scripts_models/scripts_models/dunedinpacni/dunedinpacni_feature_map_enigma.csv") %>%
  tibble::as_tibble()
    # note: cortical volume & grey/white ratio regions get mapped to the SA NAs 

#### replace NAs ----

out <- replace_from_enigma(
  target_df     = dunedinpacni,   
  enigma_df     = enigma_all, 
  id_col_target = "ID",
  id_col_enigma = "SUBJID",
  mapping       = map,
  mode          = "overwrite_if_enigma_na",
  output_path   = "dunedinpacni_enigma_aligned.csv"
)
