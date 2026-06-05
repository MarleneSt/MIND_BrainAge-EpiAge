
## BMIQ Calibration

# This script applies bmiq callibration in R (needed for subsequent AltumAge and Horvath epi clock predictions)
# https://github.com/perishky/meffonym/blob/master/vignettes/age-tutorial.md

# set working directory to the "brain_epi/epi" folder we provided 
wd = "/campaign/VB-FM5HPC-001/Vilte/Projects/Marlene-epi-clocks/brain_epi/epi" # (! update with your path to brain_epi/epi/ folder)
setwd(wd)

# set output path to bmiq folder (to store bmiq calibrated data)
output <- file.path(wd, "postprocessed/imputed/bmiq/")
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
postprocess_path = file.path(wd, "postprocessed/imputed/") # path to 'postprocessed/imputed' folder (where imputed and winsorized DNAm .rds files are stored)

# library(remotes)
# install_github("perishky/meffonym")
library(meffonym)

# define the gold standard
standard <- meffonym.horvath.standard()

# function to bmiq normalise and save methylation data (separately per timepoint and array)
process_dnam <- function(file_name, output_folder) {

  # read data
  message("Processing file: ", file_name)
  meth <- readRDS(paste0(postprocess_path, file_name))
  
  # convert to data frame if needed
  if (!is.matrix(meth)) {
    message("Converting data frame to matrix")
    meth <- as.matrix(meth)
  }
  
  # normalize using BMIQ calibration
  dnam.norm <- meffonym.bmiq.calibration(meth, standard)
  
  # save outputs
  base_name <- sub("\\.rds$", "", basename(file_name))
  message("Saving: ", paste0(base_name, ".bmiq.csv"))
  write.csv(dnam.norm, file = paste0(output_folder, base_name, ".bmiq.csv"))
  
  message("Finished processing: ", file_name)
}

# normalise each dataset and save output as .csv file
# saves output as file_name.bmiq.csv (e.g. "DNAm.imputed_wins_cord450k.bmiq.csv")
process_dnam(file_name="DNAm.imputed_wins_cord_450k.rds", output_folder=output) 
process_dnam(file_name="DNAm.imputed_wins_F7_450k.rds", output_folder=output) 
process_dnam(file_name="DNAm.imputed_wins_up15_450k.rds", output_folder=output) 
process_dnam(file_name="DNAm.imputed_wins_up15_EPIC.rds", output_folder=output)


