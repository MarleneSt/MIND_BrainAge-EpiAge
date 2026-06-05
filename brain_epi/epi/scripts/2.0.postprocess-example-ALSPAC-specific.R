## This script is specific to ALSPAC 
# it selects the relevant sample (i.e., removes twins, siblings, duplicates) 
# and saves a separate samplesheet and methylation matrix per time point and array
# for non-ALSPAC cohorts, do this step independently 

# _______________________________________________#
#                       SET UP
# _______________________________________________#
## load packages 
#library(devtools)
# install_github("MRCIEU/aries")
library(aries)
library(foreign)
library(tidyverse)

## set paths
setwd("/campaign/VB-FM5HPC-001/Vilte/Projects/Marlene-epi-clocks") 
aries_dir ="/campaign/VB-FM5HPC-001/Vilte/Projects/BioMM/woc_release/" # alspac-specific
scripts = "/campaign/VB-FM5HPC-001/Vilte/Projects/BioMM/woc_release/" # alspac-specific
output = "/campaign/VB-FM5HPC-001/Vilte/Projects/Marlene-epi-clocks/" # specify path to your output folder
pheno_path = "/campaign/VB-FM5HPC-001/Vilte/Data/ALSPAC-pheno/"

## source aries function
source(paste0(scripts, 'aries.function_3067.R')) # alspac-specific

# _______________________________________________#
#       READ IN AND SAVE METH SAMPLESHEETS
# _______________________________________________#
# first read in pheno data and create unique ID column 
df_pheno <- read.spss(paste0(pheno_path, "CIDB3067_26Aug2023.sav"), to.data.frame = TRUE)
df_pheno$cidB3067.qlet <- paste0(df_pheno$cidB3067, df_pheno$qlet)
head(df_pheno$cidB3067.qlet)

## then, read in and save meth samplsheets per each time point and array 
# first define function 
process_aries_data <- function(aries_dir, time_points = NULL, feature_sets = NULL) {
  
  # Read the data for given feature set and time point
  aries <- aries.select.mod(aries_dir, time.point = time_points, featureset = feature_sets)
  
  # samplesheet
  samplesheet = aries$samples
  samplesheet$cidB3067.qlet=paste0(samplesheet$cidB3067, samplesheet$QLET)
  
  # cell counts
  cell_counts = as.data.frame(aries$cell.counts)
  
  # add Sample_Name column
  cell_counts$Sample_Name <- rownames(cell_counts)
  
  # merge with samplesheet
  samplesheet <- merge(samplesheet, cell_counts, by = "Sample_Name")
  message("samplesheet rows: ", nrow(samplesheet))
  
  # merge in sibling column
  samplesheet <- merge(samplesheet, df_pheno[, c('cidB3067.qlet', 'mz005l')], by = "cidB3067.qlet")
  message("samplesheet rows after merging in 'mz005l' variable: ", nrow(samplesheet))
  
  return(list(samplesheet = samplesheet, aries = aries))
  
}

# 450k (not avail for F24)
samplesheet450k_cord <- process_aries_data(aries_dir, time_points = 'cord', feature_sets = '450')
samplesheet450k_F7 <- process_aries_data(aries_dir, time_points = 'F7', feature_sets = '450')
samplesheet450k_15up <- process_aries_data(aries_dir, time_points = '15up', feature_sets = '450')

# epic (not avail for cord or F7)
samplesheetEpic_15up <- process_aries_data(aries_dir, time_points = '15up', feature_sets = 'epic')
samplesheetEpic_F24 <- process_aries_data(aries_dir, time_points = 'F24', feature_sets = 'epic')

# _______________________________________________#
#      REMOVE SIBLINGS AND SAVE METH DATA
# _______________________________________________#
# combine all time points and arrays into a single results_list
results_list <- list(
  cord_450k  = samplesheet450k_cord,
  F7_450k    = samplesheet450k_F7,
  up15_450k  = samplesheet450k_15up,
  up15_epic  = samplesheetEpic_15up
)

# open connection
log_con <- file("alspac_processing_log.txt", open="wt")
sink(log_con, split=TRUE)
sink(log_con, type="message")

for (name in names(results_list)) {
  
  aries <- results_list[[name]]$aries
  samplesheet <- results_list[[name]]$samplesheet
  message("Processing: ", name)
  
  # extract methylation data
  aries$meth <- aries.methylation(aries)
  
  # STEP 1: keep only QLET == "A"
  samplesheet_unique <- samplesheet %>% filter(QLET == "A")
  message("Twins removed: ", nrow(samplesheet) - nrow(samplesheet_unique))
  
  # remove duplicates
  samplesheet_unique2 <- samplesheet_unique %>% filter(duplicate.rm == FALSE)
  message("Duplicates removed: ", nrow(samplesheet_unique) - nrow(samplesheet_unique2))
  
  # check sibling variable
  message("Sibling status:\n", capture.output(print(table(samplesheet_unique2$mz005l))))
  
  # remove siblings
  samplesheet_unique3 <- samplesheet_unique2[grepl("No, keep all", samplesheet_unique2$mz005l),]
  message("Siblings and NA removed: ", nrow(samplesheet_unique2) - nrow(samplesheet_unique3))
  
  keep_samples <- samplesheet_unique3$Sample_Name
  aries$meth <- aries$meth[, keep_samples, drop=FALSE]
  stopifnot(all(samplesheet_unique3$Sample_Name == colnames(aries$meth)))
  
  # write cleaned samplesheet to CSV
  write.csv(samplesheet_unique3, file = paste0(output, "samplesheet_", name, ".csv"), row.names = FALSE)
  message("Saved samplesheet CSV: ", paste0(output, "samplesheet_", name, ".csv"))
  
  # wrire aries meth to .rds file
  saveRDS(aries$meth, file=paste0(output, "meth_", name, ".rds"))
  message("Saved methylation data to: ", paste0(output, "meth_", name, ".rds"))
}

sink(type="message")
sink()
close(log_con)
# takes some time to run as saves DNAm for each time point and array
# after running, check log output in "alspac_processing_log.txt" 
# log file also contains file names (for samplesheets and DNAm files) to be used as input in the next script

## additional step for cord blood samplesheet
# convert bestgest to years format using this formula:
# years = 7 * (gestational weeks – 40)/365 
cord_samplesheet <- read.csv("samplesheet_cord_450k.csv")

# read in pheno data to extract bestgest 
df_pheno <- read.spss(paste0(pheno_path, "CIDB3067_26Aug2023.sav"), to.data.frame = TRUE)
df_pheno$cidB3067.qlet <- paste0(df_pheno$cidB3067, df_pheno$qlet)
head(df_pheno$cidB3067.qlet)

# merge in bestgest
cord_samplesheet <- merge(cord_samplesheet, df_pheno[, c("cidB3067.qlet", "bestgest")])
cord_samplesheet$bestgest <- as.numeric(as.character(cord_samplesheet$bestgest))
cord_samplesheet <- cord_samplesheet %>% mutate(age = (7 * (bestgest - 40)) / 365)
write.csv(cord_samplesheet, file = paste0(output, "samplesheet_cord_450k.csv"), row.names = FALSE)



