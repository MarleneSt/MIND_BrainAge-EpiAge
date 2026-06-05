# !! before running this script, make sure all output is in the same folder 
## (e.g. if you kept our folder structure, all output should already be stored in 'models' folder)

# if you have more than one time point, add subfolders within the output folder per time point and array
# e.g., "450k_cord", "450k_F7", "450k_15up", "epic_15up" 
# cd models # enter output folder 
# mkdir 450k_cord 450k_F7 450k_15up epic_15up # create 3 subfolders for the 3 timepoints (please adjust the number of folders and their names to correspond to your data)
# mv *cord* 450k_cord/ # move any output that contains 'cord' to '450k_cord'
# mv *F7* 450k_F7/ # move any output that contains 'F7' to '450k_F7'
# mv *epic_15up* epic_15up/ # same as above
# mv *450k_15up* 450k_15up/ # same as above

# load packages
library(readr)
library(openxlsx)
library(dplyr)
library(stringr)
library(purrr)

# set working directory to the "brain_epi/epi" folder we provided 
wd = "/campaign/VB-FM5HPC-001/Vilte/Projects/Marlene-epi-clocks/brain_epi/epi" # (! update with your path to brain_epi/epi/ folder)
setwd(wd)

# !! define your path to output folder (called 'models') that contains all epi age predictions
base_path <- file.path(wd, "models")
# !! define your folder names for each time point (assumes you have placed your epi predictions in separate folder per time point)
folders <- c("450k_cord", "450k_F7", "450k_15up", "epic_15up")

# save merged_data to the same folder
output_folder = file.path(wd, "models")

# helper function to standardize sample column name
standardize_sample_col <- function(df) {
  sample_cols <- c("Sample_Name", "id", "Sample")
  present <- intersect(sample_cols, names(df))
  if (length(present) == 0) stop("No sample column found!")
  if (length(present) > 1) stop("Multiple sample columns found!")
  df <- rename(df, Sample_Name = all_of(present))
  return(df)
}

# excel reading 
try_read_excel <- function(file) {
  tryCatch({
    df <- openxlsx::read.xlsx(file)
    message("Success reading with openxlsx: ", file)
    return(df)
  }, error = function(e) {
    message("Failed reading with openxlsx: ", file)
    return(NULL)
  })
}

# auto read file based on extension
read_file_auto <- function(file_path) {
  if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
    df <- readr::read_csv(file_path, show_col_types = FALSE)
    message("Success reading with readr: ", file_path)
  } else if (grepl("\\.xlsx$", file_path, ignore.case = TRUE)) {
    df <- try_read_excel(file_path)
    if (is.null(df)) stop(paste("Failed reading Excel file:", file_path))
  } else {
    stop(paste("Unsupported file type:", file_path))
  }
  df <- standardize_sample_col(df)
  return(df)
}

# main processing function for each folder
process_folder <- function(folder, base_path) {
  folder_path <- file.path(base_path, folder)
  files <- list.files(folder_path, full.names = TRUE, pattern = "\\.(csv|xlsx)$", ignore.case = TRUE)
  df_list <- map(files, read_file_auto)
  names(df_list) <- basename(files) %>% tools::file_path_sans_ext()
  
  merged <- reduce(df_list, full_join, by = "Sample_Name")
  return(merged)
}

merged_data <- map(folders, ~process_folder(.x, base_path))
names(merged_data) <- folders

# save merged_data
for (name in names(merged_data)) {
  df <- merged_data[[name]]
  write.csv(df, paste0(output_folder, "all_predictions_", name, ".csv"), row.names = FALSE)
}

# check whether NA present and whether values of clocks makes sense (replace with your folder names)
sink(file.path(output_folder, "summary_output_clocks.txt"))
cat("Summary for 450k_cord:\n")
summary(merged_data[["450k_cord"]])
cat("\nSummary for 450k_F7:\n")
summary(merged_data[["450k_F7"]])
cat("\nSummary for 450k_15up:\n")
summary(merged_data[["450k_15up"]])
cat("\nSummary for epic_15up:\n")
summary(merged_data[["epic_15up"]])
sink()



# save summary files

# function to create a summary table for a data frame
create_summary_df <- function(df) {
  data.frame(
    Column = names(df),
    NAs = colSums(is.na(df)),
    NonNAs = colSums(!is.na(df)),
    # Basic stats for numeric columns (mean, median, sd), NA for non-numeric
    Mean = sapply(df, function(x) if(is.numeric(x)) mean(x, na.rm = TRUE) else NA),
    Median = sapply(df, function(x) if(is.numeric(x)) median(x, na.rm = TRUE) else NA),
    SD = sapply(df, function(x) if(is.numeric(x)) sd(x, na.rm = TRUE) else NA)
  )
}

# !! add your path to save summary output
for (folder in folders) {
  summary_df <- create_summary_df(merged_data[[folder]])
  write.csv(summary_df, file.path(output_folder, "summary_", folder, ".csv"), row.names = FALSE)
}

## ! **if you have cord blood/birth time point**
# read in combined cord bood predictions and add in gestational clock predictions
# e.g.:
cord_all <- read.csv(file.path(wd, "models/all_predictions_450k_cord.csv")) # read in all cord blood predictions (! modify to match your file names)
gest_clocks <- read.xlsx(file.path(wd, "models/450k_cord/epi_ages_450k_cord_pred.xlsx"), sheet = 2) # read in gestational clock predictions (! modify to match your file names)
cord_merged <- merge(cord_all, gest_clocks, by.x='Sample_Name', by.y='id', all=T) # merge in gest clocks
write.csv(cord_merged, file.path(output_folder, "all_predictions_450k_cord_with_GA.csv"), row.names = FALSE)

