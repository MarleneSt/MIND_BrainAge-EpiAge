#########################################################
############# SCRIPT TO IDENTIFY OUTLIERS ###############
################ based on Euler number ###############
#########################################################
## written by Marlene Staginnus 

# Load required packages
list.of.packages <- c("dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(dplyr)

# Step 1: Read the file
aseg_data <- read.table("aseg.stats.txt", header = TRUE, sep = "\t")

# Step 2: Rename the column Measure:volume to SubjID
colnames(aseg_data)[colnames(aseg_data) == "Measure.volume"] <- "SubjID"

# Step 3: Calculate the Euler number for each hemisphere
aseg_data$euler_LH <- 2 - 2 * aseg_data$lhSurfaceHoles
aseg_data$euler_RH <- 2 - 2 * aseg_data$rhSurfaceHoles

# Step 4: Calculate the lower inner fence for outliers (Q1 - 1.5 * IQR)
calc_outliers <- function(data_column) {
  Q1 <- quantile(data_column, 0.25)
  Q3 <- quantile(data_column, 0.75)
  IQR <- IQR(data_column)
  lower_inner_fence <- Q1 - 1.5 * IQR
  return(lower_inner_fence)
}

# Get lower inner fences for both hemispheres
lower_inner_fence_LH <- calc_outliers(aseg_data$euler_LH)
lower_inner_fence_RH <- calc_outliers(aseg_data$euler_RH)

# Step 5: Identify outliers (below the lower inner fence) for each hemisphere
outliers_LH <- aseg_data$euler_LH < lower_inner_fence_LH
outliers_RH <- aseg_data$euler_RH < lower_inner_fence_RH

# Combine both conditions for participants with outlying values
outliers <- outliers_LH | outliers_RH

# Step 6: Create a list of subject IDs with outlying Euler values
outlying_subjects <- aseg_data$SubjID[outliers]

# Step 7: Output the message and the list of subject IDs
if(length(outlying_subjects) > 0) {
  cat("The following participants have an Euler number lower than Q1 - 1.5*IQR in at least one hemisphere and should be inspected further and/or excluded:\n")
  print(outlying_subjects)
  
  # Save the list of outlying subjects to a CSV file
  write.csv(outlying_subjects, "participants_outlying_euler.csv", row.names = FALSE)
} else {
  cat("No participants have an outlying Euler number.\n")
}
