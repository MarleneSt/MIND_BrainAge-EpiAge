#### Parts of this script were originally written for Brain age project by 
#### Laura Han (L.Han@ggzingeest.nl) and Esther Walton (waltonesther@gmail.com)
#### It has been adapted for the MIND brainage~epiage project by Marlene Staginnus

#### BrainAge~EpiAge project  - minimal clean version ####

# Start logging
messages = file("ENIGMA_prepare_data_for_website.log", open = "wt")
sink(messages, type = "message")
sink(messages, type = "output")

# Load only necessary libraries
cat("Prep: installing/loading libraries\n")
load.lib <- function(x) {
  for (i in x) {
    if (!require(i, character.only = TRUE)) {
      install.packages(i, dependencies = TRUE)
      library(i, character.only = TRUE)
    }
  }
}
load.lib(c("dplyr"))  # Only one actually needed

# Source only used functions
cat("Prep: sourcing functions\n")
source("get.means_MS050625.R")  
source("prepare.files_MS050625.R")

# Load imaging data
cat("Loading & averaging imaging data\n")
Thick = get.means("CorticalMeasuresENIGMA_ThickAvg.csv"); Thick$ICV = NULL
Surf  = get.means("CorticalMeasuresENIGMA_SurfAvg.csv"); Surf$ICV = NULL
Vol   = get.means("LandRvolumes.csv")

# Merge imaging data
TS  = merge(Thick, Surf, by = "row.names")
TSV = merge(TS, Vol, by.x = "Row.names", by.y = "row.names")

# Read covariates
Covs <- read.csv("Covariates.csv")
names(Covs) <- toupper(names(Covs))  # Capitalize all for standardization

# Check required columns
required_cols <- c("SUBJID", "SEX_BRAIN", "AGE_BRAIN")
missing_cols <- setdiff(required_cols, names(Covs))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns in Covariates.csv:", paste(missing_cols, collapse = ", ")))
}

# Check for duplicate subject IDs
if (anyDuplicated(Covs$SUBJID) != 0) {
  stop("Duplicate SUBJID values found in Covariates.csv. Please remove or resolve them.")
}

# Merge covariates with imaging data
data = merge(Covs, TSV, by.x = "SUBJID", by.y = "Row.names")

# Prepare sex-stratified input files
cat("Creating csv files for brainAge estimation\n")
df = prepare.files(data, names(Covs))

males = df$males
females = df$females
rm(df)

invisible(readline(prompt = "You should see 2 csv files in your working directory:\n
- females_raw.csv\n
- males_raw.csv.\n
These should be up/down-loaded as described in your documentation.\n\nPress enter again to finish the script.\n\n"))

cat("Finished with data preparation for upload!\n")
sink()
closeAllConnections()
print(readLines("ENIGMA_prepare_data_for_website.log"))
