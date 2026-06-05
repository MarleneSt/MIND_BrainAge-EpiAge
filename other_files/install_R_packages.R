# -----------------------------------------------------------------
# Run once, replacing “packages_to_install” with what you need
# Here packages were the ones called by imputation and apply R scripts for MIND Brain age models
# my_lib_path will create a new folder in your home directory (can be modified to point elsewhere
# script created by Valentine Chirokoff
# -----------------------------------------------------------------

# 1. Define the personal library path explicitly (Creates a new folder in their home directory)
my_lib_path <- file.path(Sys.getenv("HOME"), "R", "mind_project_libs")

# 2. Set CRAN mirror and library path
if (!dir.exists(my_lib_path)) {
  message("Creating personal library directory: ", my_lib_path)
  dir.create(my_lib_path, recursive = TRUE)
}
# Set the library search path to prioritize the user's writable folder
.libPaths(c(my_lib_path, .libPaths()))
# Set CRAN mirror to fix the download error
options(repos = c(CRAN = "https://cloud.r-project.org"))


# Define the list of all packages needed (This list is complete)
packages_to_install <- c(
  "skimr", 
  "psych", 
  "missForest", 
  "dplyr",
  "tibble",
  "tidyverse",
  "tidymodels", 
  "remotes", 
  "devtools",
  "data.table", 
  "readr"       
)

# Installation loop
for (pkg in packages_to_install) {
    if (!require(pkg, character.only = TRUE)) {
        message(paste("Installing", pkg, "..."))
        # Explicitly install to the user's personal library
        install.packages(pkg, dependencies = TRUE, lib = my_lib_path)
    } else {
        message(paste("Package", pkg, "already installed."))
    }
}

# Install specific xgboost version
if (!("xgboost" %in% installed.packages(lib.loc = my_lib_path)[, "Package"]) || packageVersion("xgboost") != "1.0.0.2") {
    message("Installing specific xgboost version 1.0.0.2...")
    library(remotes)
    remotes::install_version("xgboost", version = "1.0.0.2", lib = my_lib_path, upgrade = "never")
} else {
    message("xgboost 1.0.0.2 is already installed.")
}

print("All R packages installed successfully to personal library.")