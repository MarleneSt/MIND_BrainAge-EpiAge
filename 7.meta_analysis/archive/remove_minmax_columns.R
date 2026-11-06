###### remove_minmax_columns.R: one-off removal of single-person extremes ------
# Removes the min, max and whisker columns from the cohort performance files in
# data/ (e.g. AE_min, AE_max, predage_min, predage_max, PAR_lower_whisker).
# These describe a single participant's value and are not used by any script.
# All other columns are unchanged. Run once from the Reproducibilty folder.

files <- list.files("data", pattern = "_performance_(brain|epi)_.*_main\\.csv$",
                    recursive = TRUE, full.names = TRUE)

drop_pattern <- "^(AE|wAE|PAD|predage|PAR)_(min|max|lower_whisker|upper_whisker)$"

for (f in files) {
  d <- read.csv(f, check.names = FALSE)
  d <- d[, !grepl(drop_pattern, names(d)), drop = FALSE]
  write.csv(d, f, row.names = FALSE)
}

cat("Columns removed from", length(files), "performance files.\n")
