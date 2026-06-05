## ============================================================
##  Developmental Brain Age – Feature Matrix Builder
##  - Reads FS / ENIGMA files
##  - Matches features by SUBJID using feature map
##  - Keeps only subjects present in ALL data sources
##  - Outputs CSV with and without IDs
## ============================================================

cat("Starting feature matrix construction...\n\n")

## ---------- 0. Helper: simple assert function ----------

assert <- function(cond, msg) {
  if (!cond) stop(msg, call. = FALSE)
}

## ---------- 1. Read input data files ----------

cat("Reading input files...\n")

# FreeSurfer text outputs (check.names = FALSE to preserve original column names)
vol_lh <- read.table("aparc_stats_vol_lh.txt",
                     header = TRUE,
                     sep = "",
                     stringsAsFactors = FALSE,
                     check.names = FALSE)
names(vol_lh)[1] <- "SUBJID"

vol_rh <- read.table("aparc_stats_vol_rh.txt",
                     header = TRUE,
                     sep = "",
                     stringsAsFactors = FALSE,
                     check.names = FALSE)
names(vol_rh)[1] <- "SUBJID"

aseg <- read.table("aseg.stats.txt",
                   header = TRUE,
                   sep = "",
                   stringsAsFactors = FALSE,
                   check.names = FALSE)
names(aseg)[1] <- "SUBJID"


# ENIGMA CSVs (keep original names too)
SA <- read.csv("CorticalMeasuresENIGMA_SurfAvg.csv",
               stringsAsFactors = FALSE,
               check.names = FALSE)

sub <- read.csv("LandRvolumes.csv",
                stringsAsFactors = FALSE,
                check.names = FALSE)

## ---------- 2. Normalise ID column in SA and sub ----------

# Common case: ENIGMA files use "SubjID" – convert to "SUBJID"
if ("SubjID" %in% names(SA)) {
  names(SA)[names(SA) == "SubjID"] <- "SUBJID"
}
if ("SubjID" %in% names(sub)) {
  names(sub)[names(sub) == "SubjID"] <- "SUBJID"
}

## ---------- 3. Put dataframes into a named list & check IDs ----------

dfs <- list(
  vol_lh   = vol_lh,
  vol_rh   = vol_rh,
  aseg     = aseg,
  SA       = SA,
  sub      = sub
)

cat("Checking that all dataframes have SUBJID and no duplicate IDs...\n")

for (nm in names(dfs)) {
  df <- dfs[[nm]]
  
  # Ensure SUBJID exists
  assert("SUBJID" %in% names(df),
         paste("Dataframe", nm, "does not have a SUBJID column."))
  
  # Coerce to character for safe matching
  df$SUBJID <- as.character(df$SUBJID)
  
  # Stop on duplicate IDs
  if (any(duplicated(df$SUBJID))) {
    dup_ids <- unique(df$SUBJID[duplicated(df$SUBJID)])
    stop(paste("❌ Duplicated SUBJID values found in", nm, ":",
               paste(dup_ids, collapse = ", ")))
  }
  
  dfs[[nm]] <- df
}

## ---------- 4. Load and sanity-check feature map ----------

cat("Loading feature mapping file for FS below 7... \n")

feature_map <- read.delim("devbrainage_feature_map.txt",
                          sep = "\t",
                          header = TRUE,
                          stringsAsFactors = FALSE)

required_cols <- c("feature_list", "df_avail", "df_avail_feature_name")
assert(all(required_cols %in% names(feature_map)),
       paste("feature_map must contain columns:",
             paste(required_cols, collapse = ", ")))

if (any(duplicated(feature_map$feature_list))) {
  warning("Duplicate entries found in feature_map$feature_list; ",
          "later duplicates will overwrite earlier ones in the output.")
}

## ---------- 5. Build master SUBJID list across ALL files (intersection) ----------

cat("Building master subject ID list (intersection across all files)...\n")

id_lists <- lapply(dfs, function(d) unique(d$SUBJID))
master_ids <- Reduce(intersect, id_lists)
master_ids <- sort(master_ids)

cat("Number of subjects present in ALL input datasets (only these will be taken forward):", length(master_ids), "\n")

if (length(master_ids) == 0) {
  stop("❌ No common SUBJIDs across all input datasets. Check your IDs.")
}

## ---------- 6. Build feature matrix using feature map ----------

cat("Constructing model feature matrix using feature_map...\n")

feature_data <- list()
feature_data[["SUBJID"]] <- master_ids   # ID column (first)

for (i in seq_len(nrow(feature_map))) {
  feature_name <- feature_map$feature_list[i]
  df_name      <- feature_map$df_avail[i]
  orig_col     <- feature_map$df_avail_feature_name[i]
  
  # Basic sanity: df_name must be non-empty and present in dfs
  if (is.na(df_name) || df_name == "") {
    warning("Row ", i, ": df_avail is empty/NA; filling ", feature_name, " with NA.")
    feature_data[[feature_name]] <- rep(NA, length(master_ids))
    next
  }
  
  if (!df_name %in% names(dfs)) {
    warning("Row ", i, ": dataframe '", df_name,
            "' not found; filling ", feature_name, " with NA.")
    feature_data[[feature_name]] <- rep(NA, length(master_ids))
    next
  }
  
  df <- dfs[[df_name]]
  
  # Check the column exists
  if (!orig_col %in% names(df)) {
    warning("Row ", i, ": column '", orig_col,
            "' not found in dataframe '", df_name,
            "'; filling ", feature_name, " with NA.")
    feature_data[[feature_name]] <- rep(NA, length(master_ids))
    next
  }
  
  # Match this dataframe to the master ID list (intersection guarantees all present)
  idx <- match(master_ids, df$SUBJID)   # should never be NA if intersection is correct
  
  if (any(is.na(idx))) {
    warning("Some master_ids not found in dataframe '", df_name,
            "' when they should be. Filling those as NA for feature ", feature_name, ".")
  }
  
  vec <- df[[orig_col]][idx]
  feature_data[[feature_name]] <- vec
}

## ---------- 7. Combine into final data.frame and save ----------

cat("Combining mapped features into final dataframe...\n")

final_features_df <- as.data.frame(feature_data, stringsAsFactors = FALSE)

cat("Saving outputs...\n")

# With IDs
write.csv(final_features_df,
          "DevelopmentalBrainAge_input_withIDs.csv",
          row.names = FALSE)
cat("Saved DevelopmentalBrainAge_input_withIDs.csv\n")

# Without IDs
final_features_noID <- final_features_df
final_features_noID$SUBJID <- NULL

write.csv(final_features_noID,
          "DevelopmentalBrainAge_input.csv",
          row.names = FALSE)
cat("Saved DevelopmentalBrainAge_input.csv\n")

cat("\nDone.\n")
