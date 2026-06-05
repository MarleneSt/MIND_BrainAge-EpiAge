#### set WD to ENIGMA files ----
# likely: brain_epi/mri/extracted_features or a time-point specific subfolder

setwd("/Volumes/files/Psychology/ResearchProjects/TPFreeman/CannTeen/data/Marlene_MIND/brain_epi/mri/extracted_features")

#### creating folder for imputation ----
dir.create("ENIGMA_imputed", showWarnings = FALSE)
setwd("ENIGMA_imputed")

#### start logging ----
messages <- file("ENIGMA_imputation.log", open = "wt")
sink(messages, type = "message")
sink(messages, type = "output")

#### installing & loading libraries ----
cat("1. installing/loading libraries\n\n")

load.lib <- function(x) {
  for (i in x) {
    if (!require(i, character.only = TRUE)) {
      install.packages(i, dependencies = TRUE)
      library(i, character.only = TRUE)
    }
  }
}
load.lib(c("skimr", "psych", "missForest")) 

#### set seed ----
cat("2. setting seed\n\n")

set.seed(21021995)

#### reading data ----
cat("3. reading imaging and covariate files\n\n")
dat_imaging_miss_CT  <- read.csv("../CorticalMeasuresENIGMA_ThickAvg.csv")
dat_imaging_miss_SA  <- read.csv("../CorticalMeasuresENIGMA_SurfAvg.csv")
dat_imaging_miss_sub <- read.csv("../LandRvolumes.csv")

Covs <- read.csv("../../Covariates.csv")
names(Covs) <- toupper(names(Covs))  # Capitalize all for standardization

#### data checks ----

cat("4. Performing some data checks (availability, duplicate IDs)\n")

# checks for cov availability
required_cols <- c("SUBJID", "SEX_BRAIN", "AGE_BRAIN")
missing_cols <- setdiff(required_cols, names(Covs))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns in Covariates.csv:", paste(missing_cols, collapse = ", ")))
}

# checks for imaging availability 
required_cols <- c('SUBJID', 'L_bankssts_thickavg', 'L_caudalanteriorcingulate_thickavg', 'L_caudalmiddlefrontal_thickavg', 'L_cuneus_thickavg', 'L_entorhinal_thickavg', 'L_fusiform_thickavg', 'L_inferiorparietal_thickavg', 'L_inferiortemporal_thickavg', 'L_isthmuscingulate_thickavg', 'L_lateraloccipital_thickavg', 'L_lateralorbitofrontal_thickavg', 'L_lingual_thickavg', 'L_medialorbitofrontal_thickavg', 'L_middletemporal_thickavg', 'L_parahippocampal_thickavg', 'L_paracentral_thickavg', 'L_parsopercularis_thickavg', 'L_parsorbitalis_thickavg', 'L_parstriangularis_thickavg', 'L_pericalcarine_thickavg', 'L_postcentral_thickavg', 'L_posteriorcingulate_thickavg', 'L_precentral_thickavg', 'L_precuneus_thickavg', 'L_rostralanteriorcingulate_thickavg', 'L_rostralmiddlefrontal_thickavg', 'L_superiorfrontal_thickavg', 'L_superiorparietal_thickavg', 'L_superiortemporal_thickavg', 'L_supramarginal_thickavg', 'L_frontalpole_thickavg', 'L_temporalpole_thickavg', 'L_transversetemporal_thickavg', 'L_insula_thickavg', 'R_bankssts_thickavg', 'R_caudalanteriorcingulate_thickavg', 'R_caudalmiddlefrontal_thickavg', 'R_cuneus_thickavg', 'R_entorhinal_thickavg', 'R_fusiform_thickavg', 'R_inferiorparietal_thickavg', 'R_inferiortemporal_thickavg', 'R_isthmuscingulate_thickavg', 'R_lateraloccipital_thickavg', 'R_lateralorbitofrontal_thickavg', 'R_lingual_thickavg', 'R_medialorbitofrontal_thickavg', 'R_middletemporal_thickavg', 'R_parahippocampal_thickavg', 'R_paracentral_thickavg', 'R_parsopercularis_thickavg', 'R_parsorbitalis_thickavg', 'R_parstriangularis_thickavg', 'R_pericalcarine_thickavg', 'R_postcentral_thickavg', 'R_posteriorcingulate_thickavg', 'R_precentral_thickavg', 'R_precuneus_thickavg', 'R_rostralanteriorcingulate_thickavg', 'R_rostralmiddlefrontal_thickavg', 'R_superiorfrontal_thickavg', 'R_superiorparietal_thickavg', 'R_superiortemporal_thickavg', 'R_supramarginal_thickavg', 'R_frontalpole_thickavg', 'R_temporalpole_thickavg', 'R_transversetemporal_thickavg', 'R_insula_thickavg', 'LThickness', 'RThickness', 'LSurfArea', 'RSurfArea', 'ICV')
missing_cols <- setdiff(required_cols, names(dat_imaging_miss_CT))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns in CorticalMeasuresENIGMA_ThickAvg.csv.:", paste(missing_cols, collapse = ", ")))
}

required_cols <- c('SUBJID', 'L_bankssts_surfavg', 'L_caudalanteriorcingulate_surfavg', 'L_caudalmiddlefrontal_surfavg', 'L_cuneus_surfavg', 'L_entorhinal_surfavg', 'L_fusiform_surfavg', 'L_inferiorparietal_surfavg', 'L_inferiortemporal_surfavg', 'L_isthmuscingulate_surfavg', 'L_lateraloccipital_surfavg', 'L_lateralorbitofrontal_surfavg', 'L_lingual_surfavg', 'L_medialorbitofrontal_surfavg', 'L_middletemporal_surfavg', 'L_parahippocampal_surfavg', 'L_paracentral_surfavg', 'L_parsopercularis_surfavg', 'L_parsorbitalis_surfavg', 'L_parstriangularis_surfavg', 'L_pericalcarine_surfavg', 'L_postcentral_surfavg', 'L_posteriorcingulate_surfavg', 'L_precentral_surfavg', 'L_precuneus_surfavg', 'L_rostralanteriorcingulate_surfavg', 'L_rostralmiddlefrontal_surfavg', 'L_superiorfrontal_surfavg', 'L_superiorparietal_surfavg', 'L_superiortemporal_surfavg', 'L_supramarginal_surfavg', 'L_frontalpole_surfavg', 'L_temporalpole_surfavg', 'L_transversetemporal_surfavg', 'L_insula_surfavg', 'R_bankssts_surfavg', 'R_caudalanteriorcingulate_surfavg', 'R_caudalmiddlefrontal_surfavg', 'R_cuneus_surfavg', 'R_entorhinal_surfavg', 'R_fusiform_surfavg', 'R_inferiorparietal_surfavg', 'R_inferiortemporal_surfavg', 'R_isthmuscingulate_surfavg', 'R_lateraloccipital_surfavg', 'R_lateralorbitofrontal_surfavg', 'R_lingual_surfavg', 'R_medialorbitofrontal_surfavg', 'R_middletemporal_surfavg', 'R_parahippocampal_surfavg', 'R_paracentral_surfavg', 'R_parsopercularis_surfavg', 'R_parsorbitalis_surfavg', 'R_parstriangularis_surfavg', 'R_pericalcarine_surfavg', 'R_postcentral_surfavg', 'R_posteriorcingulate_surfavg', 'R_precentral_surfavg', 'R_precuneus_surfavg', 'R_rostralanteriorcingulate_surfavg', 'R_rostralmiddlefrontal_surfavg', 'R_superiorfrontal_surfavg', 'R_superiorparietal_surfavg', 'R_superiortemporal_surfavg', 'R_supramarginal_surfavg', 'R_frontalpole_surfavg', 'R_temporalpole_surfavg', 'R_transversetemporal_surfavg', 'R_insula_surfavg', 'LThickness', 'RThickness', 'LSurfArea', 'RSurfArea', 'ICV')
missing_cols <- setdiff(required_cols, names(dat_imaging_miss_SA))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns in CorticalMeasuresENIGMA_SurfAvg.csv.:", paste(missing_cols, collapse = ", ")))
}

required_cols <- c('SUBJID', 'LLatVent', 'RLatVent', 'Lthal', 'Rthal', 'Lcaud', 'Rcaud', 'Lput', 'Rput', 'Lpal', 'Rpal', 'Lhippo', 'Rhippo', 'Lamyg', 'Ramyg', 'Laccumb', 'Raccumb', 'ICV')
missing_cols <- setdiff(required_cols, names(dat_imaging_miss_sub))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns in LandRvolumes.csv.:", paste(missing_cols, collapse = ", ")))
}

# Check for duplicate subject IDs across files
if (anyDuplicated(Covs$SUBJID) != 0) {
  stop("Duplicate SUBJID values found in Covariates.csv. Please remove or resolve them.")
}
if (anyDuplicated(dat_imaging_miss_CT$SUBJID) != 0) {
  stop("Duplicate SUBJID values found CorticalMeasuresENIGMA_ThickAvg.csv. Please remove or resolve them.")
}
if (anyDuplicated(dat_imaging_miss_SA$SUBJID) != 0) {
  stop("Duplicate SUBJID values found CorticalMeasuresENIGMA_SurfAvg.csv. Please remove or resolve them.")
}
if (anyDuplicated(dat_imaging_miss_sub$SUBJID) != 0) {
  stop("Duplicate SUBJID values found LandRvolumes.csv. Please remove or resolve them.")
}

cat("All data checks passed\n\n")

#### merging files ----

cat("5. Merging files with global values \n\n")

Covs <- Covs[c("SUBJID", "SEX_BRAIN", "AGE_BRAIN")]

# Capture original column orders for later restoration
ct_cols_all  <- names(dat_imaging_miss_CT)
sa_cols_all  <- names(dat_imaging_miss_SA)
sub_cols_all <- names(dat_imaging_miss_sub)

# Define "global" columns  
globals <- c("LThickness","RThickness","LSurfArea","RSurfArea","ICV")
ct_globals  <- globals
sa_globals  <- globals   
sub_globals <- c("ICV")

# Columns to impute from each file (exclude SUBJID + globals)
ct_keep  <- setdiff(names(dat_imaging_miss_CT),  c("SUBJID", ct_globals))
sa_keep  <- setdiff(names(dat_imaging_miss_SA),  c("SUBJID", sa_globals))
sub_keep <- setdiff(names(dat_imaging_miss_sub), c("SUBJID", sub_globals))

# For CT, we WANT to include globals ONCE (so create an imputation set with CT globals)
ct_impute_cols <- c(ct_keep, globals)

# Merge, include CT globals in the CT slice; exclude from SA/sub to avoid duplicates
data <- merge(Covs, dat_imaging_miss_CT[c("SUBJID", ct_impute_cols)],  by = "SUBJID")
data <- merge(data, dat_imaging_miss_SA[c("SUBJID", sa_keep)],         by = "SUBJID")
data <- merge(data, dat_imaging_miss_sub[c("SUBJID", sub_keep)],       by = "SUBJID")

# define imaging columns (everything except ID/sex/age)
imaging_cols <- setdiff(names(data), c("SUBJID","SEX_BRAIN","AGE_BRAIN"))


#### setting correct var type ----

cat("6. Setting correct var types \n\n")

data$SUBJID    <- as.character(data$SUBJID)
data$AGE_BRAIN <- as.numeric(data$AGE_BRAIN)
data$SEX_BRAIN <- as.factor(data$SEX_BRAIN)

data[imaging_cols] <- lapply(data[imaging_cols], function(x)
  suppressWarnings(as.numeric(as.character(x)))
)

#### descriptive pre imputation ----
cat("7. Generating descriptive & missingness pre imputation\n")


desc_psych <- psych::describe(data[imaging_cols])
desc_psych_df <- data.frame(variable = rownames(desc_psych),
                            desc_psych, row.names = NULL)

skim_df <- as.data.frame(skimr::skim(data[imaging_cols]))
skim_keep <- unique(skim_df[c("skim_variable", "n_missing", "complete_rate")])

combined <- merge(desc_psych_df, skim_keep,
                  by.x = "variable", by.y = "skim_variable",
                  all.x = TRUE, sort = FALSE)

write.csv(combined, "ENIGMA_input_descriptives_preimpute.csv", row.names = FALSE)

cat("Descriptives saved as ENIGMA_input_descriptives_preimpute.csv\n")

# Overall missingness across imaging columns only 
n_missing <- sum(is.na(data[imaging_cols]))
perc_missing <- round(n_missing / (nrow(data) * length(imaging_cols)) * 100, 2)
cat("Overall,", n_missing, "datapoints are missing (", perc_missing, "% ).\n", sep = " ")

# stop immediately if there is NO missingness
if (n_missing == 0) {
  stop("No missingness detected in imaging variables. Imputation stopped.")
}

# Row-wise % missing over imaging columns (align numerator/denominator)
row_missing <- rowSums(is.na(data[imaging_cols]))
row_missing_percent <- row_missing / length(imaging_cols) * 100

# Participants with >30% missing (your step preserved)
high_missing <- data$SUBJID[row_missing_percent > 30]
if (length(high_missing) > 0) {
  cat("Participants with >30% missingness (consider removing?):\n\n")
  cat(high_missing, sep = "\n"); cat("\n")
} else {
  cat("No participants with >30% missingness.\n\n")
}

#### imputation ---- 

cat("8. imputation\n")

X <- data[c("SEX_BRAIN","AGE_BRAIN", imaging_cols)]

### CHANGE: defensive stop right before imputation too (in case upstream changes)
total_missing <- sum(is.na(X))
if (total_missing == 0) {
  stop("No missingness in the imputation matrix. Imputation stopped.")
} else {
  cat("Missingness detected. Proceeding with imputation.\n")
}

# running imputation 
mf <- missForest(X, verbose = TRUE)  # was: missForest(X, verbose = TRUE)
X_imp <- mf$ximp

# imputation quality 
cat("Overall OOBerror (Numeric-NRMSE / Categorical-PFC):", mf$OOBerror, "\n\n")

#### adding imputed values into DF ----

cat("9. Adding imputed values back into dataframes\n\n")

ct_out  <- dat_imaging_miss_CT
sa_out  <- dat_imaging_miss_SA
sub_out <- dat_imaging_miss_sub

# Align imputed rows by SUBJID
rownames(X_imp) <- data$SUBJID

inject_imputed <- function(df, keep_cols, Ximp) {
  common_ids <- intersect(df$SUBJID, rownames(Ximp))
  if (!length(common_ids)) return(df)
  df_idx <- match(common_ids, df$SUBJID)
  xi_idx <- match(common_ids, rownames(Ximp))
  cols_to_write <- intersect(keep_cols, colnames(Ximp))
  if (length(cols_to_write)) {
    df[df_idx, cols_to_write] <- Ximp[xi_idx, cols_to_write]
  }
  df
}

ct_out  <- inject_imputed(ct_out,  ct_keep,  X_imp)
sa_out  <- inject_imputed(sa_out,  sa_keep,  X_imp)
sub_out <- inject_imputed(sub_out, sub_keep, X_imp)

# also inject the SINGLE imputed set of globals to SA and sub
# we imputed CT globals; now copy them into SA (all 5) and sub (ICV)
ct_out  <- inject_imputed(ct_out,  globals,  X_imp)   # <-- write CT globals back too
sa_out  <- inject_imputed(sa_out,  globals,  X_imp)   # L/R Thickness, L/R SurfArea, ICV
sub_out <- inject_imputed(sub_out, "ICV",    X_imp)   # ICV only

# Restore original column order exactly
ct_out  <- ct_out[ct_cols_all]
sa_out  <- sa_out[sa_cols_all]
sub_out <- sub_out[sub_cols_all]

#### saving imputed data  ----

cat("10. Saving imputed data \n")

# Save round-tripped files (same columns, same order; globals untouched)
write.csv(ct_out,  "CorticalMeasuresENIGMA_ThickAvg.csv", row.names = FALSE)
write.csv(sa_out,  "CorticalMeasuresENIGMA_SurfAvg.csv",  row.names = FALSE)
write.csv(sub_out, "LandRvolumes.csv",                    row.names = FALSE)

cat("Imputed data saved in subfolder 'ENIGMA_imputed'\n\n")

#### post-imputation descriptives ----
cat("11. Generating descriptive post imputation\n")

desc_psych_post <- psych::describe(X_imp[imaging_cols])
desc_psych_post_df <- data.frame(variable = rownames(desc_psych_post),
                                 desc_psych_post, row.names = NULL)
skim_df_post <- as.data.frame(skimr::skim(X_imp[imaging_cols]))
skim_keep_post <- unique(skim_df_post[c("skim_variable", "n_missing", "complete_rate")])

combined_post <- merge(desc_psych_post_df, skim_keep_post,
                       by.x = "variable", by.y = "skim_variable",
                       all.x = TRUE, sort = FALSE)
write.csv(combined_post, "ENIGMA_input_descriptives_postimpute.csv", row.names = FALSE)

cat("Post imputation descriptives saved as ENIGMA_input_descriptives_preimpute.csv\n")

cat("Imputation script finished.\n")


#### close sinks ----
sink(type = "message"); sink(type = "output")
close(messages)

