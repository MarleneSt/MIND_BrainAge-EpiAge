# =============================================================================
# All-in-one missForest imputation helper for "one-table" brain-age models
# (ID + covariates + imaging features)
#
# Uses all other input features + age & sex for the imputation 
#
# What it does:
#   1) Validates the input (required columns, duplicate IDs)
#   2) Coerces types (ID→character, sex→factor, age→numeric, features→numeric)
#   3) Creates pre-imputation descriptives (psych + skimr)
#   4) Early-stops if there is no missingness
#   5) Runs missForest, prints human-readable OOB error
#   6) Creates post-imputation descriptives
#   7) Reinserts imputed feature values back into the original df (in-memory)
#   8) Saves: pre/post descriptives and (optionally) the imputed dataset
#   9) Optionally tees all printed output to a log file
#
# Key inputs:
#   - feature_cols   : which columns to impute/overwrite (default: all non-ID, non-covariates)
#   - out_dir        : base directory for outputs (e.g., "imputed")
#   - out_prefix     : file stem (e.g., "kaufmann/kaufmann_input")
#   - log_path       : path to a .log file (captures cat/print output)
#   - save_imputed   : write the imputed dataset CSV when TRUE
#
# Returns a list with imputed matrix, OOB error, pre/post descriptives,
# the in-memory imputed df, and the paths of saved files.
# =============================================================================

impute_brain_df <- function(
    df,                          # data.frame with SUBJID + covariates + imaging features
    id_col  = "SUBJID",          # ID column name
    sex_col = "SEX_BRAIN",       # sex covariate (categorical)
    age_col = "AGE_BRAIN",       # age covariate (numeric)
    feature_cols = NULL,         # imaging variables to impute (NULL -> auto pick)
    seed = 21021995,             # seed for reproducibility
    stop_if_no_missing = TRUE,   # stop early if nothing to impute
    verbose = TRUE,              # pass-through verbosity to missForest
    out_prefix = NULL,           # e.g., "kaufmann/kaufmann_input" (no extension)
    out_dir = ".",               # output root directory (will create subfolders)
    log_path = NULL,             # write a log file with console prints if not NULL
    save_imputed = TRUE          # save "{out_prefix}_features_imputed.csv" if TRUE
) {
  # -------- start logging everything we cat()/print() to a file -----
  # We only "sink" the *output* stream here, not the "message" stream.
  if (!is.null(log_path)) {
    dir.create(dirname(log_path), showWarnings = FALSE, recursive = TRUE)
    log_con <- file(log_path, open = "wt")         # open log file for writing
    sink(log_con, type = "output")                 # redirect console output to file
    on.exit({
      while (sink.number() > 0) sink()             # restore console output
      if (isOpen(log_con)) close(log_con)          # close log on exit
    }, add = TRUE)
  }
  
  cat("== Imputation & descriptives script started ==\n\n")
  
  # -------- Decide which feature columns to impute/overwrite ------------------
  # If not specified, treat *all* columns except ID/sex/age as features.
  if (is.null(feature_cols)) {
    feature_cols <- setdiff(names(df), c(id_col, sex_col, age_col))
  }
  
  # -------- Sanity checks: columns exist -------------------------------------
  
  cat("Running checks on required columns & duplicated IDs\n\n")
  
  
  # We require ID, sex, age, and every feature column to be present.
  required <- c(id_col, sex_col, age_col, feature_cols)
  miss <- setdiff(required, names(df))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))
  
  # -------- Sanity checks: duplicate IDs -------------------------------------
  # Imputation & reinjection require one row per ID.
  if (anyDuplicated(df[[id_col]]) > 0) {
    bad <- unique(df[[id_col]][duplicated(df[[id_col]])])
    stop("Duplicate IDs in ", id_col, ": ", paste(bad, collapse = ", "))
  }
  
  # -------- Type hygiene (by *name*, not by position) ------------------------
  # - Keep IDs as character to avoid factor/level issues
  # - Ensure SEX is factor (categorical)
  # - Ensure AGE is numeric (continuous)
  # - Ensure all features are numeric (missForest can handle mixed types,
  #   but our use-case here is almost entirely numeric imaging variables)
  
  cat("Setting correct var types\n\n")

  df[[id_col]]  <- as.character(df[[id_col]])
  df[[sex_col]] <- as.factor(df[[sex_col]])
  df[[age_col]] <- suppressWarnings(as.numeric(df[[age_col]]))
  df[feature_cols] <- lapply(
    df[feature_cols],
    function(x) suppressWarnings(as.numeric(as.character(x)))
  )
  
  # -------- Build the imputation matrix X ------------------------------------
  # X = covariates (predictors) + features (targets to impute).
  # Covariates are used to *predict* missing feature values but are not overwritten
  X <- df[c(sex_col, age_col, feature_cols)]
  
  # -------- Pre-imputation descriptives (features only) ----------------------
  # psych::describe + skimr::skim (n_missing, complete_rate)
  pre_desc <- .describe_df(df[feature_cols])
  
  # If a prefix is provided, write the CSV to disk.
  if (!is.null(out_prefix)) {
    dir.create(file.path(out_dir, dirname(out_prefix)), showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(
      pre_desc,
      paste0(out_prefix, "_descr_preimputation.csv"),
      row.names = FALSE
    )
    cat("Saved pre-imputation descriptives as ",
        file.path(paste0(out_prefix, "_descr_preimputation.csv")), "\n\n", sep = "")
  }
  
  # -------- Early stop if nothing to impute ----------------------------------
  # We inspect X (covariates + features). Because AGE/SEX are usually complete,
  # the missingness will typically come from the features.
  
  cat("Check if there is missingness\n\n")
  
  tot_miss <- sum(is.na(X))
  total_cells <- nrow(X) * ncol(X)
  perc_miss <- round((tot_miss / total_cells) * 100, 2)
  
  cat("Missing values among imaging features : ", tot_miss, " (", perc_miss, "% of total data)\n\n", sep = "")
  
  if (tot_miss == 0 && stop_if_no_missing) {
    cat("No missingness detected. Skipping imputation.\n\n")
    return(list(
      X_imp = X,
      oob_error = NA,
      id_col = id_col, sex_col = sex_col, age_col = age_col,
      feature_cols = feature_cols,
      pre_desc = pre_desc,
      post_desc = pre_desc,
      df_imputed = df,
      saved_files = if (!is.null(out_prefix))
        list(pre = paste0(out_prefix, "_descr_preimputation.csv")) else list()
    ))
  } else if (tot_miss > 0) {
    cat("Missingness detected — continuing with imputation using missForest.\n\n")
  }
  
  # -------- Run missForest ----------------------------------------------------
  # Note: set the seed for reproducibility; verbose prints iteration details.
  set.seed(seed)
  mf <- missForest::missForest(X, verbose = verbose)
  X_imp <- mf$ximp
  rownames(X_imp) <- df[[id_col]]  # tag rows by SUBJID for reinjection later
  
  # -------- Human-readable OOB error print -----------------------------------
  # missForest returns a 2-element numeric vector:
  #   err[1] = NRMSE (numeric features)     -> lower is better
  #   err[2] = PFC   (categorical features) -> lower is better
  # If SEX_BRAIN has no missingness, PFC is often 0.
  err <- mf$OOBerror
  if (is.numeric(err) && length(err) == 2) {
    cat(sprintf("Overall OOB error:\n  Numeric (NRMSE): %.4f\n  Categorical (PFC): %.4f\n\n",
                err[1], err[2]))
    # Soft warning if NRMSE is high
    if (!is.na(err[1]) && err[1] > 0.10) {
      cat("Warning: NRMSE > 0.10; imputations may be less accurate.\n\n")
    }
  } else {
    cat("Overall OOB error: ", paste(err, collapse = ", "), "\n\n", sep = "")
  }
  
  # -------- Post-imputation descriptives (features only) ---------------------
  post_desc <- .describe_df(X_imp[feature_cols])
  saved <- list()  # will collect written file paths
  if (!is.null(out_prefix)) {
    utils::write.csv(
      post_desc,
      paste0(out_prefix, "_descr_postimputation.csv"),
      row.names = FALSE
    )
    cat("Saved post-imputation descriptives as ",
        file.path(paste0(out_prefix, "_descr_postimputation.csv")), "\n\n", sep = "")
    saved$pre  <- file.path(paste0(out_prefix, "_descr_preimputation.csv"))
    saved$post <- file.path(paste0(out_prefix, "_descr_postimputation.csv"))
  }
  
  # -------- Reinstate imputed feature values into the original df ------------
  # We overwrite only feature_cols (not AGE/SEX) so covariates remain unchanged.
  df_imp <- reinstate_imputed(
    df,
    list(X_imp = X_imp, id_col = id_col, feature_cols = feature_cols)
  )
  
  # -------- Save the imputed dataset  ------------------------------
  # (ID + features only; covariates excluded
  # now containing imputed values for feature_cols.
  if (save_imputed && !is.null(out_prefix)) {
    export_df <- df_imp[c(id_col, feature_cols)]  # df_imp already has imputed features
    out_csv <- file.path(out_dir, paste0(out_prefix, "_features_imputed.csv"))
    utils::write.csv(export_df, out_csv, row.names = FALSE)
    cat("Saved imputed dataset (ID + features only) to ", out_csv, "\n\n", sep = "")
    saved$imputed <- out_csv
  }
  
  cat("== Imputation & descriptives script complete ==\n\n")
  
  # -------- Return everything useful to the caller ---------------------------
  list(
    X_imp = X_imp,                 # imputed covariates+features matrix (rownames = IDs)
    oob_error = mf$OOBerror,       # NRMSE (numeric), PFC (categorical)
    id_col = id_col,
    sex_col = sex_col,
    age_col = age_col,
    feature_cols = feature_cols,
    pre_desc = pre_desc,
    post_desc = post_desc,
    df_imputed = df_imp,           # original df with features overwritten by imputed values
    saved_files = saved            # list of file paths we wrote (if any)
  )
}


# -----------------------------------------------------------------------------
# Reinjection helper: copy imputed values back into df by matching IDs
# - df: original data.frame (must have id_col)
# - fit: list with X_imp, id_col, feature_cols (as created just above)
# - write_cols: which columns to overwrite; default = feature_cols
# -----------------------------------------------------------------------------
reinstate_imputed <- function(df, fit, write_cols = NULL) {
  id_col <- fit$id_col
  if (is.null(write_cols)) write_cols <- fit$feature_cols
  
  # Find IDs present in both df and the imputed matrix
  common <- intersect(as.character(df[[id_col]]), rownames(fit$X_imp))
  if (!length(common)) return(df)
  
  # Row index mapping between df and X_imp for the common IDs
  i_df <- match(common, df[[id_col]])
  i_x  <- match(common, rownames(fit$X_imp))
  
  # Only write columns that actually exist in X_imp
  cols <- intersect(write_cols, colnames(fit$X_imp))
  if (!length(cols)) return(df)
  
  # Overwrite df's values with the imputed ones
  df[i_df, cols] <- fit$X_imp[i_x, cols]
  df
}


# -----------------------------------------------------------------------------
# Descriptives helper:
#   - psych::describe (means, sds, etc.; one row per variable)
#   - skimr::skim (we keep n_missing & complete_rate)
# -----------------------------------------------------------------------------
.describe_df <- function(num_df) {
  desc_psych <- psych::describe(num_df)
  desc_psych_df <- data.frame(
    variable = rownames(desc_psych),
    desc_psych,
    row.names = NULL
  )
  skim_df <- as.data.frame(skimr::skim(num_df))
  skim_keep <- unique(skim_df[c("skim_variable", "n_missing", "complete_rate")])
  merge(
    desc_psych_df,
    skim_keep,
    by.x = "variable",
    by.y = "skim_variable",
    all.x = TRUE,
    sort = FALSE
  )
}