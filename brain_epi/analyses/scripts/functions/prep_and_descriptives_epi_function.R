# ============================================
#' prep_and_descr_timepoint_epi
# ============================================
#' Preprocess and generate descriptives for epigenetic age data at a given timepoint.
#'
#' - Loads data from the configured path.
#' - Coerces variable types (e.g., sets SEX to factor).
#' - Excludes fully NA records for key covariates.
#' - Saves raw data.
#' - Identifies variables to winsorise (based on config exclusions and missingness).
#' - Performs winsorisation (3*IQR rule).
#' - Outputs descriptive summaries pre- and post-winsorisation (CSV and HTML).
#' - Saves winsorised dataset.
# ============================================
#' @param timepoint A string specifying the timepoint.
#' @param epi_vars A character vector of epigenetic clock variables to summarise and winsorise, defined in the main R script
#' @param config_list A named list of configuration details (data paths, exclusions).
#' @param log_dir Directory where logs and output summaries will be saved, currently provided as "results/descriptives"
# ============================================

prep_and_descr_timepoint_epi <- function(timepoint, epi_vars, config_list, log_dir = "results/descriptives") {
  
  config <- config_list[[timepoint]]
  data_path <- config$file
  exclude_vars <- config$exclude_vars
  
  data_name <- paste0("epi_", timepoint)
  log_file <- file.path(log_dir, paste0("log_", data_name, ".txt"))
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
  
  data_dir <- "data" #to save rds files 
  dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Open log connection
  log_con <- file(log_file, open = "wt")
  
  # Ensure it gets closed on exit
  on.exit(close(log_con), add = TRUE)
  
  # Simple helper that logs AND prints to console
  log_msg <- function(...) {
    txt <- paste0(..., collapse = "")
    # write to log
    writeLines(txt, con = log_con)
    # echo to console
    message(txt)
  }
  
  log_msg("🟢 Starting processing for timepoint: ", timepoint)
  
  # Load data
  df <- read.csv(data_path)
  
  log_msg("Loaded epi data from: ", data_path)
  log_msg("Initial dimensions: ", nrow(df), " rows x ", ncol(df), " columns")
  log_msg("Column names: ", paste(names(df), collapse = ", "))
  
  # Check for some key columns
  required_cols <- c("SUBJID", "SEX_epi", "AGE_epi")
  missing_required <- setdiff(required_cols, names(df))
  if (length(missing_required) > 0) {
    stop("Missing required covariates in epi data: ",
            paste(missing_required, collapse = ", "))
  }
  
  # Remove ALSPAC-specific column if present
  if ("Sample_Name" %in% names(df)) df$Sample_Name <- NULL
  
  # Convert variables to appropriate types for analysis
  # If a sample includes additional categorical variables, please ensure that they are included here and correctly set to factor in this section
  df$SUBJID <- as.character(df$SUBJID)
  df$SEX_epi <- as.factor(df$SEX_epi) #not included in potential categorical covariates as it has to be present
  
  # List of potential categorical variables in the brain data
  cat_candidates <- c("batch") #add categorical variables here, remove batch if not numeric
  
  # Only keep those that actually exist in this cohort's data
  factor_vars <- intersect(cat_candidates, names(df))
  
  # Convert all of them to factors
  df[factor_vars] <- lapply(df[factor_vars], as.factor)
  
  # Protect ID + all factor vars from numeric conversion
  protected_cols <- c("SUBJID", "SEX_epi", factor_vars)
  
  # Convert other columns to numeric
  cols_to_convert <- setdiff(names(df), protected_cols)

  # Track NAs before conversion
  na_before <- sapply(df[cols_to_convert], function(x) sum(is.na(x)))
  
  df[cols_to_convert] <- lapply(df[cols_to_convert], function(x) as.numeric(as.character(x)))
  
  # Track NAs after conversion
  na_after <- sapply(df[cols_to_convert], function(x) sum(is.na(x)))
  
  # Report any columns where NAs increased (likely due to coercion issues)
  na_increase <- na_after - na_before
  problem_cols <- names(na_increase[na_increase > 0])
  
  if (length(problem_cols) > 0) {
    log_msg("⚠️  Numeric conversion introduced additional NAs in columns: ",
            paste(problem_cols, collapse = ", "))
    conv_details <- data.frame(
      variable = problem_cols,
      NAs_before = na_before[problem_cols],
      NAs_after = na_after[problem_cols],
      stringsAsFactors = FALSE
    )
    print(conv_details)
  }
  
  log_msg("Variable classes after type conversion:")
  var_classes <- sapply(df, function(x) paste(class(x), collapse = ", "))
  var_lines <- paste(names(var_classes), var_classes, sep = ": ") #format nicely 
  # Write to log AND console
  for (ln in var_lines) {
    log_msg("  ", ln)
  }
  
  # Reorder vars for clarity
  reorder_cols <- intersect(c("SUBJID", "SEX_epi", "AGE_epi", "AGE_gest"), names(df))
  df <- df[, c(reorder_cols, setdiff(names(df), reorder_cols))] #only required variables are considered
  
  # Track missing data removal
  n_total <- nrow(df)
  n_missing_age <- sum(is.na(df$AGE_epi))
  n_missing_sex <- sum(is.na(df$SEX_epi))
    
  # Filter out missing AGE_epi and SEX_epi
  df <- df[!is.na(df$AGE_epi) & !is.na(df$SEX_epi), ]
  n_remaining <- nrow(df)
  
  # Log how many were removed
  log_msg(sprintf("Removed %d participants with missing AGE_epi", n_missing_age))
  log_msg(sprintf("Removed %d participants with missing SEX_epi", n_missing_sex))
  log_msg(sprintf("Remaining participants after after removing those with missing SEX or AGE: %d (from original %d)", n_remaining, n_total))
  
  # Save raw data
  saveRDS(df, file = paste0(data_dir, "/", data_name, "_raw.rds"))
  
  # Warn about missing variables
  missing_vars <- setdiff(setdiff(epi_vars, config$exclude_vars), names(df))
  if (length(missing_vars) > 0) {
    warning("The following expected variables are missing from the dataset: ", 
            paste(missing_vars, collapse = ", "))
  }
  
  # Determine variables to winsorise
  vars_in_data <- intersect(epi_vars, names(df))
  vars_to_check <- setdiff(vars_in_data, exclude_vars)
  fully_missing_vars <- vars_to_check[colSums(!is.na(df[, vars_to_check, drop = FALSE])) == 0]
  vars_to_winsor <- setdiff(vars_to_check, fully_missing_vars)
  
  log_msg("Excluding fully missing variables: ", paste(fully_missing_vars, collapse = ", "))
  log_msg("Excluding configured variables: ", paste(exclude_vars, collapse = ", "))
  log_msg("Variables to winsorise: ", paste(vars_to_winsor, collapse = ", "))
  
  # ==== Descriptive stats function ====
 
  dir.create("results/descriptives", showWarnings = FALSE, recursive = TRUE)
  
  # Function for summarizing numeric variables
  summarize_numeric <- function(x, name) {
    tibble(
      variable = name,
      N = sum(!is.na(x)),
      N_missing = sum(is.na(x)),
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      min = min(x, na.rm = TRUE),
      max = max(x, na.rm = TRUE)
    )
  }
  
  # Function for summarizing categorical variables
  summarize_categorical <- function(x, name) {
    freq_tbl <- table(x, useNA = "ifany")
    prop_tbl <- prop.table(freq_tbl)
    
    tibble(
      variable = name,
      level = names(freq_tbl),
      count = as.integer(freq_tbl),
      percent = round(100 * as.numeric(prop_tbl), 1)
    )
  }
  
   generate_descriptives <- function(df_input, suffix_label) {
    df_summ <- df_input[, setdiff(names(df_input), "SUBJID")]
    
    # Summarytools HTML summary
    summary_overall <- dfSummary(df_summ)
    summarytools::view(summary_overall,
                       file = paste0("results/descriptives/summarytools_", data_name, "_", suffix_label, ".html"))
    
    # Numeric and categorical
    data_numeric <- df_input %>% select(where(is.numeric), -SUBJID)
    data_categorical <- df_input %>% select(where(~ is.factor(.) || is.character(.)), -SUBJID)
    
    summary_numeric <- imap_dfr(data_numeric, summarize_numeric)
    summary_categorical <- imap_dfr(data_categorical, summarize_categorical)
    
    write.csv(summary_numeric, paste0("results/descriptives/descriptives_numeric_", data_name, "_", suffix_label, ".csv"), row.names = FALSE)
    write.csv(summary_categorical, paste0("results/descriptives/descriptives_categorical_", data_name, "_", suffix_label, ".csv"), row.names = FALSE)
  }
  
  # ==== Descriptives: pre-winsorisation ====
  generate_descriptives(df, "raw")
  
  # ==== Winsorisation ====
   # Winsorise: Cap extreme values beyond 3*IQR (based on Q1 & Q3) for selected epi clocks
   # Adapted from Vilte Baltramonaityte's methylation scripts
   
  log_msg("Winsorizing extreme values column-wise...")
  df_wins <- df
  num_values_replaced <- 0
  winsorized_vars <- character()
  
  for (var in vars_to_winsor) {
    x <- df_wins[[var]]
    if (!is.numeric(x) || all(is.na(x))) next
    
    q1 <- quantile(x, 0.25, na.rm = TRUE)
    q3 <- quantile(x, 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    lower <- q1 - 3 * iqr
    upper <- q3 + 3 * iqr
    
    winsorized <- pmin(pmax(x, lower), upper)
    changed <- sum(winsorized != x, na.rm = TRUE)
    
    if (changed > 0) {
      df_wins[[var]] <- winsorized
      winsorized_vars <- c(winsorized_vars, var)
      num_values_replaced <- num_values_replaced + changed
    }
  }
  
  total_values <- length(vars_to_winsor) * nrow(df_wins)
  if (total_values > 0) {
    pct_values_replaced <- 100 * num_values_replaced / total_values
  } else {
    pct_values_replaced <- NA_real_ #ensures that we do not divide by 0
  }
  
  log_msg(sprintf("Winsorised %d values (%.4f%% of selected data)", num_values_replaced, pct_values_replaced))
  if (length(winsorized_vars) > 0) {
    log_msg("Winsorised variables: ", paste(winsorized_vars, collapse = ", "))
  } else {
    log_msg("No variables required winsorisation.")
  }
  
  # ==== Descriptives: post-winsorisation ====
  generate_descriptives(df_wins, "wins")
  
  # Save winsorised dataset
  saveRDS(df_wins, file = paste0(data_dir, "/", data_name, "_wins.rds"))
  
  log_msg("✅ Finished processing timepoint: ", timepoint)
}
