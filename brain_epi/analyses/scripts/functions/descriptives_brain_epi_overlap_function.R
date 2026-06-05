# ============================================
#' descr_overlap
# ============================================
#' Generate joint descriptive statistics for timepoints with brain and epi data (overlap subsample)
#'
#' Merges pre-saved brain and epi age data (raw and winsorised) at a given timepoint,
#' and outputs summary stats and HTML reports.
# ============================================
#' @param timepoint Name of the overlapping timepoint, string (e.g., '15up').
#' @param log_dir Directory to save log files and summary outputs, defined as "results/descriptives"
# ============================================

descr_overlap <- function(timepoint, log_dir = "results/descriptives") {
  suffixes <- c("raw", "wins")
  
  #logging setup
  log_file <- file.path(log_dir, paste0("log_overlap_brain_epi_", timepoint, ".txt"))
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
  
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
  
  for (suffix in suffixes) {
    brain_file <- file.path("data", paste0("brain_", timepoint, "_", suffix, ".rds"))
    epi_file <- file.path("data", paste0("epi_", timepoint, "_", suffix, ".rds"))
    
    # Check if expected data files exist before attempting merge
    if (!file.exists(brain_file) || !file.exists(epi_file)) {
      warning(paste("Missing files for timepoint:", timepoint, "with suffix:", suffix))
      next
    }
    
    # Load RDS files for brain and epi data
    df_brain <- readRDS(brain_file)
    df_epi <- readRDS(epi_file)
    
    log_msg("Loaded overlap brain file: ", brain_file,
            " [", nrow(df_brain), " rows, ", ncol(df_brain), " cols]")
    log_msg("Loaded overlap epi file: ", epi_file,
            " [", nrow(df_epi), " rows, ", ncol(df_epi), " cols]")
    
    # Check for duplicated IDs
    dups_brain <- sum(duplicated(df_brain$SUBJID))
    dups_epi <- sum(duplicated(df_epi$SUBJID))
    if (dups_brain > 0) log_msg("⚠️  Brain data has ", dups_brain, " duplicated SUBJID values.")
    if (dups_epi > 0)   log_msg("⚠️  Epi data has ", dups_epi, " duplicated SUBJID values.")
    
    # Merge brain and epi data by SUBJID
    # Adds _brain and _epi suffixes to duplicated column names
    df_merged <- merge(df_brain, df_epi, by = "SUBJID", suffixes = c("_brain", "_epi"))
    
    # Difference in years: DNAm assessment age minus brain assessment age
    df_merged$AGE_diff_brain_minus_epi <- df_merged$AGE_brain - df_merged$AGE_epi
    df_merged$AGE_diff_brain_minus_epi <- as.numeric(df_merged$AGE_diff_brain_minus_epi)
    
    # Absolute age difference between MRI and DNAm assessments in years
    df_merged$AGE_diff_abs <- abs(df_merged$AGE_diff_brain_minus_epi)
    
    log_msg("Merged dataset dimensions: ", nrow(df_merged), " rows x ", ncol(df_merged), " cols")
    
    data_name <- paste0("overlap_brain_epi_", timepoint, "_", suffix)
    
    log_msg("🧠🧬 Running joint descriptives for: ", timepoint, " [", suffix, "]")
    
    # ==== Descriptives ====
    # Exclude ID column from summary
    df_summ <- df_merged[, setdiff(names(df_merged), "SUBJID")]
    
    # Generate HTML summary using summarytools
    summary_overall <- dfSummary(df_summ)
    summarytools::view(summary_overall,
                       file = paste0("results/descriptives/summarytools_", data_name, ".html"))
    
    # Separate numeric and categorical data for separate summaries
    data_numeric <- df_merged %>% select(where(is.numeric), -SUBJID)
    data_categorical <- df_merged %>% select(where(~ is.factor(.) || is.character(.)), -SUBJID)
    
    # numeric summary function
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
    
    # categorical summary function
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
    
    # Generate and save descriptive summaries
    summary_numeric <- imap_dfr(data_numeric, summarize_numeric)
    summary_categorical <- imap_dfr(data_categorical, summarize_categorical)
    
    write.csv(summary_numeric, paste0("results/descriptives/descriptives_numeric_", data_name, ".csv"), row.names = FALSE)
    write.csv(summary_categorical, paste0("results/descriptives/descriptives_categorical_", data_name, ".csv"), row.names = FALSE)
    
    # Save merged dataset for downstream usage 
    saveRDS(df_merged, file = file.path("data", paste0("overlap_brain_epi_fordescr_", timepoint, "_", suffix, ".rds")))
    
    log_msg("✅ Completed joint descriptives for: ", timepoint, " [", suffix, "]")
  }
}
