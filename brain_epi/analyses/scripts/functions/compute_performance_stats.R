# ============================================
# compute_performance_stats.R
# Module 2: Performance Evaluation Metrics with Bootstrap SEs and Boxplot Statistics
# ============================================
#' Compute model performance metrics (main and stratified)
#'
#' - Metrics include MAE, RMSE, wMAE, nRMSE, R2, Pearson, Spearman
#' - Bootstrapped SEs and  CIs for key metrics
#' - Boxplot-style summary stats for prediction error metrics (AE, PAD, PAR)
# ============================================
#' @param df Data frame containing age variables and model predictions
#' @param age_var Column name of the chronological age variable
#' @param domain Either "brain" or "epi"
#' @param timepoint Timepoint string (used in file naming)
#' @param model_vars Optional vector of prediction model names (default: all _predage columns)
#' @param output_dir Output directory for saving plots
#' @param stratify_vars Optional vector of categorical variables for stratified results
#' @param seed set seed
# ============================================
#' 10/12/2025 - MS removed BabyPy special case from function, some simplifications
#' ============================================

compute_performance_stats <- function(
    df,
    age_var,
    model_vars,
    timepoint,
    domain,
    output_dir = "results/performance",
    stratify_vars = NULL,
    seed = 2102
) {
  
  if (!is.null(seed)) set.seed(seed) #set seed 
  
  
  # Create output directory if needed
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  message("\u25b6\ufe0f Starting performance evaluation for domain: ", domain, ", timepoint: ", timepoint)
  
  # Automatically select model variables if not provided
  if (missing(model_vars) || is.null(model_vars)) {
    model_vars <- grep("_predage$", names(df), value = TRUE)
  }
  
  # Include special models if present
  special_models <- c("DunedinPACE", "DunedinPACNI")
  model_vars <- unique(c(model_vars, intersect(special_models, names(df))))
  
  # Helper:Compute bootstrapped standard error and percentile CI for a metric
  boot_se <- function(metric_fn, data, R = 1000) {
    boot_obj <- boot::boot(data, statistic = function(d, i) metric_fn(d[i, , drop = FALSE]), R = R)
    ci <- tryCatch({
      boot::boot.ci(boot_obj, type = "perc")
    }, error = function(e) NULL)
    list(
      se = sd(boot_obj$t, na.rm = TRUE),
      ci_lower = if (!is.null(ci)) ci$percent[4] else NA_real_,
      ci_upper = if (!is.null(ci)) ci$percent[5] else NA_real_
    )
  }
  
  # Helper: Compute boxplot-style statistics, used across AE, PAD, PAR
  # get_box_stats <- function(x) {
  #   x <- x[!is.na(x)]
  #   q1 <- quantile(x, 0.25)
  #   q3 <- quantile(x, 0.75)
  #   iqr <- q3 - q1
  #   # Whisker bounds (theoretical)
  #   lw_bound <- q1 - 1.5 * iqr
  #   uw_bound <- q3 + 1.5 * iqr
  #   # Whiskers based on actual data within bounds
  #   x_lower <- x[x >= lw_bound]
  #   x_upper <- x[x <= uw_bound]
  #   data.frame(
  #     mean = mean(x),
  #     sd = sd(x),
  #     median = median(x),
  #     Q1 = q1,
  #     Q3 = q3,
  #     min = min(x),
  #     max = max(x),
  #     lower_whisker <- if (length(x_lower) > 0) min(x_lower) else NA,
  #     upper_whisker <- if (length(x_upper) > 0) max(x_upper) else NA,
  #     n = length(x)
  #   )
  # }
  
  get_box_stats <- function(x, coef = 1.5) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
      return(data.frame(
        mean = NA_real_, sd = NA_real_,
        median = NA_real_, Q1 = NA_real_, Q3 = NA_real_, IQR = NA_real_,
        min = NA_real_, max = NA_real_,
        lower_whisker = NA_real_, upper_whisker = NA_real_,
        n = 0L
      ))
    }
    
    b <- boxplot.stats(x, coef = coef)
    
    data.frame(
      mean = mean(x),
      sd = sd(x),
      median = b$stats[3],
      Q1 = b$stats[2],
      Q3 = b$stats[4],
      IQR = b$stats[4] - b$stats[2],
      min = min(x),                      # raw min (may be outlier)
      max = max(x),                      # raw max (may be outlier)
      lower_whisker = b$stats[1],        # extreme non-outlier
      upper_whisker = b$stats[5],        # extreme non-outlier
      n = b$n
    )
  }
  
  
  # Main function to compute all metrics for a given model
  get_metrics <- function(data, model, age_col) {
    ref_age_col <-  age_col
    if (!all(c(model, ref_age_col) %in% names(data))) return(NULL)
    
    valid <- complete.cases(data[[model]], data[[ref_age_col]])
    df_model <- data[valid, ]
    N <- nrow(df_model)
    if (N == 0) {
      message("Skipping model ", model, " - no cases with predicted & chronological age.")
      return(NULL)
    }
    
    # === predage, PAR, and PAD ===
    base_model <- sub("_predage$", "", model)
    
    # predage
    predage_col <- if (base_model %in% c("DunedinPACE", "DunedinPACNI")) base_model else paste0(base_model, "_predage")
    predage_stats <- if (predage_col %in% names(df_model)) get_box_stats(df_model[[predage_col]]) else get_box_stats(rep(NA_real_, nrow(df_model)))
    
    # PAR
    par_col <- paste0(base_model, "_PAR")
    par_stats <- if (par_col %in% names(df_model) && !all(is.na(df_model[[par_col]]))) {
      get_box_stats(df_model[[par_col]])
    } else {
      get_box_stats(rep(NA_real_, nrow(df_model)))
    }
    
    # PAD
    pad_col <- paste0(base_model, "_PAD")
    pad <- if (pad_col %in% names(df_model)) df_model[[pad_col]] else rep(NA_real_, nrow(df_model))
    pad_stats <- get_box_stats(pad)
    
    # === AE ===
    ae <- abs(pad)
    ae_stats <- get_box_stats(ae)
    
    # === wAE ===
    age_range_col <- age_col
    range_age <- if (age_range_col %in% names(df_model)) {
      max(df_model[[age_range_col]], na.rm = TRUE) - min(df_model[[age_range_col]], na.rm = TRUE)
    } else {
      NA_real_
    }
    wae <- ae / range_age
    wae_stats <- get_box_stats(wae)
    
    # === Standard metric functions ===
    mae_fn <- function(x) mean(abs(x[[model]] - x[[ref_age_col]]))
    rmse_fn <- function(x) sqrt(mean((x[[model]] - x[[ref_age_col]])^2))
    wmae_fn <- function(x) {
      denom <- max(x[[ref_age_col]], na.rm = TRUE) - min(x[[ref_age_col]], na.rm = TRUE)
      if (denom == 0) return(NA_real_)
      mae_fn(x) / denom 
    } #more complex wmae_fn function to ensure no errors even if max-min = 0 in resamples 
    nrmse_fn <- function(x) {
      denom <- max(x[[ref_age_col]], na.rm = TRUE) - min(x[[ref_age_col]], na.rm = TRUE)
      if (denom == 0) return(NA_real_)
      rmse_fn(x) / denom
    } #more complex nrmse_fn function to ensure no errors even if max-min = 0 in resamples
    rae_fn <- function(x) mean(abs(x[[model]] - x[[ref_age_col]])) / mean(abs(mean(x[[ref_age_col]]) - x[[ref_age_col]]))
    
    # === Bootstrapped SEs ===
    mae_boot <- boot_se(mae_fn, df_model)
    rmse_boot <- boot_se(rmse_fn, df_model)
    wmae_boot <- boot_se(wmae_fn, df_model)
    nrmse_boot <- boot_se(nrmse_fn, df_model)
    rae_boot <- boot_se(rae_fn, df_model)
    
    # === Standard metrics ===
    df_out <- data.frame(
      model = model,
      N = N,
      MAE = round(mae_fn(df_model), 5),
      MAE_SE_boot = round(mae_boot$se, 5),
      RMSE = round(rmse_fn(df_model), 5),
      RMSE_SE_boot = round(rmse_boot$se, 5),
      wMAE_test = round(wmae_fn(df_model), 5),
      wMAE_SE_boot = round(wmae_boot$se, 5),
      nRMSE = round(nrmse_fn(df_model), 5),
      nRMSE_SE_boot = round(nrmse_boot$se, 5),
      RAE = round(rae_fn(df_model), 5),
      RAE_SE_boot = round(rae_boot$se, 5),
      R2 = round(caret::R2(df_model[[model]], df_model[[ref_age_col]]), 5),
      Pearson = round(cor(df_model[[model]], df_model[[ref_age_col]], method = "pearson"), 5),
      Spearman = round(cor(df_model[[model]], df_model[[ref_age_col]], method = "spearman"), 5)
    )
    
    # === Combine all stats ===
    stats_list <- list(
      setNames(ae_stats, paste0("AE_", names(ae_stats))),
      setNames(wae_stats, paste0("wAE_", names(wae_stats))),
      setNames(pad_stats, paste0("PAD_", names(pad_stats))),
      setNames(predage_stats, paste0("predage_", names(predage_stats))),
      setNames(par_stats, paste0("PAR_", names(par_stats)))
    )
    
    do.call(cbind, c(list(df_out), stats_list))
  }
  
  # === MAIN METRICS ===
  message("\ud83d\udcca Evaluating ", length(model_vars), " models...")
  main_results <- do.call(rbind, lapply(model_vars, function(m) get_metrics(df, m, age_var)))
  write.csv(
    main_results,
    file.path(output_dir, paste0("performance_", domain, "_", timepoint, "_main.csv")),
    row.names = FALSE
  )
  
  # === STRATIFIED METRICS ===
  if (!is.null(stratify_vars)) {
    for (strat in stratify_vars) {
      message("\ud83d\udcc1 Stratifying by: ", strat)
      stratified_list <- split(df, df[[strat]])
      stratified_results <- do.call(rbind, lapply(names(stratified_list), function(level) {
        df_level <- stratified_list[[level]]
        df_metrics <- do.call(rbind, lapply(model_vars, function(m) get_metrics(df_level, m, age_var)))
        if (!is.null(df_metrics)) {
          df_metrics$stratified_by <- strat
          df_metrics$level <- level
          df_metrics[, c("stratified_by", "level", names(df_metrics)[1:(ncol(df_metrics) - 2)])]
        } else {
          NULL
        }
      }))
      
      write.csv(
        stratified_results,
        file.path(output_dir, paste0("performance_", domain, "_", timepoint, "_stratified_by_", strat, ".csv")),
        row.names = FALSE
      )
    }
  }
  
  message("\u2705 Performance stats computed for ", domain, " at timepoint: ", timepoint)
}
