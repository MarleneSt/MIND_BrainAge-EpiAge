# ============================================
# association_helpers.R
# ============================================
# includes a range of helper functions for the associaton analyses
#   - check_model_availability: checks model availability & produces log 
#   - safe_scale: function to safely scale numeric model terms 
#   - add_unsuffixed: add Dunedin models that don't have the expected suffixes 
#   - create_model_formula: create formula for association models
#   - run_lm_model_with_full_stats: run lm models 
#   - confint_rlm: function to get CIs for RLM estimates
#   - run_rlm_model: run rlm models
#   - save_combined_lm_diagnostics: save lm model diagnostics from plot()
# ============================================

# function to check which biological age models are available 
check_model_availability_agediff <- function(data, site_config) {
  suffixes <- c("_predage", "_PAR", "_PAD")
  
  brain_models <- c("Centile2", "DevBrainAge", "DBN", "ENIGMA", "Kaufmann", 
                    "PyBrainAge", "Pyment", "DunedinPACNI")
  epi_models <- c("PCBrainAge", "PCGrimAge", "CorticalClock", "cAge", "AltumAge", "Horvath2013",
                  "DamAge", "AdaptAge", "DunedinPACE", "Hannum", "PhenoAge", "skinHorvath",
                  "PedBE", "Wu", "DNAmTL", "ZhangBLUP", "ZhangEN", "Bohlin", "EPIC","Knight")
  
  # Only keep columns that are not entirely NA
  non_na_vars <- names(data)[colSums(!is.na(data)) > 0]
  
  # Clean setup
  log_lines <- c()
  log_file <- paste0("results/associations/log_modelavailability_agediff_", site_config$sample_name, "_", site_config$timepoint_var, ".txt")
  dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  
  # Find all base names present in the data (ignoring suffix)
  brain_present <- intersect(brain_models, gsub("(_predage|_PAR|_PAD)$", "", non_na_vars))
  epi_present   <- intersect(epi_models,   gsub("(_predage|_PAR|_PAD)$", "", non_na_vars))
  
  # Identify completely missing models
  missing_brain_models <- setdiff(brain_models, brain_present)
  missing_epi_models   <- setdiff(epi_models, epi_present)
  
  log_lines <- c(log_lines, paste0("📌 Sample: ", site_config$sample_name, " | Timepoint: ", site_config$timepoint_var), "")
  log_lines <- c(log_lines, sprintf("✅ %d brain age models detected:", length(brain_present)))
  log_lines <- c(log_lines, paste0(" - ", sort(brain_present)), "")
  log_lines <- c(log_lines, sprintf("✅ %d epigenetic age models detected:", length(epi_present)))
  log_lines <- c(log_lines, paste0(" - ", sort(epi_present)), "")
  
  # Log completely missing models
  if (length(missing_brain_models) > 0) {
    log_lines <- c(log_lines, "❌ The following brain models are completely missing (no suffixes found):")
    log_lines <- c(log_lines, paste0(" - ", sort(missing_brain_models)), "")
  }
  
  if (length(missing_epi_models) > 0) {
    log_lines <- c(log_lines, "❌ The following epigenetic models are completely missing (no suffixes found):")
    log_lines <- c(log_lines, paste0(" - ", sort(missing_epi_models)), "")
    log_lines <- c(log_lines, "If this is a non-birth time point, Bohlin, EPIC and Knight should be missing.\n")
  }
  
  # Function to check if all suffixes are available
  has_all_suffixes <- function(model, vars) {
    all(paste0(model, suffixes) %in% vars)
  }
  
  # Keep only models with all suffixes (except for DunedinPAC(N)I)
  brain_full <- unlist(lapply(brain_present, function(model) {
    if (model == "DunedinPACNI") return(model)
    if (has_all_suffixes(model, non_na_vars)) return(model)
    return(NULL)
  }))
  
  epi_full <- unlist(lapply(epi_present, function(model) {
    if (model == "DunedinPACE") return(model)
    if (has_all_suffixes(model, non_na_vars)) return(model)
    return(NULL)
  }))
  
  # Log what's missing
  log_lines <- c(log_lines, "🔍 Checking for missing suffix variants amongst available models:")
  missing_brain <- setdiff(brain_present, brain_full)
  missing_epi <- setdiff(epi_present, epi_full)
  
  if (length(missing_brain) > 0) {
    log_lines <- c(log_lines, "⚠️ Brain models missing full suffix set:")
    log_lines <- c(log_lines, paste0(" - ", missing_brain), "")
    log_lines <- c(log_lines, "❗Please check these models as if a model is available, all versions (predage, PAR, PAD) should be available.")
  } else {
    log_lines <- c(log_lines, "✅ All  available brain age models have all suffixes (_predage, _PAR, _PAD).", "")
  }
  
  if (length(missing_epi) > 0) {
    log_lines <- c(log_lines, "⚠️ Epigenetic models missing full suffix set:")
    log_lines <- c(log_lines, paste0(" - ", missing_epi), "")
    log_lines <- c(log_lines, "❗ Please check these models as if a model is available, all versions (predage, PAR, PAD) should be available.")
  } else {
    log_lines <- c(log_lines, "✅ All  available epigenetic age models have all suffixes (_predage, _PAR, _PAD).", "")
  }
  
  # outputting age overlap information
  
  # Requires AGE_brain and AGE_epi to exist in the merged overlap dataset.
  if (all(c("AGE_brain", "AGE_epi") %in% names(data))) {
    
    valid_age <- complete.cases(data$AGE_brain, data$AGE_epi)
    n_pair <- sum(valid_age)
    
    if (n_pair > 0) {
      age_diff <- data$AGE_brain[valid_age] - data$AGE_epi[valid_age]  # brain minus epi
      
      mean_diff <- mean(age_diff)
      sd_diff   <- sd(age_diff)
      mean_abs  <- mean(abs(age_diff))
      sd_abs    <- sd(abs(age_diff))
      
      r_pearson  <- suppressWarnings(cor(data$AGE_brain[valid_age], data$AGE_epi[valid_age], method = "pearson"))
      r_spearman <- suppressWarnings(cor(data$AGE_brain[valid_age], data$AGE_epi[valid_age], method = "spearman"))
      
      log_lines <- c(log_lines,
                     "⏱️ MRI vs DNAm assessment age alignment (overlap sample):",
                     paste0(" - N with both AGE_brain & AGE_epi: ", n_pair),
                     paste0(" - AGE_brain - AGE_epi: mean = ", round(mean_diff, 3), ", SD = ", round(sd_diff, 3)),
                     paste0(" - |AGE_brain - AGE_epi|: mean = ", round(mean_abs, 3), ", SD = ", round(sd_abs, 3)),
                     paste0(" - Correlation (Pearson):  r = ", round(r_pearson, 3)),
                     paste0(" - Correlation (Spearman): ρ = ", round(r_spearman, 3)),
                     ""
      )
    } else {
      log_lines <- c(log_lines,
                     "⚠️ MRI vs DNAm assessment age alignment (overlap sample):",
                     " - AGE_brain and AGE_epi exist, but there are 0 complete pairs (all missing) - please check and ensure these variables are available.",
                     ""
      )
    }
    
  } else {
    log_lines <- c(log_lines,
                   "⚠️️ MRI vs DNAm assessment age alignment (overlap sample):",
                   " - Skipped as AGE_brain and/or AGE_epi not found in dataset - please check and ensure these variables are available.",
                   ""
    )
  }
  
  # Save log
  cat(paste(log_lines, collapse = "\n"))
  writeLines(log_lines, log_file)
  
  # Return actual variables for ages 
  brain_vars <- unlist(lapply(brain_full, function(model) {
    if (model == "DunedinPACNI") return("DunedinPACNI")
    paste0(model, suffixes)
  }))
  
  epi_vars <- unlist(lapply(epi_full, function(model) {
    if (model == "DunedinPACE") return("DunedinPACE")
    paste0(model, suffixes)
  }))

  
  message("\n📝 Model availability and age difference log written to: ", log_file)
  return(list(
    brain_vars = brain_vars,
    epi_vars = epi_vars
  ))
}

# helper to scale all numeric variables in a df 
safe_scale <- function(x) {
  if (is.numeric(x)) return(scale(x)) else return(x)
}

# helper to add the Dunedin models which might otherwise be excluded when selecting based on suffix
add_unsuffixed <- function(vars, type, data, is_brain = FALSE) {
  extra <- if (is_brain) "DunedinPACNI" else "DunedinPACE"
  
  if (type %in% c("_predage", "_PAR", "_PAD") && extra %in% names(data)) {
    return(unique(c(vars, extra)))
  } else {
    return(vars)
  }
}

# helper to create model formula (for models 1-3)
create_model_formula <- function(outcome, predictor, model_num, config, data, var_type) {
  age_brain <- "AGE_brain"
  age_epi <- "AGE_epi"
  sex <- "SEX_brain"
  covars <- c()
  
  if (model_num >= 2) {
    if (cor(data[[age_brain]], data[[age_epi]], use="complete.obs") >= 0.8) {
      covars <- c(covars, age_brain)
    } else {
      covars <- c(covars, age_brain, age_epi)
    }
    covars <- c(covars, sex, config$batch_vars)
    
    # Only add extra covariates for model 2+
    if (!is.null(config$extra_covariates)) {
      covars <- c(covars, config$extra_covariates)
    }
  }
  
  if (model_num == 3) {
    covars <- c(covars, config$cell_type_vars)
  }
  
  rhs <- paste(c(predictor, covars), collapse = " + ")
  formula(paste(outcome, "~", rhs))
}

# function to run linear model - lm() - non-scaled & scaled 
run_lm_model_with_full_stats <- function(formula, data, model_string) {
  outcome_var <- as.character(formula[[2]])
  fit <- lm(formula, data = data)
  summ <- summary(fit)
  coefs_raw <- summ$coefficients
  coefs_hc3 <- coeftest(fit, vcov = vcovHC(fit, type = "HC3"))
  ci_raw <- confint(fit)
  r2 <- summ$r.squared
  adj_r2 <- summ$adj.r.squared
  n <- nobs(fit)
  # complicated way of getting VIF: 
      #1) Only calculate VIF if model has 2+ predictors (excluding intercept)
      #2) need to ensure that we're using VIF unless there's a categorical 
          # var with 3+ levels --> GVIF^(1/(2*Df)
          # GVIF^(1/(2*Df) has to be squared for (similar interpretation) 
          # see e.g. https://stacyderuiter.github.io/s245-notes-bookdown/collinearity-and-multicollinearity.html
          # ~1 is ideal, >5 or >10 suggests collinearity
  X <- model.matrix(fit)
  if (ncol(X) > 2) {
    vif_vals_full <- vif(fit)
    
    # Extract appropriate VIFs
    if (is.matrix(vif_vals_full)) {
      vif_vals_full <- (vif_vals_full[, "GVIF^(1/(2*Df))"])^2
    }
    
    # Pre-fill with NAs
    vif_safe <- rep(NA_real_, nrow(coefs_raw))
    coef_names <- rownames(coefs_raw)
    
    # Match VIF values to coefficient names
    for (i in seq_along(coef_names)) {
      name <- coef_names[i]
      
      # Exact match first
      if (name %in% names(vif_vals_full)) {
        vif_safe[i] <- vif_vals_full[name]
      } else {
        # Try partial match for dummy variables
        match_idx <- which(sapply(names(vif_vals_full), function(v) startsWith(name, v)))
        if (length(match_idx) == 1) {
          vif_safe[i] <- vif_vals_full[match_idx]
        }
      }
    }
  } else {
    vif_safe <- rep(NA_real_, nrow(coefs_raw))
  }
  
  #scaled model
  data_scaled <- data %>% mutate(across(where(is.numeric), scale))
  fit_scaled <- lm(formula, data = data_scaled)
  summ_scaled <- summary(fit_scaled)
  coefs_scaled<- summ_scaled$coefficients
  coefs_hc3_scaled <- coeftest(fit_scaled, vcov = vcovHC(fit_scaled, type = "HC3"))
  ci_scaled <- confint(fit_scaled)

  results <- tibble(
    Predictor = rownames(coefs_raw),
    B = coefs_raw[, 1],
    SE = coefs_raw[, 2],
    CI_lower = ci_raw[, 1],
    CI_upper = ci_raw[, 2],
    T_value = coefs_raw[, 3],
    P_value = coefs_raw[, 4],
    HC3_SE = coefs_hc3[, 2],
    HC3_t = coefs_hc3[, 3],
    HC3_p = coefs_hc3[, 4],
    beta_standardised = coefs_scaled[, 1],
    SE_standardised = coefs_scaled[, 2],
    beta_CI_lower = ci_scaled[, 1],
    beta_CI_upper = ci_scaled[, 2],
    T_value_standardised = coefs_scaled[, 3],
    P_value_standardised = coefs_scaled[, 4],
    HC3_SE_standardised = coefs_hc3_scaled[, 2],
    HC3_t_standardised = coefs_hc3_scaled[, 3],
    HC3_p_standardised = coefs_hc3_scaled[, 4],
    R_squared = r2,
    Adj_R_squared = adj_r2,
    N = n,
    VIF_GVIF = vif_safe #will be VIF if no categorical predictor with more than 2 levels is included, otherwise squared GVIF^(1/(2*Df))
  )

  return(list(results = results, fit = fit))
}

# function to estimate CIs for RLM
confint_rlm <- function(coefs, level = 0.95) {
  z <- qnorm((1 + level) / 2)
  est <- coefs[, 1]
  se  <- coefs[, 2]
  
  lower <- est - z * se
  upper <- est + z * se
  
  out <- cbind(lower, upper)
  rownames(out) <- rownames(coefs)
  colnames(out) <- c("CI_lower", "CI_upper")
  return(out)
}

# function to run rlm models - rlm() based on MASS package, scaled & non-scaled
run_rlm_model <- function(formula, data) {
  
  # model without scaling 
  fit <- rlm(formula, data = data, maxit = 200)
  
  #not scaled: rlm natural canonical SEs & CIs (but rlm does not naturally produce p-values)
  summ <- summary(fit)
  coefs_raw <- summ$coefficients
  ci_raw <- confint_rlm(coefs_raw) #compute CIs
  
  #not scaled: coeftest to produce HC consistent SEs
  coefs <- coeftest(fit, vcovHC(fit, type = "HC3"))
  ci <- confint_rlm(coefs) #compute CIs
  
  # model with scaling 
  data_scaled <- data %>% mutate(across(where(is.numeric), scale))
  fit_scaled <- rlm(formula, data = data_scaled, maxit = 200)
  
  #scaled: rlm natural canonical SEs & CIs (but rlm does not naturally produce p-values)
  summ_scaled <- summary(fit_scaled)
  coefs_raw_scaled <- summ_scaled$coefficients
  ci_raw_scaled <- confint_rlm(coefs_raw_scaled) #compute CIs
  
  #scaled: coeftest to produce HC consistent SEs
  coefs_scaled <- coeftest(fit_scaled, vcovHC(fit_scaled, type = "HC3"))
  ci_scaled <- confint_rlm(coefs_scaled) #compute CIs
  
  results <- tibble(
    Predictor = rownames(coefs),
    
    #rlm unscaled data 
    RLM_Estimate = coefs[, 1], #estimate is the same in coefs_raw & coefs
    RLM_SE_raw = coefs_raw[, 2],
    RLM_CI_lower_raw = ci_raw[, 1],
    RLM_CI_upper_raw = ci_raw[, 2],
    RLM_t_raw = coefs_raw[, 3],
    RLM_SE_coeftest_HC3 = coefs[, 2],
    RLM_CI_lower_coeftest = ci[, 1],
    RLM_CI_upper_coeftest = ci[, 2],
    RLM_z_coeftest = coefs[, 3],
    RLM_p_coeftest = coefs[, 4],
    
    #rlm scaled data 
    RLM_Estimate_scaled = coefs_scaled[, 1], #estimate is the same in coefs_raw & coefs
    RLM_SE_raw_scaled = coefs_raw_scaled[, 2],
    RLM_CI_lower_raw_scaled = ci_raw_scaled[, 1],
    RLM_CI_upper_raw_scaled = ci_raw_scaled[, 2],
    RLM_t_raw_scaled = coefs_raw_scaled[, 3],
    RLM_SE_coeftest_HC3_scaled = coefs_scaled[, 2],
    RLM_CI_lower_coeftest_scaled = ci_scaled[, 1],
    RLM_CI_upper_coeftest_scaled = ci_scaled[, 2],
    RLM_z_coeftest_scaled = coefs_scaled[, 3],
    RLM_p_coeftest_scaled = coefs_scaled[, 4],
    
    #n
    RLM_N = nobs(fit)
  )

  return(results)
}

# lm diagnostics disabled - would produce diagnostic plots for the lm model
# save_combined_lm_diagnostics <- function(lm_models, filename) {
#   if (length(lm_models) == 0) return(NULL)  # early exit if empty
#   
#   dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
#   pdf(filename)
#   
#   for (m in lm_models) {
#     par(mfrow = c(2, 2))
#     plot(m$fit, main = paste(m$brain, "~", m$epi))
#   }
#   
#   dev.off()
# }
