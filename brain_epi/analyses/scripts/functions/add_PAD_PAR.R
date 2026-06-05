# ============================================
# add_PAD_PAR.R
# Module 1: Create PAD and PAR variables 
# ============================================
#' Add PAD and PAR metrics to prediction models
#'
#' - PAD = Predicted Age Difference (predicted - chronological)
#' - PAR = Prediction residual (from linear model)
#'
# ============================================
#' @param df Data frame containing age variables and model predictions
#' @param age_var Column name of the chronological age variable
# ============================================
#' 10/12/2025 - MS removed BabyPy special case from function
#' ============================================

add_PAD_PAR <- function(df, age_var) {
  predage_vars <- grep("_predage$", names(df), value = TRUE)
  
  for (var in predage_vars) {
    base_var <- gsub("_predage$", "", var)
    
    # Define reference age
    ref_age <- age_var #removed special cases for BabyPy (MS 10/12/2025)
    
    # Check required columns exist
    if (!all(c(var, ref_age) %in% names(df))) {
      warning(paste("Missing columns:", var, "or", ref_age, "- skipping."))
      next
    }
    
    # Check if there are any non-NA values
    valid_rows <- complete.cases(df[[var]], df[[ref_age]])
    if (sum(valid_rows) == 0) {
      message(paste("Skipping model", var, "- all values are NA."))
      next
    }
    
    # Compute PAD
    pad_name <- paste0(base_var, "_PAD")
    df[[pad_name]] <- df[[var]] - df[[ref_age]]
    
    # Compute PAR (residuals of lm)
    par_name <- paste0(base_var, "_PAR")
    df[[par_name]] <- NA
    df[[par_name]][valid_rows] <- residuals(lm(df[[var]][valid_rows] ~ df[[ref_age]][valid_rows]))
  }
  
  return(df)
}
