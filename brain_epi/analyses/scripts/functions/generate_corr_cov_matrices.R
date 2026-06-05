# ============================================
# generate_corr_cov_matrices.R
# Module 3: Correlation & Covariance Matrices
# ============================================
#' Generate correlation and covariance matrices
#'
#' - Produces Pearson and Spearman correlation matrices (wide + long format)
#' - Produces covariance matrices
#' - Handles predicted ages, PAD, and PAR sets separately
# ============================================
#' @param df Data frame with model predictions and PAD/PAR metrics
#' @param age_var Column name of the chronological age variable
#' @param domain Either "brain" or "epi"
#' @param timepoint Timepoint string (used in file naming)
#' @param output_dir Output directory for saving plots
#' # ============================================
#' 10/12/2025 - MS removed BabyPy special case from function, some simplifications
#' ============================================

# ==== Main Function ====
generate_corr_cov_matrices <- function(
    df,
    age_var,
    domain,
    timepoint,
    output_dir = "results/performance/corr_cov"
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  message("▶️ Generating correlation/covariance matrices for domain: ", domain, ", timepoint: ", timepoint)
  
  # Filter to relevant columns
  special_models <- c("DunedinPACE", "DunedinPACNI") #they will be included in predage, PAD, and PAR matrices 
  
  pred_vars <- unique(c(grep("_predage$", names(df), value = TRUE), intersect(special_models, names(df))))
  pad_vars  <- unique(c(grep("_PAD$", names(df), value = TRUE), intersect(special_models, names(df))))
  par_vars  <- unique(c(grep("_PAR$", names(df), value = TRUE), intersect(special_models, names(df))))
  
  age_related <- age_var
  
  # Define matrix groupings: predicted age, PAD, PAR (incl. special models if present)
  matrix_sets <- list(
    predage = c(age_related, pred_vars),
    PAD = c(age_related, pad_vars),
    PAR = c(age_related, par_vars)
  )
  
  # Compute Pearson or Spearman correlation matrix with p-values and N
  compute_corrs <- function(data, vars, method = "pearson") {
    valid <- data[, vars[vars %in% names(data)], drop = FALSE]
    valid <- valid[, colSums(!is.na(valid)) > 1]
    cor_res <- Hmisc::rcorr(as.matrix(valid), type = method)
    list(cor = cor_res$r, p = cor_res$P, n = cor_res$n)
  }
  
  compute_covs <- function(data, vars) {
    valid <- data[, vars[vars %in% names(data)], drop = FALSE]
    valid <- valid[, colSums(!is.na(valid)) > 1]
    cov(valid, use = "pairwise.complete.obs")
  }
  
  for (setname in names(matrix_sets)) {
    vars <- matrix_sets[[setname]]
    
    pearson <- compute_corrs(df, vars, method = "pearson")
    spearman <- compute_corrs(df, vars, method = "spearman")
    covar <- compute_covs(df, vars)
    
    # Write wide-format correlation matrices
    write.csv(round(pearson$cor, 10), file.path(output_dir, paste0(domain, "_", timepoint, "_correlations_", setname, "_pearson.csv")))
    write.csv(round(spearman$cor, 10), file.path(output_dir, paste0(domain, "_", timepoint, "_correlations_", setname, "_spearman.csv")))
    write.csv(round(covar, 10), file.path(output_dir, paste0(domain, "_", timepoint, "_covariance_", setname, ".csv")))
    
    # Long-format outputs for Pearson
    cor_df <- as.data.frame(as.table(pearson$cor))
    names(cor_df) <- c("Variable1", "Variable2", "Correlation")
    pval_df <- as.data.frame(as.table(pearson$p))
    names(pval_df) <- c("Variable1", "Variable2", "P_value")
    n_df <- as.data.frame(as.table(pearson$n))
    names(n_df) <- c("Variable1", "Variable2", "N")
    
    long_merged <- dplyr::left_join(cor_df, pval_df, by = c("Variable1", "Variable2")) %>%
      dplyr::left_join(n_df, by = c("Variable1", "Variable2")) %>%
      dplyr::mutate(P_value = round(P_value, 20), Correlation = round(Correlation, 10))
    
    write.csv(
      long_merged,
      file.path(output_dir, paste0(domain, "_", timepoint, "_correlations_", setname, "_pearson_long.csv")),
      row.names = FALSE
    )
    
    # Long-format outputs for Spearman
    cor_df_s <- as.data.frame(as.table(spearman$cor))
    names(cor_df_s) <- c("Variable1", "Variable2", "Correlation")
    pval_df_s <- as.data.frame(as.table(spearman$p))
    names(pval_df_s) <- c("Variable1", "Variable2", "P_value")
    n_df_s <- as.data.frame(as.table(spearman$n))
    names(n_df_s) <- c("Variable1", "Variable2", "N")
    
    long_merged_s <- dplyr::left_join(cor_df_s, pval_df_s, by = c("Variable1", "Variable2")) %>%
      dplyr::left_join(n_df_s, by = c("Variable1", "Variable2")) %>%
      dplyr::mutate(P_value = round(P_value, 20), Correlation = round(Correlation, 5))
    
    write.csv(
      long_merged_s,
      file.path(output_dir, paste0(domain, "_", timepoint, "_correlations_", setname, "_spearman_long.csv")),
      row.names = FALSE
    )
  }
  
  message("✅ Correlation and covariance matrices saved to ", output_dir)
}
