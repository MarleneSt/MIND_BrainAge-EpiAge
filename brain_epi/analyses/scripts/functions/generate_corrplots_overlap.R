# ============================================
# generate_corrplots_overlap
# adapation of Module 4: generate_model_plots focused on corrplots
# ============================================
#' Generate correlation plots for combined brain + epi dataset
#'
#' - Produces clustered and unclustered correlation plots
#' - Includes PAD, PAR, and predicted age values
#' - Supports categorical covariates
# ============================================
#' @param df Merged data frame with brain and epi predictions
#' @param age_var Vector of age vars: c("AGE_brain", "AGE_epi")
#' @param domain "brain_epi_overlap" (default), used in file naming
#' @param timepoint Timepoint string (used in file naming)
#' @param output_dir Output directory for saving plots, "results/performance/overlap/plots"
#' @param include_categorical List of categorical variables to include in additional correlation plots
#' @param sample_name Optional string included in plot titles
# ============================================
#' 10/12/2025 - MS removed BabyPy special case from function
#' ============================================

generate_corrplots_overlap <- function(
    df,
    age_var = c("AGE_brain", "AGE_epi"),
    domain = "brain_epi_overlap",
    timepoint,
    output_dir = "results/performance/overlap/plots",
    include_categorical = NULL,
    sample_name = NULL
) {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  message("\U0001F4CA Generating correlation plots for domain: ", domain, ", timepoint: ", timepoint)
  
  # Define model types to extract
  special_models <- c("DunedinPACE", "DunedinPACNI")
  pred_vars <- unique(c(grep("_predage$", names(df), value = TRUE), intersect(special_models, names(df))))
  pad_vars  <- unique(c(grep("_PAD$", names(df), value = TRUE), intersect(special_models, names(df))))
  par_vars  <- unique(c(grep("_PAR$", names(df), value = TRUE), intersect(special_models, names(df))))
  
  # Drop entirely missing columns
  pred_vars <- pred_vars[sapply(df[pred_vars], function(x) !all(is.na(x)))]
  pad_vars  <- pad_vars[sapply(df[pad_vars], function(x) !all(is.na(x)))]
  par_vars  <- par_vars[sapply(df[par_vars], function(x) !all(is.na(x)))]
  
  # Clock generation mapping (used to order models in plots)
  generation_map <- list(
    brain = list(
      "Centile2" = "1st Gen", "DevBrainAge" = "1st Gen",
      "DBN" = "1st Gen", "ENIGMA" = "1st Gen", "Kaufmann" = "1st Gen",
      "PyBrainAge" = "1st Gen", "Pyment" = "1st Gen", "DunedinPACNI" = "3rd Gen"
    ),
    epi = list(
      "Bohlin_conv" = "Gest", "EPIC_conv" = "Gest", "Knight_conv" = "Gest",
      "Horvath2013" = "1st Gen", "skinHorvath" = "1st Gen", "Hannum" = "1st Gen",
      "Wu" = "1st Gen", "PedBE" = "1st Gen", "PCBrainAge" = "1st Gen",
      "CorticalClock" = "1st Gen", "AltumAge" = "1st Gen", "cAge" = "1st Gen",
      "ZhangEN" = "1st Gen", "ZhangBLUP" = "1st Gen",
      "PhenoAge" = "2nd Gen", "PCGrimAge" = "2nd Gen", "DNAmTL" = "2nd Gen",
      "DunedinPACE" = "3rd Gen", "AdaptAge" = "4th Gen", "DamAge" = "4th Gen"
    )
  )
  
  title_prefix <- paste0(
    if (!is.null(sample_name)) paste0("Sample: ", sample_name, " | ") else "",
    "Timepoint: ", timepoint
  )
  
  # Define categorical variables (optional)
  cat_vars <- intersect(include_categorical, names(df))
  
  # === Internal correlation plot generation ===
  # Internal helper to create corrplots for each model type (predage, PAD, PAR)
  generate_corrplots <- function(vars, suffix, include_cats = FALSE) {
    if (length(vars) < 1) return()
    
    vars_with_age <- unique(c(age_var[age_var %in% names(df)], vars))
    clean_names <- gsub("_(predage|PAD|PAR)$", "", vars_with_age)
    
    # === NEW Combined Domain + Generation Ordering ===
    
    # Flatten the generation maps and tag with domain
    gen_df <- dplyr::bind_rows(
      data.frame(model = names(generation_map$brain), 
                 gen = unlist(generation_map$brain), 
                 domain = "brain"),
      data.frame(model = names(generation_map$epi), 
                 gen = unlist(generation_map$epi), 
                 domain = "epi")
    )
    
    # Create ordered generation levels
    generation_levels <- c("Gest", "1st Gen", "2nd Gen", "3rd Gen", "4th Gen")
    gen_df$gen <- factor(gen_df$gen, levels = generation_levels)
    
    # Sort first by domain, then generation
    gen_df <- gen_df[order(gen_df$domain, gen_df$gen), ]
    
    # Make clean names (match what’s used in your df column names)
    clean_names <- gsub("_(predage|PAD|PAR)$", "", vars_with_age)
    
    # Filter the ordering to include only the models actually present in data
    ordered_models <- gen_df$model[gen_df$model %in% clean_names]
    
    # Map back to full variable names in df
    ordered_vars <- vars_with_age[order(match(clean_names, ordered_models))]
    
    df_corr <- df[, ordered_vars, drop = FALSE]
    
    if (include_cats && length(cat_vars) > 0) {
      df_corr <- cbind(df_corr, df[, cat_vars, drop = FALSE])
      df_corr[cat_vars] <- lapply(df_corr[cat_vars], function(x) as.integer(as.factor(x)))
    }
    
    colnames(df_corr) <- gsub("_(predage|PAD|PAR)$", "", colnames(df_corr))
    
    corrmat <- cor(df_corr, use = "pairwise.complete.obs")
    
    plot_suffix <- if (include_cats) paste0(suffix, "_with_categorical") else suffix
    
    save_corrplot <- function(mat, filename, cluster = FALSE) {
      png(filename, width = 5000, height = 5000, res = 300)
      corrplot(mat, type = "lower", method = "color",
               order = if (cluster) "hclust" else "original",
               addCoef.col = "black", tl.col = "black", tl.cex = 1,
               main = paste0("Correlation Matrix (", 
                             if (cluster) "Clustered" else "Unclustered", 
                             if (include_cats) ", + Categorical" else "",
                             ")\n", title_prefix),
               cex.main = 1.7,
               mar = c(0, 0, 3, 0))
      dev.off()
    }
    
    save_corrplot(corrmat, file.path(output_dir, paste0(domain, "_", timepoint, "_corrplot_", plot_suffix, "_uncluster.png")), cluster = FALSE)
    save_corrplot(corrmat, file.path(output_dir, paste0(domain, "_", timepoint, "_corrplot_", plot_suffix, "_clustered.png")), cluster = TRUE)
  }
  
  # Generate correlation plots
  generate_corrplots(pred_vars, "pred", include_cats = FALSE)
  generate_corrplots(pred_vars, "pred", include_cats = TRUE)
  generate_corrplots(pad_vars,  "pad",  include_cats = FALSE)
  generate_corrplots(pad_vars,  "pad",  include_cats = TRUE)
  generate_corrplots(par_vars,  "par",  include_cats = FALSE)
  generate_corrplots(par_vars,  "par",  include_cats = TRUE)
  
  message("\u2705 Correlation plots saved to ", output_dir)
}
