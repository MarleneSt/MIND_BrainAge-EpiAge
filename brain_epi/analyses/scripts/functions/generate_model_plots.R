# ============================================
# generate_model_plots.R
# Module 4: generate_model_plots 
# ============================================
#' Generate visual performance summaries
#'
#' - Density plots for predicted age vs actual (not super useful)
#' - Boxplots of predicted age, PAD, PAR, MAE, wMAE
#' - Scatterplots per model with R and R²
#' - Optional scatterplot matrix and corrplots
# ============================================
#' @param df Data frame containing age variables and model predictions
#' @param age_var Column name of the chronological age variable
#' @param domain Either "brain" or "epi"
#' @param timepoint Timepoint string (used in file naming)
#' @param output_dir Output directory for saving plots
#' @param include_categorical List of categorical variables to include in correlation plots
#' @param sample_name Optional string included in plot titles
# ============================================
#' 10/12/2025 - MS removed BabyPy special case from function, some simplifications
#' ============================================

generate_model_plots <- function(
    df,
    age_var,
    domain,
    timepoint,
    output_dir = "results/performance",
    include_categorical = NULL,
    sample_name = NULL
) {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  message("\U0001F4CA Generating plots for domain: ", domain, ", timepoint: ", timepoint)
  
  # Filter to relevant columns
  special_models <- c("DunedinPACE", "DunedinPACNI") #they will be included in predage, PAD, and PAR matrices 
  
  pred_vars <- unique(c(grep("_predage$", names(df), value = TRUE), intersect(special_models, names(df))))
  pad_vars  <- unique(c(grep("_PAD$", names(df), value = TRUE), intersect(special_models, names(df))))
  par_vars  <- unique(c(grep("_PAR$", names(df), value = TRUE), intersect(special_models, names(df))))
  
  # Drop entirely missing variables to avoid plot issues
  pred_vars <- pred_vars[sapply(df[pred_vars], function(x) !all(is.na(x)))]
  pad_vars  <- pad_vars[sapply(df[pad_vars],  function(x) !all(is.na(x)))]
  par_vars  <- par_vars[sapply(df[par_vars],  function(x) !all(is.na(x)))]
  
  all_pred <- grep("_predage$", names(df), value = TRUE)
  excluded_pred <- setdiff(all_pred, pred_vars)
  if (length(excluded_pred) > 0) message("⚠️Excluded predicted ages due to all NA: ", paste(excluded_pred, collapse = ", "))
  
  # Exclude Dunedin models from boxplots and density plots
  # Exclude DunedinPACE and DunedinPACNI from boxplots and density plots
  exclude_from_plots <- c("DunedinPACE", "DunedinPACNI")
  pred_vars_plot <- setdiff(grep("_predage$", pred_vars, value = TRUE), exclude_from_plots)
  pad_vars_plot  <- setdiff(grep("_PAD$", pad_vars, value = TRUE), paste0(exclude_from_plots, "_PAD"))
  par_vars_plot  <- setdiff(grep("_PAR$", par_vars, value = TRUE), paste0(exclude_from_plots, "_PAR"))
  
  # Clock generation mapping for ordering in plots 
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
  
  #define title prefix for figure headings
  title_prefix <- paste0(
    if (!is.null(sample_name)) paste0("Sample: ", sample_name, " | ") else "",
    "Timepoint: ", timepoint
  )
  
  # === Density plot ===
  age_related <- age_var
  all_density_vars <- c(age_related, pred_vars_plot)
  df_long <- df %>%
    select(all_of(all_density_vars)) %>%
    pivot_longer(cols = everything(), names_to = "model", values_to = "value") %>%
    filter(!is.na(value))
  df_long$model_clean <- gsub("_(predage)$", "", df_long$model)
  
  p_density <- ggplot(df_long, aes(x = value, color = model_clean, fill = model_clean)) +
    geom_density(alpha = 0.3) +
    theme_classic() +
    scale_color_viridis(discrete = TRUE, option = "D") +
    scale_fill_viridis(discrete = TRUE, option = "D") +
    labs(x = "Age", y = "Density", title = paste("Age Density Plot -", title_prefix)) +
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5), text = element_text(size = 12)) + 
    xlim(min(df_long$value, na.rm = TRUE), max(df_long$value, na.rm = TRUE))
  
  ggsave(file.path(output_dir, paste0(domain, "_", timepoint, "_agedensity_all.png")),
         plot = p_density, width = 12, height = 10, dpi = 300)
  
  # === BOXPLOTS ===
  for (type in c("Predicted Age", "PAD", "PAR", "MAE", "WMAE")) {
    vars <- switch(type, 
                   "Predicted Age" = pred_vars_plot, 
                   "PAD" = pad_vars_plot, 
                   "PAR" = par_vars_plot, 
                   "MAE" = pred_vars_plot, 
                   "WMAE" = pred_vars_plot)
    if (length(vars) == 0) next
    
    if (type %in% c("MAE", "WMAE")) {
      df_long <- lapply(vars, function(var) {
        
        actual <- df[[age_var]]
        
        individual_mae <- abs(df[[var]] - actual) #df[[var]] = predicted age 
        denom <- max(actual, na.rm = TRUE) - min(actual, na.rm = TRUE)
        denom <- ifelse(denom == 0, NA_real_, denom)  # safeguard
        individual_wmae <- individual_mae / denom
        
        data.frame(model = gsub("_(predage|PAD|PAR)$", "", var),
          value = if (type == "MAE") individual_mae else individual_wmae
        )
      }) %>% dplyr::bind_rows()
      df_long$model_clean <- df_long$model
    } else {
      df_long <- df[, c(vars)] %>%
        tidyr::pivot_longer(cols = everything(), names_to = "model", values_to = "value") %>%
        dplyr::mutate(model_clean = gsub("_(predage|PAD|PAR)$", "", model), type = type)
    }
    
    # Order models by generation
    model_order <- names(generation_map[[domain]])
    model_order <- c(
      model_order[grep("Gest", generation_map[[domain]][model_order])],
      model_order[!generation_map[[domain]][model_order] %in% "Gest"]
    )
    df_long$model_clean <- factor(df_long$model_clean, levels = model_order)
    
    ref_line <- ifelse(type == "Predicted Age", mean(df[[age_var]], na.rm = TRUE), 0)
    y_label <- switch(type,
                      "Predicted Age" = "Predicted Age",
                      "PAD" = "PAD",
                      "PAR" = "PAR",
                      "MAE" = "Absolute Error",
                      "WMAE" = "Weighted Absolute Error")
    
    p <- ggplot(df_long, aes(x = model_clean, y = value, fill = model_clean)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.3, size = 0.6) +
      stat_summary(fun = "mean", geom = "point", color = "red", size = 3, shape = 18) +
      labs(title = paste("Boxplot", type, "|", title_prefix), x = "Model", y = y_label) +
      theme_classic(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none",
            plot.title = element_text(hjust = 0.5))
    
    if (type %in% c("Predicted Age")) {
      p <- p + geom_hline(yintercept = mean(df[[age_var]], na.rm = TRUE), linetype = "dashed", color = "red")
    }
    if (type %in% c("PAD", "PAR")) {
      p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "red")
    }
    
    ggsave(file.path(output_dir, paste0(domain, "_", timepoint, "_boxplot_", tolower(type), ".png")),
           plot = p, width = 12, height = 8, dpi = 300)
  }
  
  # === SCATTERPLOTS ===
  for (model in pred_vars) {
    model_clean <- gsub("_predage$", "", model)
  
    df_model <- df[, c(age_var, model)]
    names(df_model) <- c("age", "pred")
    
    df_model <- df_model[complete.cases(df_model), ]
    if (nrow(df_model) == 0) next
    
    r_val <- cor(df_model$age, df_model$pred)
    r2_val <- summary(lm(pred ~ age, data = df_model))$r.squared
    
    p <- ggplot(df_model, aes(x = age, y = pred)) +
      geom_point(alpha = 0.5) +
      geom_smooth(method = "lm", se = FALSE, color = "blue") +
      annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1.5,
               label = paste0("r = ", round(r_val, 3), ", R² = ", round(r2_val, 3)), size = 4) +
      labs(title = paste("Scatterplot:", model_clean, "|", title_prefix), x = "Chronological Age", y = "Predicted Age") +
      theme_classic(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5))
    
    ggsave(file.path(output_dir, paste0(domain, "_", timepoint, "_scatter_", model_clean, ".png")),
           plot = p, width = 7, height = 6, dpi = 300)
  }
  
  if (domain == "brain") {
    df_pairs <- df[, c(age_var, pred_vars)]
    df_pairs <- df_pairs[, colSums(!is.na(df_pairs)) > 1]
    p_matrix <- GGally::ggpairs(df_pairs, upper = list(continuous = wrap("cor", size = 3)),
                                lower = list(continuous = wrap("smooth", alpha = 0.3, size = 0.2)))
    ggsave(filename = file.path(output_dir, paste0(domain, "_", timepoint, "_scattermatrix.png")),
           plot = p_matrix, width = 12, height = 10, dpi = 300)
  }
  
  # === CORRPLOTS ===
  generate_corrplots <- function(vars, suffix, include_cats = FALSE) {
    if (length(vars) < 1) return()
    
    vars_with_age <- unique(c(age_var, vars))
    clean_names <- gsub("_(predage|PAD|PAR)$", "", vars_with_age)
    
    ordered_models <- names(generation_map[[domain]])
    ordered_models <- c(
      ordered_models[grep("Gest", generation_map[[domain]][ordered_models])],
      ordered_models[!generation_map[[domain]][ordered_models] %in% "Gest"]
    )
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
      png(filename, width = 4000, height = 4000, res = 300)
      corrplot(mat, type = "lower", method = "color",
               order = if (cluster) "hclust" else "original",
               addCoef.col = "black", tl.col = "black", tl.cex = 1,
               main = paste0("Correlation Matrix (", 
                             if (cluster) "Clustered" else "Unclustered", 
                             if (include_cats) ", + Categorical" else "",
                             ")\n", title_prefix),
               mar = c(0, 0, 3, 0))
      dev.off()
    }
    
    save_corrplot(corrmat, file.path(output_dir, paste0(domain, "_", timepoint, "_corrplot_", plot_suffix, "_uncluster.png")), cluster = FALSE)
    save_corrplot(corrmat, file.path(output_dir, paste0(domain, "_", timepoint, "_corrplot_", plot_suffix, "_clustered.png")), cluster = TRUE)
  }
  
  # Define categorical variables
  cat_vars <- intersect(include_categorical, names(df))
  
  # Generate plots
  generate_corrplots(pred_vars, "pred", include_cats = FALSE)
  generate_corrplots(pred_vars, "pred", include_cats = TRUE)
  generate_corrplots(pad_vars,  "pad",  include_cats = FALSE)
  generate_corrplots(pad_vars,  "pad",  include_cats = TRUE)
  generate_corrplots(par_vars,  "par",  include_cats = FALSE)
  generate_corrplots(par_vars,  "par",  include_cats = TRUE)
  
  message("\u2705 Module 4 plots saved to ", output_dir)
}
  