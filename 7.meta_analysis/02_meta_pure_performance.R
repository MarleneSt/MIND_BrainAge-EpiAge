#### META-ANALYSES SCRIPT - PURE performance analyses
#### MIND BrainAge~EpiAge
#### last updated: 06/07/2026

#### library ----
library(tidyverse)
library(stringr)
library(metafor)
library(ggplot2)
library(pheatmap)
library(forcats)
library(tidyr)
library(purrr)
library(tibble)


#### define important folders ----
home <- "data" # [REPRODUCIBILITY] relative path
data_dir <- file.path(home)
results_dir <- "results" # [REPRODUCIBILITY] relative path

#### [REPRODUCIBILITY] helper and data loading added so this script runs on its own
#### (previously used objects still in memory from 1.main_meta_wMAE.R). Copied from
#### 01_main_meta_wMAE.R ("some helpers" and "load data" sections).

# Helper: I2 for multilevel metafor models
# code from Viechtbauer's metafor website: https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate
calc_I2_mv <- function(model, component_names = NULL) {

  W <- diag(1 / model$vi)
  X <- model.matrix(model)

  P <- W - W %*% X %*%
    solve(t(X) %*% W %*% X) %*%
    t(X) %*% W

  sampling_var <- (model$k - model$p) / sum(diag(P))

  I2_total <- 100 * sum(model$sigma2) /
    (sum(model$sigma2) + sampling_var)

  I2_components <- 100 * model$sigma2 /
    (sum(model$sigma2) + sampling_var)

  if (!is.null(component_names)) {
    names(I2_components) <- component_names
  }

  list(
    I2_total = I2_total,
    I2_components = I2_components
  )
}

perf_brain <- readRDS(file.path(results_dir, "perf_brain.rds"))

perf_epi <- readRDS(file.path(results_dir, "perf_epi.rds"))

perf_brain <- perf_brain %>%
  mutate(
    cohort = factor(cohort),
    timepoint = factor(timepoint),
    model = factor(model)
  )

perf_epi <- perf_epi %>%
  mutate(
    cohort = factor(cohort),
    timepoint = factor(timepoint),
    model = factor(model)
  )

##### pure performance datasets -----

##### brain pure sample -----

perf_brain_pure <- perf_brain %>%
  filter(!model %in% c("DunedinPACNI")) %>%
  mutate(
    overlap_lower_brain = pmax(min_age, training_age_min_brain),
    overlap_upper_brain = pmin(max_age, training_age_max_brain),
    
    overlap_brain = if_else(
      overlap_upper_brain >= overlap_lower_brain,
      "yes",
      "no"
    ),
    
    overlap_width_brain = if_else(
      overlap_brain == "yes",
      overlap_upper_brain - overlap_lower_brain,
      0
    ),
    
    study_width_brain = max_age - min_age,
    
    overlap_perc_brain = 100 * overlap_width_brain / study_width_brain
  ) %>%
  filter(overlap_perc_brain == 100) %>%
  select(
    -overlap_lower_brain,
    -overlap_upper_brain,
    -overlap_width_brain,
    -study_width_brain
  ) %>%
  droplevels()

##### epi pure sample -----

gestational_models <- c("Knight", "Bohlin", "EPIC")


perf_epi_pure <- perf_epi %>%
  mutate(
    model = recode(
      model,
      "Bohlin_conv" = "Bohlin",
      "Knight_conv" = "Knight",
      "EPIC_conv" = "EPIC"
    )
  ) %>%
  filter(!model %in% c("DunedinPACE", "DNAmTL")) %>%
  mutate(
    is_gestational_model = model %in% gestational_models,
    
    overlap_lower_epi = pmax(min_age, training_age_min_epi),
    overlap_upper_epi = pmin(max_age, training_age_max_epi),
    
    overlap_epi = case_when(
      is_gestational_model ~ "gestational_no_training_range",
      overlap_upper_epi >= overlap_lower_epi ~ "yes",
      TRUE ~ "no"
    ),
    
    overlap_width_epi = case_when(
      overlap_epi == "yes" ~ overlap_upper_epi - overlap_lower_epi,
      is_gestational_model ~ NA_real_,
      TRUE ~ 0
    ),
    
    study_width_epi = max_age - min_age,
    
    overlap_perc_epi = case_when(
      overlap_epi == "yes" ~ 100 * overlap_width_epi / study_width_epi,
      is_gestational_model ~ NA_real_,
      TRUE ~ 0
    ),
    
    tissue_overlap = case_when(
      tissue_training == "multi" ~ "yes",
      tissue_training == tissue_test ~ "yes",
      TRUE ~ "no"
    ),
    
    tissue_overlap = factor(tissue_overlap, levels = c("no", "yes"))
  ) %>%
  filter(
    tissue_overlap == "yes",
    overlap_perc_epi == 100 | is_gestational_model
  ) %>%
  mutate(
    log_wMAE = log(wMAE_test),
    log_wMAE_var = (wMAE_SE_boot / wMAE_test)^2
  ) %>%
  select(
    -overlap_lower_epi,
    -overlap_upper_epi,
    -overlap_width_epi,
    -study_width_epi
  ) %>%
  droplevels()

perf_brain_sub <- perf_brain_pure
perf_epi_sub <- perf_epi_pure

##### wMAE no DunedinPACNI & DNAmTL PURE -----

###### epi wMAE global pure ------


# meta-analysis did not converge due to differences in sampling variances which span ~7 orders of magnitude from ~0.00002 to ~34,800 when squared (also tried out scaling based on median but didn't work)

# take the log 
perf_epi_sub <- perf_epi_sub |>
  mutate(
    log_wMAE = log(wMAE_test),
    # Delta method: Var(log(X)) ≈ (SE_X / X)^2
    log_wMAE_var = (wMAE_SE_boot / wMAE_test)^2
  )

# --- set up folder ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi <- file.path(results_dir, paste0("perf_epi_wMAE_PURE_noDunedinDNAmTL_global_", run_id))

dir.create(results_dir_epi, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_epi, "analysis_log.txt")

# --- descriptives ---

write.csv(perf_epi_sub, file.path(results_dir_epi, "perf_epi_sub_descr_detailed.csv"), row.names=FALSE)

perf_epi_sub_desc <- perf_epi_sub %>%
  group_by(model, domain) %>%
  summarise(
    
    k = n(),  # number of cohort/timepoint observations per model
    
    # MAE
    MAE_mean   = mean(MAE, na.rm = TRUE),
    MAE_sd     = sd(MAE, na.rm = TRUE),
    MAE_min    = min(MAE, na.rm = TRUE),
    MAE_max    = max(MAE, na.rm = TRUE),
    
    wMAE_mean   = mean(wMAE_test, na.rm = TRUE),
    wMAE_sd     = sd(wMAE_test, na.rm = TRUE),
    wMAE_min    = min(wMAE_test, na.rm = TRUE),
    wMAE_max    = max(wMAE_test, na.rm = TRUE),
    
    # RMSE
    RMSE_mean  = mean(RMSE, na.rm = TRUE),
    RMSE_sd    = sd(RMSE, na.rm = TRUE),
    RMSE_min   = min(RMSE, na.rm = TRUE),
    RMSE_max   = max(RMSE, na.rm = TRUE),
    
    # R2 and correlations
    R2_mean      = mean(R2, na.rm = TRUE),
    Pearson_mean = mean(Pearson, na.rm = TRUE),
    Spearman_mean = mean(Spearman, na.rm = TRUE),
    
    .groups = "drop"
  )

write.csv(perf_epi_sub_desc, file.path(results_dir_epi, "model_perf_epi_desc.csv"), row.names = FALSE)

# model retention 
epi_retention <- tibble(
  metric = c("Rows retained", "Models retained", "Cohort-timepoints retained"),
  value_percent = c(
    nrow(perf_epi_pure) / nrow(perf_epi %>% filter(!model %in% c("DunedinPACE", "DNAmTL"))),
    n_distinct(perf_epi_pure$model) / n_distinct((perf_epi %>% filter(!model %in% c("DunedinPACE", "DNAmTL")))$model),
    n_distinct(paste(perf_epi_pure$cohort, perf_epi_pure$timepoint)) /
      n_distinct(paste(perf_epi$cohort, perf_epi$timepoint))
  ) * 100
)

# --- meta ---
epi_wMAE_mv_global <- rma.mv(
  yi = log_wMAE,
  V  = log_wMAE_var,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data   = perf_epi_sub,
  method = "REML"
)

summary(epi_wMAE_mv_global)

# Back-transform global pooled estimate
overall_est <- data.frame(
  model      = "Pooled",
  wMAE       = exp(epi_wMAE_mv_global$beta[1, 1]),
  ci.lb      = exp(epi_wMAE_mv_global$ci.lb[1]),
  ci.ub      = exp(epi_wMAE_mv_global$ci.ub[1]),
  pval = epi_wMAE_mv_global$pval[1],  
  generation = "Pooled"
)

overall_est

saveRDS(
  epi_wMAE_mv_global,
  file.path(results_dir_epi, "model_epi_wMAE_mv_global.rds")
)

capture.output(
  summary(epi_wMAE_mv_global),
  file = file.path(results_dir_epi, "model_summary.txt")
)

write.csv(
  overall_est,
  file.path(results_dir_epi, "overall_wMAE_estimate.csv"),
  row.names = FALSE
)

# --- heterogeneity estimates ---
I2 <- calc_I2_mv(
  epi_wMAE_mv_global,
  component_names = c("cohort", "timepoint_within_cohort", "model")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_epi, "heterogeneity_I2.csv"), row.names = FALSE)

# --- LOO cohort (for log-transformed model) ---

cohorts <- unique(perf_epi_sub$cohort)

loo_cohort <- lapply(cohorts, function(cc) {
  
  tmp <- subset(perf_epi_sub, cohort != cc)
  
  fit <- rma.mv(
    yi = log_wMAE,
    V  = log_wMAE_var,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model
    ),
    data = tmp,
    method = "REML"
  )
  
  log_est <- coef(fit)[1]
  log_se  <- sqrt(vcov(fit)[1, 1])
  
  data.frame(
    cohort_removed = cc,
    
    # log scale
    log_wMAE = log_est,
    log_se   = log_se,
    log_ci.lb = log_est - 1.96 * log_se,
    log_ci.ub = log_est + 1.96 * log_se,
    
    # back-transformed scale
    wMAE = exp(log_est),
    ci.lb = exp(log_est - 1.96 * log_se),
    ci.ub = exp(log_est + 1.96 * log_se)
  )
})

loo_cohort <- do.call(rbind, loo_cohort) %>%
  as_tibble()

loo_summary <- tibble(
  min_wMAE   = min(loo_cohort$wMAE, na.rm = TRUE),
  max_wMAE   = max(loo_cohort$wMAE, na.rm = TRUE),
  range_wMAE = max_wMAE - min_wMAE,
  sd_wMAE    = sd(loo_cohort$wMAE, na.rm = TRUE),
  
  min_log_wMAE   = min(loo_cohort$log_wMAE, na.rm = TRUE),
  max_log_wMAE   = max(loo_cohort$log_wMAE, na.rm = TRUE),
  range_log_wMAE = max_log_wMAE - min_log_wMAE,
  sd_log_wMAE    = sd(loo_cohort$log_wMAE, na.rm = TRUE)
)

write.csv(
  loo_cohort,
  file.path(results_dir_epi, "leave_one_cohort_out.csv"),
  row.names = FALSE
)

write.csv(
  loo_summary,
  file.path(results_dir_epi, "leave_one_cohort_out_summary.csv"),
  row.names = FALSE
)

full_wMAE <- exp(coef(epi_wMAE_mv_global)[1])

png(
  filename = file.path(results_dir_epi, "plot_global_wMAE_allbD_loo.png"),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_cohort,
  aes(
    x = reorder(cohort_removed, wMAE),
    y = wMAE
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = ci.lb,
      ymax = ci.ub
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept = full_wMAE,
    linetype = "dashed"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "Pooled wMAE",
    title = "Leave-one-cohort-out sensitivity of pooled wMAE"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# full analysis log

sink(log_file)

cat("Epi wMAE multilevel meta-analysis (log-transformed)\n")
cat("===================================\n\n")
cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_epi, "\n\n")

cat("Model retention\n")
cat("=============\n")
print(epi_retention)

cat("\n\nModel summary\n")
cat("=============\n")
print(summary(epi_wMAE_mv_global))

cat("\n\nModel summary - back transformed \n")
cat("=============\n")
print(overall_est)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_cohort, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_summary, n = Inf, width = Inf)

cat("\n\nDescriptives\n")
cat("============\n")
print(perf_epi_sub_desc, n = Inf, width = Inf)

sink()

###### epi wMAE modelwise PURE ------

# --- set up folder ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi_modelwise <- file.path(results_dir, paste0("perf_epi_wMAE_PURE_noDunedinDNAmTL_modelwise_", run_id))

dir.create(results_dir_epi_modelwise, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_epi_modelwise, "analysis_modelwise_log.txt")

# --- meta ---
epi_wMAE_mv_modelwise <- rma.mv(
  yi = log_wMAE,
  V  = log_wMAE_var,
  mods = ~ model - 1,
  random = list(
    ~ 1 | cohort/timepoint
  ),
  data   = perf_epi_sub,
  method = "REML"
)

summary(epi_wMAE_mv_modelwise)



# Back-transform global pooled estimate
epi_wMAE_est <- data.frame(
  model = gsub("^model", "", rownames(epi_wMAE_mv_modelwise$beta)),
  log_wMAE = epi_wMAE_mv_modelwise$beta[, 1],
  wMAE = exp(epi_wMAE_mv_modelwise$beta[, 1]),
  ci.lb = exp(epi_wMAE_mv_modelwise$ci.lb),
  ci.ub = exp(epi_wMAE_mv_modelwise$ci.ub),
  pval = epi_wMAE_mv_modelwise$pval
)

X <- model.matrix(epi_wMAE_mv_modelwise)

model_cols <- colnames(X)[
  str_detect(colnames(X), "^model") &
    !str_detect(colnames(X), ":")
]

effect_sizes_per_model <- tibble(
  model = str_remove(model_cols, "^model"),
  k = colSums(X[, model_cols, drop = FALSE] != 0)
)


epi_wMAE_est <- epi_wMAE_est %>%
  left_join(effect_sizes_per_model, by = "model")

epi_wMAE_est

saveRDS(
  epi_wMAE_mv_modelwise,
  file.path(results_dir_epi_modelwise, "model_epi_wMAE_mv_modelwise.rds")
)

capture.output(
  summary(epi_wMAE_mv_modelwise),
  file = file.path(results_dir_epi_modelwise, "model_summary.txt")
)

write.csv(
  epi_wMAE_est,
  file.path(results_dir_epi_modelwise, "modelwise_wMAE_estimates.csv"),
  row.names = FALSE
)

# --- heterogeneity ---
I2 <- calc_I2_mv(
  epi_wMAE_mv_modelwise,
  component_names = c("cohort", "timepoint_within_cohort")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_epi_modelwise, "heterogeneity_I2.csv"), row.names = FALSE)

# --- LOO cohort (for log-transformed model) ---

cohorts <- unique(perf_epi_sub$cohort)

loo_cohort_modelwise <- lapply(cohorts, function(cc) {
  
  tmp <- subset(perf_epi_sub, cohort != cc)
  
  fit <- rma.mv(
    yi = log_wMAE,
    V  = log_wMAE_var,
    mods = ~ model - 1,
    random = list(
      ~ 1 | cohort/timepoint
    ),
    data = tmp,
    method = "REML"
  )
  
  tibble(
    cohort_removed = cc,
    term = rownames(fit$b),
    
    # log-scale estimates
    log_wMAE = as.numeric(fit$b),
    log_se = fit$se,
    log_ci.lb = fit$ci.lb,
    log_ci.ub = fit$ci.ub,
    pval = fit$pval
  ) %>%
    mutate(
      model = str_remove(term, "^model"),
      
      # back-transformed estimates: geometric mean wMAE
      wMAE = exp(log_wMAE),
      ci.lb = exp(log_ci.lb),
      ci.ub = exp(log_ci.ub)
    ) %>%
    select(
      cohort_removed,
      model,
      log_wMAE,
      log_se,
      log_ci.lb,
      log_ci.ub,
      wMAE,
      ci.lb,
      ci.ub,
      pval
    )
}) %>%
  bind_rows()

write.csv(
  loo_cohort_modelwise,
  file.path(results_dir_epi_modelwise, "modelwise_leave_one_cohort_out.csv"),
  row.names = FALSE
)

modelwise_loo_summary <- loo_cohort_modelwise %>%
  group_by(model) %>%
  summarise(
    min_wMAE = min(wMAE, na.rm = TRUE),
    max_wMAE = max(wMAE, na.rm = TRUE),
    range_wMAE = max_wMAE - min_wMAE,
    sd_wMAE = sd(wMAE, na.rm = TRUE),
    
    min_log_wMAE = min(log_wMAE, na.rm = TRUE),
    max_log_wMAE = max(log_wMAE, na.rm = TRUE),
    range_log_wMAE = max_log_wMAE - min_log_wMAE,
    sd_log_wMAE = sd(log_wMAE, na.rm = TRUE),
    
    .groups = "drop"
  )

png(
  filename = file.path(results_dir_epi_modelwise, "plot_modelwise_wMAE_loo.png"),
  width = 2200,
  height = 1400,
  res = 300
)

ggplot(
  loo_cohort_modelwise,
  aes(
    x = reorder(cohort_removed, wMAE),
    y = wMAE
  )
) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = ci.lb, ymax = ci.ub),
    width = 0.15
  ) +
  facet_wrap(~ model, scales = "free_y") +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "Pooled wMAE",
    title = "Modelwise leave-one-cohort-out sensitivity of pooled wMAE"
  ) +
  theme_minimal(base_size = 14)

dev.off()

write.csv(
  modelwise_loo_summary,
  file.path(results_dir_epi_modelwise, "modelwise_leave_one_cohort_out_summary.csv"),
  row.names = FALSE
)

# full analysis log

sink(log_file)

cat("Epi wMAE multilevel meta-analysis modelwise (log-transformed)\n")
cat("===================================\n\n")
cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_epi_modelwise, "\n\n")

cat("Model summary\n")
cat("=============\n")
print(summary(epi_wMAE_mv_modelwise))

cat("Model summary - back transformed \n")
cat("=============\n")
print(epi_wMAE_est)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_cohort_modelwise, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(modelwise_loo_summary, n = Inf, width = Inf)

sink()

###### epi wMAE age global PURE ------
# all but DunedinPACNI & DNAmTL

# --- set up ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi_global_ageint <- file.path(results_dir, paste0("perf_epi_wMAE_PURE_noDunedinDNAmTL_global_ageint_", run_id))

dir.create(results_dir_epi_global_ageint, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_epi_global_ageint, "analysis_log_global_ageint.txt")

# --- data prep ---

# remove gestational clocks as it doesn't make sense to look at age here 
perf_epi_sub_nogest   <- perf_epi_sub   %>% filter(!model %in% c("Bohlin", "Knight", "EPIC")) %>% droplevels()

# center age 
perf_epi_sub_nogest <- perf_epi_sub_nogest %>%
  mutate(mean_age_c = mean_age - mean(mean_age, na.rm = TRUE))

# --- meta ---

epi_wMAE_mv_global_ageint <- rma.mv(
  yi = log_wMAE,
  V  = log_wMAE_var,
  mods = ~ mean_age_c,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data   = perf_epi_sub_nogest,
  method = "REML"
)

summary(epi_wMAE_mv_global_ageint)

# --- heterogeneity estimates ---
I2 <- calc_I2_mv(
  epi_wMAE_mv_global_ageint,
  component_names = c("cohort", "timepoint_within_cohort", "model")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_epi_global_ageint, "heterogeneity_I2.csv"), row.names = FALSE)

# back transform 
# --- back-transform main model results ---

age_model_est <- tibble(
  term = names(coef(epi_wMAE_mv_global_ageint)),
  estimate_log = as.numeric(coef(epi_wMAE_mv_global_ageint)),
  se_log = epi_wMAE_mv_global_ageint$se,
  ci.lb_log = epi_wMAE_mv_global_ageint$ci.lb,
  ci.ub_log = epi_wMAE_mv_global_ageint$ci.ub,
  pval = epi_wMAE_mv_global_ageint$pval
) %>%
  mutate(
    estimate_exp = exp(estimate_log),
    ci.lb_exp = exp(ci.lb_log),
    ci.ub_exp = exp(ci.ub_log),
    percent_change = (exp(estimate_log) - 1) * 100,
    percent_ci.lb = (exp(ci.lb_log) - 1) * 100,
    percent_ci.ub = (exp(ci.ub_log) - 1) * 100
  )

age_model_est

write.csv(
  age_model_est,
  file.path(results_dir_epi_global_ageint, "age_moderator_backtransformed_results.csv"),
  row.names = FALSE
)

saveRDS(
  epi_wMAE_mv_global_ageint,
  file.path(results_dir_epi_global_ageint, "model_epi_wMAE_mv_global_ageint_logscale.rds")
)

capture.output(
  summary(epi_wMAE_mv_global_ageint),
  file = file.path(results_dir_epi_global_ageint, "model_summary_log_scale.txt")
)


# --- raw data plot  ---

png(filename = file.path(results_dir_epi_global_ageint, "plot_meanage_wMAE_allbD_rawdata.png"), width = 1800,height = 1200,res = 300)
ggplot(perf_epi_sub_nogest,
       aes(mean_age, wMAE_test)) +
  geom_point(alpha = .7) +
  geom_smooth(method = "lm") +
  labs(
    x = "Mean sample age",
    y = "Weighted MAE (MAE ÷ age range)",
    title = "Raw association between mean age and wMAE"
  ) + 
  theme_minimal()
dev.off()

# if we wanted to plot the meta results (raw datapoints, superimposed cohort-timepoint means & regression line based on meta results )

# --- prediction grid ---

mean_age_mean <- mean(perf_epi_sub_nogest$mean_age, na.rm = TRUE)

newdat <- data.frame(
  mean_age = seq(
    min(perf_epi_sub_nogest$mean_age, na.rm = TRUE),
    max(perf_epi_sub_nogest$mean_age, na.rm = TRUE),
    length.out = 100
  )
) %>%
  mutate(
    mean_age_c = mean_age - mean_age_mean
  )

# --- model predictions on log scale ---

pred <- predict(
  epi_wMAE_mv_global_ageint,
  newmods = newdat$mean_age_c
)

# --- back-transform predictions to wMAE scale ---

newdat <- newdat %>%
  mutate(
    # log scale
    log_pred  = pred$pred,
    log_ci.lb = pred$ci.lb,
    log_ci.ub = pred$ci.ub,
    
    # back-transformed scale
    pred  = exp(log_pred),
    ci.lb = exp(log_ci.lb),
    ci.ub = exp(log_ci.ub)
  )

write.csv(
  newdat,
  file.path(results_dir_epi_global_ageint, "age_prediction_grid.csv"),
  row.names = FALSE
)

# --- optional: cohort-timepoint means for cleaner overlay ---

cohort_time_means <- perf_epi_sub_nogest %>%
  group_by(cohort, timepoint) %>%
  summarise(
    mean_age = mean(mean_age, na.rm = TRUE),
    mean_wMAE = mean(wMAE_test, na.rm = TRUE),
    .groups = "drop"
  )

# --- plot ---

png(
  filename = file.path(
    results_dir_epi_global_ageint,
    "plot_meanage_wMAE_meta_regression_backtransformed.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot() +
  
  # raw observations
  geom_point(
    data = perf_epi_sub_nogest,
    aes(
      x = mean_age,
      y = wMAE_test,
      colour = cohort
    ),
    alpha = .25,
    size = 2
  ) +
  
  # cohort-timepoint means
  geom_point(
    data = cohort_time_means,
    aes(
      x = mean_age,
      y = mean_wMAE,
      colour = cohort
    ),
    size = 4
  ) +
  
  # meta-regression line
  geom_line(
    data = newdat,
    aes(
      x = mean_age,
      y = pred
    ),
    linewidth = 1.2,
    colour = "black"
  ) +
  
  # confidence interval ribbon
  geom_ribbon(
    data = newdat,
    aes(
      x = mean_age,
      ymin = ci.lb,
      ymax = ci.ub
    ),
    alpha = .15
  ) +
  
  labs(
    title = "Meta-regression of epigenetic clock performance on mean age",
    x = "Mean sample age (years)",
    y = "Geometric mean weighted MAE",
    colour = "Cohort"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "right"
  )

dev.off()

# --- LOO ---

cohorts <- unique(perf_epi_sub_nogest$cohort)

loo_age <- lapply(cohorts, function(cc) {
  
  tmp <- perf_epi_sub_nogest %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = log_wMAE,
    V  = log_wMAE_var,
    mods = ~ mean_age_c,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model
    ),
    data = tmp,
    method = "REML"
  )
  
  log_intercept <- coef(fit)[1]
  log_slope     <- coef(fit)[2]
  
  log_intercept_se <- sqrt(vcov(fit)[1, 1])
  log_slope_se     <- sqrt(vcov(fit)[2, 2])
  
  data.frame(
    cohort_removed = cc,
    
    log_intercept = log_intercept,
    log_slope = log_slope,
    log_intercept_se = log_intercept_se,
    log_slope_se = log_slope_se,
    
    intercept_wMAE = exp(log_intercept),
    intercept_ci.lb = exp(log_intercept - 1.96 * log_intercept_se),
    intercept_ci.ub = exp(log_intercept + 1.96 * log_intercept_se),
    
    slope_ratio_per_year = exp(log_slope),
    slope_percent_per_year = (exp(log_slope) - 1) * 100,
    
    slope_percent_ci.lb = (exp(log_slope - 1.96 * log_slope_se) - 1) * 100,
    slope_percent_ci.ub = (exp(log_slope + 1.96 * log_slope_se) - 1) * 100
  )
})

loo_age_df <- do.call(rbind, loo_age) %>%
  as_tibble()

loo_age_summary <- tibble(
  min_intercept_wMAE = min(loo_age_df$intercept_wMAE, na.rm = TRUE),
  max_intercept_wMAE = max(loo_age_df$intercept_wMAE, na.rm = TRUE),
  range_intercept_wMAE = max_intercept_wMAE - min_intercept_wMAE,
  sd_intercept_wMAE = sd(loo_age_df$intercept_wMAE, na.rm = TRUE),
  
  min_percent_per_year = min(loo_age_df$slope_percent_per_year, na.rm = TRUE),
  max_percent_per_year = max(loo_age_df$slope_percent_per_year, na.rm = TRUE),
  range_percent_per_year = max_percent_per_year - min_percent_per_year,
  sd_percent_per_year = sd(loo_age_df$slope_percent_per_year, na.rm = TRUE)
)

write.csv(
  loo_age_df,
  file.path(results_dir_epi_global_ageint, "leave_one_cohort_out_age_model.csv"),
  row.names = FALSE
)

write.csv(
  loo_age_summary,
  file.path(results_dir_epi_global_ageint, "leave_one_cohort_out_age_model_summary.csv"),
  row.names = FALSE
)

full_percent_slope <- (exp(coef(epi_wMAE_mv_global_ageint)[2]) - 1) * 100

png(
  filename = file.path(results_dir_epi_global_ageint, "plot_meanage_wMAE_loo_percent_change.png"),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_age_df,
  aes(
    x = reorder(cohort_removed, slope_percent_per_year),
    y = slope_percent_per_year
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = slope_percent_ci.lb,
      ymax = slope_percent_ci.ub
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept = full_percent_slope,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 0,
    colour = "grey60"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "% change in geometric mean wMAE per year",
    title = "Leave-one-cohort-out sensitivity of age effect"
  ) +
  theme_minimal(base_size = 14)

dev.off()
# --- log ---

sink(log_file)

cat("Epi wMAE multilevel meta-analysis - global age moderator, log-transformed\n")
cat("==========================================================================\n\n")
cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_epi_global_ageint, "\n\n")

cat("Model summary, log scale\n")
cat("========================\n")
print(summary(epi_wMAE_mv_global_ageint))

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nBack-transformed model results\n")
cat("==============================\n")
print(age_model_est, n = Inf, width = Inf)

cat("\n\nEffect sizes per model\n")
cat("======================\n")
print(effect_sizes_per_model, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_age_df, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_age_summary, n = Inf, width = Inf)

sink()


##### brain pure sample -----

###### brain wMAE global PURE  ------

# --- set up folder ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_brain <- file.path(results_dir, paste0("perf_brain_wMAE_PURE_noDunedin_global_", run_id))
dir.create(results_dir_brain, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_brain, "analysis_log.txt")

# --- descriptives ---

write.csv(perf_brain_sub, file.path(results_dir_brain, "perf_brain_sub_descr_detailed.csv"), row.names=FALSE)

perf_brain_sub_desc <- perf_brain_sub %>%
  group_by(model, domain) %>%
  summarise(
    
    k = n(),  # number of cohort/timepoint observations per model
    
    # MAE
    MAE_mean   = mean(MAE, na.rm = TRUE),
    MAE_sd     = sd(MAE, na.rm = TRUE),
    MAE_min    = min(MAE, na.rm = TRUE),
    MAE_max    = max(MAE, na.rm = TRUE),
    
    wMAE_mean   = mean(wMAE_test, na.rm = TRUE),
    wMAE_sd     = sd(wMAE_test, na.rm = TRUE),
    wMAE_min    = min(wMAE_test, na.rm = TRUE),
    wMAE_max    = max(wMAE_test, na.rm = TRUE),
    
    # RMSE
    RMSE_mean  = mean(RMSE, na.rm = TRUE),
    RMSE_sd    = sd(RMSE, na.rm = TRUE),
    RMSE_min   = min(RMSE, na.rm = TRUE),
    RMSE_max   = max(RMSE, na.rm = TRUE),
    
    # R2 and correlations
    R2_mean      = mean(R2, na.rm = TRUE),
    Pearson_mean = mean(Pearson, na.rm = TRUE),
    Spearman_mean = mean(Spearman, na.rm = TRUE),
    
    .groups = "drop"
  )

write.csv(perf_brain_sub_desc, file.path(results_dir_brain, "model_perf_brain_desc.csv"), row.names = FALSE)


# model retention
brain_retention <- tibble(
  metric = c("Rows retained", "Models retained", "Cohort-timepoints retained"),
  value_percent = c(
    nrow(perf_brain_pure) / nrow(perf_brain %>% filter(!model %in% c("DunedinPACNI"))),
    n_distinct(perf_brain_pure$model) / n_distinct((perf_brain %>% filter(!model %in% c("DunedinPACNI")))$model),
    n_distinct(paste(perf_brain_pure$cohort, perf_brain_pure$timepoint)) /
      n_distinct(paste(perf_brain$cohort, perf_brain$timepoint))
  ) * 100
)


# --- meta ---
brain_wMAE_mv_global <- rma.mv(
  yi = wMAE_test,
  V  = wMAE_SE_boot^2,
  random = list(
    ~ 1 | cohort/timepoint,  # timepoints nested within cohorts 
    ~ 1 | model              # different true means across models
  ),
  data   = perf_brain_sub,
  method = "REML",
  verbose = TRUE
)
summary(brain_wMAE_mv_global)

saveRDS(
  brain_wMAE_mv_global,
  file.path(results_dir_brain, "model_brain_wMAE_mv_global.rds")
)

capture.output(
  summary(brain_wMAE_mv_global),
  file = file.path(results_dir_brain, "model_summary.txt")
)

brain_wMAE_global_est <- data.frame(
  model = "Pooled",
  wMAE = brain_wMAE_mv_global$beta[1, 1],
  se = brain_wMAE_mv_global$se[1],
  ci.lb = brain_wMAE_mv_global$ci.lb[1],
  ci.ub = brain_wMAE_mv_global$ci.ub[1],
  pval = brain_wMAE_mv_global$pval[1]
)

write.csv(
  brain_wMAE_global_est,
  file.path(results_dir_brain, "overall_wMAE_estimate.csv"),
  row.names = FALSE
)

# --- heterogeneity estimated ---
I2 <- calc_I2_mv(
  brain_wMAE_mv_global,
  component_names = c("cohort", "timepoint_within_cohort", "model")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_brain, "heterogeneity_I2.csv"), row.names = FALSE)

# --- LOO cohort ---

cohorts <- unique(perf_brain_sub$cohort)

loo_cohort <- lapply(cohorts, function(cc) {
  tmp <- subset(perf_brain_sub, cohort != cc)
  
  fit <- rma.mv(
    yi = wMAE_test,
    V  = wMAE_SE_boot^2,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model
    ),
    data = tmp,
    method = "REML"
  )
  
  est <- coef(fit)[1]
  se <- sqrt(vcov(fit)[1, 1])
  
  data.frame(
    cohort_removed = cc,
    wMAE = est,
    se = se,
    ci.lb = est - 1.96 * se,
    ci.ub = est + 1.96 * se,
    pval = fit$pval[1]
  )
})

loo_cohort <- do.call(rbind, loo_cohort) %>%
  as_tibble()

loo_cohort
range(loo_cohort$wMAE)
sd(loo_cohort$wMAE)

loo_summary <- tibble(
  min_wMAE   = min(loo_cohort$wMAE, na.rm = TRUE),
  max_wMAE   = max(loo_cohort$wMAE, na.rm = TRUE),
  range_wMAE = max_wMAE - min_wMAE,
  sd_wMAE    = sd(loo_cohort$wMAE, na.rm = TRUE)
)

write.csv(
  loo_cohort,
  file.path(results_dir_brain, "leave_one_cohort_out.csv"),
  row.names = FALSE
)

write.csv(
  loo_summary,
  file.path(results_dir_brain, "leave_one_cohort_out_summary.csv"),
  row.names = FALSE
)

full_wMAE <- coef(brain_wMAE_mv_global)[1]

png(
  filename = file.path(results_dir_brain, "plot_global_wMAE_allbD_loo.png"),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_cohort,
  aes(
    x = reorder(cohort_removed, wMAE),
    y = wMAE
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = ci.lb,
      ymax = ci.ub
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept = full_wMAE,
    linetype = "dashed"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "Pooled wMAE",
    title = "Leave-one-cohort-out sensitivity of pooled wMAE"
  ) +
  theme_minimal(base_size = 14)

dev.off()
# --- caterpillar/forest plot? ---

# full analysis log

sink(log_file)

cat("Brain wMAE multilevel meta-analysis\n")
cat("===================================\n\n")
cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_brain, "\n\n")

cat("Model retention\n")
cat("=============\n")
print(brain_retention)

cat("\n\nModel summary\n")
cat("=============\n")
print(summary(brain_wMAE_mv_global))

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_cohort)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_summary)

cat("\n\nDescriptives\n")
cat("============\n")
print(perf_brain_sub_desc, n = Inf, width = Inf)

sink()

###### wMAE modelwise PURE ------

# --- set-up ---
run_id_modelwise <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_brain_modelwise<- file.path(results_dir, paste0("perf_brain_wMAE_PURE_noDunedin_modelwise_", run_id_modelwise))

dir.create(results_dir_brain_modelwise, recursive = TRUE, showWarnings = FALSE)

log_file_modelwise <- file.path(results_dir_brain_modelwise, "analysis_log_modelwise.txt")

# --- meta ---
brain_wMAE_mv_modelwise <- rma.mv(
  yi   = wMAE_test,
  V    = wMAE_SE_boot^2,
  mods = ~ model - 1,
  random = ~ 1 | cohort/timepoint,
  data   = perf_brain_sub,
  method = "REML"
)

summary(brain_wMAE_mv_modelwise)

saveRDS(
  brain_wMAE_mv_modelwise,
  file.path(results_dir_brain_modelwise, "model_brain_wMAE_mv_modelwise.rds")
)

capture.output(
  summary(brain_wMAE_mv_modelwise),
  file = file.path(results_dir_brain_modelwise, "modelwise_model_summary.txt")
)

# --- extract modelwise estimates ---

brain_wMAE_est <- tibble(
  term  = rownames(brain_wMAE_mv_modelwise$b),
  wMAE     = as.numeric(brain_wMAE_mv_modelwise$b),
  se    = brain_wMAE_mv_modelwise$se,
  ci.lb = brain_wMAE_mv_modelwise$ci.lb,
  ci.ub = brain_wMAE_mv_modelwise$ci.ub,
  pval  = brain_wMAE_mv_modelwise$pval
) %>%
  mutate(model = str_remove(term, "^model")) %>%
  select(-term) 


brain_wMAE_est

# --- extract effect sizes ---
X <- model.matrix(brain_wMAE_mv_modelwise)

effect_sizes_per_model <- tibble(
  model = str_remove(colnames(X), "^model"),
  k = colSums(X != 0)
)

brain_wMAE_est <- brain_wMAE_est %>%
  left_join(effect_sizes_per_model, by = "model")

write.csv(
  brain_wMAE_est,
  file.path(results_dir_brain_modelwise, "modelwise_wMAE_estimates.csv"),
  row.names = FALSE
)


# --- heterogeneity ---
I2 <- calc_I2_mv(
  brain_wMAE_mv_modelwise,
  component_names = c("cohort", "timepoint_within_cohort")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_brain_modelwise, "heterogeneity_I2.csv"), row.names = FALSE)

# --- LOO ---

cohorts <- unique(perf_brain_sub$cohort)

loo_cohort_modelwise <- lapply(cohorts, function(cc) {
  
  tmp <- subset(perf_brain_sub, cohort != cc)
  
  fit <- rma.mv(
    yi   = wMAE_test,
    V    = wMAE_SE_boot^2,
    mods = ~ model - 1,
    random = ~ 1 | cohort/timepoint,
    data = tmp,
    method = "REML"
  )
  
  tibble(
    cohort_removed = cc,
    term = rownames(fit$b),
    wMAE = as.numeric(fit$b),
    se = fit$se,
    ci.lb = fit$ci.lb,
    ci.ub = fit$ci.ub,
    pval = fit$pval
  ) %>%
    mutate(model = str_remove(term, "^model")) %>%
    select(cohort_removed, model, wMAE, se, ci.lb, ci.ub, pval)
}) %>%
  bind_rows()

write.csv(
  loo_cohort_modelwise,
  file.path(results_dir_brain_modelwise, "modelwise_leave_one_cohort_out.csv"),
  row.names = FALSE
)

modelwise_loo_summary <- loo_cohort_modelwise %>%
  group_by(model) %>%
  summarise(
    min_wMAE = min(wMAE, na.rm = TRUE),
    max_wMAE = max(wMAE, na.rm = TRUE),
    range_wMAE = max_wMAE - min_wMAE,
    sd_wMAE = sd(wMAE, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  modelwise_loo_summary,
  file.path(results_dir_brain_modelwise, "modelwise_leave_one_cohort_out_summary.csv"),
  row.names = FALSE
)

# --- loo plot ---
png(
  filename = file.path(results_dir_brain_modelwise, "plot_modelwise_wMAE_loo.png"),
  width = 2200,
  height = 1400,
  res = 300
)

ggplot(
  loo_cohort_modelwise,
  aes(
    x = reorder(cohort_removed, wMAE),
    y = wMAE
  )
) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = ci.lb, ymax = ci.ub),
    width = 0.15
  ) +
  facet_wrap(~ model, scales = "free_y") +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "Pooled wMAE",
    title = "Modelwise leave-one-cohort-out sensitivity of pooled brain wMAE"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- log ---

sink(log_file_modelwise)

cat("Brain wMAE modelwise multilevel meta-analysis\n")
cat("============================================\n\n")
cat("Run ID:", run_id_modelwise, "\n")
cat("Output directory:", results_dir_brain_modelwise, "\n\n")

cat("Modelwise model summary\n")
cat("=======================\n")
print(summary(brain_wMAE_mv_modelwise))

cat("\n\nModelwise estimates\n")
cat("===================\n")
print(brain_wMAE_est, n = Inf, width = Inf)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nModelwise leave-one-cohort-out results\n")
cat("=====================================\n")
print(loo_cohort_modelwise, n = Inf, width = Inf)

cat("\n\nModelwise leave-one-cohort-out summary\n")
cat("=====================================\n")
print(modelwise_loo_summary, n = Inf, width = Inf)

cat("\n\nSession info\n")
cat("============\n")
print(sessionInfo())

sink()

###### wMAE age global PURE ------
# all but DunedinPACNI

# --- set up ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_brain_global_ageint <- file.path(results_dir, paste0("perf_brain_wMAE_PURE_noDunedin_global_ageint_", run_id))

dir.create(results_dir_brain_global_ageint, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_brain_global_ageint, "analysis_log_global_ageint.txt")

# --- data prep ---

# center age 
perf_brain_sub <- perf_brain_sub %>%
  mutate(mean_age_c = mean_age - mean(mean_age, na.rm = TRUE))

# --- meta ---

brain_wMAE_mv_global_ageint <- rma.mv(
  yi = wMAE_test,
  V  = wMAE_SE_boot^2,
  mods = ~ mean_age_c,
  random = list(
    ~ 1 | cohort/timepoint,  # timepoints nested within cohorts as soon as more data is included change to cohort/timepoint
    ~ 1 | model  
  ),
  data   = perf_brain_sub,
  method = "REML",
  verbose = TRUE
)

summary(brain_wMAE_mv_global_ageint)

saveRDS(
  brain_wMAE_mv_global_ageint,
  file.path(results_dir_brain_global_ageint, "model_brain_wMAE_mv_global_ageint.rds")
)

capture.output(
  summary(brain_wMAE_mv_global_ageint),
  file = file.path(results_dir_brain_global_ageint, "model_summary.txt")
)

# --- predicted values ---

brain_age_est <- tibble(
  term = names(coef(brain_wMAE_mv_global_ageint)),
  estimate = as.numeric(coef(brain_wMAE_mv_global_ageint)),
  se = brain_wMAE_mv_global_ageint$se,
  ci.lb = brain_wMAE_mv_global_ageint$ci.lb,
  ci.ub = brain_wMAE_mv_global_ageint$ci.ub,
  pval = brain_wMAE_mv_global_ageint$pval
)

write.csv(
  brain_age_est,
  file.path(results_dir_brain_global_ageint, "age_moderator_results.csv"),
  row.names = FALSE
)

# --- heterogeneity ---
I2 <- calc_I2_mv(
  brain_wMAE_mv_global_ageint,
  component_names = c("cohort", "timepoint_within_cohort", "model")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_brain_global_ageint, "heterogeneity_I2.csv"), row.names = FALSE)

# --- plotting (should be adapted) ---

png(filename = file.path(results_dir_brain_global_ageint, "plot_meanage_wMAE_allbD_rawdata.png"), width = 1800,height = 1200,res = 300)
ggplot(perf_brain_sub,
       aes(mean_age, wMAE_test)) +
  geom_point(alpha = .7) +
  geom_smooth(method = "lm") +
  labs(
    x = "Mean sample age",
    y = "Weighted MAE (MAE ÷ age range)"
  ) + 
  theme_minimal()
dev.off()

# if we wanted to plot the meta results (raw datapoints, superimposed cohort-timepoint means & regression line based on meta results )
newdat <- data.frame(
  mean_age = seq(
    min(perf_brain_sub$mean_age),
    max(perf_brain_sub$mean_age),
    length.out = 100
  )
) %>%
  mutate(
    mean_age_c = mean_age - mean(perf_brain_sub$mean_age, na.rm = TRUE)
  )

pred <- predict(
  brain_wMAE_mv_global_ageint,
  newmods = newdat$mean_age_c
)

newdat$pred  <- pred$pred
newdat$ci.lb <- pred$ci.lb
newdat$ci.ub <- pred$ci.ub

cohort_time_means <- perf_brain_sub %>%
  group_by(cohort, timepoint) %>%
  summarise(
    mean_age = mean(mean_age, na.rm = TRUE),
    mean_wMAE = mean(wMAE_test),
    .groups = "drop"
  )

png(filename = file.path(results_dir_brain_global_ageint, "plot_meanage_wMAE_allbD_metadata.png"), width = 1800,height = 1200,res = 300)
ggplot() +
  
  geom_point(
    data = perf_brain_sub,
    aes(mean_age, wMAE_test),
    alpha = .25
  ) +
  
  geom_point(
    data = cohort_time_means,
    aes(mean_age, mean_wMAE),
    size = 3
  ) +
  
  geom_line(
    data = newdat,
    aes(mean_age, pred),
    linewidth = 1
  ) +
  
  geom_ribbon(
    data = newdat,
    aes(mean_age, ymin = ci.lb, ymax = ci.ub),
    alpha = .15
  ) +
  
  labs(
    x = "Mean sample age",
    y = "Weighted MAE",
    title = "Meta-regression of model performance on mean age"
  ) + 
  theme_minimal()
dev.off()

# --- LOO ---

cohorts <- unique(perf_brain_sub$cohort)

loo_age <- lapply(cohorts, function(cc) {
  
  tmp <- perf_brain_sub %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = wMAE_test,
    V  = wMAE_SE_boot^2,
    mods = ~ mean_age_c,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model
    ),
    data = tmp,
    method = "REML"
  )
  
  data.frame(
    cohort_removed = cc,
    intercept = coef(fit)[1],
    slope = coef(fit)[2],
    intercept_se = sqrt(vcov(fit)[1, 1]),
    slope_se = sqrt(vcov(fit)[2, 2]),
    intercept_ci.lb = coef(fit)[1] - 1.96 * sqrt(vcov(fit)[1, 1]),
    intercept_ci.ub = coef(fit)[1] + 1.96 * sqrt(vcov(fit)[1, 1]),
    slope_ci.lb = coef(fit)[2] - 1.96 * sqrt(vcov(fit)[2, 2]),
    slope_ci.ub = coef(fit)[2] + 1.96 * sqrt(vcov(fit)[2, 2]),
    intercept_pval = fit$pval[1],
    slope_pval = fit$pval[2]
  )
})

loo_age_df <- do.call(rbind, loo_age) %>%
  as_tibble()

loo_age_summary <- tibble(
  min_intercept   = min(loo_age_df$intercept, na.rm = TRUE),
  max_intercept   = max(loo_age_df$intercept, na.rm = TRUE),
  range_intercept = max_intercept - min_intercept,
  sd_intercept    = sd(loo_age_df$intercept, na.rm = TRUE),
  
  min_slope       = min(loo_age_df$slope, na.rm = TRUE),
  max_slope       = max(loo_age_df$slope, na.rm = TRUE),
  range_slope     = max_slope - min_slope,
  sd_slope        = sd(loo_age_df$slope, na.rm = TRUE)
)

write.csv(
  loo_age_df,
  file.path(results_dir_brain_global_ageint, "leave_one_cohort_out_age_model.csv"),
  row.names = FALSE
)

write.csv(
  loo_age_summary,
  file.path(results_dir_brain_global_ageint, "leave_one_cohort_out_age_model_summary.csv"),
  row.names = FALSE
)

png(filename = file.path(results_dir_brain_global_ageint, "plot_meanage_wMAE_allbD_loo.png"), width = 1800,height = 1200,res = 300)
ggplot(loo_age_df,
       aes(x = reorder(cohort_removed, slope),
           y = slope)) +
  geom_point(size = 3) +
  geom_hline(
    yintercept = coef(brain_wMAE_mv_global_ageint)[2],
    linetype = "dashed"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "Estimated age slope",
    title = "Leave-one-cohort-out sensitivity of age effect"
  ) + 
  theme_minimal()
dev.off()

# --- log ---

sink(log_file)

cat("Brain wMAE multilevel meta-analysis - global age interaction\n")
cat("===================================\n\n")
cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_brain_global_ageint, "\n\n")

cat("Model summary\n")
cat("=============\n")
print(summary(brain_wMAE_mv_global_ageint))

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_age_df)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_age_summary)

cat("\n\nDescriptives\n")
cat("============\n")
print(perf_brain_sub_desc, n = Inf, width = Inf)

sink()
