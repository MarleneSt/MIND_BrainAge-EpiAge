#### META-ANALYSES SCRIPT PEARSON CORRELATION
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


#### define important folders ----
home <- "data" # [REPRODUCIBILITY] relative path
data_dir <- file.path(home)
results_dir <- "results" # [REPRODUCIBILITY] relative path


#### Correlation analyses ----

# load data 
perf_brain <- readRDS(file.path(results_dir, "perf_brain.rds"))
perf_epi <- readRDS(file.path(results_dir, "perf_epi.rds"))

# set var type 
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

##### epi pearson (all clocks) -----
###### epi pearson global ------

# --- prep ---
# exclude specific models if needed 
#perf_epi   <- perf_epi   %>% filter(!model %in% c("DunedinPACE", "DNAmTL")) %>% droplevels()

# rename gestational clocks to not have _conv as suffix 
perf_epi <- perf_epi %>%
  mutate(
    model = recode(
      model,
      "Bohlin_conv" = "Bohlin",
      "Knight_conv" = "Knight",
      "EPIC_conv" = "EPIC"
    )
  )

# --- set up folder ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi_r <- file.path(results_dir, paste0("perf_epi_pearson_global_", run_id))

dir.create(results_dir_epi_r, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_epi_r, "analysis_log.txt")

# --- descriptives ---

write.csv(perf_epi, file.path(results_dir_epi_r, "perf_epi_descr_detailed_withDunedin_DNAmTL.csv"), row.names=FALSE)

perf_epi_desc <- perf_epi %>%
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

write.csv(perf_epi_desc, file.path(results_dir_epi_r, "model_perf_epi_desc_withDunedin_DNAmTL.csv.csv"), row.names = FALSE)

# --- meta ---
epi_pearson_mv_global <- rma.mv(
  yi = pearson_z,
  V  = pearson_var,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data   = perf_epi,
  method = "REML"
)

summary(epi_pearson_mv_global)

# heterogeneity estimates

W <- diag(1/epi_pearson_mv_global$vi)
X <- model.matrix(epi_pearson_mv_global)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_total <- 100 * sum(epi_pearson_mv_global$sigma2) / (sum(epi_pearson_mv_global$sigma2) + (epi_pearson_mv_global$k-epi_pearson_mv_global$p)/sum(diag(P)))


I2_components <- 100 * epi_pearson_mv_global$sigma2 / (
  sum(epi_pearson_mv_global$sigma2) +
    (epi_pearson_mv_global$k - epi_pearson_mv_global$p) /
    sum(diag(P))
)

names(I2_components) <- c(
  "cohort",
  "timepoint_within_cohort",
  "model"
)

I2_total
I2_components

# Back-transform global pooled estimate
overall_est <- data.frame(
  model = "Pooled",
  pearson_z = epi_pearson_mv_global$beta[1, 1],
  pearson_r = transf.ztor(epi_pearson_mv_global$beta[1, 1]),
  ci.lb = transf.ztor(epi_pearson_mv_global$ci.lb[1]),
  ci.ub = transf.ztor(epi_pearson_mv_global$ci.ub[1]),
  pval = epi_pearson_mv_global$pval[1]
)

write.csv(
  overall_est,
  file.path(results_dir_epi_r, "overall_pearson_estimate.csv"),
  row.names = FALSE
)

saveRDS(
  epi_pearson_mv_global,
  file.path(results_dir_epi_r, "model_epi_pearson_mv_global.rds")
)

capture.output(
  summary(epi_pearson_mv_global),
  file = file.path(results_dir_epi_r, "model_summary.txt")
)

# --- LOO cohort ---

cohorts <- unique(perf_epi$cohort)

loo_cohort <- lapply(cohorts, function(cc) {
  
  tmp <- subset(perf_epi, cohort != cc)
  
  fit <- rma.mv(
    yi = pearson_z,
    V  = pearson_var,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model
    ),
    data = tmp,
    method = "REML"
  )
  
  z_est <- coef(fit)[1]
  z_se  <- sqrt(vcov(fit)[1, 1])
  
  data.frame(
    cohort_removed = cc,
    
    pearson_z = z_est,
    pearson_z_se = z_se,
    z_ci.lb = z_est - 1.96 * z_se,
    z_ci.ub = z_est + 1.96 * z_se,
    
    pearson_r = transf.ztor(z_est),
    ci.lb = transf.ztor(z_est - 1.96 * z_se),
    ci.ub = transf.ztor(z_est + 1.96 * z_se)
  )
}) %>%
  bind_rows()

loo_summary <- tibble(
  min_r = min(loo_cohort$pearson_r, na.rm = TRUE),
  max_r = max(loo_cohort$pearson_r, na.rm = TRUE),
  range_r = max_r - min_r,
  sd_r = sd(loo_cohort$pearson_r, na.rm = TRUE),
  
  min_z = min(loo_cohort$pearson_z, na.rm = TRUE),
  max_z = max(loo_cohort$pearson_z, na.rm = TRUE),
  range_z = max_z - min_z,
  sd_z = sd(loo_cohort$pearson_z, na.rm = TRUE)
)

write.csv(
  loo_cohort,
  file.path(results_dir_epi_r, "leave_one_cohort_out.csv"),
  row.names = FALSE
)

write.csv(
  loo_summary,
  file.path(results_dir_epi_r, "leave_one_cohort_out_summary.csv"),
  row.names = FALSE
)

# --- global LOO plot: Pearson r scale ---

full_r <- transf.ztor(coef(epi_pearson_mv_global)[1])

png(
  filename = file.path(results_dir_epi_r, "plot_global_pearson_loo.png"),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_cohort,
  aes(
    x = reorder(cohort_removed, pearson_r),
    y = pearson_r
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = ci.lb, ymax = ci.ub),
    width = 0.15
  ) +
  geom_hline(
    yintercept = full_r,
    linetype = "dashed"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "Pooled Pearson r",
    title = "Leave-one-cohort-out sensitivity of pooled Pearson correlation"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# full analysis log

sink(log_file)

cat("Epi Pearson multilevel meta-analysis\n")
cat("====================================\n\n")
cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_epi_r, "\n\n")

cat("Model summary, Fisher-z scale\n")
cat("=============================\n")
print(summary(epi_pearson_mv_global))

cat("\n\nModel summary - back-transformed to Pearson r\n")
cat("=============================================\n")
print(overall_est)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")

cat("Total I2 (%):\n")
print(I2_total)

cat("\nVariance decomposition (% of total heterogeneity):\n")
print(I2_components)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_cohort)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_summary)

cat("\n\nDescriptives\n")
cat("============\n")
print(perf_epi_desc, n = Inf, width = Inf)

sink()

###### epi Pearson modelwise ------

# --- set up folder ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi_pearson_modelwise <- file.path(
  results_dir,
  paste0("perf_epi_pearson_modelwise_", run_id)
)

dir.create(results_dir_epi_pearson_modelwise, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_epi_pearson_modelwise, "analysis_modelwise_log.txt")

# --- meta ---

epi_pearson_mv_modelwise <- rma.mv(
  yi = pearson_z,
  V  = pearson_var,
  mods = ~ model - 1,
  random = list(
    ~ 1 | cohort/timepoint
  ),
  data = perf_epi,
  method = "REML"
)

summary(epi_pearson_mv_modelwise)


# heterogeneity estimates

W <- diag(1/epi_pearson_mv_modelwise$vi)
X <- model.matrix(epi_pearson_mv_modelwise)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_total <- 100 * sum(epi_pearson_mv_modelwise$sigma2) / (sum(epi_pearson_mv_modelwise$sigma2) + (epi_pearson_mv_modelwise$k-epi_pearson_mv_modelwise$p)/sum(diag(P)))


I2_components <- 100 * epi_pearson_mv_modelwise$sigma2 / (
  sum(epi_pearson_mv_modelwise$sigma2) +
    (epi_pearson_mv_modelwise$k - epi_pearson_mv_modelwise$p) /
    sum(diag(P))
)

names(I2_components) <- c(
  "cohort",
  "timepoint_within_cohort"
)

I2_total
I2_components

# --- back-transform modelwise pooled estimates ---

epi_pearson_est <- data.frame(
  model = gsub("^model", "", rownames(epi_pearson_mv_modelwise$beta)),
  pearson_z = epi_pearson_mv_modelwise$beta[, 1],
  pearson_r = transf.ztor(epi_pearson_mv_modelwise$beta[, 1]),
  ci.lb = transf.ztor(epi_pearson_mv_modelwise$ci.lb),
  ci.ub = transf.ztor(epi_pearson_mv_modelwise$ci.ub),
  pval = epi_pearson_mv_modelwise$pval
)

X <- model.matrix(epi_pearson_mv_modelwise)

model_cols <- colnames(X)[
  str_detect(colnames(X), "^model") &
    !str_detect(colnames(X), ":")
]

effect_sizes_per_model <- tibble(
  model = str_remove(model_cols, "^model"),
  k = colSums(X[, model_cols, drop = FALSE] != 0)
)


epi_pearson_est <- epi_pearson_est %>%
  left_join(effect_sizes_per_model, by = "model")

epi_pearson_est

write.csv(
  epi_pearson_est,
  file.path(results_dir_epi_pearson_modelwise, "modelwise_pearson_estimates.csv"),
  row.names = FALSE
)

saveRDS(
  epi_pearson_mv_modelwise,
  file.path(results_dir_epi_pearson_modelwise, "model_epi_pearson_mv_modelwise.rds")
)

capture.output(
  summary(epi_pearson_mv_modelwise),
  file = file.path(results_dir_epi_pearson_modelwise, "model_summary.txt")
)

# --- LOO cohort for Fisher-z Pearson model ---

cohorts <- unique(perf_epi$cohort)

loo_cohort_modelwise <- lapply(cohorts, function(cc) {
  
  tmp <- subset(perf_epi, cohort != cc)
  
  fit <- rma.mv(
    yi = pearson_z,
    V  = pearson_var,
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
    
    # Fisher-z scale estimates
    pearson_z = as.numeric(fit$b),
    pearson_z_se = fit$se,
    z_ci.lb = fit$ci.lb,
    z_ci.ub = fit$ci.ub,
    pval = fit$pval
  ) %>%
    mutate(
      model = str_remove(term, "^model"),
      
      # back-transformed estimates: Pearson r
      pearson_r = transf.ztor(pearson_z),
      ci.lb = transf.ztor(z_ci.lb),
      ci.ub = transf.ztor(z_ci.ub)
    ) %>%
    select(
      cohort_removed,
      model,
      pearson_z,
      pearson_z_se,
      z_ci.lb,
      z_ci.ub,
      pearson_r,
      ci.lb,
      ci.ub,
      pval
    )
}) %>%
  bind_rows()

write.csv(
  loo_cohort_modelwise,
  file.path(results_dir_epi_pearson_modelwise, "modelwise_leave_one_cohort_out.csv"),
  row.names = FALSE
)

modelwise_loo_summary <- loo_cohort_modelwise %>%
  group_by(model) %>%
  summarise(
    min_r = min(pearson_r, na.rm = TRUE),
    max_r = max(pearson_r, na.rm = TRUE),
    range_r = max_r - min_r,
    sd_r = sd(pearson_r, na.rm = TRUE),
    
    min_z = min(pearson_z, na.rm = TRUE),
    max_z = max(pearson_z, na.rm = TRUE),
    range_z = max_z - min_z,
    sd_z = sd(pearson_z, na.rm = TRUE),
    
    .groups = "drop"
  )

write.csv(
  modelwise_loo_summary,
  file.path(results_dir_epi_pearson_modelwise, "modelwise_leave_one_cohort_out_summary.csv"),
  row.names = FALSE
)

# --- optional LOO plot by model ---

png(
  filename = file.path(results_dir_epi_pearson_modelwise, "plot_modelwise_pearson_loo.png"),
  width = 2200,
  height = 1400,
  res = 300
)

ggplot(
  loo_cohort_modelwise,
  aes(
    x = reorder(cohort_removed, pearson_r),
    y = pearson_r
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
    y = "Pooled Pearson r",
    title = "Modelwise leave-one-cohort-out sensitivity of pooled Pearson correlations"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- full analysis log ---

sink(log_file)

cat("Epi Pearson multilevel meta-analysis modelwise\n")
cat("==============================================\n\n")
cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_epi_pearson_modelwise, "\n\n")

cat("Model summary, Fisher-z scale\n")
cat("=============================\n")
print(summary(epi_pearson_mv_modelwise))

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")

cat("Total I2 (%):\n")
print(I2_total)

cat("\nVariance decomposition (% of total heterogeneity):\n")
print(I2_components)

cat("\n\nModel summary - back-transformed to Pearson r\n")
cat("=============================================\n")
print(epi_pearson_est)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_cohort_modelwise, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(modelwise_loo_summary, n = Inf, width = Inf)

sink()

###### epi Pearson age global ------

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi_pearson_age <- file.path(
  results_dir,
  paste0("perf_epi_pearson_global_ageint_", run_id)
)

dir.create(results_dir_epi_pearson_age, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_epi_pearson_age, "analysis_log_global_ageint.txt")

# --- data prep ---

perf_epi <- perf_epi %>%
  mutate(mean_age_c = mean_age - mean(mean_age, na.rm = TRUE))

# remove gestational clocks for age moderation
perf_epi_nogest <- perf_epi %>%
  filter(!model %in% c("Bohlin", "Knight", "EPIC")) %>%
  droplevels()

# --- meta-regression ---

epi_pearson_mv_global_ageint <- rma.mv(
  yi = pearson_z,
  V  = pearson_var,
  mods = ~ mean_age_c,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = perf_epi_nogest,
  method = "REML"
)

summary(epi_pearson_mv_global_ageint)

# --- effect sizes per model ---
effect_sizes_per_model <- perf_epi_nogest %>%
  count(model, name = "k")

write.csv(
  effect_sizes_per_model,
  file.path(results_dir_epi_pearson_age, "effect_sizes_per_model.csv"),
  row.names = FALSE
)


# --- I2 calculation ---

W <- diag(1 / epi_pearson_mv_global_ageint$vi)
X <- model.matrix(epi_pearson_mv_global_ageint)

P <- W - W %*% X %*%
  solve(t(X) %*% W %*% X) %*%
  t(X) %*% W

I2_total_age <- 100 * sum(epi_pearson_mv_global_ageint$sigma2) / (
  sum(epi_pearson_mv_global_ageint$sigma2) +
    (epi_pearson_mv_global_ageint$k - epi_pearson_mv_global_ageint$p) /
    sum(diag(P))
)

I2_components_age <- 100 * epi_pearson_mv_global_ageint$sigma2 / (
  sum(epi_pearson_mv_global_ageint$sigma2) +
    (epi_pearson_mv_global_ageint$k - epi_pearson_mv_global_ageint$p) /
    sum(diag(P))
)

names(I2_components_age) <- c(
  "cohort",
  "timepoint_within_cohort",
  "model"
)

heterogeneity_partition_age <- 100 * epi_pearson_mv_global_ageint$sigma2 /
  sum(epi_pearson_mv_global_ageint$sigma2)

names(heterogeneity_partition_age) <- c(
  "cohort",
  "timepoint_within_cohort",
  "model"
)

# --- extract model results ---

age_model_est <- tibble(
  term = names(coef(epi_pearson_mv_global_ageint)),
  estimate_z = as.numeric(coef(epi_pearson_mv_global_ageint)),
  se_z = epi_pearson_mv_global_ageint$se,
  ci.lb_z = epi_pearson_mv_global_ageint$ci.lb,
  ci.ub_z = epi_pearson_mv_global_ageint$ci.ub,
  pval = epi_pearson_mv_global_ageint$pval
) %>%
  mutate(
    estimate_r = if_else(term == "intrcpt", transf.ztor(estimate_z), NA_real_),
    ci.lb_r = if_else(term == "intrcpt", transf.ztor(ci.lb_z), NA_real_),
    ci.ub_r = if_else(term == "intrcpt", transf.ztor(ci.ub_z), NA_real_)
  )

write.csv(
  age_model_est,
  file.path(results_dir_epi_pearson_age, "age_moderator_results.csv"),
  row.names = FALSE
)

saveRDS(
  epi_pearson_mv_global_ageint,
  file.path(results_dir_epi_pearson_age, "model_epi_pearson_mv_global_ageint.rds")
)

# --- prediction grid ---

mean_age_mean <- mean(perf_epi_nogest$mean_age, na.rm = TRUE)

newdat <- data.frame(
  mean_age = seq(
    min(perf_epi_nogest$mean_age, na.rm = TRUE),
    max(perf_epi_nogest$mean_age, na.rm = TRUE),
    length.out = 100
  )
) %>%
  mutate(mean_age_c = mean_age - mean_age_mean)

pred <- predict(
  epi_pearson_mv_global_ageint,
  newmods = newdat$mean_age_c
)

newdat <- newdat %>%
  mutate(
    pred_z = pred$pred,
    z_ci.lb = pred$ci.lb,
    z_ci.ub = pred$ci.ub,
    pred_r = transf.ztor(pred_z),
    ci.lb = transf.ztor(z_ci.lb),
    ci.ub = transf.ztor(z_ci.ub)
  )

png(
  filename = file.path(results_dir_epi_pearson_age, "plot_meanage_pearson_meta_regression.png"),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot() +
  geom_point(
    data = perf_epi_nogest,
    aes(x = mean_age, y = Pearson, colour = cohort),
    alpha = .25,
    size = 2
  ) +
  geom_line(
    data = newdat,
    aes(x = mean_age, y = pred_r),
    linewidth = 1.2,
    colour = "black"
  ) +
  geom_ribbon(
    data = newdat,
    aes(x = mean_age, ymin = ci.lb, ymax = ci.ub),
    alpha = .15
  ) +
  labs(
    title = "Meta-regression of Pearson correlation on mean age",
    x = "Mean sample age",
    y = "Predicted Pearson r",
    colour = "Cohort"
  ) +
  theme_minimal(base_size = 16)

dev.off()

# --- LOO age moderation ---

cohorts <- unique(perf_epi_nogest$cohort)

loo_age <- lapply(cohorts, function(cc) {
  
  tmp <- perf_epi_nogest %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = pearson_z,
    V  = pearson_var,
    mods = ~ mean_age_c,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model
    ),
    data = tmp,
    method = "REML"
  )
  
  z_intercept <- coef(fit)[1]
  z_slope <- coef(fit)[2]
  
  z_intercept_se <- sqrt(vcov(fit)[1, 1])
  z_slope_se <- sqrt(vcov(fit)[2, 2])
  
  data.frame(
    cohort_removed = cc,
    
    intercept_z = z_intercept,
    slope_z_per_year = z_slope,
    intercept_z_se = z_intercept_se,
    slope_z_se = z_slope_se,
    
    intercept_r = transf.ztor(z_intercept),
    intercept_ci.lb = transf.ztor(z_intercept - 1.96 * z_intercept_se),
    intercept_ci.ub = transf.ztor(z_intercept + 1.96 * z_intercept_se),
    
    slope_z_ci.lb = z_slope - 1.96 * z_slope_se,
    slope_z_ci.ub = z_slope + 1.96 * z_slope_se
  )
}) %>%
  bind_rows()

loo_age_summary <- tibble(
  min_intercept_r = min(loo_age$intercept_r, na.rm = TRUE),
  max_intercept_r = max(loo_age$intercept_r, na.rm = TRUE),
  range_intercept_r = max_intercept_r - min_intercept_r,
  sd_intercept_r = sd(loo_age$intercept_r, na.rm = TRUE),
  
  min_slope_z = min(loo_age$slope_z_per_year, na.rm = TRUE),
  max_slope_z = max(loo_age$slope_z_per_year, na.rm = TRUE),
  range_slope_z = max_slope_z - min_slope_z,
  sd_slope_z = sd(loo_age$slope_z_per_year, na.rm = TRUE)
)

write.csv(
  loo_age,
  file.path(results_dir_epi_pearson_age, "leave_one_cohort_out_age_model.csv"),
  row.names = FALSE
)

write.csv(
  loo_age_summary,
  file.path(results_dir_epi_pearson_age, "leave_one_cohort_out_age_model_summary.csv"),
  row.names = FALSE
)

full_slope_z <- coef(epi_pearson_mv_global_ageint)[2]

png(
  filename = file.path(results_dir_epi_pearson_age, "plot_meanage_pearson_loo_slope.png"),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_age,
  aes(
    x = reorder(cohort_removed, slope_z_per_year),
    y = slope_z_per_year
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = slope_z_ci.lb, ymax = slope_z_ci.ub),
    width = 0.15
  ) +
  geom_hline(yintercept = full_slope_z, linetype = "dashed") +
  geom_hline(yintercept = 0, colour = "grey60") +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "Change in Fisher-z Pearson correlation per year",
    title = "Leave-one-cohort-out sensitivity of age effect"
  ) +
  theme_minimal(base_size = 14)

dev.off()

sink(log_file)

cat("Epi Pearson multilevel meta-analysis - global age moderator\n")
cat("============================================================\n\n")
cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_epi_pearson_age, "\n\n")

cat("Model summary, Fisher-z scale\n")
cat("=============================\n")
print(summary(epi_pearson_mv_global_ageint))

cat("\n\nModel results\n")
cat("=============\n")
print(age_model_est, n = Inf, width = Inf)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")

cat("Total I2 (%):\n")
print(I2_total_age)

cat("\nI2 components (% of total observed variance):\n")
print(I2_components_age)
cat("\n=each random effect component as % of heterogeneity+sampling variance\n, components sum to I^2 but not to 100%\n")

cat("\nHeterogeneity partition (% of random-effect heterogeneity only):\n")
print(heterogeneity_partition_age)
cat("\n=each random effect component as % of heterogeneity only\n, components sum to 100%\n")

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_age)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_age_summary, n = Inf, width = Inf)

sink()

###### epi Pearson age modelwise ------

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi_pearson_modelwise_age <- file.path(
  results_dir,
  paste0("perf_epi_pearson_modelwise_ageint_", run_id)
)

dir.create(results_dir_epi_pearson_modelwise_age, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  results_dir_epi_pearson_modelwise_age,
  "analysis_log_modelwise_ageint.txt"
)

# --- data prep ---

perf_epi_nogest <- perf_epi %>%
  mutate(mean_age_c = mean_age - mean(mean_age, na.rm = TRUE)) %>%
  filter(!model %in% c("Bohlin", "Knight", "EPIC")) %>%
  droplevels()

# --- meta: model-specific intercepts and model-specific age slopes ---

epi_pearson_mv_modelwise_ageint <- rma.mv(
  yi = pearson_z,
  V  = pearson_var,
  mods = ~ model + model:mean_age_c - 1,
  random = list(
    ~ 1 | cohort/timepoint
  ),
  data = perf_epi_nogest,
  method = "REML"
)

summary(epi_pearson_mv_modelwise_ageint)

# compare common age slope model vs model-specific age slopes
anova(
  epi_pearson_mv_global_ageint,
  epi_pearson_mv_modelwise_ageint,
  refit = TRUE
)

# --- I2 calculation ---

W <- diag(1 / epi_pearson_mv_modelwise_ageint$vi)
X <- model.matrix(epi_pearson_mv_modelwise_ageint)

P <- W - W %*% X %*%
  solve(t(X) %*% W %*% X) %*%
  t(X) %*% W

I2_total_modelwise_age <- 100 * sum(epi_pearson_mv_modelwise_ageint$sigma2) / (
  sum(epi_pearson_mv_modelwise_ageint$sigma2) +
    (epi_pearson_mv_modelwise_ageint$k - epi_pearson_mv_modelwise_ageint$p) /
    sum(diag(P))
)

I2_components_modelwise_age <- 100 * epi_pearson_mv_modelwise_ageint$sigma2 / (
  sum(epi_pearson_mv_modelwise_ageint$sigma2) +
    (epi_pearson_mv_modelwise_ageint$k - epi_pearson_mv_modelwise_ageint$p) /
    sum(diag(P))
)

names(I2_components_modelwise_age) <- c(
  "cohort",
  "timepoint_within_cohort"
)

heterogeneity_partition_modelwise_age <- 100 * epi_pearson_mv_modelwise_ageint$sigma2 /
  sum(epi_pearson_mv_modelwise_ageint$sigma2)

names(heterogeneity_partition_modelwise_age) <- c(
  "cohort",
  "timepoint_within_cohort"
)

# --- extract modelwise results ---

modelwise_age_est <- tibble(
  term = rownames(epi_pearson_mv_modelwise_ageint$b),
  estimate_z = as.numeric(epi_pearson_mv_modelwise_ageint$b),
  se_z = epi_pearson_mv_modelwise_ageint$se,
  ci.lb_z = epi_pearson_mv_modelwise_ageint$ci.lb,
  ci.ub_z = epi_pearson_mv_modelwise_ageint$ci.ub,
  pval = epi_pearson_mv_modelwise_ageint$pval
) %>%
  mutate(
    type = if_else(str_detect(term, "mean_age_c"), "age_slope", "intercept"),
    model = term %>%
      str_remove("^model") %>%
      str_remove(":mean_age_c") %>%
      str_remove("\\.mean_age_c"),
    
    # intercepts can be back-transformed to Pearson r
    estimate_r = if_else(type == "intercept", transf.ztor(estimate_z), NA_real_),
    ci.lb_r = if_else(type == "intercept", transf.ztor(ci.lb_z), NA_real_),
    ci.ub_r = if_else(type == "intercept", transf.ztor(ci.ub_z), NA_real_)
  )

X <- model.matrix(epi_pearson_mv_modelwise_ageint)

model_cols <- colnames(X)[
  str_detect(colnames(X), "^model") &
    !str_detect(colnames(X), ":")
]

effect_sizes_per_model <- tibble(
  model = str_remove(model_cols, "^model"),
  k = colSums(X[, model_cols, drop = FALSE] != 0)
)

modelwise_age_est <- modelwise_age_est %>%
  left_join(effect_sizes_per_model, by = "model")

write.csv(
  modelwise_age_est,
  file.path(
    results_dir_epi_pearson_modelwise_age,
    "modelwise_age_interaction_results.csv"
  ),
  row.names = FALSE
)

saveRDS(
  epi_pearson_mv_modelwise_ageint,
  file.path(
    results_dir_epi_pearson_modelwise_age,
    "model_epi_pearson_mv_modelwise_ageint.rds"
  )
)

# --- plot model-specific age effects (r scale) ---

mean_age_mean <- mean(perf_epi_nogest$mean_age, na.rm = TRUE)

newdat <- expand.grid(
  mean_age_c = seq(
    min(perf_epi_nogest$mean_age_c, na.rm = TRUE),
    max(perf_epi_nogest$mean_age_c, na.rm = TRUE),
    length.out = 100
  ),
  model = unique(perf_epi_nogest$model)
)

pred <- predict(
  epi_pearson_mv_modelwise_ageint,
  newmods = model.matrix(~ model + model:mean_age_c - 1, newdat)
)

newdat <- newdat %>%
  mutate(
    mean_age = mean_age_c + mean_age_mean,
    
    pred_z = pred$pred,
    z_ci.lb = pred$ci.lb,
    z_ci.ub = pred$ci.ub,
    
    pred_r = transf.ztor(pred_z),
    ci.lb = transf.ztor(z_ci.lb),
    ci.ub = transf.ztor(z_ci.ub)
  )

png(
  filename = file.path(
    results_dir_epi_pearson_modelwise_age,
    "plot_epi_meanage_pearson_modelwise.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(newdat, aes(mean_age, pred_r, colour = model)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = ci.lb, ymax = ci.ub, fill = model),
    alpha = .12,
    colour = NA
  ) +
  labs(
    x = "Mean sample age",
    y = "Predicted Pearson r",
    title = "Model-specific age effects on Pearson correlation"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- LOO modelwise age slopes ---
cohorts <- unique(perf_epi_nogest$cohort)

loo_modelwise_age <- lapply(cohorts, function(cc) {
  
  tmp <- perf_epi_nogest %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = pearson_z,
    V  = pearson_var,
    mods = ~ model + model:mean_age_c - 1,
    random = list(
      ~ 1 | cohort/timepoint
    ),
    data = tmp,
    method = "REML"
  )
  
  tibble(
    cohort_removed = cc,
    term = rownames(fit$b),
    estimate_z = as.numeric(fit$b),
    se_z = fit$se,
    ci.lb_z = fit$ci.lb,
    ci.ub_z = fit$ci.ub,
    pval = fit$pval
  ) %>%
    mutate(
      type = if_else(str_detect(term, "mean_age_c"), "age_slope", "intercept"),
      model = term %>%
        str_remove("^model") %>%
        str_remove(":mean_age_c") %>%
        str_remove("\\.mean_age_c"),
      
      estimate_r = if_else(type == "intercept", transf.ztor(estimate_z), NA_real_),
      ci.lb_r = if_else(type == "intercept", transf.ztor(ci.lb_z), NA_real_),
      ci.ub_r = if_else(type == "intercept", transf.ztor(ci.ub_z), NA_real_)
    )
}) %>%
  bind_rows()

write.csv(
  loo_modelwise_age,
  file.path(
    results_dir_epi_pearson_modelwise_age,
    "modelwise_leave_one_cohort_out_age_interactions.csv"
  ),
  row.names = FALSE
)

loo_modelwise_age_summary <- loo_modelwise_age %>%
  filter(type == "age_slope") %>%
  group_by(model) %>%
  summarise(
    min_slope_z = min(estimate_z, na.rm = TRUE),
    max_slope_z = max(estimate_z, na.rm = TRUE),
    range_slope_z = max_slope_z - min_slope_z,
    sd_slope_z = sd(estimate_z, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  loo_modelwise_age_summary,
  file.path(
    results_dir_epi_pearson_modelwise_age,
    "modelwise_leave_one_cohort_out_age_interactions_summary.csv"
  ),
  row.names = FALSE
)

# --- modelwise LOO plot for age slopes ---

# full-model age slopes
full_modelwise_slopes <- modelwise_age_est %>%
  filter(type == "age_slope") %>%
  select(
    model,
    full_slope_z = estimate_z
  )

# LOO age slopes only
plot_dat_slopes <- loo_modelwise_age %>%
  filter(type == "age_slope") %>%
  left_join(full_modelwise_slopes, by = "model") %>%
  mutate(
    model = factor(model),
    cohort_removed = factor(cohort_removed)
  )

full_modelwise_slopes_plot <- plot_dat_slopes %>%
  distinct(model, full_slope_z)

png(
  filename = file.path(
    results_dir_epi_pearson_modelwise_age,
    "plot_modelwise_meanage_loo_slopes.png"
  ),
  width = 3000,
  height = 2200,
  res = 300
)

p <- ggplot(
  plot_dat_slopes,
  aes(
    x = estimate_z,
    y = reorder(cohort_removed, estimate_z)
  )
) +
  # full-model age slope per clock
  geom_vline(
    data = full_modelwise_slopes_plot,
    aes(xintercept = full_slope_z),
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  
  # null effect
  geom_vline(
    xintercept = 0,
    colour = "grey60"
  ) +
  
  # LOO point estimates
  geom_point(size = 1.8) +
  
  # confidence intervals
  geom_segment(
    aes(
      x = ci.lb_z,
      xend = ci.ub_z,
      yend = reorder(cohort_removed, estimate_z)
    ),
    linewidth = 0.5
  ) +
  
  facet_wrap(~ model, scales = "free_x") +
  
  labs(
    x = "Age slope on Fisher-z scale\n(change in correlation per year)",
    y = "Cohort removed",
    title = "Modelwise leave-one-cohort-out sensitivity of mean-age effects"
  ) +
  
  theme_minimal(base_size = 11)

print(p)

dev.off()

# --- log ---
sink(log_file)

cat("Epi Pearson multilevel meta-analysis - modelwise age interactions\n")
cat("==================================================================\n\n")
cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_epi_pearson_modelwise_age, "\n\n")

cat("Model summary, Fisher-z scale\n")
cat("=============================\n")
print(summary(epi_pearson_mv_modelwise_ageint))

cat("\n\nModel comparison: common age slope vs model-specific age slopes\n")
cat("===============================================================\n")
print(anova(
  epi_pearson_mv_global_ageint,
  epi_pearson_mv_modelwise_ageint,
  refit = TRUE
))

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")

cat("Total I2 (%):\n")
print(I2_total_modelwise_age)

cat("\nI2 components (% of total observed variance):\n")
print(I2_components_modelwise_age)

cat("\nHeterogeneity partition (% of random-effect heterogeneity only):\n")
print(heterogeneity_partition_modelwise_age)

cat("\n\nModelwise age interaction results\n")
cat("=================================\n")
print(modelwise_age_est, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_modelwise_age, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary for age slopes\n")
cat("===========================================\n")
print(loo_modelwise_age_summary, n = Inf, width = Inf)

sink()

###### epi Pearson overlap_epi global (excluding gestational clocks) ------

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi_overlap_binary <- file.path(
  results_dir,
  paste0("perf_epi_pearson_noGest_global_ageoverlap_binary_", run_id)
)

dir.create(
  results_dir_epi_overlap_binary,
  recursive = TRUE,
  showWarnings = FALSE
)

log_file <- file.path(
  results_dir_epi_overlap_binary,
  "analysis_log_ageoverlap_binary.txt"
)

# --- data prep ---

# make overlap var 
perf_epi <- perf_epi %>%
  mutate(
    
    overlap_lower = pmax(min_age, training_age_min_epi),
    overlap_upper = pmin(max_age, training_age_max_epi),
    
    overlap_epi = if_else(
      overlap_upper >= overlap_lower,
      "yes",
      "no"
    ),
    
    # inclusive overlap width
    overlap_width = if_else(
      overlap_epi == "yes",
      overlap_upper - overlap_lower,
      0
    ),
    
    # inclusive study range width
    study_width = max_age - min_age,
    
    overlap_perc_epi =
      100 * overlap_width / study_width
  ) %>%
  select(-overlap_lower, -overlap_upper,
         -overlap_width, -study_width)

# exclude gestational clocks ---

perf_epi_nogest <- perf_epi %>%
  filter(
    !model %in% c("Bohlin", "Knight", "EPIC")
  ) %>%
  mutate(
    overlap_epi = factor(
      overlap_epi,
      levels = c("no", "yes")
    ),
    
    cohort = factor(cohort),
    timepoint = factor(timepoint),
    model = factor(model)
  ) %>%
  droplevels()

# --- meta ---

epi_pearson_mv_overlap_binary <- rma.mv(
  yi = pearson_z,
  V  = pearson_var,
  mods = ~ overlap_epi,
  
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  
  data = perf_epi_nogest,
  method = "REML"
)

summary(epi_pearson_mv_overlap_binary)

# --- I2 calculation ---

W <- diag(1 / epi_pearson_mv_overlap_binary$vi)

X <- model.matrix(epi_pearson_mv_overlap_binary)

P <- W - W %*% X %*%
  solve(t(X) %*% W %*% X) %*%
  t(X) %*% W

I2_total <- 100 * sum(epi_pearson_mv_overlap_binary$sigma2) / (
  sum(epi_pearson_mv_overlap_binary$sigma2) +
    (epi_pearson_mv_overlap_binary$k -
       epi_pearson_mv_overlap_binary$p) /
    sum(diag(P))
)

I2_components <- 100 *
  epi_pearson_mv_overlap_binary$sigma2 / (
    sum(epi_pearson_mv_overlap_binary$sigma2) +
      (epi_pearson_mv_overlap_binary$k -
         epi_pearson_mv_overlap_binary$p) /
      sum(diag(P))
  )

names(I2_components) <- c(
  "cohort",
  "timepoint_within_cohort",
  "model"
)

heterogeneity_partition <- 100 *
  epi_pearson_mv_overlap_binary$sigma2 /
  sum(epi_pearson_mv_overlap_binary$sigma2)

names(heterogeneity_partition) <- c(
  "cohort",
  "timepoint_within_cohort",
  "model"
)

# --- model results ---

overlap_binary_est <- tibble(
  term = names(coef(epi_pearson_mv_overlap_binary)),
  estimate_z = as.numeric(coef(epi_pearson_mv_overlap_binary)),
  se_z = epi_pearson_mv_overlap_binary$se,
  ci.lb_z = epi_pearson_mv_overlap_binary$ci.lb,
  ci.ub_z = epi_pearson_mv_overlap_binary$ci.ub,
  pval = epi_pearson_mv_overlap_binary$pval
)

write.csv(
  overlap_binary_est,
  file.path(
    results_dir_epi_overlap_binary,
    "ageoverlap_binary_results.csv"
  ),
  row.names = FALSE
)

# --- predicted values ---

newdat_overlap_binary <- data.frame(
  overlap_epi = levels(perf_epi_nogest$overlap_epi)
) %>%
  mutate(
    overlap_epi = factor(
      overlap_epi,
      levels = levels(perf_epi_nogest$overlap_epi)
    )
  )

X_new_overlap_binary <- model.matrix(
  ~ overlap_epi,
  data = newdat_overlap_binary
)

pred_overlap_binary <- predict(
  epi_pearson_mv_overlap_binary,
  newmods = X_new_overlap_binary[, -1, drop = FALSE]
)

newdat_overlap_binary <- newdat_overlap_binary %>%
  mutate(
    pred_z = pred_overlap_binary$pred,
    z_ci.lb = pred_overlap_binary$ci.lb,
    z_ci.ub = pred_overlap_binary$ci.ub,
    
    pred_r = transf.ztor(pred_z),
    ci.lb = transf.ztor(z_ci.lb),
    ci.ub = transf.ztor(z_ci.ub)
  )

write.csv(
  newdat_overlap_binary,
  file.path(
    results_dir_epi_overlap_binary,
    "ageoverlap_binary_predicted_values.csv"
  ),
  row.names = FALSE
)

# --- plot ---

png(
  filename = file.path(
    results_dir_epi_overlap_binary,
    "plot_ageoverlap_binary_pearson_meta.png"
  ),
  width = 1600,
  height = 1200,
  res = 300
)

ggplot() +
  
  geom_jitter(
    data = perf_epi_nogest,
    aes(
      x = overlap_epi,
      y = Pearson
    ),
    width = 0.12,
    alpha = 0.25
  ) +
  
  geom_point(
    data = newdat_overlap_binary,
    aes(
      x = overlap_epi,
      y = pred_r
    ),
    size = 4
  ) +
  
  geom_errorbar(
    data = newdat_overlap_binary,
    aes(
      x = overlap_epi,
      ymin = ci.lb,
      ymax = ci.ub
    ),
    width = 0.15
  ) +
  
  scale_y_continuous(
    limits = c(-1, 1)
  ) +
  
  labs(
    x = "Epigenetic clock training age overlap",
    y = "Predicted Pearson r",
    title = "Meta-regression of epigenetic age correlation by training age overlap"
  ) +
  
  theme_minimal(base_size = 14)

dev.off()

# --- save outputs ---

saveRDS(
  epi_pearson_mv_overlap_binary,
  file.path(
    results_dir_epi_overlap_binary,
    "model_epi_pearson_mv_ageoverlap_binary.rds"
  )
)

capture.output(
  summary(epi_pearson_mv_overlap_binary),
  file = file.path(
    results_dir_epi_overlap_binary,
    "model_summary.txt"
  )
)

# --- LOO cohort ---

cohorts <- unique(perf_epi_nogest$cohort)

loo_overlap_binary <- lapply(cohorts, function(cc) {
  
  tmp <- perf_epi_nogest %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = pearson_z,
    V  = pearson_var,
    mods = ~ overlap_epi,
    
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model
    ),
    
    data = tmp,
    method = "REML"
  )
  
  tibble(
    cohort_removed = cc,
    
    intercept_z = coef(fit)[1],
    intercept_z_se = sqrt(vcov(fit)[1, 1]),
    
    overlap_effect_z = coef(fit)[2],
    overlap_effect_z_se = sqrt(vcov(fit)[2, 2]),
    
    overlap_ci.lb_z =
      overlap_effect_z -
      1.96 * overlap_effect_z_se,
    
    overlap_ci.ub_z =
      overlap_effect_z +
      1.96 * overlap_effect_z_se
  )
  
}) %>%
  bind_rows()

loo_overlap_binary_summary <- tibble(
  
  min_overlap_effect_z =
    min(loo_overlap_binary$overlap_effect_z, na.rm = TRUE),
  
  max_overlap_effect_z =
    max(loo_overlap_binary$overlap_effect_z, na.rm = TRUE),
  
  range_overlap_effect_z =
    max_overlap_effect_z -
    min_overlap_effect_z,
  
  sd_overlap_effect_z =
    sd(loo_overlap_binary$overlap_effect_z, na.rm = TRUE)
)

write.csv(
  loo_overlap_binary,
  file.path(
    results_dir_epi_overlap_binary,
    "leave_one_cohort_out_ageoverlap_binary.csv"
  ),
  row.names = FALSE
)

write.csv(
  loo_overlap_binary_summary,
  file.path(
    results_dir_epi_overlap_binary,
    "leave_one_cohort_out_ageoverlap_binary_summary.csv"
  ),
  row.names = FALSE
)

# --- LOO plot ---

png(
  filename = file.path(
    results_dir_epi_overlap_binary,
    "plot_ageoverlap_binary_loo.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_overlap_binary,
  aes(
    x = reorder(cohort_removed, overlap_effect_z),
    y = overlap_effect_z
  )
) +
  
  geom_point(size = 3) +
  
  geom_errorbar(
    aes(
      ymin = overlap_ci.lb_z,
      ymax = overlap_ci.ub_z
    ),
    width = 0.15
  ) +
  
  geom_hline(
    yintercept = coef(epi_pearson_mv_overlap_binary)[2],
    linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = 0,
    colour = "grey60"
  ) +
  
  coord_flip() +
  
  labs(
    x = "Cohort removed",
    y = "Overlap effect on Fisher-z correlation",
    title = "Leave-one-cohort-out sensitivity of age-overlap effect"
  ) +
  
  theme_minimal(base_size = 14)

dev.off()

# --- full analysis log ---

sink(log_file)

cat("Epi Pearson multilevel meta-analysis - age overlap binary moderator\n")
cat("===================================================================\n\n")

cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_epi_overlap_binary, "\n\n")

cat("Model summary, Fisher-z scale\n")
cat("=============================\n")
print(summary(epi_pearson_mv_overlap_binary))

cat("\n\nModel results\n")
cat("=============\n")
print(overlap_binary_est, n = Inf, width = Inf)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")

cat("Total I2 (%):\n")
print(I2_total)

cat("\nI2 components (% of total observed variance):\n")
print(I2_components)

cat("\nHeterogeneity partition (% of random-effect heterogeneity only):\n")
print(heterogeneity_partition)

cat("\n\nPredicted values, back-transformed to Pearson r\n")
cat("===============================================\n")
print(newdat_overlap_binary)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_overlap_binary, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_overlap_binary_summary, n = Inf, width = Inf)

cat("\n\nSession info\n")
cat("============\n")
print(sessionInfo())

sink()

###### epi Pearson overlap_perc_epi global (excluding gestational clocks) ------

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi_overlap_perc <- file.path(
  results_dir,
  paste0("perf_epi_pearson_noGest_global_ageoverlap_perc_", run_id)
)

dir.create(
  results_dir_epi_overlap_perc,
  recursive = TRUE,
  showWarnings = FALSE
)

log_file <- file.path(
  results_dir_epi_overlap_perc,
  "analysis_log_ageoverlap_perc.txt"
)

# --- meta ---

epi_pearson_mv_overlap_perc <- rma.mv(
  yi = pearson_z,
  V  = pearson_var,
  mods = ~ overlap_perc_epi,
  
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  
  data = perf_epi_nogest,
  method = "REML"
)

summary(epi_pearson_mv_overlap_perc)

# --- I2 calculation ---

W <- diag(1 / epi_pearson_mv_overlap_perc$vi)
X <- model.matrix(epi_pearson_mv_overlap_perc)

P <- W - W %*% X %*%
  solve(t(X) %*% W %*% X) %*%
  t(X) %*% W

I2_total <- 100 * sum(epi_pearson_mv_overlap_perc$sigma2) / (
  sum(epi_pearson_mv_overlap_perc$sigma2) +
    (epi_pearson_mv_overlap_perc$k - epi_pearson_mv_overlap_perc$p) /
    sum(diag(P))
)

I2_components <- 100 * epi_pearson_mv_overlap_perc$sigma2 / (
  sum(epi_pearson_mv_overlap_perc$sigma2) +
    (epi_pearson_mv_overlap_perc$k - epi_pearson_mv_overlap_perc$p) /
    sum(diag(P))
)

names(I2_components) <- c(
  "cohort",
  "timepoint_within_cohort",
  "model"
)

heterogeneity_partition <- 100 * epi_pearson_mv_overlap_perc$sigma2 /
  sum(epi_pearson_mv_overlap_perc$sigma2)

names(heterogeneity_partition) <- c(
  "cohort",
  "timepoint_within_cohort",
  "model"
)

# --- model results ---

overlap_perc_est <- tibble(
  term = names(coef(epi_pearson_mv_overlap_perc)),
  estimate_z = as.numeric(coef(epi_pearson_mv_overlap_perc)),
  se_z = epi_pearson_mv_overlap_perc$se,
  ci.lb_z = epi_pearson_mv_overlap_perc$ci.lb,
  ci.ub_z = epi_pearson_mv_overlap_perc$ci.ub,
  pval = epi_pearson_mv_overlap_perc$pval
)

write.csv(
  overlap_perc_est,
  file.path(
    results_dir_epi_overlap_perc,
    "ageoverlap_perc_results.csv"
  ),
  row.names = FALSE
)

# --- predicted values ---

newdat_overlap_perc <- data.frame(
  overlap_perc_epi = seq(
    min(perf_epi_nogest$overlap_perc_epi, na.rm = TRUE),
    max(perf_epi_nogest$overlap_perc_epi, na.rm = TRUE),
    length.out = 100
  )
)

pred_overlap_perc <- predict(
  epi_pearson_mv_overlap_perc,
  newmods = newdat_overlap_perc$overlap_perc_epi
)

newdat_overlap_perc <- newdat_overlap_perc %>%
  mutate(
    pred_z = pred_overlap_perc$pred,
    z_ci.lb = pred_overlap_perc$ci.lb,
    z_ci.ub = pred_overlap_perc$ci.ub,
    
    pred_r = transf.ztor(pred_z),
    ci.lb = transf.ztor(z_ci.lb),
    ci.ub = transf.ztor(z_ci.ub)
  )

write.csv(
  newdat_overlap_perc,
  file.path(
    results_dir_epi_overlap_perc,
    "ageoverlap_perc_predicted_values.csv"
  ),
  row.names = FALSE
)

# --- plot ---

png(
  filename = file.path(
    results_dir_epi_overlap_perc,
    "plot_ageoverlap_perc_pearson_meta.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot() +
  geom_point(
    data = perf_epi_nogest,
    aes(
      x = overlap_perc_epi,
      y = Pearson
    ),
    alpha = 0.25
  ) +
  geom_line(
    data = newdat_overlap_perc,
    aes(
      x = overlap_perc_epi,
      y = pred_r
    ),
    linewidth = 1
  ) +
  geom_ribbon(
    data = newdat_overlap_perc,
    aes(
      x = overlap_perc_epi,
      ymin = ci.lb,
      ymax = ci.ub
    ),
    alpha = 0.15
  ) +
  scale_y_continuous(
    limits = c(-1, 1)
  ) +
  labs(
    x = "Epigenetic clock training age overlap (%)",
    y = "Predicted Pearson r",
    title = "Meta-regression of epigenetic age correlation by proportional training age overlap"
  ) +
  theme_minimal(base_size = 14)

dev.off

# --- save outputs ---

saveRDS(
  epi_pearson_mv_overlap_perc,
  file.path(
    results_dir_epi_overlap_perc,
    "model_epi_pearson_mv_ageoverlap_perc.rds"
  )
)

capture.output(
  summary(epi_pearson_mv_overlap_perc),
  file = file.path(
    results_dir_epi_overlap_perc,
    "model_summary.txt"
  )
)

# --- LOO cohort ---

cohorts <- unique(perf_epi_nogest$cohort)

loo_overlap_perc <- lapply(cohorts, function(cc) {
  
  tmp <- perf_epi_nogest %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = pearson_z,
    V  = pearson_var,
    mods = ~ overlap_perc_epi,
    
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model
    ),
    
    data = tmp,
    method = "REML"
  )
  
  tibble(
    cohort_removed = cc,
    
    intercept_z = coef(fit)[1],
    intercept_z_se = sqrt(vcov(fit)[1, 1]),
    
    overlap_perc_effect_z = coef(fit)[2],
    overlap_perc_effect_z_se = sqrt(vcov(fit)[2, 2]),
    
    overlap_perc_ci.lb_z =
      overlap_perc_effect_z -
      1.96 * overlap_perc_effect_z_se,
    
    overlap_perc_ci.ub_z =
      overlap_perc_effect_z +
      1.96 * overlap_perc_effect_z_se
  )
  
}) %>%
  bind_rows()

loo_overlap_perc_summary <- tibble(
  
  min_overlap_perc_effect_z =
    min(loo_overlap_perc$overlap_perc_effect_z, na.rm = TRUE),
  
  max_overlap_perc_effect_z =
    max(loo_overlap_perc$overlap_perc_effect_z, na.rm = TRUE),
  
  range_overlap_perc_effect_z =
    max_overlap_perc_effect_z -
    min_overlap_perc_effect_z,
  
  sd_overlap_perc_effect_z =
    sd(loo_overlap_perc$overlap_perc_effect_z, na.rm = TRUE)
)

write.csv(
  loo_overlap_perc,
  file.path(
    results_dir_epi_overlap_perc,
    "leave_one_cohort_out_ageoverlap_perc.csv"
  ),
  row.names = FALSE
)

write.csv(
  loo_overlap_perc_summary,
  file.path(
    results_dir_epi_overlap_perc,
    "leave_one_cohort_out_ageoverlap_perc_summary.csv"
  ),
  row.names = FALSE
)

# --- LOO plot ---

png(
  filename = file.path(
    results_dir_epi_overlap_perc,
    "plot_ageoverlap_perc_loo.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_overlap_perc,
  aes(
    x = reorder(cohort_removed, overlap_perc_effect_z),
    y = overlap_perc_effect_z
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = overlap_perc_ci.lb_z,
      ymax = overlap_perc_ci.ub_z
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept = coef(epi_pearson_mv_overlap_perc)[2],
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 0,
    colour = "grey60"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "Overlap percentage effect on Fisher-z correlation",
    title = "Leave-one-cohort-out sensitivity of proportional age-overlap effect"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- full analysis log ---

sink(log_file)

cat("Epi Pearson multilevel meta-analysis - age overlap percentage moderator\n")
cat("=======================================================================\n\n")

cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_epi_overlap_perc, "\n\n")

cat("Model summary, Fisher-z scale\n")
cat("=============================\n")
print(summary(epi_pearson_mv_overlap_perc))

cat("\n\nModel results\n")
cat("=============\n")
print(overlap_perc_est, n = Inf, width = Inf)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")

cat("Total I2 (%):\n")
print(I2_total)

cat("\nI2 components (% of total observed variance):\n")
print(I2_components)

cat("\nHeterogeneity partition (% of random-effect heterogeneity only):\n")
print(heterogeneity_partition)

cat("\n\nPredicted values, back-transformed to Pearson r\n")
cat("===============================================\n")
print(newdat_overlap_perc)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_overlap_perc, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_overlap_perc_summary, n = Inf, width = Inf)

cat("\n\nSession info\n")
cat("============\n")
print(sessionInfo())

sink()

###### epi Pearson generation_epi global ------

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi_generation_epi <- file.path(
  results_dir,
  paste0("perf_epi_pearson_global_generation_epi_", run_id)
)

dir.create(results_dir_epi_generation_epi, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_epi_generation_epi, "analysis_log_generation_epi.txt")

# --- data prep ---

perf_epi_gen <- perf_epi %>%
  mutate(
    generation_epi = factor(generation_epi),
    cohort = factor(cohort),
    timepoint = factor(timepoint),
    model = factor(model)
  ) %>%
  droplevels()

# --- meta ---

epi_pearson_mv_generation_epi <- rma.mv(
  yi = pearson_z,
  V  = pearson_var,
  mods = ~ generation_epi - 1,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = perf_epi_gen,
  method = "REML"
)

summary(epi_pearson_mv_generation_epi)

# --- model results ---

generation_epi_est <- tibble(
  term = names(coef(epi_pearson_mv_generation_epi)),
  generation_epi = str_remove(term, "^generation_epi"),
  estimate_z = as.numeric(coef(epi_pearson_mv_generation_epi)),
  se_z = epi_pearson_mv_generation_epi$se,
  ci.lb_z = epi_pearson_mv_generation_epi$ci.lb,
  ci.ub_z = epi_pearson_mv_generation_epi$ci.ub,
  pval = epi_pearson_mv_generation_epi$pval
) %>%
  mutate(
    estimate_r = transf.ztor(estimate_z),
    ci.lb = transf.ztor(ci.lb_z),
    ci.ub = transf.ztor(ci.ub_z)
  )

write.csv(
  generation_epi_est,
  file.path(results_dir_epi_generation_epi, "generation_epi_results.csv"),
  row.names = FALSE
)

# --- predicted values ---

newdat_generation_epi <- data.frame(
  generation_epi = levels(perf_epi_gen$generation_epi)
) %>%
  mutate(
    generation_epi = factor(
      generation_epi,
      levels = levels(perf_epi_gen$generation_epi)
    )
  )

X_new_generation_epi <- model.matrix(
  ~ generation_epi - 1,
  data = newdat_generation_epi
)

pred_generation_epi <- predict(
  epi_pearson_mv_generation_epi,
  newmods = X_new_generation_epi
)

newdat_generation_epi <- newdat_generation_epi %>%
  mutate(
    pred_z = pred_generation_epi$pred,
    z_ci.lb = pred_generation_epi$ci.lb,
    z_ci.ub = pred_generation_epi$ci.ub,
    pred_r = transf.ztor(pred_z),
    ci.lb = transf.ztor(z_ci.lb),
    ci.ub = transf.ztor(z_ci.ub)
  )

write.csv(
  newdat_generation_epi,
  file.path(results_dir_epi_generation_epi, "generation_epi_predicted_values.csv"),
  row.names = FALSE
)

# --- plot ---

png(
  filename = file.path(
    results_dir_epi_generation_epi,
    "plot_generation_epi_pearson_meta.png"
  ),
  width = 1600,
  height = 1200,
  res = 300
)

ggplot() +
  geom_jitter(
    data = perf_epi_gen,
    aes(x = generation_epi, y = Pearson),
    width = 0.12,
    alpha = 0.25
  ) +
  geom_point(
    data = newdat_generation_epi,
    aes(x = generation_epi, y = pred_r),
    size = 4
  ) +
  geom_errorbar(
    data = newdat_generation_epi,
    aes(x = generation_epi, ymin = ci.lb, ymax = ci.ub),
    width = 0.15
  ) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(
    x = "Epi clock generation",
    y = "Predicted Pearson r",
    title = "Meta-regression of epigenetic age correlation by clock generation"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- save outputs ---

saveRDS(
  epi_pearson_mv_generation_epi,
  file.path(results_dir_epi_generation_epi, "model_epi_pearson_mv_generation_epi.rds")
)

capture.output(
  summary(epi_pearson_mv_generation_epi),
  file = file.path(results_dir_epi_generation_epi, "model_summary.txt")
)

# --- LOO cohort ---

cohorts <- unique(perf_epi_gen$cohort)

loo_generation_epi <- lapply(cohorts, function(cc) {
  
  tmp <- perf_epi_gen %>%
    filter(cohort != cc) %>%
    droplevels()
  
  fit <- rma.mv(
    yi = pearson_z,
    V  = pearson_var,
    mods = ~ generation_epi - 1,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model
    ),
    data = tmp,
    method = "REML"
  )
  
  tibble(
    cohort_removed = cc,
    term = names(coef(fit)),
    estimate_z = as.numeric(coef(fit)),
    se_z = fit$se,
    ci.lb_z = fit$ci.lb,
    ci.ub_z = fit$ci.ub
  ) %>%
    mutate(
      generation_epi = str_remove(term, "^generation_epi"),
      estimate_r = transf.ztor(estimate_z),
      ci.lb = transf.ztor(ci.lb_z),
      ci.ub = transf.ztor(ci.ub_z)
    )
}) %>%
  bind_rows()

loo_generation_epi_summary <- loo_generation_epi %>%
  group_by(generation_epi) %>%
  summarise(
    min_r = min(estimate_r, na.rm = TRUE),
    max_r = max(estimate_r, na.rm = TRUE),
    range_r = max_r - min_r,
    sd_r = sd(estimate_r, na.rm = TRUE),
    min_z = min(estimate_z, na.rm = TRUE),
    max_z = max(estimate_z, na.rm = TRUE),
    range_z = max_z - min_z,
    sd_z = sd(estimate_z, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  loo_generation_epi,
  file.path(results_dir_epi_generation_epi, "leave_one_cohort_out_generation_epi.csv"),
  row.names = FALSE
)

write.csv(
  loo_generation_epi_summary,
  file.path(results_dir_epi_generation_epi, "leave_one_cohort_out_generation_epi_summary.csv"),
  row.names = FALSE
)

# --- LOO plot ---

png(
  filename = file.path(results_dir_epi_generation_epi, "plot_generation_epi_loo.png"),
  width = 2000,
  height = 1400,
  res = 300
)

ggplot(
  loo_generation_epi,
  aes(
    x = reorder(cohort_removed, estimate_r),
    y = estimate_r
  )
) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = ci.lb, ymax = ci.ub),
    width = 0.15
  ) +
  coord_flip() +
  facet_wrap(~ generation_epi, scales = "free_x") +
  labs(
    x = "Cohort removed",
    y = "Predicted Pearson r",
    title = "Leave-one-cohort-out sensitivity by epi clock generation"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- log ---

sink(log_file)

cat("Epi Pearson multilevel meta-analysis - generation_epi\n")
cat("=====================================================\n\n")
cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_epi_generation_epi, "\n\n")

cat("Model summary, Fisher-z scale\n")
cat("=============================\n")
print(summary(epi_pearson_mv_generation_epi))

cat("\n\nModel results, back-transformed to Pearson r\n")
cat("============================================\n")
print(generation_epi_est, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_generation_epi, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_generation_epi_summary, n = Inf, width = Inf)

cat("\n\nSession info\n")
cat("============\n")
print(sessionInfo())

sink()

###### epi Pearson age generation_epi ------

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi_generation_epi_ageint <- file.path(
  results_dir,
  paste0("perf_epi_pearson_generation_epi_ageint_", run_id)
)

dir.create(results_dir_epi_generation_epi_ageint, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  results_dir_epi_generation_epi_ageint,
  "analysis_log_generation_epi_ageint.txt"
)

# --- data prep ---

perf_epi_gen_nogest <- perf_epi %>%
  filter(generation_epi != "1_gest") %>%
  droplevels() %>%
  mutate(
    generation_epi = factor(generation_epi),
    cohort = factor(cohort),
    timepoint = factor(timepoint),
    model = factor(model),
    mean_age_c = mean_age - mean(mean_age, na.rm = TRUE)
  )

# --- meta ---

epi_pearson_mv_generation_epi_ageint <- rma.mv(
  yi = pearson_z,
  V  = pearson_var,
  mods = ~ generation_epi * mean_age_c,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = perf_epi_gen_nogest,
  method = "REML"
)

summary(epi_pearson_mv_generation_epi_ageint)

# --- model results ---

generation_epi_age_est <- tibble(
  term = names(coef(epi_pearson_mv_generation_epi_ageint)),
  estimate_z = as.numeric(coef(epi_pearson_mv_generation_epi_ageint)),
  se_z = epi_pearson_mv_generation_epi_ageint$se,
  ci.lb_z = epi_pearson_mv_generation_epi_ageint$ci.lb,
  ci.ub_z = epi_pearson_mv_generation_epi_ageint$ci.ub,
  pval = epi_pearson_mv_generation_epi_ageint$pval
) %>%
  mutate(
    estimate_r = transf.ztor(estimate_z),
    ci.lb = transf.ztor(ci.lb_z),
    ci.ub = transf.ztor(ci.ub_z),
    type = case_when(
      term == "intrcpt" | term == "(Intercept)" ~ "reference_intercept",
      term == "mean_age_c" ~ "reference_age_slope",
      str_detect(term, ":mean_age_c") ~ "age_slope_difference",
      TRUE ~ "intercept_difference"
    )
  )

write.csv(
  generation_epi_age_est,
  file.path(results_dir_epi_generation_epi_ageint, "generation_epi_age_interaction_results.csv"),
  row.names = FALSE
)

saveRDS(
  epi_pearson_mv_generation_epi_ageint,
  file.path(results_dir_epi_generation_epi_ageint, "model_epi_pearson_mv_generation_epi_ageint.rds")
)

capture.output(
  summary(epi_pearson_mv_generation_epi_ageint),
  file = file.path(results_dir_epi_generation_epi_ageint, "model_summary.txt")
)

# --- predicted values ---

mean_age_mean <- mean(perf_epi_gen_nogest$mean_age, na.rm = TRUE)

newdat <- expand.grid(
  mean_age_c = seq(
    min(perf_epi_gen_nogest$mean_age_c, na.rm = TRUE),
    max(perf_epi_gen_nogest$mean_age_c, na.rm = TRUE),
    length.out = 100
  ),
  generation_epi = levels(perf_epi_gen_nogest$generation_epi)
)

X_new <- model.matrix(
  ~ generation_epi * mean_age_c,
  data = newdat
)

pred <- predict(
  epi_pearson_mv_generation_epi_ageint,
  newmods = X_new[, -1, drop = FALSE]
)

newdat <- newdat %>%
  mutate(
    mean_age = mean_age_c + mean_age_mean,
    pred_z = pred$pred,
    z_ci.lb = pred$ci.lb,
    z_ci.ub = pred$ci.ub,
    pred_r = transf.ztor(pred_z),
    ci.lb = transf.ztor(z_ci.lb),
    ci.ub = transf.ztor(z_ci.ub)
  )

write.csv(
  newdat,
  file.path(results_dir_epi_generation_epi_ageint, "generation_epi_age_interaction_predicted_values.csv"),
  row.names = FALSE
)

# --- plot ---

png(
  filename = file.path(
    results_dir_epi_generation_epi_ageint,
    "plot_epi_meanage_pearson_generation_epi.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(newdat, aes(mean_age, pred_r, colour = generation_epi)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = ci.lb, ymax = ci.ub, fill = generation_epi),
    alpha = .12,
    colour = NA
  ) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(
    x = "Mean sample age",
    y = "Predicted Pearson r",
    title = "Age effects on epi clock correlation by generation"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- LOO ---

cohorts <- unique(perf_epi_gen_nogest$cohort)

loo_generation_epi_age <- lapply(cohorts, function(cc) {
  
  tmp <- perf_epi_gen_nogest %>%
    filter(cohort != cc) %>%
    droplevels()
  
  fit <- rma.mv(
    yi = pearson_z,
    V  = pearson_var,
    mods = ~ generation_epi * mean_age_c,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model
    ),
    data = tmp,
    method = "REML"
  )
  
  tibble(
    cohort_removed = cc,
    term = names(coef(fit)),
    estimate_z = as.numeric(coef(fit)),
    se_z = fit$se,
    ci.lb_z = fit$ci.lb,
    ci.ub_z = fit$ci.ub,
    pval = fit$pval
  ) %>%
    mutate(
      estimate_r = transf.ztor(estimate_z),
      ci.lb = transf.ztor(ci.lb_z),
      ci.ub = transf.ztor(ci.ub_z),
      type = case_when(
        term == "intrcpt" | term == "(Intercept)" ~ "reference_intercept",
        term == "mean_age_c" ~ "reference_age_slope",
        str_detect(term, ":mean_age_c") ~ "age_slope_difference",
        TRUE ~ "intercept_difference"
      )
    )
}) %>%
  bind_rows()

write.csv(
  loo_generation_epi_age,
  file.path(results_dir_epi_generation_epi_ageint, "generation_epi_leave_one_cohort_out_age_interactions.csv"),
  row.names = FALSE
)

loo_generation_epi_age_summary <- loo_generation_epi_age %>%
  filter(type %in% c("reference_age_slope", "age_slope_difference")) %>%
  group_by(term, type) %>%
  summarise(
    min_z = min(estimate_z, na.rm = TRUE),
    max_z = max(estimate_z, na.rm = TRUE),
    range_z = max_z - min_z,
    sd_z = sd(estimate_z, na.rm = TRUE),
    min_r = min(estimate_r, na.rm = TRUE),
    max_r = max(estimate_r, na.rm = TRUE),
    range_r = max_r - min_r,
    sd_r = sd(estimate_r, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  loo_generation_epi_age_summary,
  file.path(results_dir_epi_generation_epi_ageint, "generation_epi_leave_one_cohort_out_age_interactions_summary.csv"),
  row.names = FALSE
)

# --- log ---

sink(log_file)

cat("Epi Pearson multilevel meta-analysis - generation_epi age interactions\n")
cat("=====================================================================\n\n")
cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_epi_generation_epi_ageint, "\n\n")

cat("Model summary, Fisher-z scale\n")
cat("=============================\n")
print(summary(epi_pearson_mv_generation_epi_ageint))

cat("\n\nModel results\n")
cat("=============\n")
print(generation_epi_age_est, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_generation_epi_age, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary for age-related terms\n")
cat("==================================================\n")
print(loo_generation_epi_age_summary, n = Inf, width = Inf)

cat("\n\nSession info\n")
cat("============\n")
print(sessionInfo())

sink()


###### epi wMAE tissue  ------

# create tissue overlap variable 
perf_epi <- perf_epi %>%
  mutate(
    tissue_overlap = case_when(
      
      # multi-tissue training clocks count as overlap for all
      tissue_training == "multi" ~ "yes",
      
      # exact tissue match
      tissue_training == tissue_test ~ "yes",
      
      # otherwise no overlap
      TRUE ~ "no"
    ),
    
    tissue_overlap = factor(
      tissue_overlap,
      levels = c("no", "yes")
    )
  )

table(
  perf_epi$tissue_training,
  perf_epi$tissue_test,
  perf_epi$tissue_overlap
)

###### epi Pearson tissue_overlap global ------

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi_tissue_overlap <- file.path(
  results_dir,
  paste0("perf_epi_pearson_global_tissue_overlap_", run_id)
)

dir.create(
  results_dir_epi_tissue_overlap,
  recursive = TRUE,
  showWarnings = FALSE
)

log_file <- file.path(
  results_dir_epi_tissue_overlap,
  "analysis_log_tissue_overlap.txt"
)

# --- data prep ---

perf_epi <- perf_epi %>%
  mutate(
    tissue_overlap = factor(
      tissue_overlap,
      levels = c("no", "yes")
    ),
    
    cohort = factor(cohort),
    timepoint = factor(timepoint),
    model = factor(model)
  )

# --- meta ---

epi_pearson_mv_tissue_overlap <- rma.mv(
  yi = pearson_z,
  V  = pearson_var,
  mods = ~ tissue_overlap,
  
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  
  data = perf_epi,
  method = "REML"
)

summary(epi_pearson_mv_tissue_overlap)

# --- I2 calculation ---

W <- diag(1 / epi_pearson_mv_tissue_overlap$vi)

X <- model.matrix(epi_pearson_mv_tissue_overlap)

P <- W - W %*% X %*%
  solve(t(X) %*% W %*% X) %*%
  t(X) %*% W

I2_total <- 100 * sum(epi_pearson_mv_tissue_overlap$sigma2) / (
  sum(epi_pearson_mv_tissue_overlap$sigma2) +
    (epi_pearson_mv_tissue_overlap$k -
       epi_pearson_mv_tissue_overlap$p) /
    sum(diag(P))
)

I2_components <- 100 *
  epi_pearson_mv_tissue_overlap$sigma2 / (
    sum(epi_pearson_mv_tissue_overlap$sigma2) +
      (epi_pearson_mv_tissue_overlap$k -
         epi_pearson_mv_tissue_overlap$p) /
      sum(diag(P))
  )

names(I2_components) <- c(
  "cohort",
  "timepoint_within_cohort",
  "model"
)

heterogeneity_partition <- 100 *
  epi_pearson_mv_tissue_overlap$sigma2 /
  sum(epi_pearson_mv_tissue_overlap$sigma2)

names(heterogeneity_partition) <- c(
  "cohort",
  "timepoint_within_cohort",
  "model"
)

# --- model results ---

tissue_overlap_est <- tibble(
  term = names(coef(epi_pearson_mv_tissue_overlap)),
  estimate_z = as.numeric(coef(epi_pearson_mv_tissue_overlap)),
  se_z = epi_pearson_mv_tissue_overlap$se,
  ci.lb_z = epi_pearson_mv_tissue_overlap$ci.lb,
  ci.ub_z = epi_pearson_mv_tissue_overlap$ci.ub,
  pval = epi_pearson_mv_tissue_overlap$pval
)

write.csv(
  tissue_overlap_est,
  file.path(
    results_dir_epi_tissue_overlap,
    "tissue_overlap_results.csv"
  ),
  row.names = FALSE
)

# --- predicted values ---

newdat_tissue_overlap <- data.frame(
  tissue_overlap = levels(perf_epi$tissue_overlap)
) %>%
  mutate(
    tissue_overlap = factor(
      tissue_overlap,
      levels = levels(perf_epi$tissue_overlap)
    )
  )

X_new_tissue_overlap <- model.matrix(
  ~ tissue_overlap,
  data = newdat_tissue_overlap
)

pred_tissue_overlap <- predict(
  epi_pearson_mv_tissue_overlap,
  newmods = X_new_tissue_overlap[, -1, drop = FALSE]
)

newdat_tissue_overlap <- newdat_tissue_overlap %>%
  mutate(
    pred_z = pred_tissue_overlap$pred,
    z_ci.lb = pred_tissue_overlap$ci.lb,
    z_ci.ub = pred_tissue_overlap$ci.ub,
    
    pred_r = transf.ztor(pred_z),
    ci.lb = transf.ztor(z_ci.lb),
    ci.ub = transf.ztor(z_ci.ub)
  )

write.csv(
  newdat_tissue_overlap,
  file.path(
    results_dir_epi_tissue_overlap,
    "tissue_overlap_predicted_values.csv"
  ),
  row.names = FALSE
)

# --- plot ---

png(
  filename = file.path(
    results_dir_epi_tissue_overlap,
    "plot_tissue_overlap_pearson_meta.png"
  ),
  width = 1600,
  height = 1200,
  res = 300
)

ggplot() +
  
  geom_jitter(
    data = perf_epi,
    aes(
      x = tissue_overlap,
      y = Pearson
    ),
    width = 0.12,
    alpha = 0.25
  ) +
  
  geom_point(
    data = newdat_tissue_overlap,
    aes(
      x = tissue_overlap,
      y = pred_r
    ),
    size = 4
  ) +
  
  geom_errorbar(
    data = newdat_tissue_overlap,
    aes(
      x = tissue_overlap,
      ymin = ci.lb,
      ymax = ci.ub
    ),
    width = 0.15
  ) +
  
  scale_y_continuous(
    limits = c(-1, 1)
  ) +
  
  labs(
    x = "Training-test tissue overlap",
    y = "Predicted Pearson r",
    title = "Meta-regression of epigenetic age correlation by tissue overlap"
  ) +
  
  theme_minimal(base_size = 14)

dev.off()

# --- save outputs ---

saveRDS(
  epi_pearson_mv_tissue_overlap,
  file.path(
    results_dir_epi_tissue_overlap,
    "model_epi_pearson_mv_tissue_overlap.rds"
  )
)

capture.output(
  summary(epi_pearson_mv_tissue_overlap),
  file = file.path(
    results_dir_epi_tissue_overlap,
    "model_summary.txt"
  )
)

# --- LOO cohort ---

cohorts <- unique(perf_epi$cohort)

loo_tissue_overlap <- lapply(cohorts, function(cc) {
  
  tmp <- perf_epi %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = pearson_z,
    V  = pearson_var,
    mods = ~ tissue_overlap,
    
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model
    ),
    
    data = tmp,
    method = "REML"
  )
  
  tibble(
    cohort_removed = cc,
    
    intercept_z = coef(fit)[1],
    intercept_z_se = sqrt(vcov(fit)[1, 1]),
    
    tissue_overlap_effect_z = coef(fit)[2],
    tissue_overlap_effect_z_se = sqrt(vcov(fit)[2, 2]),
    
    tissue_overlap_ci.lb_z =
      tissue_overlap_effect_z -
      1.96 * tissue_overlap_effect_z_se,
    
    tissue_overlap_ci.ub_z =
      tissue_overlap_effect_z +
      1.96 * tissue_overlap_effect_z_se
  )
  
}) %>%
  bind_rows()

loo_tissue_overlap_summary <- tibble(
  
  min_tissue_overlap_effect_z =
    min(loo_tissue_overlap$tissue_overlap_effect_z, na.rm = TRUE),
  
  max_tissue_overlap_effect_z =
    max(loo_tissue_overlap$tissue_overlap_effect_z, na.rm = TRUE),
  
  range_tissue_overlap_effect_z =
    max_tissue_overlap_effect_z -
    min_tissue_overlap_effect_z,
  
  sd_tissue_overlap_effect_z =
    sd(loo_tissue_overlap$tissue_overlap_effect_z, na.rm = TRUE)
)

write.csv(
  loo_tissue_overlap,
  file.path(
    results_dir_epi_tissue_overlap,
    "leave_one_cohort_out_tissue_overlap.csv"
  ),
  row.names = FALSE
)

write.csv(
  loo_tissue_overlap_summary,
  file.path(
    results_dir_epi_tissue_overlap,
    "leave_one_cohort_out_tissue_overlap_summary.csv"
  ),
  row.names = FALSE
)

# --- LOO plot ---

png(
  filename = file.path(
    results_dir_epi_tissue_overlap,
    "plot_tissue_overlap_loo.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_tissue_overlap,
  aes(
    x = reorder(
      cohort_removed,
      tissue_overlap_effect_z
    ),
    y = tissue_overlap_effect_z
  )
) +
  
  geom_point(size = 3) +
  
  geom_errorbar(
    aes(
      ymin = tissue_overlap_ci.lb_z,
      ymax = tissue_overlap_ci.ub_z
    ),
    width = 0.15
  ) +
  
  geom_hline(
    yintercept = coef(epi_pearson_mv_tissue_overlap)[2],
    linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = 0,
    colour = "grey60"
  ) +
  
  coord_flip() +
  
  labs(
    x = "Cohort removed",
    y = "Tissue-overlap effect on Fisher-z correlation",
    title = "Leave-one-cohort-out sensitivity of tissue-overlap effect"
  ) +
  
  theme_minimal(base_size = 14)

dev.off()

# --- full analysis log ---

sink(log_file)

cat("Epi Pearson multilevel meta-analysis - tissue overlap moderator\n")
cat("================================================================\n\n")

cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_epi_tissue_overlap, "\n\n")

cat("Model summary, Fisher-z scale\n")
cat("=============================\n")
print(summary(epi_pearson_mv_tissue_overlap))

cat("\n\nModel results\n")
cat("=============\n")
print(tissue_overlap_est, n = Inf, width = Inf)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")

cat("Total I2 (%):\n")
print(I2_total)

cat("\nI2 components (% of total observed variance):\n")
print(I2_components)

cat("\nHeterogeneity partition (% of random-effect heterogeneity only):\n")
print(heterogeneity_partition)

cat("\n\nPredicted values, back-transformed to Pearson r\n")
cat("===============================================\n")
print(newdat_tissue_overlap)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_tissue_overlap, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_tissue_overlap_summary, n = Inf, width = Inf)

cat("\n\nSession info\n")
cat("============\n")
print(sessionInfo())

sink()

###### epi Pearson tissue_training global ------

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi_tissue_training <- file.path(
  results_dir,
  paste0("perf_epi_pearson_global_tissue_training_", run_id)
)

dir.create(results_dir_epi_tissue_training, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_epi_tissue_training, "analysis_log_tissue_training.txt")

# data prep 
perf_epi <- perf_epi %>%
  mutate(
    tissue_training = factor(tissue_training),
    cohort = factor(cohort),
    timepoint = factor(timepoint),
    model = factor(model)
  )

# meta 
epi_pearson_mv_tissue_training <- rma.mv(
  yi = pearson_z,
  V  = pearson_var,
  mods = ~ tissue_training - 1,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = perf_epi,
  method = "REML"
)

summary(epi_pearson_mv_tissue_training)

# --- I2 calculation ---

W <- diag(1 / epi_pearson_mv_tissue_training$vi)
X <- model.matrix(epi_pearson_mv_tissue_training)

P <- W - W %*% X %*%
  solve(t(X) %*% W %*% X) %*%
  t(X) %*% W

I2_total <- 100 * sum(epi_pearson_mv_tissue_training$sigma2) / (
  sum(epi_pearson_mv_tissue_training$sigma2) +
    (epi_pearson_mv_tissue_training$k - epi_pearson_mv_tissue_training$p) /
    sum(diag(P))
)

I2_components <- 100 * epi_pearson_mv_tissue_training$sigma2 / (
  sum(epi_pearson_mv_tissue_training$sigma2) +
    (epi_pearson_mv_tissue_training$k - epi_pearson_mv_tissue_training$p) /
    sum(diag(P))
)

names(I2_components) <- c("cohort", "timepoint_within_cohort", "model")

heterogeneity_partition <- 100 * epi_pearson_mv_tissue_training$sigma2 /
  sum(epi_pearson_mv_tissue_training$sigma2)

names(heterogeneity_partition) <- c("cohort", "timepoint_within_cohort", "model")

# --- model results ---

tissue_training_est <- tibble(
  term = names(coef(epi_pearson_mv_tissue_training)),
  tissue_training = str_remove(term, "^tissue_training"),
  estimate_z = as.numeric(coef(epi_pearson_mv_tissue_training)),
  se_z = epi_pearson_mv_tissue_training$se,
  ci.lb_z = epi_pearson_mv_tissue_training$ci.lb,
  ci.ub_z = epi_pearson_mv_tissue_training$ci.ub,
  pval = epi_pearson_mv_tissue_training$pval
) %>%
  mutate(
    pearson_r = transf.ztor(estimate_z),
    ci.lb = transf.ztor(ci.lb_z),
    ci.ub = transf.ztor(ci.ub_z)
  )

write.csv(
  tissue_training_est,
  file.path(results_dir_epi_tissue_training, "tissue_training_results.csv"),
  row.names = FALSE
)

# --- predicted values ---

newdat_tissue_training <- data.frame(
  tissue_training = levels(perf_epi$tissue_training)
) %>%
  mutate(
    tissue_training = factor(
      tissue_training,
      levels = levels(perf_epi$tissue_training)
    )
  )

X_new_tissue_training <- model.matrix(
  ~ tissue_training - 1,
  data = newdat_tissue_training
)

pred_tissue_training <- predict(
  epi_pearson_mv_tissue_training,
  newmods = X_new_tissue_training
)

newdat_tissue_training <- newdat_tissue_training %>%
  mutate(
    pred_z = pred_tissue_training$pred,
    z_ci.lb = pred_tissue_training$ci.lb,
    z_ci.ub = pred_tissue_training$ci.ub,
    pred_r = transf.ztor(pred_z),
    ci.lb = transf.ztor(z_ci.lb),
    ci.ub = transf.ztor(z_ci.ub)
  )

write.csv(
  newdat_tissue_training,
  file.path(results_dir_epi_tissue_training, "tissue_training_predicted_values.csv"),
  row.names = FALSE
)

# --- plot ---

png(
  filename = file.path(results_dir_epi_tissue_training, "plot_tissue_training_pearson_meta.png"),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot() +
  geom_jitter(
    data = perf_epi,
    aes(x = tissue_training, y = Pearson),
    width = 0.12,
    alpha = 0.25
  ) +
  geom_point(
    data = newdat_tissue_training,
    aes(x = tissue_training, y = pred_r),
    size = 4
  ) +
  geom_errorbar(
    data = newdat_tissue_training,
    aes(x = tissue_training, ymin = ci.lb, ymax = ci.ub),
    width = 0.15
  ) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(
    x = "Training tissue",
    y = "Predicted Pearson r",
    title = "Meta-regression of epigenetic age correlation by training tissue"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- save outputs ---

saveRDS(
  epi_pearson_mv_tissue_training,
  file.path(results_dir_epi_tissue_training, "model_epi_pearson_mv_tissue_training.rds")
)

capture.output(
  summary(epi_pearson_mv_tissue_training),
  file = file.path(results_dir_epi_tissue_training, "model_summary.txt")
)

# --- LOO cohort ---

cohorts <- unique(perf_epi$cohort)

loo_tissue_training <- lapply(cohorts, function(cc) {
  
  tmp <- perf_epi %>%
    filter(cohort != cc) %>%
    droplevels()
  
  fit <- rma.mv(
    yi = pearson_z,
    V  = pearson_var,
    mods = ~ tissue_training - 1,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model
    ),
    data = tmp,
    method = "REML"
  )
  
  tibble(
    cohort_removed = cc,
    term = names(coef(fit)),
    estimate_z = as.numeric(coef(fit)),
    se_z = fit$se,
    ci.lb_z = fit$ci.lb,
    ci.ub_z = fit$ci.ub
  ) %>%
    mutate(
      tissue_training = str_remove(term, "^tissue_training"),
      pearson_r = transf.ztor(estimate_z),
      ci.lb = transf.ztor(ci.lb_z),
      ci.ub = transf.ztor(ci.ub_z)
    )
}) %>%
  bind_rows()

loo_tissue_training_summary <- loo_tissue_training %>%
  group_by(tissue_training) %>%
  summarise(
    min_r = min(pearson_r, na.rm = TRUE),
    max_r = max(pearson_r, na.rm = TRUE),
    range_r = max_r - min_r,
    sd_r = sd(pearson_r, na.rm = TRUE),
    
    min_z = min(estimate_z, na.rm = TRUE),
    max_z = max(estimate_z, na.rm = TRUE),
    range_z = max_z - min_z,
    sd_z = sd(estimate_z, na.rm = TRUE),
    
    .groups = "drop"
  )

write.csv(
  loo_tissue_training,
  file.path(results_dir_epi_tissue_training, "leave_one_cohort_out_tissue_training.csv"),
  row.names = FALSE
)

write.csv(
  loo_tissue_training_summary,
  file.path(results_dir_epi_tissue_training, "leave_one_cohort_out_tissue_training_summary.csv"),
  row.names = FALSE
)

# --- LOO plot ---

png(
  filename = file.path(results_dir_epi_tissue_training, "plot_tissue_training_loo.png"),
  width = 2200,
  height = 1400,
  res = 300
)

ggplot(
  loo_tissue_training,
  aes(
    x = reorder(cohort_removed, pearson_r),
    y = pearson_r
  )
) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = ci.lb, ymax = ci.ub),
    width = 0.15
  ) +
  coord_flip() +
  facet_wrap(~ tissue_training, scales = "free_x") +
  labs(
    x = "Cohort removed",
    y = "Predicted Pearson r",
    title = "Leave-one-cohort-out sensitivity by training tissue"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- full analysis log ---

sink(log_file)

cat("Epi Pearson multilevel meta-analysis - tissue training\n")
cat("======================================================\n\n")

cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_epi_tissue_training, "\n\n")

cat("Model summary, Fisher-z scale\n")
cat("=============================\n")
print(summary(epi_pearson_mv_tissue_training))

cat("\n\nModel results, back-transformed to Pearson r\n")
cat("============================================\n")
print(tissue_training_est, n = Inf, width = Inf)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
cat("Total I2 (%):\n")
print(I2_total)

cat("\nI2 components (% of total observed variance):\n")
print(I2_components)

cat("\nHeterogeneity partition (% of random-effect heterogeneity only):\n")
print(heterogeneity_partition)

cat("\n\nPredicted values, back-transformed to Pearson r\n")
cat("===============================================\n")
print(newdat_tissue_training)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_tissue_training, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_tissue_training_summary, n = Inf, width = Inf)

cat("\n\nSession info\n")
cat("============\n")
print(sessionInfo())

sink()

###### epi Pearson tissue_test global ------

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi_tissue_test <- file.path(
  results_dir,
  paste0("perf_epi_pearson_global_tissue_test_", run_id)
)

dir.create(results_dir_epi_tissue_test, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_epi_tissue_test, "analysis_log_tissue_test.txt")

# data prep 

perf_epi <- perf_epi %>%
  mutate(
    tissue_test = factor(tissue_test),
    cohort = factor(cohort),
    timepoint = factor(timepoint),
    model = factor(model)
  )

# meta 
epi_pearson_mv_tissue_test <- rma.mv(
  yi = pearson_z,
  V  = pearson_var,
  mods = ~ tissue_test - 1,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = perf_epi,
  method = "REML"
)

summary(epi_pearson_mv_tissue_test)

# --- I2 calculation ---

W <- diag(1 / epi_pearson_mv_tissue_test$vi)
X <- model.matrix(epi_pearson_mv_tissue_test)

P <- W - W %*% X %*%
  solve(t(X) %*% W %*% X) %*%
  t(X) %*% W

I2_total <- 100 * sum(epi_pearson_mv_tissue_test$sigma2) / (
  sum(epi_pearson_mv_tissue_test$sigma2) +
    (epi_pearson_mv_tissue_test$k - epi_pearson_mv_tissue_test$p) /
    sum(diag(P))
)

I2_components <- 100 * epi_pearson_mv_tissue_test$sigma2 / (
  sum(epi_pearson_mv_tissue_test$sigma2) +
    (epi_pearson_mv_tissue_test$k - epi_pearson_mv_tissue_test$p) /
    sum(diag(P))
)

names(I2_components) <- c("cohort", "timepoint_within_cohort", "model")

heterogeneity_partition <- 100 * epi_pearson_mv_tissue_test$sigma2 /
  sum(epi_pearson_mv_tissue_test$sigma2)

names(heterogeneity_partition) <- c("cohort", "timepoint_within_cohort", "model")

# --- model results ---

tissue_test_est <- tibble(
  term = names(coef(epi_pearson_mv_tissue_test)),
  tissue_test = str_remove(term, "^tissue_test"),
  estimate_z = as.numeric(coef(epi_pearson_mv_tissue_test)),
  se_z = epi_pearson_mv_tissue_test$se,
  ci.lb_z = epi_pearson_mv_tissue_test$ci.lb,
  ci.ub_z = epi_pearson_mv_tissue_test$ci.ub,
  pval = epi_pearson_mv_tissue_test$pval
) %>%
  mutate(
    pearson_r = transf.ztor(estimate_z),
    ci.lb = transf.ztor(ci.lb_z),
    ci.ub = transf.ztor(ci.ub_z)
  )

write.csv(
  tissue_test_est,
  file.path(results_dir_epi_tissue_test, "tissue_test_results.csv"),
  row.names = FALSE
)

# --- predicted values ---

newdat_tissue_test <- data.frame(
  tissue_test = levels(perf_epi$tissue_test)
) %>%
  mutate(
    tissue_test = factor(
      tissue_test,
      levels = levels(perf_epi$tissue_test)
    )
  )

X_new_tissue_test <- model.matrix(
  ~ tissue_test - 1,
  data = newdat_tissue_test
)

pred_tissue_test <- predict(
  epi_pearson_mv_tissue_test,
  newmods = X_new_tissue_test
)

newdat_tissue_test <- newdat_tissue_test %>%
  mutate(
    pred_z = pred_tissue_test$pred,
    z_ci.lb = pred_tissue_test$ci.lb,
    z_ci.ub = pred_tissue_test$ci.ub,
    pred_r = transf.ztor(pred_z),
    ci.lb = transf.ztor(z_ci.lb),
    ci.ub = transf.ztor(z_ci.ub)
  )

write.csv(
  newdat_tissue_test,
  file.path(results_dir_epi_tissue_test, "tissue_test_predicted_values.csv"),
  row.names = FALSE
)

# --- plot ---

png(
  filename = file.path(results_dir_epi_tissue_test, "plot_tissue_test_pearson_meta.png"),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot() +
  geom_jitter(
    data = perf_epi,
    aes(x = tissue_test, y = Pearson),
    width = 0.12,
    alpha = 0.25
  ) +
  geom_point(
    data = newdat_tissue_test,
    aes(x = tissue_test, y = pred_r),
    size = 4
  ) +
  geom_errorbar(
    data = newdat_tissue_test,
    aes(x = tissue_test, ymin = ci.lb, ymax = ci.ub),
    width = 0.15
  ) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(
    x = "test tissue",
    y = "Predicted Pearson r",
    title = "Meta-regression of epigenetic age correlation by test tissue"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- save outputs ---

saveRDS(
  epi_pearson_mv_tissue_test,
  file.path(results_dir_epi_tissue_test, "model_epi_pearson_mv_tissue_test.rds")
)

capture.output(
  summary(epi_pearson_mv_tissue_test),
  file = file.path(results_dir_epi_tissue_test, "model_summary.txt")
)

# --- LOO cohort ---

cohorts <- unique(perf_epi$cohort)

loo_tissue_test <- lapply(cohorts, function(cc) {
  
  tmp <- perf_epi %>%
    filter(cohort != cc) %>%
    droplevels()
  
  fit <- rma.mv(
    yi = pearson_z,
    V  = pearson_var,
    mods = ~ tissue_test - 1,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model
    ),
    data = tmp,
    method = "REML"
  )
  
  tibble(
    cohort_removed = cc,
    term = names(coef(fit)),
    estimate_z = as.numeric(coef(fit)),
    se_z = fit$se,
    ci.lb_z = fit$ci.lb,
    ci.ub_z = fit$ci.ub
  ) %>%
    mutate(
      tissue_test = str_remove(term, "^tissue_test"),
      pearson_r = transf.ztor(estimate_z),
      ci.lb = transf.ztor(ci.lb_z),
      ci.ub = transf.ztor(ci.ub_z)
    )
}) %>%
  bind_rows()

loo_tissue_test_summary <- loo_tissue_test %>%
  group_by(tissue_test) %>%
  summarise(
    min_r = min(pearson_r, na.rm = TRUE),
    max_r = max(pearson_r, na.rm = TRUE),
    range_r = max_r - min_r,
    sd_r = sd(pearson_r, na.rm = TRUE),
    
    min_z = min(estimate_z, na.rm = TRUE),
    max_z = max(estimate_z, na.rm = TRUE),
    range_z = max_z - min_z,
    sd_z = sd(estimate_z, na.rm = TRUE),
    
    .groups = "drop"
  )

write.csv(
  loo_tissue_test,
  file.path(results_dir_epi_tissue_test, "leave_one_cohort_out_tissue_test.csv"),
  row.names = FALSE
)

write.csv(
  loo_tissue_test_summary,
  file.path(results_dir_epi_tissue_test, "leave_one_cohort_out_tissue_test_summary.csv"),
  row.names = FALSE
)

# --- LOO plot ---

png(
  filename = file.path(results_dir_epi_tissue_test, "plot_tissue_test_loo.png"),
  width = 2200,
  height = 1400,
  res = 300
)

ggplot(
  loo_tissue_test,
  aes(
    x = reorder(cohort_removed, pearson_r),
    y = pearson_r
  )
) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = ci.lb, ymax = ci.ub),
    width = 0.15
  ) +
  coord_flip() +
  facet_wrap(~ tissue_test, scales = "free_x") +
  labs(
    x = "Cohort removed",
    y = "Predicted Pearson r",
    title = "Leave-one-cohort-out sensitivity by test tissue"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- full analysis log ---

sink(log_file)

cat("Epi Pearson multilevel meta-analysis - tissue test\n")
cat("======================================================\n\n")

cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_epi_tissue_test, "\n\n")

cat("Model summary, Fisher-z scale\n")
cat("=============================\n")
print(summary(epi_pearson_mv_tissue_test))

cat("\n\nModel results, back-transformed to Pearson r\n")
cat("============================================\n")
print(tissue_test_est, n = Inf, width = Inf)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
cat("Total I2 (%):\n")
print(I2_total)

cat("\nI2 components (% of total observed variance):\n")
print(I2_components)

cat("\nHeterogeneity partition (% of random-effect heterogeneity only):\n")
print(heterogeneity_partition)

cat("\n\nPredicted values, back-transformed to Pearson r\n")
cat("===============================================\n")
print(newdat_tissue_test)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_tissue_test, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_tissue_test_summary, n = Inf, width = Inf)

cat("\n\nSession info\n")
cat("============\n")
print(sessionInfo())

sink()


##### brain pearson (all clocks) -----

###### brain pearson global ------

# --- prep ---
# exclude specific models if needed 
#perf_brain   <- perf_brain   %>% filter(!model %in% c("DunedinPACE", "DNAmTL")) %>% droplevels()

# --- set up folder ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_brain_r <- file.path(results_dir, paste0("perf_brain_pearson_global_", run_id))

dir.create(results_dir_brain_r, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_brain_r, "analysis_log.txt")

# --- descriptives ---

write.csv(perf_brain, file.path(results_dir_brain_r, "perf_brain_descr_detailed.csv"), row.names=FALSE)

perf_brain_desc <- perf_brain %>%
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

write.csv(perf_brain_desc, file.path(results_dir_brain_r, "model_perf_brain_desc.csv"), row.names = FALSE)

# --- meta ---
brain_pearson_mv_global <- rma.mv(
  yi = pearson_z,
  V  = pearson_var,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data   = perf_brain,
  method = "REML"
)

summary(brain_pearson_mv_global)

# heterogeneity estimates

W <- diag(1/brain_pearson_mv_global$vi)
X <- model.matrix(brain_pearson_mv_global)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_total <- 100 * sum(brain_pearson_mv_global$sigma2) / (sum(brain_pearson_mv_global$sigma2) + (brain_pearson_mv_global$k-brain_pearson_mv_global$p)/sum(diag(P)))


I2_components <- 100 * brain_pearson_mv_global$sigma2 / (
  sum(brain_pearson_mv_global$sigma2) +
    (brain_pearson_mv_global$k - brain_pearson_mv_global$p) /
    sum(diag(P))
)

names(I2_components) <- c(
  "cohort",
  "timepoint_within_cohort",
  "model"
)

I2_total
I2_components

# Back-transform global pooled estimate
overall_est <- data.frame(
  model = "Pooled",
  pearson_z = brain_pearson_mv_global$beta[1, 1],
  pearson_r = transf.ztor(brain_pearson_mv_global$beta[1, 1]),
  ci.lb = transf.ztor(brain_pearson_mv_global$ci.lb[1]),
  ci.ub = transf.ztor(brain_pearson_mv_global$ci.ub[1]),
  pval = brain_pearson_mv_global$pval[1]
  
)

write.csv(
  overall_est,
  file.path(results_dir_brain_r, "overall_pearson_estimate.csv"),
  row.names = FALSE
)

saveRDS(
  brain_pearson_mv_global,
  file.path(results_dir_brain_r, "model_brain_pearson_mv_global.rds")
)

capture.output(
  summary(brain_pearson_mv_global),
  file = file.path(results_dir_brain_r, "model_summary.txt")
)

# --- LOO cohort ---

cohorts <- unique(perf_brain$cohort)

loo_cohort <- lapply(cohorts, function(cc) {
  
  tmp <- subset(perf_brain, cohort != cc)
  
  fit <- rma.mv(
    yi = pearson_z,
    V  = pearson_var,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model
    ),
    data = tmp,
    method = "REML"
  )
  
  z_est <- coef(fit)[1]
  z_se  <- sqrt(vcov(fit)[1, 1])
  
  data.frame(
    cohort_removed = cc,
    
    pearson_z = z_est,
    pearson_z_se = z_se,
    z_ci.lb = z_est - 1.96 * z_se,
    z_ci.ub = z_est + 1.96 * z_se,
    
    pearson_r = transf.ztor(z_est),
    ci.lb = transf.ztor(z_est - 1.96 * z_se),
    ci.ub = transf.ztor(z_est + 1.96 * z_se)
  )
}) %>%
  bind_rows()

loo_summary <- tibble(
  min_r = min(loo_cohort$pearson_r, na.rm = TRUE),
  max_r = max(loo_cohort$pearson_r, na.rm = TRUE),
  range_r = max_r - min_r,
  sd_r = sd(loo_cohort$pearson_r, na.rm = TRUE),
  
  min_z = min(loo_cohort$pearson_z, na.rm = TRUE),
  max_z = max(loo_cohort$pearson_z, na.rm = TRUE),
  range_z = max_z - min_z,
  sd_z = sd(loo_cohort$pearson_z, na.rm = TRUE)
)

write.csv(
  loo_cohort,
  file.path(results_dir_brain_r, "leave_one_cohort_out.csv"),
  row.names = FALSE
)

write.csv(
  loo_summary,
  file.path(results_dir_brain_r, "leave_one_cohort_out_summary.csv"),
  row.names = FALSE
)

# --- global LOO plot: Pearson r scale ---

full_r <- transf.ztor(coef(brain_pearson_mv_global)[1])

png(
  filename = file.path(results_dir_brain_r, "plot_global_pearson_loo.png"),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_cohort,
  aes(
    x = reorder(cohort_removed, pearson_r),
    y = pearson_r
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = ci.lb, ymax = ci.ub),
    width = 0.15
  ) +
  geom_hline(
    yintercept = full_r,
    linetype = "dashed"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "Pooled Pearson r",
    title = "Leave-one-cohort-out sensitivity of pooled Pearson correlation"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# full analysis log

sink(log_file)

cat("brain Pearson multilevel meta-analysis\n")
cat("====================================\n\n")
cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_brain_r, "\n\n")

cat("Model summary, Fisher-z scale\n")
cat("=============================\n")
print(summary(brain_pearson_mv_global))

cat("\n\nModel summary - back-transformed to Pearson r\n")
cat("=============================================\n")
print(overall_est)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")

cat("Total I2 (%):\n")
print(I2_total)

cat("\nVariance decomposition (% of total heterogeneity):\n")
print(I2_components)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_cohort)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_summary)

cat("\n\nDescriptives\n")
cat("============\n")
print(perf_brain_desc, n = Inf, width = Inf)

sink()

###### brain Pearson modelwise ------

# --- set up folder ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_brain_pearson_modelwise <- file.path(
  results_dir,
  paste0("perf_brain_pearson_modelwise_", run_id)
)

dir.create(results_dir_brain_pearson_modelwise, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_brain_pearson_modelwise, "analysis_modelwise_log.txt")

# --- meta ---

brain_pearson_mv_modelwise <- rma.mv(
  yi = pearson_z,
  V  = pearson_var,
  mods = ~ model - 1,
  random = list(
    ~ 1 | cohort/timepoint
  ),
  data = perf_brain,
  method = "REML"
)

summary(brain_pearson_mv_modelwise)


# heterogeneity estimates

W <- diag(1/brain_pearson_mv_modelwise$vi)
X <- model.matrix(brain_pearson_mv_modelwise)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_total <- 100 * sum(brain_pearson_mv_modelwise$sigma2) / (sum(brain_pearson_mv_modelwise$sigma2) + (brain_pearson_mv_modelwise$k-brain_pearson_mv_modelwise$p)/sum(diag(P)))


I2_components <- 100 * brain_pearson_mv_modelwise$sigma2 / (
  sum(brain_pearson_mv_modelwise$sigma2) +
    (brain_pearson_mv_modelwise$k - brain_pearson_mv_modelwise$p) /
    sum(diag(P))
)

names(I2_components) <- c(
  "cohort",
  "timepoint_within_cohort"
)

I2_total
I2_components

# --- back-transform modelwise pooled estimates ---

brain_pearson_est <- data.frame(
  model = gsub("^model", "", rownames(brain_pearson_mv_modelwise$beta)),
  pearson_z = brain_pearson_mv_modelwise$beta[, 1],
  pearson_r = transf.ztor(brain_pearson_mv_modelwise$beta[, 1]),
  ci.lb = transf.ztor(brain_pearson_mv_modelwise$ci.lb),
  ci.ub = transf.ztor(brain_pearson_mv_modelwise$ci.ub),
  pval = brain_pearson_mv_modelwise$pval
)

X <- model.matrix(brain_pearson_mv_modelwise)

model_cols <- colnames(X)[
  str_detect(colnames(X), "^model") &
    !str_detect(colnames(X), ":")
]

effect_sizes_per_model <- tibble(
  model = str_remove(model_cols, "^model"),
  k = colSums(X[, model_cols, drop = FALSE] != 0)
)

brain_pearson_est <- brain_pearson_est %>%
  left_join(effect_sizes_per_model, by = "model")

brain_pearson_est

write.csv(
  brain_pearson_est,
  file.path(results_dir_brain_pearson_modelwise, "modelwise_pearson_estimates.csv"),
  row.names = FALSE
)

saveRDS(
  brain_pearson_mv_modelwise,
  file.path(results_dir_brain_pearson_modelwise, "model_brain_pearson_mv_modelwise.rds")
)

capture.output(
  summary(brain_pearson_mv_modelwise),
  file = file.path(results_dir_brain_pearson_modelwise, "model_summary.txt")
)

# --- LOO cohort for Fisher-z Pearson model ---

cohorts <- unique(perf_brain$cohort)

loo_cohort_modelwise <- lapply(cohorts, function(cc) {
  
  tmp <- subset(perf_brain, cohort != cc)
  
  fit <- rma.mv(
    yi = pearson_z,
    V  = pearson_var,
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
    
    # Fisher-z scale estimates
    pearson_z = as.numeric(fit$b),
    pearson_z_se = fit$se,
    z_ci.lb = fit$ci.lb,
    z_ci.ub = fit$ci.ub,
    pval = fit$pval
  ) %>%
    mutate(
      model = str_remove(term, "^model"),
      
      # back-transformed estimates: Pearson r
      pearson_r = transf.ztor(pearson_z),
      ci.lb = transf.ztor(z_ci.lb),
      ci.ub = transf.ztor(z_ci.ub)
    ) %>%
    select(
      cohort_removed,
      model,
      pearson_z,
      pearson_z_se,
      z_ci.lb,
      z_ci.ub,
      pearson_r,
      ci.lb,
      ci.ub,
      pval
    )
}) %>%
  bind_rows()

write.csv(
  loo_cohort_modelwise,
  file.path(results_dir_brain_pearson_modelwise, "modelwise_leave_one_cohort_out.csv"),
  row.names = FALSE
)

modelwise_loo_summary <- loo_cohort_modelwise %>%
  group_by(model) %>%
  summarise(
    min_r = min(pearson_r, na.rm = TRUE),
    max_r = max(pearson_r, na.rm = TRUE),
    range_r = max_r - min_r,
    sd_r = sd(pearson_r, na.rm = TRUE),
    
    min_z = min(pearson_z, na.rm = TRUE),
    max_z = max(pearson_z, na.rm = TRUE),
    range_z = max_z - min_z,
    sd_z = sd(pearson_z, na.rm = TRUE),
    
    .groups = "drop"
  )

write.csv(
  modelwise_loo_summary,
  file.path(results_dir_brain_pearson_modelwise, "modelwise_leave_one_cohort_out_summary.csv"),
  row.names = FALSE
)

# --- optional LOO plot by model ---

png(
  filename = file.path(results_dir_brain_pearson_modelwise, "plot_modelwise_pearson_loo.png"),
  width = 2200,
  height = 2200,
  res = 300
)

ggplot(
  loo_cohort_modelwise,
  aes(
    x = reorder(cohort_removed, pearson_r),
    y = pearson_r
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
    y = "Pooled Pearson r",
    title = "Modelwise leave-one-cohort-out sensitivity of pooled Pearson correlations"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- full analysis log ---

sink(log_file)

cat("brain Pearson multilevel meta-analysis modelwise\n")
cat("==============================================\n\n")
cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_brain_pearson_modelwise, "\n\n")

cat("Model summary, Fisher-z scale\n")
cat("=============================\n")
print(summary(brain_pearson_mv_modelwise))

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")

cat("Total I2 (%):\n")
print(I2_total)

cat("\nVariance decomposition (% of total heterogeneity):\n")
print(I2_components)

cat("\n\nModel summary - back-transformed to Pearson r\n")
cat("=============================================\n")
print(brain_pearson_est)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_cohort_modelwise, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(modelwise_loo_summary, n = Inf, width = Inf)

sink()

###### brain Pearson age global ------

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_brain_pearson_age <- file.path(
  results_dir,
  paste0("perf_brain_pearson_global_ageint_", run_id)
)

dir.create(results_dir_brain_pearson_age, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_brain_pearson_age, "analysis_log_global_ageint.txt")

# --- data prep ---

perf_brain <- perf_brain %>%
  mutate(mean_age_c = mean_age - mean(mean_age, na.rm = TRUE))

# remove gestational clocks for age moderation
perf_brain_nogest <- perf_brain %>%
  filter(!model %in% c("Bohlin", "Knight", "brainC")) %>%
  droplevels()

# --- meta-regression ---

brain_pearson_mv_global_ageint <- rma.mv(
  yi = pearson_z,
  V  = pearson_var,
  mods = ~ mean_age_c,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = perf_brain_nogest,
  method = "REML"
)

summary(brain_pearson_mv_global_ageint)

# --- I2 calculation ---

W <- diag(1 / brain_pearson_mv_global_ageint$vi)
X <- model.matrix(brain_pearson_mv_global_ageint)

P <- W - W %*% X %*%
  solve(t(X) %*% W %*% X) %*%
  t(X) %*% W

I2_total_age <- 100 * sum(brain_pearson_mv_global_ageint$sigma2) / (
  sum(brain_pearson_mv_global_ageint$sigma2) +
    (brain_pearson_mv_global_ageint$k - brain_pearson_mv_global_ageint$p) /
    sum(diag(P))
)

I2_components_age <- 100 * brain_pearson_mv_global_ageint$sigma2 / (
  sum(brain_pearson_mv_global_ageint$sigma2) +
    (brain_pearson_mv_global_ageint$k - brain_pearson_mv_global_ageint$p) /
    sum(diag(P))
)

names(I2_components_age) <- c(
  "cohort",
  "timepoint_within_cohort",
  "model"
)

heterogeneity_partition_age <- 100 * brain_pearson_mv_global_ageint$sigma2 /
  sum(brain_pearson_mv_global_ageint$sigma2)

names(heterogeneity_partition_age) <- c(
  "cohort",
  "timepoint_within_cohort",
  "model"
)

# --- extract model results ---

age_model_est <- tibble(
  term = names(coef(brain_pearson_mv_global_ageint)),
  estimate_z = as.numeric(coef(brain_pearson_mv_global_ageint)),
  se_z = brain_pearson_mv_global_ageint$se,
  ci.lb_z = brain_pearson_mv_global_ageint$ci.lb,
  ci.ub_z = brain_pearson_mv_global_ageint$ci.ub,
  pval = brain_pearson_mv_global_ageint$pval
) %>%
  mutate(
    estimate_r = if_else(term == "intrcpt", transf.ztor(estimate_z), NA_real_),
    ci.lb_r = if_else(term == "intrcpt", transf.ztor(ci.lb_z), NA_real_),
    ci.ub_r = if_else(term == "intrcpt", transf.ztor(ci.ub_z), NA_real_)
  )

write.csv(
  age_model_est,
  file.path(results_dir_brain_pearson_age, "age_moderator_results.csv"),
  row.names = FALSE
)

saveRDS(
  brain_pearson_mv_global_ageint,
  file.path(results_dir_brain_pearson_age, "model_brain_pearson_mv_global_ageint.rds")
)

# --- prediction grid ---

mean_age_mean <- mean(perf_brain_nogest$mean_age, na.rm = TRUE)

newdat <- data.frame(
  mean_age = seq(
    min(perf_brain_nogest$mean_age, na.rm = TRUE),
    max(perf_brain_nogest$mean_age, na.rm = TRUE),
    length.out = 100
  )
) %>%
  mutate(mean_age_c = mean_age - mean_age_mean)

pred <- predict(
  brain_pearson_mv_global_ageint,
  newmods = newdat$mean_age_c
)

newdat <- newdat %>%
  mutate(
    pred_z = pred$pred,
    z_ci.lb = pred$ci.lb,
    z_ci.ub = pred$ci.ub,
    pred_r = transf.ztor(pred_z),
    ci.lb = transf.ztor(z_ci.lb),
    ci.ub = transf.ztor(z_ci.ub)
  )

png(
  filename = file.path(results_dir_brain_pearson_age, "plot_meanage_pearson_meta_regression.png"),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot() +
  geom_point(
    data = perf_brain_nogest,
    aes(x = mean_age, y = Pearson, colour = cohort),
    alpha = .25,
    size = 2
  ) +
  geom_line(
    data = newdat,
    aes(x = mean_age, y = pred_r),
    linewidth = 1.2,
    colour = "black"
  ) +
  geom_ribbon(
    data = newdat,
    aes(x = mean_age, ymin = ci.lb, ymax = ci.ub),
    alpha = .15
  ) +
  labs(
    title = "Meta-regression of Pearson correlation on mean age",
    x = "Mean sample age",
    y = "Predicted Pearson r",
    colour = "Cohort"
  ) +
  theme_minimal(base_size = 16)

dev.off()

# --- LOO age moderation ---

cohorts <- unique(perf_brain_nogest$cohort)

loo_age <- lapply(cohorts, function(cc) {
  
  tmp <- perf_brain_nogest %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = pearson_z,
    V  = pearson_var,
    mods = ~ mean_age_c,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model
    ),
    data = tmp,
    method = "REML"
  )
  
  z_intercept <- coef(fit)[1]
  z_slope <- coef(fit)[2]
  
  z_intercept_se <- sqrt(vcov(fit)[1, 1])
  z_slope_se <- sqrt(vcov(fit)[2, 2])
  
  data.frame(
    cohort_removed = cc,
    
    intercept_z = z_intercept,
    slope_z_per_year = z_slope,
    intercept_z_se = z_intercept_se,
    slope_z_se = z_slope_se,
    
    intercept_r = transf.ztor(z_intercept),
    intercept_ci.lb = transf.ztor(z_intercept - 1.96 * z_intercept_se),
    intercept_ci.ub = transf.ztor(z_intercept + 1.96 * z_intercept_se),
    
    slope_z_ci.lb = z_slope - 1.96 * z_slope_se,
    slope_z_ci.ub = z_slope + 1.96 * z_slope_se
  )
}) %>%
  bind_rows()

loo_age_summary <- tibble(
  min_intercept_r = min(loo_age$intercept_r, na.rm = TRUE),
  max_intercept_r = max(loo_age$intercept_r, na.rm = TRUE),
  range_intercept_r = max_intercept_r - min_intercept_r,
  sd_intercept_r = sd(loo_age$intercept_r, na.rm = TRUE),
  
  min_slope_z = min(loo_age$slope_z_per_year, na.rm = TRUE),
  max_slope_z = max(loo_age$slope_z_per_year, na.rm = TRUE),
  range_slope_z = max_slope_z - min_slope_z,
  sd_slope_z = sd(loo_age$slope_z_per_year, na.rm = TRUE)
)

write.csv(
  loo_age,
  file.path(results_dir_brain_pearson_age, "leave_one_cohort_out_age_model.csv"),
  row.names = FALSE
)

write.csv(
  loo_age_summary,
  file.path(results_dir_brain_pearson_age, "leave_one_cohort_out_age_model_summary.csv"),
  row.names = FALSE
)

full_slope_z <- coef(brain_pearson_mv_global_ageint)[2]

png(
  filename = file.path(results_dir_brain_pearson_age, "plot_meanage_pearson_loo_slope.png"),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_age,
  aes(
    x = reorder(cohort_removed, slope_z_per_year),
    y = slope_z_per_year
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = slope_z_ci.lb, ymax = slope_z_ci.ub),
    width = 0.15
  ) +
  geom_hline(yintercept = full_slope_z, linetype = "dashed") +
  geom_hline(yintercept = 0, colour = "grey60") +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "Change in Fisher-z Pearson correlation per year",
    title = "Leave-one-cohort-out sensitivity of age effect"
  ) +
  theme_minimal(base_size = 14)

dev.off()

sink(log_file)

cat("brain Pearson multilevel meta-analysis - global age moderator\n")
cat("============================================================\n\n")
cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_brain_pearson_age, "\n\n")

cat("Model summary, Fisher-z scale\n")
cat("=============================\n")
print(summary(brain_pearson_mv_global_ageint))

cat("\n\nModel results\n")
cat("=============\n")
print(age_model_est, n = Inf, width = Inf)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")

cat("Total I2 (%):\n")
print(I2_total_age)

cat("\nI2 components (% of total observed variance):\n")
print(I2_components_age)
cat("\n=each random effect component as % of heterogeneity+sampling variance\n, components sum to I^2 but not to 100%\n")

cat("\nHeterogeneity partition (% of random-effect heterogeneity only):\n")
print(heterogeneity_partition_age)
cat("\n=each random effect component as % of heterogeneity only\n, components sum to 100%\n")

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_age)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_age_summary, n = Inf, width = Inf)

sink()

###### brain Pearson age modelwise ------

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_brain_pearson_modelwise_age <- file.path(
  results_dir,
  paste0("perf_brain_pearson_modelwise_ageint_", run_id)
)

dir.create(results_dir_brain_pearson_modelwise_age, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  results_dir_brain_pearson_modelwise_age,
  "analysis_log_modelwise_ageint.txt"
)

# --- data prep ---

perf_brain_nogest <- perf_brain %>%
  mutate(mean_age_c = mean_age - mean(mean_age, na.rm = TRUE)) %>%
  filter(!model %in% c("Bohlin", "Knight", "brainC")) %>%
  droplevels()

# --- meta: model-specific intercepts and model-specific age slopes ---

brain_pearson_mv_modelwise_ageint <- rma.mv(
  yi = pearson_z,
  V  = pearson_var,
  mods = ~ model + model:mean_age_c - 1,
  random = list(
    ~ 1 | cohort/timepoint
  ),
  data = perf_brain_nogest,
  method = "REML"
)

summary(brain_pearson_mv_modelwise_ageint)

# compare common age slope model vs model-specific age slopes
anova(
  brain_pearson_mv_global_ageint,
  brain_pearson_mv_modelwise_ageint,
  refit = TRUE
)


# --- I2 calculation ---

W <- diag(1 / brain_pearson_mv_modelwise_ageint$vi)
X <- model.matrix(brain_pearson_mv_modelwise_ageint)

P <- W - W %*% X %*%
  solve(t(X) %*% W %*% X) %*%
  t(X) %*% W

I2_total_modelwise_age <- 100 * sum(brain_pearson_mv_modelwise_ageint$sigma2) / (
  sum(brain_pearson_mv_modelwise_ageint$sigma2) +
    (brain_pearson_mv_modelwise_ageint$k - brain_pearson_mv_modelwise_ageint$p) /
    sum(diag(P))
)

I2_components_modelwise_age <- 100 * brain_pearson_mv_modelwise_ageint$sigma2 / (
  sum(brain_pearson_mv_modelwise_ageint$sigma2) +
    (brain_pearson_mv_modelwise_ageint$k - brain_pearson_mv_modelwise_ageint$p) /
    sum(diag(P))
)

names(I2_components_modelwise_age) <- c(
  "cohort",
  "timepoint_within_cohort"
)

heterogeneity_partition_modelwise_age <- 100 * brain_pearson_mv_modelwise_ageint$sigma2 /
  sum(brain_pearson_mv_modelwise_ageint$sigma2)

names(heterogeneity_partition_modelwise_age) <- c(
  "cohort",
  "timepoint_within_cohort"
)

# --- extract modelwise results ---

modelwise_age_est <- tibble(
  term = rownames(brain_pearson_mv_modelwise_ageint$b),
  estimate_z = as.numeric(brain_pearson_mv_modelwise_ageint$b),
  se_z = brain_pearson_mv_modelwise_ageint$se,
  ci.lb_z = brain_pearson_mv_modelwise_ageint$ci.lb,
  ci.ub_z = brain_pearson_mv_modelwise_ageint$ci.ub,
  pval = brain_pearson_mv_modelwise_ageint$pval
) %>%
  mutate(
    type = if_else(str_detect(term, "mean_age_c"), "age_slope", "intercept"),
    model = term %>%
      str_remove("^model") %>%
      str_remove(":mean_age_c") %>%
      str_remove("\\.mean_age_c"),
    
    # intercepts can be back-transformed to Pearson r
    estimate_r = if_else(type == "intercept", transf.ztor(estimate_z), NA_real_),
    ci.lb_r = if_else(type == "intercept", transf.ztor(ci.lb_z), NA_real_),
    ci.ub_r = if_else(type == "intercept", transf.ztor(ci.ub_z), NA_real_)
  )

X <- model.matrix(brain_pearson_mv_modelwise_ageint)

model_cols <- colnames(X)[
  str_detect(colnames(X), "^model") &
    !str_detect(colnames(X), ":")
]

effect_sizes_per_model <- tibble(
  model = str_remove(model_cols, "^model"),
  k = colSums(X[, model_cols, drop = FALSE] != 0)
)

modelwise_age_est <- modelwise_age_est %>%
  left_join(effect_sizes_per_model, by = "model")

write.csv(
  modelwise_age_est,
  file.path(
    results_dir_brain_pearson_modelwise_age,
    "modelwise_age_interaction_results.csv"
  ),
  row.names = FALSE
)

saveRDS(
  brain_pearson_mv_modelwise_ageint,
  file.path(
    results_dir_brain_pearson_modelwise_age,
    "model_brain_pearson_mv_modelwise_ageint.rds"
  )
)

# --- plot model-specific age effects (r scale) ---

mean_age_mean <- mean(perf_brain_nogest$mean_age, na.rm = TRUE)

newdat <- expand.grid(
  mean_age_c = seq(
    min(perf_brain_nogest$mean_age_c, na.rm = TRUE),
    max(perf_brain_nogest$mean_age_c, na.rm = TRUE),
    length.out = 100
  ),
  model = unique(perf_brain_nogest$model)
)

pred <- predict(
  brain_pearson_mv_modelwise_ageint,
  newmods = model.matrix(~ model + model:mean_age_c - 1, newdat)
)

newdat <- newdat %>%
  mutate(
    mean_age = mean_age_c + mean_age_mean,
    
    pred_z = pred$pred,
    z_ci.lb = pred$ci.lb,
    z_ci.ub = pred$ci.ub,
    
    pred_r = transf.ztor(pred_z),
    ci.lb = transf.ztor(z_ci.lb),
    ci.ub = transf.ztor(z_ci.ub)
  )

png(
  filename = file.path(
    results_dir_brain_pearson_modelwise_age,
    "plot_brain_meanage_pearson_modelwise.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(newdat, aes(mean_age, pred_r, colour = model)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = ci.lb, ymax = ci.ub, fill = model),
    alpha = .12,
    colour = NA
  ) +
  labs(
    x = "Mean sample age",
    y = "Predicted Pearson r",
    title = "Model-specific age effects on Pearson correlation"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- LOO modelwise age slopes ---
cohorts <- unique(perf_brain_nogest$cohort)

loo_modelwise_age <- lapply(cohorts, function(cc) {
  
  tmp <- perf_brain_nogest %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = pearson_z,
    V  = pearson_var,
    mods = ~ model + model:mean_age_c - 1,
    random = list(
      ~ 1 | cohort/timepoint
    ),
    data = tmp,
    method = "REML"
  )
  
  tibble(
    cohort_removed = cc,
    term = rownames(fit$b),
    estimate_z = as.numeric(fit$b),
    se_z = fit$se,
    ci.lb_z = fit$ci.lb,
    ci.ub_z = fit$ci.ub,
    pval = fit$pval
  ) %>%
    mutate(
      type = if_else(str_detect(term, "mean_age_c"), "age_slope", "intercept"),
      model = term %>%
        str_remove("^model") %>%
        str_remove(":mean_age_c") %>%
        str_remove("\\.mean_age_c"),
      
      estimate_r = if_else(type == "intercept", transf.ztor(estimate_z), NA_real_),
      ci.lb_r = if_else(type == "intercept", transf.ztor(ci.lb_z), NA_real_),
      ci.ub_r = if_else(type == "intercept", transf.ztor(ci.ub_z), NA_real_)
    )
}) %>%
  bind_rows()

write.csv(
  loo_modelwise_age,
  file.path(
    results_dir_brain_pearson_modelwise_age,
    "modelwise_leave_one_cohort_out_age_interactions.csv"
  ),
  row.names = FALSE
)

loo_modelwise_age_summary <- loo_modelwise_age %>%
  filter(type == "age_slope") %>%
  group_by(model) %>%
  summarise(
    min_slope_z = min(estimate_z, na.rm = TRUE),
    max_slope_z = max(estimate_z, na.rm = TRUE),
    range_slope_z = max_slope_z - min_slope_z,
    sd_slope_z = sd(estimate_z, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  loo_modelwise_age_summary,
  file.path(
    results_dir_brain_pearson_modelwise_age,
    "modelwise_leave_one_cohort_out_age_interactions_summary.csv"
  ),
  row.names = FALSE
)

# --- modelwise LOO plot for age slopes ---

# full-model age slopes
full_modelwise_slopes <- modelwise_age_est %>%
  filter(type == "age_slope") %>%
  select(
    model,
    full_slope_z = estimate_z
  )

# LOO age slopes only
plot_dat_slopes <- loo_modelwise_age %>%
  filter(type == "age_slope") %>%
  left_join(full_modelwise_slopes, by = "model") %>%
  mutate(
    model = factor(model),
    cohort_removed = factor(cohort_removed)
  )

full_modelwise_slopes_plot <- plot_dat_slopes %>%
  distinct(model, full_slope_z)

png(
  filename = file.path(
    results_dir_brain_pearson_modelwise_age,
    "plot_modelwise_meanage_loo_slopes.png"
  ),
  width = 3000,
  height = 2200,
  res = 300
)

p <- ggplot(
  plot_dat_slopes,
  aes(
    x = estimate_z,
    y = reorder(cohort_removed, estimate_z)
  )
) +
  # full-model age slope per clock
  geom_vline(
    data = full_modelwise_slopes_plot,
    aes(xintercept = full_slope_z),
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  
  # null effect
  geom_vline(
    xintercept = 0,
    colour = "grey60"
  ) +
  
  # LOO point estimates
  geom_point(size = 1.8) +
  
  # confidence intervals
  geom_segment(
    aes(
      x = ci.lb_z,
      xend = ci.ub_z,
      yend = reorder(cohort_removed, estimate_z)
    ),
    linewidth = 0.5
  ) +
  
  facet_wrap(~ model, scales = "free_x") +
  
  labs(
    x = "Age slope on Fisher-z scale\n(change in correlation per year)",
    y = "Cohort removed",
    title = "Modelwise leave-one-cohort-out sensitivity of mean-age effects"
  ) +
  
  theme_minimal(base_size = 11)

print(p)

dev.off()

# --- log ---
sink(log_file)

cat("brain Pearson multilevel meta-analysis - modelwise age interactions\n")
cat("==================================================================\n\n")
cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_brain_pearson_modelwise_age, "\n\n")

cat("Model summary, Fisher-z scale\n")
cat("=============================\n")
print(summary(brain_pearson_mv_modelwise_ageint))

cat("\n\nModel comparison: common age slope vs model-specific age slopes\n")
cat("===============================================================\n")
print(anova(
  brain_pearson_mv_global_ageint,
  brain_pearson_mv_modelwise_ageint,
  refit = TRUE
))

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")

cat("Total I2 (%):\n")
print(I2_total_modelwise_age)

cat("\nI2 components (% of total observed variance):\n")
print(I2_components_modelwise_age)

cat("\nHeterogeneity partition (% of random-effect heterogeneity only):\n")
print(heterogeneity_partition_modelwise_age)

cat("\n\nModelwise age interaction results\n")
cat("=================================\n")
print(modelwise_age_est, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_modelwise_age, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary for age slopes\n")
cat("===========================================\n")
print(loo_modelwise_age_summary, n = Inf, width = Inf)

sink()

###### brain Pearson overlap_brain global ------

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_brain_ageoverlap_binary <- file.path(
  results_dir,
  paste0("perf_brain_pearson_global_ageoverlap_binary_", run_id)
)

dir.create(results_dir_brain_ageoverlap_binary, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  results_dir_brain_ageoverlap_binary,
  "analysis_log_ageoverlap_binary.txt"
)

# --- data prep ---

perf_brain <- perf_brain %>%
  mutate(
    
    overlap_lower = pmax(min_age, training_age_min_brain),
    overlap_upper = pmin(max_age, training_age_max_brain),
    
    overlap_brain = if_else(
      overlap_upper >= overlap_lower,
      "yes",
      "no"
    ),
    
    # inclusive overlap width
    overlap_width = if_else(
      overlap_brain == "yes",
      overlap_upper - overlap_lower,
      0
    ),
    
    # inclusive study range width
    study_width = max_age - min_age,
    
    overlap_perc_brain =
      100 * overlap_width / study_width
  ) %>%
  select(-overlap_lower, -overlap_upper,
         -overlap_width, -study_width)

perf_brain <- perf_brain %>%
  mutate(
    overlap_brain = factor(overlap_brain)
  )

# --- meta ---

brain_pearson_mv_ageoverlap_binary <- rma.mv(
  yi = pearson_z,
  V  = pearson_var,
  mods = ~ overlap_brain,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = perf_brain,
  method = "REML"
)

summary(brain_pearson_mv_ageoverlap_binary)

# --- model results ---

ageoverlap_binary_est <- tibble(
  term = names(coef(brain_pearson_mv_ageoverlap_binary)),
  estimate_z = as.numeric(coef(brain_pearson_mv_ageoverlap_binary)),
  se_z = brain_pearson_mv_ageoverlap_binary$se,
  ci.lb_z = brain_pearson_mv_ageoverlap_binary$ci.lb,
  ci.ub_z = brain_pearson_mv_ageoverlap_binary$ci.ub,
  pval = brain_pearson_mv_ageoverlap_binary$pval
)

write.csv(
  ageoverlap_binary_est,
  file.path(
    results_dir_brain_ageoverlap_binary,
    "ageoverlap_binary_results.csv"
  ),
  row.names = FALSE
)

# --- I2 calculation ---

W <- diag(1 / brain_pearson_mv_ageoverlap_binary$vi)
X <- model.matrix(brain_pearson_mv_ageoverlap_binary)

P <- W - W %*% X %*%
  solve(t(X) %*% W %*% X) %*%
  t(X) %*% W

I2_total <- 100 * sum(brain_pearson_mv_ageoverlap_binary$sigma2) / (
  sum(brain_pearson_mv_ageoverlap_binary$sigma2) +
    (brain_pearson_mv_ageoverlap_binary$k - brain_pearson_mv_ageoverlap_binary$p) /
    sum(diag(P))
)

I2_components <- 100 * brain_pearson_mv_ageoverlap_binary$sigma2 / (
  sum(brain_pearson_mv_ageoverlap_binary$sigma2) +
    (brain_pearson_mv_ageoverlap_binary$k - brain_pearson_mv_ageoverlap_binary$p) /
    sum(diag(P))
)

names(I2_components) <- c(
  "cohort",
  "timepoint_within_cohort",
  "model"
)

heterogeneity_partition <- 100 * brain_pearson_mv_ageoverlap_binary$sigma2 /
  sum(brain_pearson_mv_ageoverlap_binary$sigma2)

names(heterogeneity_partition) <- c(
  "cohort",
  "timepoint_within_cohort",
  "model"
)

# --- predicted values ---

newdat_ageoverlap_binary <- data.frame(
  overlap_brain = levels(perf_brain$overlap_brain)
)

newdat_ageoverlap_binary <- newdat_ageoverlap_binary %>%
  mutate(
    overlap_brain = factor(
      overlap_brain,
      levels = levels(perf_brain$overlap_brain)
    )
  )

X_new_ageoverlap_binary <- model.matrix(
  ~ overlap_brain,
  data = newdat_ageoverlap_binary
)

pred_ageoverlap_binary <- predict(
  brain_pearson_mv_ageoverlap_binary,
  newmods = X_new_ageoverlap_binary[, -1, drop = FALSE]
)

newdat_ageoverlap_binary <- newdat_ageoverlap_binary %>%
  mutate(
    pred_z = pred_ageoverlap_binary$pred,
    z_ci.lb = pred_ageoverlap_binary$ci.lb,
    z_ci.ub = pred_ageoverlap_binary$ci.ub,
    
    pred_r = transf.ztor(pred_z),
    ci.lb = transf.ztor(z_ci.lb),
    ci.ub = transf.ztor(z_ci.ub)
  )

write.csv(
  newdat_ageoverlap_binary,
  file.path(
    results_dir_brain_ageoverlap_binary,
    "ageoverlap_binary_predicted_values.csv"
  ),
  row.names = FALSE
)

# --- plot predicted Pearson r ---

png(
  filename = file.path(
    results_dir_brain_ageoverlap_binary,
    "plot_ageoverlap_binary_pearson_meta.png"
  ),
  width = 1600,
  height = 1200,
  res = 300
)

ggplot() +
  geom_jitter(
    data = perf_brain,
    aes(
      x = overlap_brain,
      y = Pearson
    ),
    width = 0.12,
    alpha = 0.25
  ) +
  geom_point(
    data = newdat_ageoverlap_binary,
    aes(
      x = overlap_brain,
      y = pred_r
    ),
    size = 4
  ) +
  geom_errorbar(
    data = newdat_ageoverlap_binary,
    aes(
      x = overlap_brain,
      ymin = ci.lb,
      ymax = ci.ub
    ),
    width = 0.15
  ) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(
    x = "Brain age training overlap",
    y = "Predicted Pearson r",
    title = "Meta-regression of brain-age correlation by training age overlap"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- save model outputs ---

saveRDS(
  brain_pearson_mv_ageoverlap_binary,
  file.path(
    results_dir_brain_ageoverlap_binary,
    "model_brain_pearson_mv_ageoverlap_binary.rds"
  )
)

capture.output(
  summary(brain_pearson_mv_ageoverlap_binary),
  file = file.path(
    results_dir_brain_ageoverlap_binary,
    "model_summary.txt"
  )
)

# --- LOO cohort ---

cohorts <- unique(perf_brain$cohort)

loo_ageoverlap_binary <- lapply(cohorts, function(cc) {
  
  tmp <- perf_brain %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = pearson_z,
    V  = pearson_var,
    mods = ~ overlap_brain,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model
    ),
    data = tmp,
    method = "REML"
  )
  
  tibble(
    cohort_removed = cc,
    
    intercept_z = coef(fit)[1],
    intercept_z_se = sqrt(vcov(fit)[1, 1]),
    
    overlap_effect_z = coef(fit)[2],
    overlap_effect_z_se = sqrt(vcov(fit)[2, 2]),
    
    overlap_ci.lb_z = overlap_effect_z - 1.96 * overlap_effect_z_se,
    overlap_ci.ub_z = overlap_effect_z + 1.96 * overlap_effect_z_se
  )
}) %>%
  bind_rows()

loo_ageoverlap_binary_summary <- tibble(
  min_overlap_effect_z = min(loo_ageoverlap_binary$overlap_effect_z, na.rm = TRUE),
  max_overlap_effect_z = max(loo_ageoverlap_binary$overlap_effect_z, na.rm = TRUE),
  range_overlap_effect_z = max_overlap_effect_z - min_overlap_effect_z,
  sd_overlap_effect_z = sd(loo_ageoverlap_binary$overlap_effect_z, na.rm = TRUE)
)

write.csv(
  loo_ageoverlap_binary,
  file.path(
    results_dir_brain_ageoverlap_binary,
    "leave_one_cohort_out_ageoverlap_binary.csv"
  ),
  row.names = FALSE
)

write.csv(
  loo_ageoverlap_binary_summary,
  file.path(
    results_dir_brain_ageoverlap_binary,
    "leave_one_cohort_out_ageoverlap_binary_summary.csv"
  ),
  row.names = FALSE
)

# --- LOO plot: overlap effect on Fisher-z scale ---

png(
  filename = file.path(
    results_dir_brain_ageoverlap_binary,
    "plot_ageoverlap_binary_loo.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_ageoverlap_binary,
  aes(
    x = reorder(cohort_removed, overlap_effect_z),
    y = overlap_effect_z
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = overlap_ci.lb_z,
      ymax = overlap_ci.ub_z
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept = coef(brain_pearson_mv_ageoverlap_binary)[2],
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 0,
    colour = "grey60"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "Overlap effect on Fisher-z correlation",
    title = "Leave-one-cohort-out sensitivity of age-overlap effect"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- full analysis log ---

sink(log_file)

cat("Brain Pearson multilevel meta-analysis - age overlap binary moderator\n")
cat("=====================================================================\n\n")

cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_brain_ageoverlap_binary, "\n\n")

cat("Model summary, Fisher-z scale\n")
cat("=============================\n")
print(summary(brain_pearson_mv_ageoverlap_binary))

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")

cat("Total I2 (%):\n")
print(I2_total)

cat("\nI2 components (% of total observed variance):\n")
print(I2_components)

cat("\nHeterogeneity partition (% of random-effect heterogeneity only):\n")
print(heterogeneity_partition)

cat("\n\nPredicted values, back-transformed to Pearson r\n")
cat("===============================================\n")
print(newdat_ageoverlap_binary)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_ageoverlap_binary, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_ageoverlap_binary_summary, n = Inf, width = Inf)

cat("\n\nSession info\n")
cat("============\n")
print(sessionInfo())

sink()

###### brain Pearson overlap_perc_brain global ------

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_brain_ageoverlap_perc <- file.path(
  results_dir,
  paste0("perf_brain_pearson_global_ageoverlap_perc_", run_id)
)

dir.create(results_dir_brain_ageoverlap_perc, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  results_dir_brain_ageoverlap_perc,
  "analysis_log_ageoverlap_perc.txt"
)

# --- meta ---

brain_pearson_mv_ageoverlap_perc <- rma.mv(
  yi = pearson_z,
  V  = pearson_var,
  mods = ~ overlap_perc_brain,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = perf_brain,
  method = "REML"
)

summary(brain_pearson_mv_ageoverlap_perc)

# --- model results ---

ageoverlap_perc_est <- tibble(
  term = names(coef(brain_pearson_mv_ageoverlap_perc)),
  estimate_z = as.numeric(coef(brain_pearson_mv_ageoverlap_perc)),
  se_z = brain_pearson_mv_ageoverlap_perc$se,
  ci.lb_z = brain_pearson_mv_ageoverlap_perc$ci.lb,
  ci.ub_z = brain_pearson_mv_ageoverlap_perc$ci.ub,
  pval = brain_pearson_mv_ageoverlap_perc$pval
)

write.csv(
  ageoverlap_perc_est,
  file.path(
    results_dir_brain_ageoverlap_perc,
    "ageoverlap_perc_results.csv"
  ),
  row.names = FALSE
)

# --- I2 calculation ---

W <- diag(1 / brain_pearson_mv_ageoverlap_perc$vi)

X <- model.matrix(brain_pearson_mv_ageoverlap_perc)

P <- W - W %*% X %*%
  solve(t(X) %*% W %*% X) %*%
  t(X) %*% W

I2_total <- 100 * sum(brain_pearson_mv_ageoverlap_perc$sigma2) / (
  sum(brain_pearson_mv_ageoverlap_perc$sigma2) +
    (brain_pearson_mv_ageoverlap_perc$k -
       brain_pearson_mv_ageoverlap_perc$p) /
    sum(diag(P))
)

I2_components <- 100 *
  brain_pearson_mv_ageoverlap_perc$sigma2 / (
    sum(brain_pearson_mv_ageoverlap_perc$sigma2) +
      (brain_pearson_mv_ageoverlap_perc$k -
         brain_pearson_mv_ageoverlap_perc$p) /
      sum(diag(P))
  )

names(I2_components) <- c(
  "cohort",
  "timepoint_within_cohort",
  "model"
)

heterogeneity_partition <- 100 *
  brain_pearson_mv_ageoverlap_perc$sigma2 /
  sum(brain_pearson_mv_ageoverlap_perc$sigma2)

names(heterogeneity_partition) <- c(
  "cohort",
  "timepoint_within_cohort",
  "model"
)

# --- predicted values ---

newdat_ageoverlap_perc <- data.frame(
  overlap_perc_brain = seq(
    min(perf_brain$overlap_perc_brain, na.rm = TRUE),
    max(perf_brain$overlap_perc_brain, na.rm = TRUE),
    length.out = 100
  )
)

pred_ageoverlap_perc <- predict(
  brain_pearson_mv_ageoverlap_perc,
  newmods = newdat_ageoverlap_perc$overlap_perc_brain
)

newdat_ageoverlap_perc <- newdat_ageoverlap_perc %>%
  mutate(
    pred_z = pred_ageoverlap_perc$pred,
    z_ci.lb = pred_ageoverlap_perc$ci.lb,
    z_ci.ub = pred_ageoverlap_perc$ci.ub,
    
    pred_r = transf.ztor(pred_z),
    ci.lb = transf.ztor(z_ci.lb),
    ci.ub = transf.ztor(z_ci.ub)
  )

write.csv(
  newdat_ageoverlap_perc,
  file.path(
    results_dir_brain_ageoverlap_perc,
    "ageoverlap_perc_predicted_values.csv"
  ),
  row.names = FALSE
)

# --- cohort-timepoint means ---

cohort_time_means_ageoverlap_perc <- perf_brain %>%
  group_by(cohort, timepoint) %>%
  summarise(
    overlap_perc_brain =
      mean(overlap_perc_brain, na.rm = TRUE),
    
    mean_pearson =
      mean(Pearson, na.rm = TRUE),
    
    .groups = "drop"
  )

# --- plot predicted Pearson r ---

png(
  filename = file.path(
    results_dir_brain_ageoverlap_perc,
    "plot_ageoverlap_perc_pearson_meta.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot() +
  
  geom_point(
    data = perf_brain,
    aes(
      x = overlap_perc_brain,
      y = Pearson
    ),
    alpha = 0.25
  ) +
  
  geom_point(
    data = cohort_time_means_ageoverlap_perc,
    aes(
      x = overlap_perc_brain,
      y = mean_pearson
    ),
    size = 3
  ) +
  
  geom_line(
    data = newdat_ageoverlap_perc,
    aes(
      x = overlap_perc_brain,
      y = pred_r
    ),
    linewidth = 1
  ) +
  
  geom_ribbon(
    data = newdat_ageoverlap_perc,
    aes(
      x = overlap_perc_brain,
      ymin = ci.lb,
      ymax = ci.ub
    ),
    alpha = 0.15
  ) +
  
  scale_y_continuous(
    limits = c(-1, 1)
  ) +
  
  labs(
    x = "Brain-age training overlap (%)",
    y = "Predicted Pearson r",
    title = "Meta-regression of brain-age correlation by proportional age overlap"
  ) +
  
  theme_minimal(base_size = 14)

dev.off()

# --- save outputs ---

saveRDS(
  brain_pearson_mv_ageoverlap_perc,
  file.path(
    results_dir_brain_ageoverlap_perc,
    "model_brain_pearson_mv_ageoverlap_perc.rds"
  )
)

capture.output(
  summary(brain_pearson_mv_ageoverlap_perc),
  file = file.path(
    results_dir_brain_ageoverlap_perc,
    "model_summary.txt"
  )
)

# --- LOO cohort ---

cohorts <- unique(perf_brain$cohort)

loo_ageoverlap_perc <- lapply(cohorts, function(cc) {
  
  tmp <- perf_brain %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = pearson_z,
    V  = pearson_var,
    mods = ~ overlap_perc_brain,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model
    ),
    data = tmp,
    method = "REML"
  )
  
  tibble(
    cohort_removed = cc,
    
    intercept_z = coef(fit)[1],
    intercept_z_se = sqrt(vcov(fit)[1, 1]),
    
    overlap_perc_effect_z = coef(fit)[2],
    overlap_perc_effect_z_se = sqrt(vcov(fit)[2, 2]),
    
    overlap_perc_ci.lb_z =
      overlap_perc_effect_z -
      1.96 * overlap_perc_effect_z_se,
    
    overlap_perc_ci.ub_z =
      overlap_perc_effect_z +
      1.96 * overlap_perc_effect_z_se
  )
}) %>%
  bind_rows()

loo_ageoverlap_perc_summary <- tibble(
  
  min_overlap_perc_effect_z =
    min(loo_ageoverlap_perc$overlap_perc_effect_z, na.rm = TRUE),
  
  max_overlap_perc_effect_z =
    max(loo_ageoverlap_perc$overlap_perc_effect_z, na.rm = TRUE),
  
  range_overlap_perc_effect_z =
    max_overlap_perc_effect_z -
    min_overlap_perc_effect_z,
  
  sd_overlap_perc_effect_z =
    sd(loo_ageoverlap_perc$overlap_perc_effect_z, na.rm = TRUE)
)

write.csv(
  loo_ageoverlap_perc,
  file.path(
    results_dir_brain_ageoverlap_perc,
    "leave_one_cohort_out_ageoverlap_perc.csv"
  ),
  row.names = FALSE
)

write.csv(
  loo_ageoverlap_perc_summary,
  file.path(
    results_dir_brain_ageoverlap_perc,
    "leave_one_cohort_out_ageoverlap_perc_summary.csv"
  ),
  row.names = FALSE
)

# --- LOO plot ---

png(
  filename = file.path(
    results_dir_brain_ageoverlap_perc,
    "plot_ageoverlap_perc_loo.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_ageoverlap_perc,
  aes(
    x = reorder(cohort_removed, overlap_perc_effect_z),
    y = overlap_perc_effect_z
  )
) +
  
  geom_point(size = 3) +
  
  geom_errorbar(
    aes(
      ymin = overlap_perc_ci.lb_z,
      ymax = overlap_perc_ci.ub_z
    ),
    width = 0.15
  ) +
  
  geom_hline(
    yintercept = coef(brain_pearson_mv_ageoverlap_perc)[2],
    linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = 0,
    colour = "grey60"
  ) +
  
  coord_flip() +
  
  labs(
    x = "Cohort removed",
    y = "Overlap percentage effect on Fisher-z correlation",
    title = "Leave-one-cohort-out sensitivity of proportional age-overlap effect"
  ) +
  
  theme_minimal(base_size = 14)

dev.off()

# --- full analysis log ---

sink(log_file)

cat("Brain Pearson multilevel meta-analysis - age overlap percentage moderator\n")
cat("=========================================================================\n\n")

cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_brain_ageoverlap_perc, "\n\n")

cat("Model summary, Fisher-z scale\n")
cat("=============================\n")
print(summary(brain_pearson_mv_ageoverlap_perc))

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")

cat("Total I2 (%):\n")
print(I2_total)

cat("\nI2 components (% of total observed variance):\n")
print(I2_components)

cat("\nHeterogeneity partition (% of random-effect heterogeneity only):\n")
print(heterogeneity_partition)

cat("\n\nPredicted values, back-transformed to Pearson r\n")
cat("===============================================\n")
print(newdat_ageoverlap_perc)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_ageoverlap_perc, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_ageoverlap_perc_summary, n = Inf, width = Inf)

cat("\n\nSession info\n")
cat("============\n")
print(sessionInfo())

sink()

###### brain Pearson processing global ------

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_brain_processing <- file.path(
  results_dir,
  paste0("perf_brain_pearson_global_processing_", run_id)
)

dir.create(results_dir_brain_processing, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  results_dir_brain_processing,
  "analysis_log_processing.txt"
)

# --- data prep ---

perf_brain <- perf_brain %>%
  mutate(
    processing = factor(
      processing,
      levels = c("minimal", "extracted")
    )
  )

# --- meta ---

brain_pearson_mv_processing <- rma.mv(
  yi = pearson_z,
  V  = pearson_var,
  mods = ~ processing,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = perf_brain,
  method = "REML"
)

summary(brain_pearson_mv_processing)

# --- model results ---

processing_est <- tibble(
  term = names(coef(brain_pearson_mv_processing)),
  estimate_z = as.numeric(coef(brain_pearson_mv_processing)),
  se_z = brain_pearson_mv_processing$se,
  ci.lb_z = brain_pearson_mv_processing$ci.lb,
  ci.ub_z = brain_pearson_mv_processing$ci.ub,
  pval = brain_pearson_mv_processing$pval
)

write.csv(
  processing_est,
  file.path(
    results_dir_brain_processing,
    "processing_results.csv"
  ),
  row.names = FALSE
)

# --- I2 calculation ---

W <- diag(1 / brain_pearson_mv_processing$vi)
X <- model.matrix(brain_pearson_mv_processing)

P <- W - W %*% X %*%
  solve(t(X) %*% W %*% X) %*%
  t(X) %*% W

I2_total <- 100 * sum(brain_pearson_mv_processing$sigma2) / (
  sum(brain_pearson_mv_processing$sigma2) +
    (brain_pearson_mv_processing$k - brain_pearson_mv_processing$p) /
    sum(diag(P))
)

I2_components <- 100 * brain_pearson_mv_processing$sigma2 / (
  sum(brain_pearson_mv_processing$sigma2) +
    (brain_pearson_mv_processing$k - brain_pearson_mv_processing$p) /
    sum(diag(P))
)

names(I2_components) <- c(
  "cohort",
  "timepoint_within_cohort",
  "model"
)

heterogeneity_partition <- 100 * brain_pearson_mv_processing$sigma2 /
  sum(brain_pearson_mv_processing$sigma2)

names(heterogeneity_partition) <- c(
  "cohort",
  "timepoint_within_cohort",
  "model"
)

# --- predicted values ---

newdat_processing <- data.frame(
  processing = levels(perf_brain$processing)
)

newdat_processing <- newdat_processing %>%
  mutate(
    processing = factor(
      processing,
      levels = levels(perf_brain$processing)
    )
  )

X_new_processing <- model.matrix(
  ~ processing,
  data = newdat_processing
)

pred_processing <- predict(
  brain_pearson_mv_processing,
  newmods = X_new_processing[, -1, drop = FALSE]
)

newdat_processing <- newdat_processing %>%
  mutate(
    pred_z = pred_processing$pred,
    z_ci.lb = pred_processing$ci.lb,
    z_ci.ub = pred_processing$ci.ub,
    
    pred_r = transf.ztor(pred_z),
    ci.lb = transf.ztor(z_ci.lb),
    ci.ub = transf.ztor(z_ci.ub)
  )

write.csv(
  newdat_processing,
  file.path(
    results_dir_brain_processing,
    "processing_predicted_values.csv"
  ),
  row.names = FALSE
)

# --- plot predicted Pearson r ---

png(
  filename = file.path(
    results_dir_brain_processing,
    "plot_processing_pearson_meta.png"
  ),
  width = 1600,
  height = 1200,
  res = 300
)

ggplot() +
  geom_jitter(
    data = perf_brain,
    aes(
      x = processing,
      y = Pearson
    ),
    width = 0.12,
    alpha = 0.25
  ) +
  geom_point(
    data = newdat_processing,
    aes(
      x = processing,
      y = pred_r
    ),
    size = 4
  ) +
  geom_errorbar(
    data = newdat_processing,
    aes(
      x = processing,
      ymin = ci.lb,
      ymax = ci.ub
    ),
    width = 0.15
  ) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(
    x = "Brain age training overlap",
    y = "Predicted Pearson r",
    title = "Meta-regression of brain-age correlation by processing level"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- save model outputs ---

saveRDS(
  brain_pearson_mv_processing,
  file.path(
    results_dir_brain_processing,
    "model_brain_pearson_mv_processing.rds"
  )
)

capture.output(
  summary(brain_pearson_mv_processing),
  file = file.path(
    results_dir_brain_processing,
    "model_summary.txt"
  )
)

# --- LOO cohort ---

cohorts <- unique(perf_brain$cohort)

loo_processing <- lapply(cohorts, function(cc) {
  
  tmp <- perf_brain %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = pearson_z,
    V  = pearson_var,
    mods = ~ processing,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model
    ),
    data = tmp,
    method = "REML"
  )
  
  tibble(
    cohort_removed = cc,
    
    intercept_z = coef(fit)[1],
    intercept_z_se = sqrt(vcov(fit)[1, 1]),
    
    processing_effect_z = coef(fit)[2],
    processing_effect_z_se = sqrt(vcov(fit)[2, 2]),
    
    processing_ci.lb_z = processing_effect_z - 1.96 * processing_effect_z_se,
    processing_ci.ub_z = processing_effect_z + 1.96 * processing_effect_z_se
  )
}) %>%
  bind_rows()

loo_processing_summary <- tibble(
  min_processing_effect_z = min(loo_processing$processing_effect_z, na.rm = TRUE),
  max_processing_effect_z = max(loo_processing$processing_effect_z, na.rm = TRUE),
  range_processing_effect_z = max_processing_effect_z - min_processing_effect_z,
  sd_processing_effect_z = sd(loo_processing$processing_effect_z, na.rm = TRUE)
)

write.csv(
  loo_processing,
  file.path(
    results_dir_brain_processing,
    "leave_one_cohort_out_processing.csv"
  ),
  row.names = FALSE
)

write.csv(
  loo_processing_summary,
  file.path(
    results_dir_brain_processing,
    "leave_one_cohort_out_processing_summary.csv"
  ),
  row.names = FALSE
)

# --- LOO plot: processing effect on Fisher-z scale ---

png(
  filename = file.path(
    results_dir_brain_processing,
    "plot_processing_loo.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_processing,
  aes(
    x = reorder(cohort_removed, processing_effect_z),
    y = processing_effect_z
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = processing_ci.lb_z,
      ymax = processing_ci.ub_z
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept = coef(brain_pearson_mv_processing)[2],
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 0,
    colour = "grey60"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "processing effect on Fisher-z correlation",
    title = "Leave-one-cohort-out sensitivity of processing effect"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- full analysis log ---

sink(log_file)

cat("Brain Pearson multilevel meta-analysis - processing moderator\n")
cat("=====================================================================\n\n")

cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_brain_processing, "\n\n")

cat("Model summary, Fisher-z scale\n")
cat("=============================\n")
print(summary(brain_pearson_mv_processing))

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")

cat("Total I2 (%):\n")
print(I2_total)

cat("\nI2 components (% of total observed variance):\n")
print(I2_components)

cat("\nHeterogeneity partition (% of random-effect heterogeneity only):\n")
print(heterogeneity_partition)

cat("\n\nPredicted values, back-transformed to Pearson r\n")
cat("===============================================\n")
print(newdat_processing)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_processing, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_processing_summary, n = Inf, width = Inf)

cat("\n\nSession info\n")
cat("============\n")
print(sessionInfo())

sink()
