#### META-ANALYSES SCRIPT - unweighted MAE
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


#### unweighted MAE analyses ----

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

# --- prep ---

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

# create MAE sampling variance from bootstrapped SE
perf_epi <- perf_epi %>%
  mutate(
    MAE_var_boot = MAE_SE_boot^2
  )

#### analyses perf epi ----

###### epi MAE global ------

# --- prep ---
# exclude specific models if needed 
perf_epi_sub   <- perf_epi   %>% filter(!model %in% c("DunedinPACE", "DNAmTL")) %>% droplevels()

# rename gestational clocks to not have _conv as suffix 
perf_epi_sub <- perf_epi_sub %>%
  mutate(
    model = recode(
      model,
      "Bohlin_conv" = "Bohlin",
      "Knight_conv" = "Knight",
      "EPIC_conv" = "EPIC"
    )
  )

# meta-analysis did not converge due to differences in sampling variances which span ~7 orders of magnitude from ~0.00002 to ~34,800 when squared (also tried out scaling based on median but didn't work)

# take the log 
perf_epi_sub <- perf_epi_sub |>
  mutate(
    log_MAE = log(MAE),
    # Delta method: Var(log(X)) ≈ (SE_X / X)^2
    log_MAE_var = (MAE_SE_boot / MAE)^2
  )


# --- set up folder ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi <- file.path(results_dir, paste0("perf_epi_MAE_noDunedinDNAmTL_global_BFGSoptimiser_", run_id))

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

# --- meta ---
epi_MAE_mv_global <- rma.mv(
  yi = log_MAE,
  V  = log_MAE_var,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data   = perf_epi_sub,
  method = "REML",
  control = list(
    optimizer = "optim",
    optmethod = "BFGS",
    maxit = 10000
  )
)

summary(epi_MAE_mv_global)

# heterogeneity estimates

W <- diag(1/epi_MAE_mv_global$vi)
X <- model.matrix(epi_MAE_mv_global)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_total <- 100 * sum(epi_MAE_mv_global$sigma2) / (sum(epi_MAE_mv_global$sigma2) + (epi_MAE_mv_global$k-epi_MAE_mv_global$p)/sum(diag(P)))


I2_components <- 100 * epi_MAE_mv_global$sigma2 / (
  sum(epi_MAE_mv_global$sigma2) +
    (epi_MAE_mv_global$k - epi_MAE_mv_global$p) /
    sum(diag(P))
)

names(I2_components) <- c(
  "cohort",
  "timepoint_within_cohort",
  "model")

I2_total
I2_components

# Back-transform global pooled estimate
overall_est <- data.frame(
  model      = "Pooled",
  MAE       = exp(epi_MAE_mv_global$beta[1, 1]),
  ci.lb      = exp(epi_MAE_mv_global$ci.lb[1]),
  ci.ub      = exp(epi_MAE_mv_global$ci.ub[1]),
  pval = epi_MAE_mv_global$pval[1], 
  generation = "Pooled"
)

overall_est

saveRDS(
  epi_MAE_mv_global,
  file.path(results_dir_epi, "model_epi_MAE_mv_global.rds")
)

capture.output(
  summary(epi_MAE_mv_global),
  file = file.path(results_dir_epi, "model_summary.txt")
)

# --- LOO cohort (for log-transformed model) ---

cohorts <- unique(perf_epi_sub$cohort)

loo_cohort <- lapply(cohorts, function(cc) {
  
  tmp <- subset(perf_epi_sub, cohort != cc)
  
  fit <- rma.mv(
    yi = log_MAE,
    V  = log_MAE_var,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model
    ),
    data = tmp,
    method = "REML",
    control = list(
      optimizer = "optim",
      optmethod = "BFGS",
      maxit = 10000
    )
  )
  
  log_est <- coef(fit)[1]
  log_se  <- sqrt(vcov(fit)[1, 1])
  
  data.frame(
    cohort_removed = cc,
    
    # log scale
    log_MAE = log_est,
    log_se   = log_se,
    log_ci.lb = log_est - 1.96 * log_se,
    log_ci.ub = log_est + 1.96 * log_se,
    
    # back-transformed scale
    MAE = exp(log_est),
    ci.lb = exp(log_est - 1.96 * log_se),
    ci.ub = exp(log_est + 1.96 * log_se)
  )
})

loo_cohort <- do.call(rbind, loo_cohort) %>%
  as_tibble()

loo_summary <- tibble(
  min_MAE   = min(loo_cohort$MAE, na.rm = TRUE),
  max_MAE   = max(loo_cohort$MAE, na.rm = TRUE),
  range_MAE = max_MAE - min_MAE,
  sd_MAE    = sd(loo_cohort$MAE, na.rm = TRUE),
  
  min_log_MAE   = min(loo_cohort$log_MAE, na.rm = TRUE),
  max_log_MAE   = max(loo_cohort$log_MAE, na.rm = TRUE),
  range_log_MAE = max_log_MAE - min_log_MAE,
  sd_log_MAE    = sd(loo_cohort$log_MAE, na.rm = TRUE)
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

full_MAE <- exp(coef(epi_MAE_mv_global)[1])

png(
  filename = file.path(results_dir_epi, "plot_MAE_allbD_loo.png"),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_cohort,
  aes(
    x = reorder(cohort_removed, MAE),
    y = MAE
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
    yintercept = full_MAE,
    linetype = "dashed"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "Pooled MAE",
    title = "Leave-one-cohort-out sensitivity of pooled MAE"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# full analysis log

sink(log_file)

cat("Epi MAE multilevel meta-analysis (log-transformed)\n")
cat("===================================\n\n")
cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_epi, "\n\n")

cat("Model summary\n")
cat("=============\n")
print(summary(epi_MAE_mv_global))

cat("Model summary - back transformed \n")
cat("=============\n")
print(overall_est)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")

cat("Total I2 (%):\n")
print(I2_total)

cat("\nVariance decomposition (% of total heterogeneity):\n")
print(I2_components)

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

###### epi MAE modelwise ------
# --- set up folder ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi_modelwise <- file.path(
  results_dir,
  paste0("perf_epi_MAE_noDunedinDNAmTL_modelwise_BFGSoptimiser_", run_id)
)

dir.create(results_dir_epi_modelwise, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_epi_modelwise, "analysis_modelwise_log.txt")

# --- meta ---

epi_MAE_mv_modelwise <- rma.mv(
  yi = log_MAE,
  V  = log_MAE_var,
  mods = ~ model - 1,
  random = list(
    ~ 1 | cohort/timepoint
  ),
  data = perf_epi_sub,
  method = "REML",
  control = list(
    optimizer = "optim",
    optmethod = "BFGS",
    maxit = 10000
  )
)

summary(epi_MAE_mv_modelwise)

# --- heterogeneity estimates ---

W <- diag(1 / epi_MAE_mv_modelwise$vi)
X <- model.matrix(epi_MAE_mv_modelwise)

P <- W - W %*% X %*%
  solve(t(X) %*% W %*% X) %*%
  t(X) %*% W

I2_total <- 100 * sum(epi_MAE_mv_modelwise$sigma2) /
  (
    sum(epi_MAE_mv_modelwise$sigma2) +
      (epi_MAE_mv_modelwise$k - epi_MAE_mv_modelwise$p) / sum(diag(P))
  )

I2_components <- 100 * epi_MAE_mv_modelwise$sigma2 /
  (
    sum(epi_MAE_mv_modelwise$sigma2) +
      (epi_MAE_mv_modelwise$k - epi_MAE_mv_modelwise$p) / sum(diag(P))
  )

names(I2_components) <- c(
  "cohort",
  "timepoint_within_cohort"
)

I2_total
I2_components

# --- back-transform modelwise pooled estimates ---

epi_MAE_est <- data.frame(
  model = gsub("^model", "", rownames(epi_MAE_mv_modelwise$beta)),
  log_MAE = epi_MAE_mv_modelwise$beta[, 1],
  log_se = epi_MAE_mv_modelwise$se,
  log_ci.lb = epi_MAE_mv_modelwise$ci.lb,
  log_ci.ub = epi_MAE_mv_modelwise$ci.ub,
  MAE = exp(epi_MAE_mv_modelwise$beta[, 1]),
  ci.lb = exp(epi_MAE_mv_modelwise$ci.lb),
  ci.ub = exp(epi_MAE_mv_modelwise$ci.ub),
  pval = epi_MAE_mv_modelwise$pval
) 

X <- model.matrix(epi_MAE_mv_modelwise)

model_cols <- colnames(X)[
  str_detect(colnames(X), "^model") &
    !str_detect(colnames(X), ":")
]

effect_sizes_per_model <- tibble(
  model = str_remove(model_cols, "^model"),
  k = colSums(X[, model_cols, drop = FALSE] != 0)
)


epi_MAE_est <- epi_MAE_est %>%
  left_join(effect_sizes_per_model, by = "model")


epi_MAE_est

write.csv(
  epi_MAE_est,
  file.path(results_dir_epi_modelwise, "modelwise_MAE_estimates.csv"),
  row.names = FALSE
)

saveRDS(
  epi_MAE_mv_modelwise,
  file.path(results_dir_epi_modelwise, "model_epi_MAE_mv_modelwise.rds")
)

capture.output(
  summary(epi_MAE_mv_modelwise),
  file = file.path(results_dir_epi_modelwise, "model_summary.txt")
)

# --- LOO cohort modelwise ---

cohorts <- unique(perf_epi_sub$cohort)

loo_cohort_modelwise <- lapply(cohorts, function(cc) {
  
  tmp <- subset(perf_epi_sub, cohort != cc)
  
  fit <- rma.mv(
    yi = log_MAE,
    V  = log_MAE_var,
    mods = ~ model - 1,
    random = list(
      ~ 1 | cohort/timepoint
    ),
    data = tmp,
    method = "REML",
    control = list(
      optimizer = "optim",
      optmethod = "BFGS",
      maxit = 10000
    )
  )
  
  tibble(
    cohort_removed = cc,
    term = rownames(fit$b),
    
    log_MAE = as.numeric(fit$b),
    log_se = fit$se,
    log_ci.lb = fit$ci.lb,
    log_ci.ub = fit$ci.ub,
    pval = fit$pval
  ) %>%
    mutate(
      model = str_remove(term, "^model"),
      
      MAE = exp(log_MAE),
      ci.lb = exp(log_ci.lb),
      ci.ub = exp(log_ci.ub)
    ) %>%
    select(
      cohort_removed,
      model,
      log_MAE,
      log_se,
      log_ci.lb,
      log_ci.ub,
      MAE,
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
    min_MAE = min(MAE, na.rm = TRUE),
    max_MAE = max(MAE, na.rm = TRUE),
    range_MAE = max_MAE - min_MAE,
    sd_MAE = sd(MAE, na.rm = TRUE),
    
    min_log_MAE = min(log_MAE, na.rm = TRUE),
    max_log_MAE = max(log_MAE, na.rm = TRUE),
    range_log_MAE = max_log_MAE - min_log_MAE,
    sd_log_MAE = sd(log_MAE, na.rm = TRUE),
    
    .groups = "drop"
  )

write.csv(
  modelwise_loo_summary,
  file.path(results_dir_epi_modelwise, "modelwise_leave_one_cohort_out_summary.csv"),
  row.names = FALSE
)

# --- optional LOO plot by model ---

png(
  filename = file.path(results_dir_epi_modelwise, "plot_modelwise_MAE_loo.png"),
  width = 2200,
  height = 1400,
  res = 300
)

ggplot(
  loo_cohort_modelwise,
  aes(
    x = reorder(cohort_removed, MAE),
    y = MAE
  )
) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(
      ymin = ci.lb,
      ymax = ci.ub
    ),
    width = 0.15
  ) +
  facet_wrap(~ model, scales = "free_y") +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "Pooled MAE",
    title = "Modelwise leave-one-cohort-out sensitivity of pooled MAE"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- full analysis log ---

sink(log_file)

cat("Epi MAE multilevel meta-analysis modelwise (log-transformed)\n")
cat("===========================================================\n\n")
cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_epi_modelwise, "\n\n")

cat("Model summary\n")
cat("=============\n")
print(summary(epi_MAE_mv_modelwise))

cat("\n\nModel summary - back transformed\n")
cat("================================\n")
print(epi_MAE_est)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
cat("Total I2 (%):\n")
print(I2_total)

cat("\nVariance decomposition (% of total heterogeneity):\n")
print(I2_components)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_cohort_modelwise, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(modelwise_loo_summary, n = Inf, width = Inf)

sink()

###### epi MAE age global (not run) ------

# --- set up ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi_global_ageint <- file.path(
  results_dir,
  paste0("perf_epi_MAE_noDunedinDNAmTL_global_ageint_BFGSoptimiser_", run_id)
)

dir.create(results_dir_epi_global_ageint, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_epi_global_ageint, "analysis_log_global_ageint.txt")

# --- data prep ---

perf_epi_sub <- perf_epi_sub %>%
  mutate(
    mean_age_c = mean_age - mean(mean_age, na.rm = TRUE)
  )

# remove gestational clocks for age moderation
perf_epi_sub_nogest <- perf_epi_sub %>%
  filter(!model %in% c("Bohlin", "Knight", "EPIC")) %>%
  droplevels()

# --- meta ---

epi_MAE_mv_global_ageint <- rma.mv(
  yi = log_MAE,
  V  = log_MAE_var,
  mods = ~ mean_age_c,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = perf_epi_sub_nogest,
  method = "REML",
  control = list(
    optimizer = "optim",
    optmethod = "BFGS",
    maxit = 10000
  )
)

summary(epi_MAE_mv_global_ageint)

# --- heterogeneity estimates ---

W <- diag(1 / epi_MAE_mv_global_ageint$vi)
X <- model.matrix(epi_MAE_mv_global_ageint)

P <- W - W %*% X %*%
  solve(t(X) %*% W %*% X) %*%
  t(X) %*% W

I2_total <- 100 * sum(epi_MAE_mv_global_ageint$sigma2) /
  (
    sum(epi_MAE_mv_global_ageint$sigma2) +
      (epi_MAE_mv_global_ageint$k - epi_MAE_mv_global_ageint$p) / sum(diag(P))
  )

I2_components <- 100 * epi_MAE_mv_global_ageint$sigma2 /
  (
    sum(epi_MAE_mv_global_ageint$sigma2) +
      (epi_MAE_mv_global_ageint$k - epi_MAE_mv_global_ageint$p) / sum(diag(P))
  )

names(I2_components) <- c(
  "cohort",
  "timepoint_within_cohort",
  "model"
)

I2_total
I2_components

# --- back-transform main model results ---

age_model_est <- tibble(
  term = names(coef(epi_MAE_mv_global_ageint)),
  estimate_log = as.numeric(coef(epi_MAE_mv_global_ageint)),
  se_log = epi_MAE_mv_global_ageint$se,
  ci.lb_log = epi_MAE_mv_global_ageint$ci.lb,
  ci.ub_log = epi_MAE_mv_global_ageint$ci.ub,
  pval = epi_MAE_mv_global_ageint$pval
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
  epi_MAE_mv_global_ageint,
  file.path(results_dir_epi_global_ageint, "model_epi_MAE_mv_global_ageint_logscale.rds")
)

capture.output(
  summary(epi_MAE_mv_global_ageint),
  file = file.path(results_dir_epi_global_ageint, "model_summary_log_scale.txt")
)

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

pred <- predict(
  epi_MAE_mv_global_ageint,
  newmods = newdat$mean_age_c
)

newdat <- newdat %>%
  mutate(
    log_pred = pred$pred,
    log_ci.lb = pred$ci.lb,
    log_ci.ub = pred$ci.ub,
    pred = exp(log_pred),
    ci.lb = exp(log_ci.lb),
    ci.ub = exp(log_ci.ub)
  )

cohort_time_means <- perf_epi_sub_nogest %>%
  group_by(cohort, timepoint) %>%
  summarise(
    mean_age = mean(mean_age, na.rm = TRUE),
    mean_MAE = mean(MAE, na.rm = TRUE),
    .groups = "drop"
  )

png(
  filename = file.path(
    results_dir_epi_global_ageint,
    "plot_meanage_MAE_meta_regression_backtransformed.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot() +
  geom_point(
    data = perf_epi_sub_nogest,
    aes(
      x = mean_age,
      y = MAE,
      colour = cohort
    ),
    alpha = .25,
    size = 2
  ) +
  geom_point(
    data = cohort_time_means,
    aes(
      x = mean_age,
      y = mean_MAE,
      colour = cohort
    ),
    size = 4
  ) +
  geom_line(
    data = newdat,
    aes(
      x = mean_age,
      y = pred
    ),
    linewidth = 1.2,
    colour = "black"
  ) +
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
    title = "Meta-regression of epigenetic clock MAE on mean age",
    x = "Mean sample age",
    y = "Geometric mean MAE",
    colour = "Cohort"
  ) +
  theme_minimal(base_size = 16)

dev.off()

# --- LOO cohort ---

cohorts <- unique(perf_epi_sub_nogest$cohort)

loo_age <- lapply(cohorts, function(cc) {
  
  tmp <- perf_epi_sub_nogest %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = log_MAE,
    V  = log_MAE_var,
    mods = ~ mean_age_c,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model
    ),
    data = tmp,
    method = "REML",
    control = list(
      optimizer = "optim",
      optmethod = "BFGS",
      maxit = 10000
    )
  )
  
  log_intercept <- coef(fit)[1]
  log_slope <- coef(fit)[2]
  
  log_intercept_se <- sqrt(vcov(fit)[1, 1])
  log_slope_se <- sqrt(vcov(fit)[2, 2])
  
  data.frame(
    cohort_removed = cc,
    
    log_intercept = log_intercept,
    log_slope = log_slope,
    log_intercept_se = log_intercept_se,
    log_slope_se = log_slope_se,
    
    intercept_MAE = exp(log_intercept),
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
  min_intercept_MAE = min(loo_age_df$intercept_MAE, na.rm = TRUE),
  max_intercept_MAE = max(loo_age_df$intercept_MAE, na.rm = TRUE),
  range_intercept_MAE = max_intercept_MAE - min_intercept_MAE,
  sd_intercept_MAE = sd(loo_age_df$intercept_MAE, na.rm = TRUE),
  
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

full_percent_slope <- (exp(coef(epi_MAE_mv_global_ageint)[2]) - 1) * 100

png(
  filename = file.path(results_dir_epi_global_ageint, "plot_meanage_MAE_loo_percent_change.png"),
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
    y = "% change in geometric mean MAE per year",
    title = "Leave-one-cohort-out sensitivity of age effect"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- log ---

sink(log_file)

cat("Epi MAE multilevel meta-analysis - global age moderator, log-transformed\n")
cat("=======================================================================\n\n")
cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_epi_global_ageint, "\n\n")

cat("Model summary, log scale\n")
cat("========================\n")
print(summary(epi_MAE_mv_global_ageint))

cat("\n\nBack-transformed model results\n")
cat("==============================\n")
print(age_model_est, n = Inf, width = Inf)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
cat("Total I2 (%):\n")
print(I2_total)

cat("\nVariance decomposition (% of total heterogeneity):\n")
print(I2_components)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_age_df, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_age_summary, n = Inf, width = Inf)

sink()

###### epi MAE age modelwise (not run) ------

# --- set up ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi_modelwise_ageint <- file.path(
  results_dir,
  paste0("perf_epi_MAE_noDunedinDNAmTL_modelwise_ageint_BFGSoptimiser_", run_id)
)

dir.create(results_dir_epi_modelwise_ageint, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_epi_modelwise_ageint, "analysis_log_modelwise_ageint.txt")

# --- data prep ---

perf_epi_sub_nogest <- perf_epi_sub_nogest %>%
  mutate(
    mean_age_c = mean_age - mean(mean_age, na.rm = TRUE)
  )

# --- meta ---

epi_MAE_mv_modelwise_ageint <- rma.mv(
  yi = log_MAE,
  V  = log_MAE_var,
  mods = ~ model + model:mean_age_c - 1,
  random = list(
    ~ 1 | cohort/timepoint
  ),
  data = perf_epi_sub_nogest,
  method = "REML",
  control = list(
    optimizer = "optim",
    optmethod = "BFGS",
    maxit = 10000
  )
)

summary(epi_MAE_mv_modelwise_ageint)

# --- heterogeneity estimates ---

W <- diag(1 / epi_MAE_mv_modelwise_ageint$vi)
X <- model.matrix(epi_MAE_mv_modelwise_ageint)

P <- W - W %*% X %*%
  solve(t(X) %*% W %*% X) %*%
  t(X) %*% W

I2_total <- 100 * sum(epi_MAE_mv_modelwise_ageint$sigma2) /
  (
    sum(epi_MAE_mv_modelwise_ageint$sigma2) +
      (epi_MAE_mv_modelwise_ageint$k - epi_MAE_mv_modelwise_ageint$p) / sum(diag(P))
  )

I2_components <- 100 * epi_MAE_mv_modelwise_ageint$sigma2 /
  (
    sum(epi_MAE_mv_modelwise_ageint$sigma2) +
      (epi_MAE_mv_modelwise_ageint$k - epi_MAE_mv_modelwise_ageint$p) / sum(diag(P))
  )

names(I2_components) <- c(
  "cohort",
  "timepoint_within_cohort"
)

I2_total
I2_components

# compare common age slope model with model-specific age slopes
anova(
  epi_MAE_mv_global_ageint,
  epi_MAE_mv_modelwise_ageint,
  refit = TRUE
)

modelwise_age_est <- tibble(
  term = rownames(epi_MAE_mv_modelwise_ageint$b),
  estimate_log = as.numeric(epi_MAE_mv_modelwise_ageint$b),
  se_log = epi_MAE_mv_modelwise_ageint$se,
  ci.lb_log = epi_MAE_mv_modelwise_ageint$ci.lb,
  ci.ub_log = epi_MAE_mv_modelwise_ageint$ci.ub,
  pval = epi_MAE_mv_modelwise_ageint$pval
) %>%
  mutate(
    type = if_else(str_detect(term, "mean_age_c"), "age_slope", "intercept"),
    model = term %>%
      str_remove("^model") %>%
      str_remove(":mean_age_c") %>%
      str_remove("\\.mean_age_c")
  ) %>%
  mutate(
    estimate_exp = exp(estimate_log),
    ci.lb_exp = exp(ci.lb_log),
    ci.ub_exp = exp(ci.ub_log),
    percent_change = (exp(estimate_log) - 1) * 100,
    percent_ci.lb = (exp(ci.lb_log) - 1) * 100,
    percent_ci.ub = (exp(ci.ub_log) - 1) * 100
  )


X <- model.matrix(epi_MAE_mv_modelwise_ageint)

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
  file.path(results_dir_epi_modelwise_ageint, "modelwise_age_interaction_backtransformed_results.csv"),
  row.names = FALSE
)

saveRDS(
  epi_MAE_mv_modelwise_ageint,
  file.path(results_dir_epi_modelwise_ageint, "model_epi_MAE_mv_modelwise_ageint.rds")
)

capture.output(
  summary(epi_MAE_mv_modelwise_ageint),
  file = file.path(results_dir_epi_modelwise_ageint, "model_summary_log_scale.txt")
)

# --- model-specific age-effect plot ---

mean_age_mean <- mean(perf_epi_sub_nogest$mean_age, na.rm = TRUE)

newdat <- expand.grid(
  mean_age_c = seq(
    min(perf_epi_sub_nogest$mean_age_c, na.rm = TRUE),
    max(perf_epi_sub_nogest$mean_age_c, na.rm = TRUE),
    length.out = 100
  ),
  model = unique(perf_epi_sub_nogest$model)
)

pred <- predict(
  epi_MAE_mv_modelwise_ageint,
  newmods = model.matrix(~ model + model:mean_age_c - 1, newdat)
)

newdat <- newdat %>%
  mutate(
    mean_age = mean_age_c + mean_age_mean,
    log_pred = pred$pred,
    log_ci.lb = pred$ci.lb,
    log_ci.ub = pred$ci.ub,
    pred = exp(log_pred),
    ci.lb = exp(log_ci.lb),
    ci.ub = exp(log_ci.ub)
  )

png(
  filename = file.path(
    results_dir_epi_modelwise_ageint,
    "plot_epi_meanage_MAE_modelwise_backtransformed.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(newdat, aes(mean_age, pred, colour = model)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = ci.lb, ymax = ci.ub, fill = model),
    alpha = .12,
    colour = NA
  ) +
  labs(
    x = "Mean sample age",
    y = "Predicted geometric mean MAE",
    title = "Model-specific age effects"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- LOO cohort ---

cohorts <- unique(perf_epi_sub_nogest$cohort)

loo_modelwise_age <- lapply(cohorts, function(cc) {
  
  tmp <- perf_epi_sub_nogest %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = log_MAE,
    V  = log_MAE_var,
    mods = ~ model + model:mean_age_c - 1,
    random = list(
      ~ 1 | cohort/timepoint
    ),
    data = tmp,
    method = "REML",
    control = list(
      optimizer = "optim",
      optmethod = "BFGS",
      maxit = 10000
    )
  )
  
  tibble(
    cohort_removed = cc,
    term = rownames(fit$b),
    estimate_log = as.numeric(fit$b),
    se_log = fit$se,
    ci.lb_log = fit$ci.lb,
    ci.ub_log = fit$ci.ub,
    pval = fit$pval
  ) %>%
    mutate(
      type = if_else(str_detect(term, "mean_age_c"), "age_slope", "intercept"),
      model = term %>%
        str_remove("^model") %>%
        str_remove(":mean_age_c") %>%
        str_remove("\\.mean_age_c"),
      estimate_exp = exp(estimate_log),
      ci.lb_exp = exp(ci.lb_log),
      ci.ub_exp = exp(ci.ub_log),
      percent_change = (exp(estimate_log) - 1) * 100,
      percent_ci.lb = (exp(ci.lb_log) - 1) * 100,
      percent_ci.ub = (exp(ci.ub_log) - 1) * 100
    )
}) %>%
  bind_rows()

write.csv(
  loo_modelwise_age,
  file.path(results_dir_epi_modelwise_ageint, "modelwise_leave_one_cohort_out_age_interactions.csv"),
  row.names = FALSE
)

loo_modelwise_age_summary <- loo_modelwise_age %>%
  filter(type == "age_slope") %>%
  group_by(model) %>%
  summarise(
    min_percent_per_year = min(percent_change, na.rm = TRUE),
    max_percent_per_year = max(percent_change, na.rm = TRUE),
    range_percent_per_year = max_percent_per_year - min_percent_per_year,
    sd_percent_per_year = sd(percent_change, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  loo_modelwise_age_summary,
  file.path(results_dir_epi_modelwise_ageint, "modelwise_leave_one_cohort_out_age_interactions_summary.csv"),
  row.names = FALSE
)

# --- log ---

sink(log_file)

cat("Epi MAE multilevel meta-analysis - modelwise age interactions, log-transformed\n")
cat("=============================================================================\n\n")
cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_epi_modelwise_ageint, "\n\n")

cat("Model summary, log scale\n")
cat("========================\n")
print(summary(epi_MAE_mv_modelwise_ageint))

cat("\n\nBack-transformed modelwise age interaction results\n")
cat("==================================================\n")
print(modelwise_age_est, n = Inf, width = Inf)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
cat("Total I2 (%):\n")
print(I2_total)

cat("\nVariance decomposition (% of total heterogeneity):\n")
print(I2_components)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_modelwise_age, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary for age slopes\n")
cat("===========================================\n")
print(loo_modelwise_age_summary, n = Inf, width = Inf)

sink()

#### analyses perf brain: absolute MAE ----

##### brain MAE no DunedinPACNI -----

perf_brain_sub <- perf_brain %>%
  filter(!model %in% c("DunedinPACNI")) %>%
  droplevels()


###### brain MAE global ------

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_brain <- file.path(
  results_dir,
  paste0("perf_brain_MAE_noDunedin_global_", run_id)
)

dir.create(results_dir_brain, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_brain, "analysis_log.txt")

# --- descriptives ---

write.csv(
  perf_brain_sub,
  file.path(results_dir_brain, "perf_brain_sub_descr_detailed.csv"),
  row.names = FALSE
)

perf_brain_sub_desc <- perf_brain_sub %>%
  group_by(model, domain) %>%
  summarise(
    k = n(),
    
    MAE_mean = mean(MAE, na.rm = TRUE),
    MAE_sd   = sd(MAE, na.rm = TRUE),
    MAE_min  = min(MAE, na.rm = TRUE),
    MAE_max  = max(MAE, na.rm = TRUE),
    
    wMAE_mean = mean(wMAE_test, na.rm = TRUE),
    wMAE_sd   = sd(wMAE_test, na.rm = TRUE),
    wMAE_min  = min(wMAE_test, na.rm = TRUE),
    wMAE_max  = max(wMAE_test, na.rm = TRUE),
    
    RMSE_mean = mean(RMSE, na.rm = TRUE),
    RMSE_sd   = sd(RMSE, na.rm = TRUE),
    RMSE_min  = min(RMSE, na.rm = TRUE),
    RMSE_max  = max(RMSE, na.rm = TRUE),
    
    R2_mean       = mean(R2, na.rm = TRUE),
    Pearson_mean  = mean(Pearson, na.rm = TRUE),
    Spearman_mean = mean(Spearman, na.rm = TRUE),
    
    .groups = "drop"
  )

write.csv(
  perf_brain_sub_desc,
  file.path(results_dir_brain, "model_perf_brain_desc.csv"),
  row.names = FALSE
)

# --- meta ---

brain_MAE_mv_global <- rma.mv(
  yi = MAE,
  V  = MAE_SE_boot^2,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = perf_brain_sub,
  method = "REML"
)

summary(brain_MAE_mv_global)

# --- heterogeneity estimates ---

W <- diag(1 / brain_MAE_mv_global$vi)
X <- model.matrix(brain_MAE_mv_global)

P <- W - W %*% X %*%
  solve(t(X) %*% W %*% X) %*%
  t(X) %*% W

I2_total <- 100 * sum(brain_MAE_mv_global$sigma2) /
  (
    sum(brain_MAE_mv_global$sigma2) +
      (brain_MAE_mv_global$k - brain_MAE_mv_global$p) / sum(diag(P))
  )

I2_components <- 100 * brain_MAE_mv_global$sigma2 /
  (
    sum(brain_MAE_mv_global$sigma2) +
      (brain_MAE_mv_global$k - brain_MAE_mv_global$p) / sum(diag(P))
  )

names(I2_components) <- c(
  "cohort",
  "timepoint_within_cohort",
  "model"
)

overall_est <- data.frame(
  model = "Pooled",
  MAE = brain_MAE_mv_global$beta[1, 1],
  ci.lb = brain_MAE_mv_global$ci.lb[1],
  ci.ub = brain_MAE_mv_global$ci.ub[1],
  pval = brain_MAE_mv_global$pval[1],
  generation = "Pooled"
)

overall_est

saveRDS(
  brain_MAE_mv_global,
  file.path(results_dir_brain, "model_brain_MAE_mv_global.rds")
)

capture.output(
  summary(brain_MAE_mv_global),
  file = file.path(results_dir_brain, "model_summary.txt")
)

# --- LOO cohort ---

cohorts <- unique(perf_brain_sub$cohort)

loo_cohort <- lapply(cohorts, function(cc) {
  
  tmp <- subset(perf_brain_sub, cohort != cc)
  
  fit <- rma.mv(
    yi = MAE,
    V  = MAE_SE_boot^2,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model
    ),
    data = tmp,
    method = "REML"
  )
  
  mae_est <- coef(fit)[1]
  mae_se <- sqrt(vcov(fit)[1, 1])
  
  data.frame(
    cohort_removed = cc,
    MAE = mae_est,
    se = mae_se,
    ci.lb = mae_est - 1.96 * mae_se,
    ci.ub = mae_est + 1.96 * mae_se
  )
}) %>%
  bind_rows()

loo_summary <- tibble(
  min_MAE = min(loo_cohort$MAE, na.rm = TRUE),
  max_MAE = max(loo_cohort$MAE, na.rm = TRUE),
  range_MAE = max_MAE - min_MAE,
  sd_MAE = sd(loo_cohort$MAE, na.rm = TRUE)
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

full_MAE <- coef(brain_MAE_mv_global)[1]

png(
  filename = file.path(results_dir_brain, "plot_MAE_loo.png"),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_cohort,
  aes(
    x = reorder(cohort_removed, MAE),
    y = MAE
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = ci.lb, ymax = ci.ub),
    width = 0.15
  ) +
  geom_hline(
    yintercept = full_MAE,
    linetype = "dashed"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "Pooled MAE",
    title = "Leave-one-cohort-out sensitivity of pooled brain-age MAE"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- log ---

sink(log_file)

cat("Brain MAE multilevel meta-analysis\n")
cat("==================================\n\n")
cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_brain, "\n\n")

cat("Model summary\n")
cat("=============\n")
print(summary(brain_MAE_mv_global))

cat("\n\nPooled MAE estimate\n")
cat("===================\n")
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
print(loo_summary, n = Inf, width = Inf)

cat("\n\nDescriptives\n")
cat("============\n")
print(perf_brain_sub_desc, n = Inf, width = Inf)

sink()

###### brain MAE modelwise ------

run_id_modelwise <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_brain_modelwise <- file.path(
  results_dir,
  paste0("perf_brain_MAE_noDunedin_modelwise_", run_id_modelwise)
)

dir.create(results_dir_brain_modelwise, recursive = TRUE, showWarnings = FALSE)

log_file_modelwise <- file.path(
  results_dir_brain_modelwise,
  "analysis_log_modelwise.txt"
)

# --- meta ---

brain_MAE_mv_modelwise <- rma.mv(
  yi = MAE,
  V  = MAE_SE_boot^2,
  mods = ~ model - 1,
  random = list(
    ~ 1 | cohort/timepoint
  ),
  data = perf_brain_sub,
  method = "REML"
)

summary(brain_MAE_mv_modelwise)

# --- heterogeneity estimates ---

W <- diag(1 / brain_MAE_mv_modelwise$vi)
X <- model.matrix(brain_MAE_mv_modelwise)

P <- W - W %*% X %*%
  solve(t(X) %*% W %*% X) %*%
  t(X) %*% W

I2_total <- 100 * sum(brain_MAE_mv_modelwise$sigma2) /
  (
    sum(brain_MAE_mv_modelwise$sigma2) +
      (brain_MAE_mv_modelwise$k - brain_MAE_mv_modelwise$p) / sum(diag(P))
  )

I2_components <- 100 * brain_MAE_mv_modelwise$sigma2 /
  (
    sum(brain_MAE_mv_modelwise$sigma2) +
      (brain_MAE_mv_modelwise$k - brain_MAE_mv_modelwise$p) / sum(diag(P))
  )

names(I2_components) <- c(
  "cohort",
  "timepoint_within_cohort"
)

# --- extract modelwise estimates ---

brain_MAE_est <- tibble(
  term = rownames(brain_MAE_mv_modelwise$b),
  MAE = as.numeric(brain_MAE_mv_modelwise$b),
  se = brain_MAE_mv_modelwise$se,
  ci.lb = brain_MAE_mv_modelwise$ci.lb,
  ci.ub = brain_MAE_mv_modelwise$ci.ub,
  pval = brain_MAE_mv_modelwise$pval
) %>%
  mutate(model = str_remove(term, "^model")) %>%
  select(-term) %>%
  arrange(MAE)

X <- model.matrix(brain_MAE_mv_modelwise)

model_cols <- colnames(X)[
  str_detect(colnames(X), "^model") &
    !str_detect(colnames(X), ":")
]

effect_sizes_per_model <- tibble(
  model = str_remove(model_cols, "^model"),
  k = colSums(X[, model_cols, drop = FALSE] != 0)
)

brain_MAE_est <- brain_MAE_est %>%
  left_join(effect_sizes_per_model, by = "model")

brain_MAE_est

write.csv(
  brain_MAE_est,
  file.path(results_dir_brain_modelwise, "modelwise_MAE_estimates.csv"),
  row.names = FALSE
)

saveRDS(
  brain_MAE_mv_modelwise,
  file.path(results_dir_brain_modelwise, "model_brain_MAE_mv_modelwise.rds")
)

capture.output(
  summary(brain_MAE_mv_modelwise),
  file = file.path(results_dir_brain_modelwise, "modelwise_model_summary.txt")
)

# --- LOO cohort modelwise ---

cohorts <- unique(perf_brain_sub$cohort)

loo_cohort_modelwise <- lapply(cohorts, function(cc) {
  
  tmp <- subset(perf_brain_sub, cohort != cc)
  
  fit <- rma.mv(
    yi = MAE,
    V  = MAE_SE_boot^2,
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
    MAE = as.numeric(fit$b),
    se = fit$se,
    ci.lb = fit$ci.lb,
    ci.ub = fit$ci.ub,
    pval = fit$pval
  ) %>%
    mutate(model = str_remove(term, "^model")) %>%
    select(cohort_removed, model, MAE, se, ci.lb, ci.ub, pval)
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
    min_MAE = min(MAE, na.rm = TRUE),
    max_MAE = max(MAE, na.rm = TRUE),
    range_MAE = max_MAE - min_MAE,
    sd_MAE = sd(MAE, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  modelwise_loo_summary,
  file.path(results_dir_brain_modelwise, "modelwise_leave_one_cohort_out_summary.csv"),
  row.names = FALSE
)

# --- optional LOO plot by model ---

png(
  filename = file.path(results_dir_brain_modelwise, "plot_modelwise_MAE_loo.png"),
  width = 2200,
  height = 1400,
  res = 300
)

ggplot(
  loo_cohort_modelwise,
  aes(
    x = reorder(cohort_removed, MAE),
    y = MAE
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
    y = "Pooled MAE",
    title = "Modelwise leave-one-cohort-out sensitivity of brain-age MAE"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- log ---

sink(log_file_modelwise)

cat("Brain MAE modelwise multilevel meta-analysis\n")
cat("===========================================\n\n")
cat("Run ID:", run_id_modelwise, "\n")
cat("Output directory:", results_dir_brain_modelwise, "\n\n")

cat("Modelwise model summary\n")
cat("=======================\n")
print(summary(brain_MAE_mv_modelwise))

cat("\n\nModelwise estimates\n")
cat("===================\n")
print(brain_MAE_est, n = Inf, width = Inf)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
cat("Total I2 (%):\n")
print(I2_total)

cat("\nVariance decomposition (% of total heterogeneity):\n")
print(I2_components)

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

###### brain MAE global age moderator (not run) ------

# --- set up folder ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_brain_global_ageint <- file.path(
  results_dir,
  paste0("perf_brain_MAE_noDunedin_global_ageint_", run_id)
)

dir.create(results_dir_brain_global_ageint, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  results_dir_brain_global_ageint,
  "analysis_log_global_ageint.txt"
)

# --- prep ---

perf_brain_sub <- perf_brain_sub %>%
  mutate(
    mean_age_c = mean_age - mean(mean_age, na.rm = TRUE)
  )

# --- meta ---

brain_MAE_mv_global_ageint <- rma.mv(
  yi = MAE,
  V  = MAE_SE_boot^2,
  mods = ~ mean_age_c,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = perf_brain_sub,
  method = "REML"
)

summary(brain_MAE_mv_global_ageint)

# --- heterogeneity estimates ---

W <- diag(1 / brain_MAE_mv_global_ageint$vi)
X <- model.matrix(brain_MAE_mv_global_ageint)

P <- W - W %*% X %*%
  solve(t(X) %*% W %*% X) %*%
  t(X) %*% W

I2_total <- 100 * sum(brain_MAE_mv_global_ageint$sigma2) /
  (
    sum(brain_MAE_mv_global_ageint$sigma2) +
      (brain_MAE_mv_global_ageint$k - brain_MAE_mv_global_ageint$p) /
      sum(diag(P))
  )

I2_components <- 100 * brain_MAE_mv_global_ageint$sigma2 /
  (
    sum(brain_MAE_mv_global_ageint$sigma2) +
      (brain_MAE_mv_global_ageint$k - brain_MAE_mv_global_ageint$p) /
      sum(diag(P))
  )

names(I2_components) <- c(
  "cohort",
  "timepoint_within_cohort",
  "model"
)

# --- extract estimates ---

age_model_est <- tibble(
  term = names(coef(brain_MAE_mv_global_ageint)),
  estimate = as.numeric(coef(brain_MAE_mv_global_ageint)),
  se = brain_MAE_mv_global_ageint$se,
  ci.lb = brain_MAE_mv_global_ageint$ci.lb,
  ci.ub = brain_MAE_mv_global_ageint$ci.ub,
  pval = brain_MAE_mv_global_ageint$pval
)

write.csv(
  age_model_est,
  file.path(results_dir_brain_global_ageint, "age_moderator_results.csv"),
  row.names = FALSE
)

saveRDS(
  brain_MAE_mv_global_ageint,
  file.path(results_dir_brain_global_ageint, "model_brain_MAE_mv_global_ageint.rds")
)

capture.output(
  summary(brain_MAE_mv_global_ageint),
  file = file.path(results_dir_brain_global_ageint, "model_summary.txt")
)

# --- prediction grid ---

mean_age_mean <- mean(perf_brain_sub$mean_age, na.rm = TRUE)

newdat <- data.frame(
  mean_age = seq(
    min(perf_brain_sub$mean_age, na.rm = TRUE),
    max(perf_brain_sub$mean_age, na.rm = TRUE),
    length.out = 100
  )
) %>%
  mutate(
    mean_age_c = mean_age - mean_age_mean
  )

pred_obj <- predict(
  brain_MAE_mv_global_ageint,
  newmods = newdat$mean_age_c
)

newdat <- newdat %>%
  mutate(
    pred_MAE = pred_obj$pred,
    ci.lb = pred_obj$ci.lb,
    ci.ub = pred_obj$ci.ub
  )

cohort_time_means <- perf_brain_sub %>%
  group_by(cohort, timepoint) %>%
  summarise(
    mean_age = mean(mean_age, na.rm = TRUE),
    mean_MAE = mean(MAE, na.rm = TRUE),
    .groups = "drop"
  )

png(
  filename = file.path(
    results_dir_brain_global_ageint,
    "plot_meanage_MAE_meta_regression.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot() +
  geom_point(
    data = perf_brain_sub,
    aes(
      x = mean_age,
      y = MAE,
      colour = cohort
    ),
    alpha = .25,
    size = 2
  ) +
  geom_point(
    data = cohort_time_means,
    aes(
      x = mean_age,
      y = mean_MAE,
      colour = cohort
    ),
    size = 4
  ) +
  geom_line(
    data = newdat,
    aes(
      x = mean_age,
      y = pred_MAE
    ),
    linewidth = 1.2,
    colour = "black"
  ) +
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
    title = "Meta-regression of brain-age MAE on mean age",
    x = "Mean sample age",
    y = "Predicted MAE",
    colour = "Cohort"
  ) +
  theme_minimal(base_size = 16)

dev.off()

# --- LOO cohort ---

cohorts <- unique(perf_brain_sub$cohort)

loo_age <- lapply(cohorts, function(cc) {
  
  tmp <- perf_brain_sub %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = MAE,
    V  = MAE_SE_boot^2,
    mods = ~ mean_age_c,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model
    ),
    data = tmp,
    method = "REML"
  )
  
  est <- coef(fit)[1]
  slope <- coef(fit)[2]
  
  est_se <- sqrt(vcov(fit)[1, 1])
  slope_se <- sqrt(vcov(fit)[2, 2])
  
  data.frame(
    cohort_removed = cc,
    
    intercept = est,
    intercept_se = est_se,
    intercept_ci.lb = est - 1.96 * est_se,
    intercept_ci.ub = est + 1.96 * est_se,
    
    slope = slope,
    slope_se = slope_se,
    slope_ci.lb = slope - 1.96 * slope_se,
    slope_ci.ub = slope + 1.96 * slope_se
  )
})

loo_age_df <- do.call(rbind, loo_age) %>%
  as_tibble()

loo_age_summary <- tibble(
  min_intercept = min(loo_age_df$intercept, na.rm = TRUE),
  max_intercept = max(loo_age_df$intercept, na.rm = TRUE),
  range_intercept = max_intercept - min_intercept,
  sd_intercept = sd(loo_age_df$intercept, na.rm = TRUE),
  
  min_slope = min(loo_age_df$slope, na.rm = TRUE),
  max_slope = max(loo_age_df$slope, na.rm = TRUE),
  range_slope = max_slope - min_slope,
  sd_slope = sd(loo_age_df$slope, na.rm = TRUE)
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

full_slope <- coef(brain_MAE_mv_global_ageint)[2]

png(
  filename = file.path(
    results_dir_brain_global_ageint,
    "plot_meanage_MAE_loo_slope.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_age_df,
  aes(
    x = reorder(cohort_removed, slope),
    y = slope
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = slope_ci.lb,
      ymax = slope_ci.ub
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept = full_slope,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 0,
    colour = "grey60"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "MAE change per year",
    title = "Leave-one-cohort-out sensitivity of age effect"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- log ---

sink(log_file)

cat("Brain MAE multilevel meta-analysis - global age moderator\n")
cat("=========================================================\n\n")

cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_brain_global_ageint, "\n\n")

cat("Model summary\n")
cat("=============\n")
print(summary(brain_MAE_mv_global_ageint))

cat("\n\nAge moderator results\n")
cat("=====================\n")
print(age_model_est, n = Inf, width = Inf)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")

cat("Total I2 (%):\n")
print(I2_total)

cat("\nVariance decomposition (% of total heterogeneity):\n")
print(I2_components)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_age_df, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_age_summary, n = Inf, width = Inf)

sink()

###### brain MAE modelwise age interactions (not run) ------

# --- set up folder ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_brain_modelwise_ageint <- file.path(
  results_dir,
  paste0("perf_brain_MAE_modelwise_ageint_", run_id)
)

dir.create(results_dir_brain_modelwise_ageint, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  results_dir_brain_modelwise_ageint,
  "analysis_log_modelwise_ageint.txt"
)

# --- prep ---

perf_brain_sub <- perf_brain_sub %>%
  mutate(
    mean_age_c = mean_age - mean(mean_age, na.rm = TRUE)
  )

# --- meta ---

brain_MAE_mv_modelwise_ageint <- rma.mv(
  yi = MAE,
  V  = MAE_SE_boot^2,
  mods = ~ model + model:mean_age_c - 1,
  random = list(
    ~ 1 | cohort/timepoint
  ),
  data = perf_brain_sub,
  method = "REML"
)

summary(brain_MAE_mv_modelwise_ageint)

# --- heterogeneity estimates ---

W <- diag(1 / brain_MAE_mv_modelwise_ageint$vi)
X <- model.matrix(brain_MAE_mv_modelwise_ageint)

P <- W - W %*% X %*%
  solve(t(X) %*% W %*% X) %*%
  t(X) %*% W

I2_total <- 100 * sum(brain_MAE_mv_modelwise_ageint$sigma2) /
  (
    sum(brain_MAE_mv_modelwise_ageint$sigma2) +
      (brain_MAE_mv_modelwise_ageint$k - brain_MAE_mv_modelwise_ageint$p) /
      sum(diag(P))
  )

I2_components <- 100 * brain_MAE_mv_modelwise_ageint$sigma2 /
  (
    sum(brain_MAE_mv_modelwise_ageint$sigma2) +
      (brain_MAE_mv_modelwise_ageint$k - brain_MAE_mv_modelwise_ageint$p) /
      sum(diag(P))
  )

names(I2_components) <- c(
  "cohort",
  "timepoint_within_cohort"
)

# --- compare with common slope model ---

anova(
  brain_MAE_mv_global_ageint,
  brain_MAE_mv_modelwise_ageint,
  refit = TRUE
)

# --- extract estimates ---

modelwise_age_est <- tibble(
  term = rownames(brain_MAE_mv_modelwise_ageint$b),
  estimate = as.numeric(brain_MAE_mv_modelwise_ageint$b),
  se = brain_MAE_mv_modelwise_ageint$se,
  ci.lb = brain_MAE_mv_modelwise_ageint$ci.lb,
  ci.ub = brain_MAE_mv_modelwise_ageint$ci.ub,
  pval = brain_MAE_mv_modelwise_ageint$pval
) %>%
  mutate(
    type = if_else(
      str_detect(term, "mean_age_c"),
      "age_slope",
      "intercept"
    ),
    model = term %>%
      str_remove("^model") %>%
      str_remove(":mean_age_c") %>%
      str_remove("\\.mean_age_c")
  )

X <- model.matrix(brain_MAE_mv_modelwise_ageint)

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
    results_dir_brain_modelwise_ageint,
    "modelwise_age_interaction_results.csv"
  ),
  row.names = FALSE
)

saveRDS(
  brain_MAE_mv_modelwise_ageint,
  file.path(
    results_dir_brain_modelwise_ageint,
    "model_brain_MAE_mv_modelwise_ageint.rds"
  )
)

capture.output(
  summary(brain_MAE_mv_modelwise_ageint),
  file = file.path(
    results_dir_brain_modelwise_ageint,
    "model_summary.txt"
  )
)

# --- prediction grid ---

mean_age_mean <- mean(perf_brain_sub$mean_age, na.rm = TRUE)

newdat <- expand.grid(
  mean_age_c = seq(
    min(perf_brain_sub$mean_age_c, na.rm = TRUE),
    max(perf_brain_sub$mean_age_c, na.rm = TRUE),
    length.out = 100
  ),
  model = unique(perf_brain_sub$model)
)

pred_obj <- predict(
  brain_MAE_mv_modelwise_ageint,
  newmods = model.matrix(
    ~ model + model:mean_age_c - 1,
    newdat
  )
)

newdat <- newdat %>%
  mutate(
    mean_age = mean_age_c + mean_age_mean,
    pred_MAE = pred_obj$pred,
    ci.lb = pred_obj$ci.lb,
    ci.ub = pred_obj$ci.ub
  )

png(
  filename = file.path(
    results_dir_brain_modelwise_ageint,
    "plot_brain_meanage_MAE_modelwise.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  newdat,
  aes(mean_age, pred_MAE, colour = model)
) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(
      ymin = ci.lb,
      ymax = ci.ub,
      fill = model
    ),
    alpha = .12,
    colour = NA
  ) +
  labs(
    x = "Mean sample age",
    y = "Predicted MAE",
    title = "Model-specific age effects"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- LOO cohort ---

cohorts <- unique(perf_brain_sub$cohort)

loo_modelwise_age <- lapply(cohorts, function(cc) {
  
  tmp <- perf_brain_sub %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = MAE,
    V  = MAE_SE_boot^2,
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
    estimate = as.numeric(fit$b),
    se = fit$se,
    ci.lb = fit$ci.lb,
    ci.ub = fit$ci.ub,
    pval = fit$pval
  ) %>%
    mutate(
      type = if_else(
        str_detect(term, "mean_age_c"),
        "age_slope",
        "intercept"
      ),
      model = term %>%
        str_remove("^model") %>%
        str_remove(":mean_age_c") %>%
        str_remove("\\.mean_age_c")
    )
}) %>%
  bind_rows()

write.csv(
  loo_modelwise_age,
  file.path(
    results_dir_brain_modelwise_ageint,
    "modelwise_leave_one_cohort_out_age_interactions.csv"
  ),
  row.names = FALSE
)

loo_modelwise_age_summary <- loo_modelwise_age %>%
  filter(type == "age_slope") %>%
  group_by(model) %>%
  summarise(
    min_slope = min(estimate, na.rm = TRUE),
    max_slope = max(estimate, na.rm = TRUE),
    range_slope = max_slope - min_slope,
    sd_slope = sd(estimate, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  loo_modelwise_age_summary,
  file.path(
    results_dir_brain_modelwise_ageint,
    "modelwise_leave_one_cohort_out_age_interactions_summary.csv"
  ),
  row.names = FALSE
)

# --- log ---

sink(log_file)

cat("Brain MAE multilevel meta-analysis - modelwise age interactions\n")
cat("================================================================\n\n")

cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_brain_modelwise_ageint, "\n\n")

cat("Model summary\n")
cat("=============\n")
print(summary(brain_MAE_mv_modelwise_ageint))

cat("\n\nModelwise age interaction results\n")
cat("=================================\n")
print(modelwise_age_est, n = Inf, width = Inf)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")

cat("Total I2 (%):\n")
print(I2_total)

cat("\nVariance decomposition (% of total heterogeneity):\n")
print(I2_components)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_modelwise_age, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("===========================================\n")
print(loo_modelwise_age_summary, n = Inf, width = Inf)

sink()