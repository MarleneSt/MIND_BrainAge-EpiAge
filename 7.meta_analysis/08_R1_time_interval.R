#### R1: time interval between DNAm and MRI assessment (assoc moderator) ----
# [REPRODUCIBILITY] header added so this script runs on its own 
# Library, helper and data loading copied from 01_main_meta_wMAE.R.

#### [REPRODUCIBILITY] library ----
library(tidyverse)
library(metafor)

#### [REPRODUCIBILITY] define important folders ----
results_dir <- "results"

#### [REPRODUCIBILITY] helper (copied from 01_main_meta_wMAE.R) ----
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

#### [REPRODUCIBILITY] load data (copied from 01_main_meta_wMAE.R) ----
mod2_eff2 <- readRDS(file.path(results_dir, "mod2_eff2.rds"))

mod2_eff2 <- mod2_eff2 %>%
  mutate(
    cohort = factor(cohort),
    timepoint = factor(timepoint),
    model_combi = factor(model_combi)
  )


###### assoc age diff global ------
# --- set up ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_assoc_global_agediff <- file.path(results_dir, paste0("assoc_model2_PAR_global_agediff_", run_id))

dir.create(results_dir_assoc_global_agediff, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_assoc_global_agediff, "analysis_log_global_agediff.txt")

# --- data prep ---

# center age difference
mod2_eff2 <- mod2_eff2 %>%
  mutate(
    mean_age_diff_c = age_diff_brain_minus_epi_mean - mean(age_diff_brain_minus_epi_mean, na.rm = TRUE)
  )

# save age diff by cohort and timepoint 


age_diff_table <- mod2_eff2 %>%
  group_by(cohort, timepoint) %>%
  summarise(
    n = sum(!is.na(age_diff_brain_minus_epi_mean)),
    mean = round(mean(age_diff_brain_minus_epi_mean, na.rm = TRUE), 2),
    sd = round(sd(age_diff_brain_minus_epi_mean, na.rm = TRUE), 2),
    median = round(median(age_diff_brain_minus_epi_mean, na.rm = TRUE), 2),
    min = round(min(age_diff_brain_minus_epi_mean, na.rm = TRUE), 2),
    max = round(max(age_diff_brain_minus_epi_mean, na.rm = TRUE), 2),
    .groups = "drop"
  )

write.csv(
  age_diff_table,
  file.path(results_dir_assoc_global_agediff, "age_diff_by_cohort_timepoint.csv"),
  row.names = FALSE
)

# --- meta ---

mod2_res_mv_global_agediff <- rma.mv(
  yi = RLM_Estimate_scaled,
  V  = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ mean_age_diff_c,
  random = list(
    ~ 1 | cohort/timepoint, 
    ~ 1 | model_combi  
  ),
  data   = mod2_eff2,
  method = "REML",
  verbose = TRUE
)

summary(mod2_res_mv_global_agediff)

saveRDS(
  mod2_res_mv_global_agediff,
  file.path(results_dir_assoc_global_agediff, "model_assoc_mv_global_agediff.rds")
)

capture.output(
  summary(mod2_res_mv_global_agediff),
  file = file.path(results_dir_assoc_global_agediff, "global_agediff_model_summary.txt")
)

# --- heterogeneity ---
I2 <- calc_I2_mv(
  mod2_res_mv_global_agediff,
  component_names = c("cohort", "timepoint_within_cohort", "model_combi")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_assoc_global_agediff, "heterogeneity_I2.csv"), row.names = FALSE)

# --- extract estimates ---
assoc_global_agediff_est <- tibble(
  term = names(coef(mod2_res_mv_global_agediff)),
  estimate = as.numeric(coef(mod2_res_mv_global_agediff)),
  se = mod2_res_mv_global_agediff$se,
  ci.lb = mod2_res_mv_global_agediff$ci.lb,
  ci.ub = mod2_res_mv_global_agediff$ci.ub,
  pval = mod2_res_mv_global_agediff$pval
)

write.csv(
  assoc_global_agediff_est,
  file.path(results_dir_assoc_global_agediff, "age_moderator_results.csv"),
  row.names = FALSE
)

# --- plotting (should be adapted) ---

# plotting 

newdat <- data.frame(
  mean_age_diff_c = seq(
    min(mod2_eff2$mean_age_diff_c, na.rm = TRUE),
    max(mod2_eff2$mean_age_diff_c, na.rm = TRUE),
    length.out = 100
  )
)

pred <- predict(mod2_res_mv_global_agediff, newmods = newdat$mean_age_diff_c)

newdat$pred  <- pred$pred
newdat$ci.lb <- pred$ci.lb
newdat$ci.ub <- pred$ci.ub

newdat$mean_age_brain <- newdat$mean_age_diff_c + mean(mod2_eff2$mean_age_brain, na.rm = TRUE)

png(
  filename = file.path(results_dir_assoc_global_agediff, "association_global_meanage_meta_regression.png"),
  width = 1800, height = 1200, res = 300
)

ggplot() +
  geom_point(
    data = mod2_eff2,
    aes(x = mean_age_brain, y = RLM_Estimate_scaled, colour = cohort),
    alpha = .5
  ) +
  geom_line(
    data = newdat,
    aes(x = mean_age_brain, y = pred),
    colour = "black",
    linewidth = 1
  ) +
  geom_ribbon(
    data = newdat,
    aes(x = mean_age_brain, ymin = ci.lb, ymax = ci.ub),
    alpha = .15
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Brain-PAR - Epi-PAR association by mean age",
    x = "Mean age",
    y = "Standardised beta (RLM)"
  )

dev.off()

# --- LOO ---

# --- LOO cohort ---

cohorts <- unique(mod2_eff2$cohort)

loo_assoc_ageint <- lapply(cohorts, function(cc) {
  
  tmp <- subset(mod2_eff2, cohort != cc)
  
  fit <- rma.mv(
    yi = RLM_Estimate_scaled,
    V  = RLM_SE_coeftest_HC3_scaled^2,
    mods = ~ mean_age_diff_c,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model_combi
    ),
    data = tmp,
    method = "REML"
  )
  
  tibble(
    cohort_removed = cc,
    
    intercept = coef(fit)[1],
    intercept_se = sqrt(vcov(fit)[1,1]),
    
    age_slope = coef(fit)[2],
    age_slope_se = sqrt(vcov(fit)[2,2])
  )
})

loo_assoc_ageint <- bind_rows(loo_assoc_ageint)


loo_assoc_agediff_summary <- loo_assoc_ageint %>%
  summarise(
    min_intercept = min(intercept, na.rm = TRUE),
    max_intercept = max(intercept, na.rm = TRUE),
    range_intercept = max_intercept - min_intercept,
    sd_intercept = sd(intercept, na.rm = TRUE),
    
    min_age_slope = min(age_slope, na.rm = TRUE),
    max_age_slope = max(age_slope, na.rm = TRUE),
    range_age_slope = max_age_slope - min_age_slope,
    sd_age_slope = sd(age_slope, na.rm = TRUE)
  )

write.csv(
  loo_assoc_ageint,
  file.path(
    results_dir_assoc_global_agediff,
    "leave_one_cohort_out_age_interaction.csv"
  ),
  row.names = FALSE
)

png(
  filename = file.path(
    results_dir_assoc_global_agediff,
    "plot_assoc_meanage_brain_loo.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_assoc_ageint,
  aes(
    x = reorder(cohort_removed, age_slope),
    y = age_slope
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = age_slope - 1.96 * age_slope_se,
      ymax = age_slope + 1.96 * age_slope_se
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept = coef(mod2_res_mv_global_agediff)[2],
    linetype = "dashed"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "Estimated brain age slope",
    title = "Leave-one-cohort-out sensitivity of brain age effect on associations"
  ) +
  theme_minimal()

dev.off()


# --- full analysis log ---

sink(log_file)

cat("Association global age interaction meta-analysis: brain mean age\n")
cat("===============================================================\n\n")

cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_assoc_global_agediff, "\n\n")

cat("Model summary\n")
cat("=============\n")
print(summary(mod2_res_mv_global_agediff))

cat("\n\nGlobal age interaction estimates\n")
cat("================================\n")
print(assoc_global_agediff_est)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_assoc_ageint, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_assoc_agediff_summary, n = Inf, width = Inf)

cat("\n\nSession info\n")
cat("============\n")
print(sessionInfo())

sink()
