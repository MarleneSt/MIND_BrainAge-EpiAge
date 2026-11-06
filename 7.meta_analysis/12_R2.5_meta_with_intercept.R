#### META-ANALYSES SCRIPT
#### MIND BrainAge~EpiAge
#### last updated: 10/07/2026

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

#### some helpers ----

# Helper: I2 for multilevel metafor models
# code from Viechtbauer's metafor website: https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate

# note: I2 was not included from the beginning so there can be some inconsistencies in how it's applied throughout


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

#### load data ----
# load data if the above script is not run in full (recommended)

perf_brain <- readRDS(file.path(results_dir, "perf_brain.rds"))

perf_epi <- readRDS(file.path(results_dir, "perf_epi.rds"))

mod2_eff2 <- readRDS(file.path(results_dir, "mod2_eff2.rds"))

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

mod2_eff2 <- mod2_eff2 %>%
  mutate(
    cohort = factor(cohort),
    timepoint = factor(timepoint),
    model_combi = factor(model_combi)
  )


# [REPRODUCIBILITY] GenR dedup step removed: the duplicates came from an identical
# copy of the GenR association files (assoc_old folder), now removed at the source.
nrow(mod2_eff2)   # check: should be 2448


###### assoc modelwise ------

# --- set up folder ---

run_id_modelwise <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_assoc_modelwise <- file.path("results/revision/R2.5", paste0("assoc_model2_PAR_modelwise_", run_id_modelwise)) # [REPRODUCIBILITY] relative path
dir.create(results_dir_assoc_modelwise, recursive = TRUE, showWarnings = FALSE)

# results_dir_assoc_modelwise <- file.path(
#   results_dir,
#   paste0("assoc_model2_PAR_modelwise_", run_id_modelwise)
# )
# 
# dir.create(results_dir_assoc_modelwise, recursive = TRUE, showWarnings = FALSE)
# 
# log_file_modelwise <- file.path(results_dir_assoc_modelwise, "analysis_log_modelwise.txt")



###### R2.5
# --- meta ---

# original Marlene ran without the intercept
mod2_res_mv_modelwise <- rma.mv(
  yi   = RLM_Estimate_scaled,
  V    = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ model_combi - 1,
  random = ~ 1 | cohort/timepoint,
  data   = mod2_eff2,
  method = "REML"
)

# updated one with the intercept
mod2_res_mv_modelwise_with_intercept <- rma.mv(
  yi   = RLM_Estimate_scaled,
  V    = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ model_combi,
  random = ~ 1 | cohort/timepoint,
  data   = mod2_eff2,
  method = "REML"
)

# original without intercept
summary(mod2_res_mv_modelwise)
range(coef(mod2_res_mv_modelwise))

# updated with intercept
summary(mod2_res_mv_modelwise_with_intercept)
range(coef(mod2_res_mv_modelwise_with_intercept))

# save original (without intercept)
saveRDS(
  mod2_res_mv_modelwise,
  file.path(results_dir_assoc_modelwise, "model_assoc_mv_modelwise_original_R2.5.rds")
)

capture.output(
  summary(mod2_res_mv_modelwise),
  file = file.path(results_dir_assoc_modelwise, "modelwise_model_summary_original_R2.5.txt")
)

# save updated (with intercept)
saveRDS(
  mod2_res_mv_modelwise_with_intercept,
  file.path(results_dir_assoc_modelwise, "model_assoc_mv_modelwise_with_intercept_R2.5.rds")
)

capture.output(
  summary(mod2_res_mv_modelwise_with_intercept),
  file = file.path(results_dir_assoc_modelwise, "modelwise_model_summary_with_intercept_R2.5.txt")
)



# --- sensitivity: correlated sampling errors within cohort/timepoint ---> NOT USED FOR NOW

mod2_eff2 <- mod2_eff2 %>%
  mutate(vi = RLM_SE_coeftest_HC3_scaled^2,
         cl = interaction(cohort, timepoint, drop = TRUE))

sens_rho <- map_dfr(c(0.3, 0.5, 0.7), function(rho) {
  V <- vcalc(vi, cluster = cl, obs = model_combi, rho = rho, data = mod2_eff2)
  fit <- rma.mv(yi = RLM_Estimate_scaled, V = V,
                mods = ~ model_combi,
                random = ~ 1 | cohort/timepoint,
                data = mod2_eff2, method = "REML")
  tibble(rho = rho, QM = fit$QM, df = fit$QMdf[1], p = fit$QMp)
})

sens_rho

write_csv(sens_rho, file.path(results_dir_assoc_modelwise, "modelwise_QM_sensitivity_rho_R2.5.csv"))


## I^2 for the text edits after removing duplicates

# --- Overall (global) model: average association across all pairs ---
# Same structure as the original script (line 2995): random effects for
# timepoints nested in cohorts, plus a separate random effect for model pair.
# Re-run on the deduplicated data (GenR duplicates removed).
mod2_res_mv_global <- rma.mv(
  yi = RLM_Estimate_scaled,
  V  = RLM_SE_coeftest_HC3_scaled^2,
  random = list(~ 1 | cohort/timepoint, ~ 1 | model_combi),
  data = mod2_eff2,
  method = "REML"
)

# --- Heterogeneity: I² ---
# Share of variation that is real (not sampling error), in total and
# split into cohort, timepoint and model pair. Replaces I² = 14.64%.
calc_I2_mv(mod2_res_mv_global, c("cohort", "timepoint", "model_combi"))

# --- Heterogeneity: Q test p-value ---
# Is there more variation than expected by chance? Replaces p = 0.578.
mod2_res_mv_global$QEp

# --- Overall association with cluster-robust inference ---
# Adjusts SEs for the fact that estimates from the same cohort share people
# (12 cohorts -> t(11)). Gives β, 95% CI, t(11) and p for the first sentence.
robust(mod2_res_mv_global, cluster = mod2_eff2$cohort)


###### omnibus tests of differences between levels (main text) ------

# The original modelwise/generation models were fit without an intercept
# (mods = ~ X - 1), so their QM(k) tested "all levels = 0". Refitting with
# an intercept (mods = ~ X) makes QM test "all levels equal" (k - 1 df).
# Per-level estimates for tables/figures still come from the original
# no-intercept models (unchanged for performance).

# --- brain wMAE between models (was QM(7) = 7699.08; original script line 1071) ---
perf_brain_sub <- perf_brain %>% filter(!model %in% c("DunedinPACNI")) %>% droplevels()

brain_wMAE_mv_modelwise_with_intercept <- rma.mv(
  yi   = wMAE_test,
  V    = wMAE_SE_boot^2,
  mods = ~ model,
  random = ~ 1 | cohort/timepoint,
  data   = perf_brain_sub,
  method = "REML"
)

# --- epi wMAE between models (was QM(18) = 225166.02; original script lines 1812-1833, 2075) ---
perf_epi_sub <- perf_epi %>%
  filter(!model %in% c("DunedinPACE", "DNAmTL")) %>%
  droplevels() %>%
  mutate(
    model = recode(model, "Bohlin_conv" = "Bohlin", "Knight_conv" = "Knight", "EPIC_conv" = "EPIC"),
    log_wMAE = log(wMAE_test),
    log_wMAE_var = (wMAE_SE_boot / wMAE_test)^2
  )

epi_wMAE_mv_modelwise_with_intercept <- rma.mv(
  yi   = log_wMAE,
  V    = log_wMAE_var,
  mods = ~ model,
  random = ~ 1 | cohort/timepoint,
  data   = perf_epi_sub,
  method = "REML"
)

# --- assoc between epi clock generations (was QM(4) = 20.11; original script line 8243) ---
# uses the deduplicated mod2_eff2
mod2_eff2 <- mod2_eff2 %>% mutate(generation_epi = factor(generation_epi))

mod2_res_mv_generation_epi_with_intercept <- rma.mv(
  yi   = RLM_Estimate_scaled,
  V    = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ generation_epi,
  random = list(~ 1 | cohort/timepoint, ~ 1 | model_combi),
  data   = mod2_eff2,
  method = "REML"
)

# --- collect and save ---
qm_row <- function(fit) tibble(QM = fit$QM, df = fit$QMdf[1], p = fit$QMp)

omnibus_diff <- bind_rows(
  assoc_pairs       = qm_row(mod2_res_mv_modelwise_with_intercept),   # QM(135) = 192.0
  brain_wMAE_models = qm_row(brain_wMAE_mv_modelwise_with_intercept), # df should be 6
  epi_wMAE_models   = qm_row(epi_wMAE_mv_modelwise_with_intercept),   # df should be 17
  assoc_generation  = qm_row(mod2_res_mv_generation_epi_with_intercept), # df should be 3
  .id = "test"
)

omnibus_diff

write_csv(omnibus_diff, file.path(results_dir_assoc_modelwise, "omnibus_differences_main_text_R2.5.csv"))

