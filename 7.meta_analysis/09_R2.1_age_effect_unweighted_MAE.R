###### R2.1: age effect with unweighted MAE, with test-sample age range as moderator ------
#
# Two models (performance data only; no GenR duplicates there):
#   (1) epi   unweighted MAE ~ age range + mean age   (coefficient of interest: mean age)
#   (2) brain unweighted MAE ~ age range + mean age   (coefficient of interest: mean age)
#
# Labels:
#   [ORIGINAL]  copied from the original scripts (file and line given), unchanged
#   [ORIGINAL, CHANGED]  copied from the original, with the change noted
#   [FROM tmp_agerange.R] copied from tmp_agerange.R (R2.1 folder)
#   [NEW]       added for R2.1


#### [NEW] libraries and folders ----
library(tidyverse)
library(metafor)

results_dir <- "results"  # perf_*.rds # [REPRODUCIBILITY] relative path
out_dir     <- "results/revision/R2.1" # [REPRODUCIBILITY] relative path
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE) # [REPRODUCIBILITY] create folder on a fresh run


#### [NEW] load data ----
perf_epi   <- readRDS(file.path(results_dir, "perf_epi.rds"))
perf_brain <- readRDS(file.path(results_dir, "perf_brain.rds"))


#### data prep ----
# [ORIGINAL] epi subset and model names (4.meta_unweightedMAE.R, lines 68-80)
perf_epi_sub   <- perf_epi   %>% filter(!model %in% c("DunedinPACE", "DNAmTL")) %>% droplevels()

perf_epi_sub <- perf_epi_sub %>%
  mutate(
    model = recode(
      model,
      "Bohlin_conv" = "Bohlin",
      "Knight_conv" = "Knight",
      "EPIC_conv" = "EPIC"
    )
  )

# [ORIGINAL] log unweighted MAE (4.meta_unweightedMAE.R, lines 84-89)
perf_epi_sub <- perf_epi_sub |>
  mutate(
    log_MAE = log(MAE),
    # Delta method: Var(log(X)) ≈ (SE_X / X)^2
    log_MAE_var = (MAE_SE_boot / MAE)^2
  )

# [ORIGINAL] centre age, then remove gestational clocks (4.meta_unweightedMAE.R, lines 633-642)
perf_epi_sub <- perf_epi_sub %>%
  mutate(
    mean_age_c = mean_age - mean(mean_age, na.rm = TRUE)
  )

# remove gestational clocks for age moderation
perf_epi_sub_nogest <- perf_epi_sub %>%
  filter(!model %in% c("Bohlin", "Knight", "EPIC")) %>%
  droplevels()

# [ORIGINAL] brain subset (4.meta_unweightedMAE.R, lines 1285-1287)
perf_brain_sub <- perf_brain %>%
  filter(!model %in% c("DunedinPACNI")) %>%
  droplevels()

# [ORIGINAL] centre age (4.meta_unweightedMAE.R, lines 1782-1785)
perf_brain_sub <- perf_brain_sub %>%
  mutate(
    mean_age_c = mean_age - mean(mean_age, na.rm = TRUE)
  )

# [FROM tmp_agerange.R] centred test-sample age range
perf_epi_sub_nogest <- perf_epi_sub_nogest %>% mutate(age_range_c = age_range - mean(age_range, na.rm = TRUE))
perf_brain_sub      <- perf_brain_sub      %>% mutate(age_range_c = age_range - mean(age_range, na.rm = TRUE))


#### meta-analyses ----

###### (1) Epigenetic clocks: unweighted MAE ~ age range + mean age ------
# outcome: log MAE (years) | data: perf_epi_sub_nogest | random: cohort/timepoint + model
# [FROM tmp_agerange.R, CHANGED] epi_MAE_mv_global with mods = ~ age_range_c (tmp_agerange.R, line 64).
# CHANGED: + mean_age_c added (coefficient of interest); data = perf_epi_sub_nogest
# (gestational clocks excluded, as in all original age models); object name.
epi_MAE_mv_global_age_agerange <- rma.mv(
  yi = log_MAE,
  V  = log_MAE_var,
  mods = ~ age_range_c + mean_age_c,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data   = perf_epi_sub_nogest,
  method = "REML",
  control = list(
    optimizer = "optim",
    optmethod = "BFGS",
    maxit = 10000
  )
)
summary(epi_MAE_mv_global_age_agerange)

###### (2) Brain age models: unweighted MAE ~ age range + mean age ------
# outcome: MAE (years, raw scale) | data: perf_brain_sub | random: cohort/timepoint + model
# [ORIGINAL, CHANGED] brain_MAE_mv_global_ageint (4.meta_unweightedMAE.R, line 1789).
# CHANGED: mods = ~ mean_age_c  ->  ~ age_range_c + mean_age_c; object name.
brain_MAE_mv_global_age_agerange <- rma.mv(
  yi = MAE,
  V  = MAE_SE_boot^2,
  mods = ~ age_range_c + mean_age_c,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = perf_brain_sub,
  method = "REML"
)
summary(brain_MAE_mv_global_age_agerange)


#### [NEW] results table ----
# coefficients: [2] = age range, [3] = mean age
# epi: log scale (% change per year = (exp(b) - 1) * 100); brain: raw scale (years)
res_1_age <- data.frame(
  model = "(1) epi MAE ~ age range + age", term = "mean age",
  k = epi_MAE_mv_global_age_agerange$k,
  b = epi_MAE_mv_global_age_agerange$b["mean_age_c", 1],
  ci.lb = epi_MAE_mv_global_age_agerange$ci.lb[3], ci.ub = epi_MAE_mv_global_age_agerange$ci.ub[3],
  p = epi_MAE_mv_global_age_agerange$pval[3]
)
res_1_range <- data.frame(
  model = "(1) epi MAE ~ age range + age", term = "age range",
  k = epi_MAE_mv_global_age_agerange$k,
  b = epi_MAE_mv_global_age_agerange$b["age_range_c", 1],
  ci.lb = epi_MAE_mv_global_age_agerange$ci.lb[2], ci.ub = epi_MAE_mv_global_age_agerange$ci.ub[2],
  p = epi_MAE_mv_global_age_agerange$pval[2]
)
res_2_age <- data.frame(
  model = "(2) brain MAE ~ age range + age", term = "mean age",
  k = brain_MAE_mv_global_age_agerange$k,
  b = brain_MAE_mv_global_age_agerange$b["mean_age_c", 1],
  ci.lb = brain_MAE_mv_global_age_agerange$ci.lb[3], ci.ub = brain_MAE_mv_global_age_agerange$ci.ub[3],
  p = brain_MAE_mv_global_age_agerange$pval[3]
)
res_2_range <- data.frame(
  model = "(2) brain MAE ~ age range + age", term = "age range",
  k = brain_MAE_mv_global_age_agerange$k,
  b = brain_MAE_mv_global_age_agerange$b["age_range_c", 1],
  ci.lb = brain_MAE_mv_global_age_agerange$ci.lb[2], ci.ub = brain_MAE_mv_global_age_agerange$ci.ub[2],
  p = brain_MAE_mv_global_age_agerange$pval[2]
)

age_effects <- rbind(res_1_age, res_1_range, res_2_age, res_2_range)

# [NEW] % change per year, back-transformed from log scale (as in SM Table S9);
# epi only (log scale), brain is on the raw scale (years per year)
is_epi <- startsWith(age_effects$model, "(1) epi")
age_effects$pct       <- ifelse(is_epi, (exp(age_effects$b) - 1) * 100, NA)
age_effects$pct_ci.lb <- ifelse(is_epi, (exp(age_effects$ci.lb) - 1) * 100, NA)
age_effects$pct_ci.ub <- ifelse(is_epi, (exp(age_effects$ci.ub) - 1) * 100, NA)

age_effects

write.csv(age_effects, file.path(out_dir, "R2.1_age_effects.csv"), row.names = FALSE)
