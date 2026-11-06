###### R2.1 (follow-up): per-clock age effects, MAEw vs unweighted MAE ------
#
# Four models (performance data only; no GenR duplicates there):
#   (1) epi   MAEw per clock ~ mean age                         (as in paper, unchanged)
#   (2) epi   unweighted MAE per clock ~ mean age + age range
#   (3) brain MAEw per clock ~ mean age                         (as in paper, unchanged)
#   (4) brain unweighted MAE per clock ~ mean age + age range
#
# Labels:
#   [ORIGINAL]  copied from the original scripts (file and line given), unchanged
#   [ORIGINAL, CHANGED]  copied from the original, with the change noted
#   [FROM R2.1] copied from R2.1_age_effect_unweighted_MAE_age_range.R
#   [NEW]       added for this script


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

# [ORIGINAL] log MAEw (1.main_meta_wMAE.R, lines 1828-1833)
perf_epi_sub <- perf_epi_sub |>
  mutate(
    log_wMAE = log(wMAE_test),
    # Delta method: Var(log(X)) ≈ (SE_X / X)^2
    log_wMAE_var = (wMAE_SE_boot / wMAE_test)^2
  )

# [ORIGINAL] remove gestational clocks for age moderation (1.main_meta_wMAE.R, line 2318)
perf_epi_sub_nogest <- perf_epi_sub %>% filter(!model %in% c("Bohlin", "Knight", "EPIC")) %>% droplevels()

# [ORIGINAL] centre age (1.main_meta_wMAE.R, lines 2702-2703; 4.meta_unweightedMAE.R, lines 999-1002)
perf_epi_sub_nogest <- perf_epi_sub_nogest %>%
  mutate(mean_age_c = mean_age - mean(mean_age, na.rm = TRUE))

# [ORIGINAL] brain subset (1.main_meta_wMAE.R, line 842)
perf_brain_sub <- perf_brain %>% filter(!model %in% c("DunedinPACNI")) %>% droplevels()

# [ORIGINAL] centre age (1.main_meta_wMAE.R, lines 1555-1560)
perf_brain_sub <- perf_brain_sub %>%
  mutate(
    mean_age_c = mean_age - mean(mean_age, na.rm = TRUE),
    model = droplevels(factor(model))
  )

# [FROM R2.1] centred test-sample age range
perf_epi_sub_nogest <- perf_epi_sub_nogest %>% mutate(age_range_c = age_range - mean(age_range, na.rm = TRUE))
perf_brain_sub      <- perf_brain_sub      %>% mutate(age_range_c = age_range - mean(age_range, na.rm = TRUE))


#### meta-analyses ----

###### (1) Epigenetic clocks: MAEw per clock ~ mean age (as in paper) ------
# outcome: log MAEw | data: perf_epi_sub_nogest | random: cohort/timepoint
# [ORIGINAL] epi_wMAE_mv_modelwise_ageint (1.main_meta_wMAE.R, line 2707)
epi_wMAE_mv_modelwise_ageint <- rma.mv(
  yi = log_wMAE,
  V  = log_wMAE_var,
  mods = ~ model + model:mean_age_c - 1,
  random = list(
    ~ 1 | cohort/timepoint
  ),
  data = perf_epi_sub_nogest,
  method = "REML"
)
summary(epi_wMAE_mv_modelwise_ageint)

###### (2) Epigenetic clocks: unweighted MAE per clock ~ mean age + age range ------
# outcome: log MAE (years) | data: perf_epi_sub_nogest | random: cohort/timepoint
# [ORIGINAL, CHANGED] epi_MAE_mv_modelwise_ageint (4.meta_unweightedMAE.R, line 1006; marked "not run" there).
# CHANGED: + age_range_c added to mods; object name.
epi_MAE_mv_modelwise_ageint_agerange <- rma.mv(
  yi = log_MAE,
  V  = log_MAE_var,
  mods = ~ model + model:mean_age_c + age_range_c - 1,
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
summary(epi_MAE_mv_modelwise_ageint_agerange)

###### (3) Brain age models: MAEw per model ~ mean age (as in paper) ------
# outcome: MAEw (raw scale) | data: perf_brain_sub | random: cohort/timepoint
# [ORIGINAL] brain_wMAE_mv_modelwise_ageint (1.main_meta_wMAE.R, line 1564)
brain_wMAE_mv_modelwise_ageint <- rma.mv(
  yi = wMAE_test,
  V  = wMAE_SE_boot^2,
  mods = ~ model + model:mean_age_c - 1,
  random = list(
    ~ 1 | cohort/timepoint
  ),
  data = perf_brain_sub,
  method = "REML"
)
summary(brain_wMAE_mv_modelwise_ageint)

###### (4) Brain age models: unweighted MAE per model ~ mean age + age range ------
# outcome: MAE (years, raw scale) | data: perf_brain_sub | random: cohort/timepoint
# [ORIGINAL, CHANGED] brain_MAE_mv_modelwise_ageint (4.meta_unweightedMAE.R, line 2130; marked "not run" there).
# CHANGED: + age_range_c added to mods; object name.
brain_MAE_mv_modelwise_ageint_agerange <- rma.mv(
  yi = MAE,
  V  = MAE_SE_boot^2,
  mods = ~ model + model:mean_age_c + age_range_c - 1,
  random = list(
    ~ 1 | cohort/timepoint
  ),
  data = perf_brain_sub,
  method = "REML"
)
summary(brain_MAE_mv_modelwise_ageint_agerange)


#### [NEW] results table: per-clock age effect, MAEw vs unweighted MAE ----
# rows with ":mean_age_c" are the per-clock age slopes
# epi: log scale (% change per year = (exp(b) - 1) * 100); brain: raw scale

# --- epigenetic clocks ---
i_w <- grep(":mean_age_c", rownames(epi_wMAE_mv_modelwise_ageint$b))
i_u <- grep(":mean_age_c", rownames(epi_MAE_mv_modelwise_ageint_agerange$b))

epi_per_clock <- data.frame(
  domain        = "epi",
  model         = sub("^model(.*):mean_age_c$", "\\1", rownames(epi_wMAE_mv_modelwise_ageint$b)[i_w]),
  MAEw_b        = epi_wMAE_mv_modelwise_ageint$b[i_w],
  MAEw_ci.lb    = epi_wMAE_mv_modelwise_ageint$ci.lb[i_w],
  MAEw_ci.ub    = epi_wMAE_mv_modelwise_ageint$ci.ub[i_w],
  MAEw_p        = epi_wMAE_mv_modelwise_ageint$pval[i_w],
  MAE_adj_b     = epi_MAE_mv_modelwise_ageint_agerange$b[i_u],
  MAE_adj_ci.lb = epi_MAE_mv_modelwise_ageint_agerange$ci.lb[i_u],
  MAE_adj_ci.ub = epi_MAE_mv_modelwise_ageint_agerange$ci.ub[i_u],
  MAE_adj_p     = epi_MAE_mv_modelwise_ageint_agerange$pval[i_u]
)

# --- brain age models ---
j_w <- grep(":mean_age_c", rownames(brain_wMAE_mv_modelwise_ageint$b))
j_u <- grep(":mean_age_c", rownames(brain_MAE_mv_modelwise_ageint_agerange$b))

brain_per_clock <- data.frame(
  domain        = "brain",
  model         = sub("^model(.*):mean_age_c$", "\\1", rownames(brain_wMAE_mv_modelwise_ageint$b)[j_w]),
  MAEw_b        = brain_wMAE_mv_modelwise_ageint$b[j_w],
  MAEw_ci.lb    = brain_wMAE_mv_modelwise_ageint$ci.lb[j_w],
  MAEw_ci.ub    = brain_wMAE_mv_modelwise_ageint$ci.ub[j_w],
  MAEw_p        = brain_wMAE_mv_modelwise_ageint$pval[j_w],
  MAE_adj_b     = brain_MAE_mv_modelwise_ageint_agerange$b[j_u],
  MAE_adj_ci.lb = brain_MAE_mv_modelwise_ageint_agerange$ci.lb[j_u],
  MAE_adj_ci.ub = brain_MAE_mv_modelwise_ageint_agerange$ci.ub[j_u],
  MAE_adj_p     = brain_MAE_mv_modelwise_ageint_agerange$pval[j_u]
)

# check: clock order is the same in both models of each domain
stopifnot(all(epi_per_clock$model ==
                sub("^model(.*):mean_age_c$", "\\1", rownames(epi_MAE_mv_modelwise_ageint_agerange$b)[i_u])))
stopifnot(all(brain_per_clock$model ==
                sub("^model(.*):mean_age_c$", "\\1", rownames(brain_MAE_mv_modelwise_ageint_agerange$b)[j_u])))

per_clock_age <- rbind(epi_per_clock, brain_per_clock)

# [NEW] FDR correction for unweighted MAE age effects across models within each domain
# (as in 8.fdr.R: "Modelwise results were corrected across the number of models")
per_clock_age$MAE_adj_pFDR <- ave(per_clock_age$MAE_adj_p, per_clock_age$domain,
                                  FUN = function(p) p.adjust(p, method = "fdr"))

# [NEW] % change per year, back-transformed from log scale (as in SM Table S9); epi only,
# brain is on the raw scale (years per year)
per_clock_age <- per_clock_age %>%
  mutate(
    MAE_adj_pct       = if_else(domain == "epi", (exp(MAE_adj_b) - 1) * 100, NA_real_),
    MAE_adj_pct_ci.lb = if_else(domain == "epi", (exp(MAE_adj_ci.lb) - 1) * 100, NA_real_),
    MAE_adj_pct_ci.ub = if_else(domain == "epi", (exp(MAE_adj_ci.ub) - 1) * 100, NA_real_)
  )

per_clock_age

write.csv(per_clock_age, file.path(out_dir, "R2.1_per_clock_age_effects.csv"), row.names = FALSE)
