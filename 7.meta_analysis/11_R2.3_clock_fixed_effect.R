## Data prep and check-1 models are taken from the original script (1.main_meta_wMAE.R).
## New: (1) moderator models refit with clock as fixed effect, (2) sensitivity
## analyses excluding gestational clocks and adjusting for mean age and age range.
## Labels: [ORIGINAL] = taken from 1.main_meta_wMAE.R; [NEW] = added for R2.3.

###### R2.3: age range overlap / tissue match vs. clock (generation) ------

# [NEW] --- libraries, folders and data (so the script runs on its own) ---
library(tidyverse)
library(metafor)

results_dir <- "results" # [REPRODUCIBILITY] relative path
results_dir_R2.3 <- "results/revision/R2.3" # [REPRODUCIBILITY] relative path
dir.create(results_dir_R2.3, recursive = TRUE, showWarnings = FALSE) # [REPRODUCIBILITY] create folder on a fresh run

perf_epi <- readRDS(file.path(results_dir, "perf_epi.rds")) %>%
  mutate(
    cohort = factor(cohort),
    timepoint = factor(timepoint),
    model = factor(model)
  )


# [ORIGINAL] --- data prep: same epi subset and variables as original script ---
# (1.main_meta_wMAE.R, lines 1812-1833, 5900-6000 and 6625-6650; same logic,
# condensed into one block)
perf_epi_sub <- perf_epi %>%
  # same exclusions as original epi MAEw analyses
  filter(!model %in% c("DunedinPACE", "DNAmTL")) %>%
  droplevels() %>%
  mutate(
    # drop "_conv" suffix from gestational clocks
    model = recode(model, "Bohlin_conv" = "Bohlin", "Knight_conv" = "Knight", "EPIC_conv" = "EPIC"),
    # MAEw on log scale (original models did not converge on raw scale);
    # variance via delta method: Var(log X) ≈ (SE_X / X)^2
    log_wMAE = log(wMAE_test),
    log_wMAE_var = (wMAE_SE_boot / wMAE_test)^2,
    # PhenoAge training age range corrected as in original script
    training_age_min_epi = if_else(model == "PhenoAge", 20, training_age_min_epi),
    training_age_max_epi = if_else(model == "PhenoAge", 85, training_age_max_epi),
    # age overlap between training sample and test sample (original definitions)
    overlap_lower = pmax(min_age, training_age_min_epi),
    overlap_upper = pmin(max_age, training_age_max_epi),
    # binary: any overlap (boundaries count as overlap)
    overlap_epi = factor(if_else(overlap_upper >= overlap_lower, "yes", "no"), levels = c("no", "yes")),
    # proportional: % of the test sample's age range covered by training range
    overlap_perc_epi = 100 * if_else(overlap_upper >= overlap_lower, overlap_upper - overlap_lower, 0) /
      (max_age - min_age),
    # tissue match: multi-tissue clocks always count as matched
    tissue_overlap = factor(case_when(
      tissue_training == "multi" ~ "yes",
      tissue_training == tissue_test ~ "yes",
      TRUE ~ "no"), levels = c("no", "yes"))
  )

# [ORIGINAL] age-overlap analyses exclude gestational clocks
perf_epi_sub_nogest <- perf_epi_sub %>%
  filter(!model %in% c("Bohlin", "Knight", "EPIC")) %>%
  droplevels()


# --- check 1: original models (clock as random effect) reproduce preprint? ---
# [ORIGINAL] binary overlap model (original line ~5975; gave b = -0.353)
rma.mv(log_wMAE, log_wMAE_var, mods = ~ overlap_epi,
       random = list(~ 1 | cohort/timepoint, ~ 1 | model),
       data = perf_epi_sub_nogest, method = "REML")

# [ORIGINAL] proportional overlap model (original line ~6305; preprint: b_log = -0.01)
rma.mv(log_wMAE, log_wMAE_var, mods = ~ overlap_perc_epi,
       random = list(~ 1 | cohort/timepoint, ~ 1 | model),
       data = perf_epi_sub_nogest, method = "REML")

# [ORIGINAL] tissue match model (original line ~6670; preprint: b_log = -0.22)
rma.mv(log_wMAE, log_wMAE_var, mods = ~ tissue_overlap,
       random = list(~ 1 | cohort/timepoint, ~ 1 | model),
       data = perf_epi_sub, method = "REML")


# [NEW] --- check 2: which clocks vary in overlap / tissue match across samples? ---
# Only clocks with both "no" and "yes" (or a range of %) can inform the
# within-clock estimates below; clocks that never or always overlap/match
# (e.g. AdaptAge, DamAge, PCGrimAge; multi-tissue clocks) drop out.
# “no”: the cohort’s age range lies entirely outside the ages the clock was trained on 
# “yes”: the age ranges overlap at least partly
table(perf_epi_sub_nogest$model, perf_epi_sub_nogest$overlap_epi)
table(perf_epi_sub$model, perf_epi_sub$tissue_overlap) # same as above but for tissue
perf_epi_sub_nogest %>%
  group_by(model) %>%
  summarise(min_perc = min(overlap_perc_epi), max_perc = max(overlap_perc_epi))


# [NEW] --- joint models: moderator + clock as fixed effect ---
# Clock as fixed effect absorbs all differences between clocks (incl.
# generation), so the moderator effect is estimated only from variation
# WITHIN clocks (across samples). ~1|model is dropped because clock is now
# a fixed effect.


# CHANGED vs original: "+ model" added to mods; "~ 1 | model" removed from random
# binary age overlap (SM only)
epi_wMAE_mv_ageoverlap_binary_clockfixed <- rma.mv(
  yi = log_wMAE,
  V  = log_wMAE_var,
  mods = ~ overlap_epi + model, # in original script, "+ model" was not included
  random = ~ 1 | cohort/timepoint, # in original script was "list(~ 1 | cohort/timepoint, ~ 1 | model)". The new version keeps only the cohort and timepoint random effects. Clock is no longer here because it’s now in mods (+ model) as a fixed effect.
  data = perf_epi_sub_nogest,
  method = "REML"
)
summary(epi_wMAE_mv_ageoverlap_binary_clockfixed)

# proportional age overlap (main text)
epi_wMAE_mv_ageoverlap_perc_clockfixed <- rma.mv(
  yi = log_wMAE,
  V  = log_wMAE_var,
  mods = ~ overlap_perc_epi + model,
  random = ~ 1 | cohort/timepoint,
  data = perf_epi_sub_nogest,
  method = "REML"
)
summary(epi_wMAE_mv_ageoverlap_perc_clockfixed)

# tissue match, incl. gestational clocks (as original)
epi_wMAE_mv_tissue_overlap_clockfixed <- rma.mv(
  yi = log_wMAE,
  V  = log_wMAE_var,
  mods = ~ tissue_overlap + model,
  random = ~ 1 | cohort/timepoint,
  data = perf_epi_sub,
  method = "REML"
)
summary(epi_wMAE_mv_tissue_overlap_clockfixed)

# [NEW] sensitivity: tissue match, excluding gestational clocks
epi_wMAE_mv_tissue_overlap_clockfixed_nogest <- rma.mv(
  yi = log_wMAE,
  V  = log_wMAE_var,
  mods = ~ tissue_overlap + model,
  random = ~ 1 | cohort/timepoint,
  data = perf_epi_sub_nogest,
  method = "REML"
)
summary(epi_wMAE_mv_tissue_overlap_clockfixed_nogest)


## combine above results 
# results for the moderator (= 2nd coefficient, after the intercept) of each model

res_binary <- tibble(
  moderator = "overlap_epi (binary)",
  b     = coef(epi_wMAE_mv_ageoverlap_binary_clockfixed)[2],
  ci.lb = epi_wMAE_mv_ageoverlap_binary_clockfixed$ci.lb[2],
  ci.ub = epi_wMAE_mv_ageoverlap_binary_clockfixed$ci.ub[2],
  p     = epi_wMAE_mv_ageoverlap_binary_clockfixed$pval[2]
)

res_perc <- tibble(
  moderator = "overlap_perc_epi (%)",
  b     = coef(epi_wMAE_mv_ageoverlap_perc_clockfixed)[2],
  ci.lb = epi_wMAE_mv_ageoverlap_perc_clockfixed$ci.lb[2],
  ci.ub = epi_wMAE_mv_ageoverlap_perc_clockfixed$ci.ub[2],
  p     = epi_wMAE_mv_ageoverlap_perc_clockfixed$pval[2]
)

res_tissue <- tibble(
  moderator = "tissue_overlap",
  b     = coef(epi_wMAE_mv_tissue_overlap_clockfixed)[2],
  ci.lb = epi_wMAE_mv_tissue_overlap_clockfixed$ci.lb[2],
  ci.ub = epi_wMAE_mv_tissue_overlap_clockfixed$ci.ub[2],
  p     = epi_wMAE_mv_tissue_overlap_clockfixed$pval[2]
)

res_tissue_nogest <- tibble(
  moderator = "tissue_overlap (no gestational clocks)",
  b     = coef(epi_wMAE_mv_tissue_overlap_clockfixed_nogest)[2],
  ci.lb = epi_wMAE_mv_tissue_overlap_clockfixed_nogest$ci.lb[2],
  ci.ub = epi_wMAE_mv_tissue_overlap_clockfixed_nogest$ci.ub[2],
  p     = epi_wMAE_mv_tissue_overlap_clockfixed_nogest$pval[2]
)

# combine into one table and add % change in MAEw (back-transformed from log scale)
joint_res <- bind_rows(res_binary, res_perc, res_tissue, res_tissue_nogest) %>%
  mutate(pct_change = (exp(b) - 1) * 100)

joint_res


# [NEW] --- sensitivity: additionally adjust for sample mean age and age range ---
# Within a clock, non-overlapping samples tend to be the youngest/narrowest;
# MAEw also depends on age range. Adjusting for both checks that the
# overlap/tissue effects are not driven by sample age characteristics.

epi_wMAE_mv_ageoverlap_binary_clockfixed_adj <- rma.mv(
  yi = log_wMAE,
  V  = log_wMAE_var,
  mods = ~ overlap_epi + model + mean_age + log(age_range),
  random = ~ 1 | cohort/timepoint,
  data = perf_epi_sub_nogest,
  method = "REML"
)

epi_wMAE_mv_ageoverlap_perc_clockfixed_adj <- rma.mv(
  yi = log_wMAE,
  V  = log_wMAE_var,
  mods = ~ overlap_perc_epi + model + mean_age + log(age_range),
  random = ~ 1 | cohort/timepoint,
  data = perf_epi_sub_nogest,
  method = "REML"
)

epi_wMAE_mv_tissue_overlap_clockfixed_adj <- rma.mv(
  yi = log_wMAE,
  V  = log_wMAE_var,
  mods = ~ tissue_overlap + model + mean_age + log(age_range),
  random = ~ 1 | cohort/timepoint,
  data = perf_epi_sub,
  method = "REML"
)

# results for the moderator (= 2nd coefficient) of each model
res_binary_adj <- tibble(
  moderator = "overlap_epi (binary)",
  b     = coef(epi_wMAE_mv_ageoverlap_binary_clockfixed_adj)[2],
  ci.lb = epi_wMAE_mv_ageoverlap_binary_clockfixed_adj$ci.lb[2],
  ci.ub = epi_wMAE_mv_ageoverlap_binary_clockfixed_adj$ci.ub[2],
  p     = epi_wMAE_mv_ageoverlap_binary_clockfixed_adj$pval[2]
)

res_perc_adj <- tibble(
  moderator = "overlap_perc_epi (%)",
  b     = coef(epi_wMAE_mv_ageoverlap_perc_clockfixed_adj)[2],
  ci.lb = epi_wMAE_mv_ageoverlap_perc_clockfixed_adj$ci.lb[2],
  ci.ub = epi_wMAE_mv_ageoverlap_perc_clockfixed_adj$ci.ub[2],
  p     = epi_wMAE_mv_ageoverlap_perc_clockfixed_adj$pval[2]
)

res_tissue_adj <- tibble(
  moderator = "tissue_overlap",
  b     = coef(epi_wMAE_mv_tissue_overlap_clockfixed_adj)[2],
  ci.lb = epi_wMAE_mv_tissue_overlap_clockfixed_adj$ci.lb[2],
  ci.ub = epi_wMAE_mv_tissue_overlap_clockfixed_adj$ci.ub[2],
  p     = epi_wMAE_mv_tissue_overlap_clockfixed_adj$pval[2]
)

joint_res_adj <- bind_rows(res_binary_adj, res_perc_adj, res_tissue_adj) %>%
  mutate(adjusted_for = "mean_age + log(age_range)",
         pct_change = (exp(b) - 1) * 100)

joint_res_adj

# [NEW] --- save for SM table ---
write_csv(joint_res,     file.path(results_dir_R2.3, "overlap_tissue_clock_fixed_R2.3.csv"))
write_csv(joint_res_adj, file.path(results_dir_R2.3, "overlap_tissue_clock_fixed_adj_age_R2.3.csv"))
