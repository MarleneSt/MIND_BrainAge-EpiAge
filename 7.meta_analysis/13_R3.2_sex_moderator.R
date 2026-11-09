###### R3.2: proportion of females as a moderator ------
# Reviewer: "While information about the % of males and females in each cohort
# is provided, data on sex differences in epigenetic aging and brain aging are
# missing."
#
# Cohorts only returned results pooled across sexes, so sex differences cannot
# be estimated directly. As a sample-level check, the proportion of females per
# sample (% female) is added as a moderator to the three main global models,
# using the same structure as the original moderator analyses
# (1.main_meta_wMAE.R): epi MAEw, brain MAEw, and epigenetic PAR - brain PAR
# associations (model 2). Sensitivity analysis: excluding TAG (100% female).
# Note: this is an ecological (sample-level) analysis, not a sex difference.
#
# Labels:
#   [ORIGINAL]  copied from 1.main_meta_wMAE.R (line numbers given), unchanged
#   [ORIGINAL, CHANGED]  copied from the original, with the change noted
#   [FROM R2.5] copied from R2.5_meta_with_intercept.R
#   [NEW]       added for R3.2


#### libraries and folders ----
library(tidyverse)
library(metafor)

results_dir <- "results"  # mod2_eff2.rds # [REPRODUCIBILITY] relative path
out_dir     <- "results/revision/R3.2" # [REPRODUCIBILITY] relative path
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE) # [REPRODUCIBILITY] create folder on a fresh run


#### [NEW] load data ----
perf_epi   <- readRDS(file.path(results_dir, "perf_epi.rds"))
perf_brain <- readRDS(file.path(results_dir, "perf_brain.rds"))
mod2_eff2  <- readRDS(file.path(results_dir, "mod2_eff2.rds"))

# [FROM R2.5] remove GenR duplicates
mod2_eff2 <- dplyr::distinct(mod2_eff2)
nrow(mod2_eff2)   # should be 2448


#### data prep ----
# [ORIGINAL, CHANGED] epi subset, model names and log transform (lines 1812-1833).
# CHANGED: the three original steps are combined into one pipe; same result.
perf_epi_sub <- perf_epi %>%
  filter(!model %in% c("DunedinPACE", "DNAmTL")) %>%
  droplevels() %>%
  mutate(
    model = recode(model, "Bohlin_conv" = "Bohlin", "Knight_conv" = "Knight", "EPIC_conv" = "EPIC"),
    log_wMAE     = log(wMAE_test),
    log_wMAE_var = (wMAE_SE_boot / wMAE_test)^2
  )

# [ORIGINAL] brain subset (line 842)
perf_brain_sub <- perf_brain %>% filter(!model %in% c("DunedinPACNI")) %>% droplevels()

# [NEW] % female per sample
# (association: overlap sample with both measures = n_male.x / n_female.x)
perf_epi_sub   <- perf_epi_sub   %>% mutate(perc_female = 100 * n_female / (n_male + n_female))
perf_brain_sub <- perf_brain_sub %>% mutate(perc_female = 100 * n_female / (n_male + n_female))
mod2_eff2      <- mod2_eff2      %>% mutate(perc_female = 100 * n_female.x / (n_male.x + n_female.x))

# [NEW] check: range of % female
summary(perf_epi_sub$perc_female)
summary(perf_brain_sub$perc_female)
summary(mod2_eff2$perc_female)


#### meta-analyses: % female as moderator ----

###### (1) Epigenetic clocks: MAEw ~ % female (all samples) ------
# outcome: log MAEw | data: perf_epi_sub (k = 667) | random: cohort/timepoint + model
# [ORIGINAL, CHANGED] epi MAEw moderator model, copied from epi_wMAE_mv_generation_epi (line 7647).
# CHANGED: mods = ~ generation_epi - 1  ->  mods = ~ perc_female; object name.
epi_wMAE_mv_perc_female <- rma.mv(
  yi = log_wMAE,
  V  = log_wMAE_var,
  mods = ~ perc_female,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = perf_epi_sub,
  method = "REML"
)
summary(epi_wMAE_mv_perc_female)

###### (2) Brain age models: MAEw ~ % female (all samples) ------
# outcome: MAEw (raw scale) | data: perf_brain_sub (k = 184) | random: cohort/timepoint + model
# [ORIGINAL, CHANGED] brain MAEw moderator model, copied from brain_wMAE_mv_ageoverlap_binary (line 4993).
# CHANGED: mods = ~ overlap_brain  ->  mods = ~ perc_female; object name.
brain_wMAE_mv_perc_female <- rma.mv(
  yi = wMAE_test,
  V  = wMAE_SE_boot^2,
  mods = ~ perc_female,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = perf_brain_sub,
  method = "REML"
)
summary(brain_wMAE_mv_perc_female)

###### (3) Epigenetic PAR - brain PAR associations ~ % female (all samples) ------
# outcome: standardised beta (model 2) | data: mod2_eff2, deduplicated (k = 2448) | random: cohort/timepoint + model_combi
# [ORIGINAL, CHANGED] association moderator model, copied from mod2_res_mv_tissue_cortical (line 10389).
# CHANGED: mods = ~ tissue_cortical  ->  mods = ~ perc_female; object name.
mod2_res_mv_perc_female <- rma.mv(
  yi = RLM_Estimate_scaled,
  V  = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ perc_female,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model_combi
  ),
  data = mod2_eff2,
  method = "REML"
)
summary(mod2_res_mv_perc_female)


#### [NEW] sensitivity: excluding TAG (100% female) ----
# Same three models as above; CHANGED: data filtered to cohort != "TAG"; object name.

###### (4) = (1) excluding TAG: Epigenetic clocks: MAEw ~ % female ------
# outcome: log MAEw | data: perf_epi_sub without TAG (k = 652) | random: cohort/timepoint + model
epi_wMAE_mv_perc_female_noTAG <- rma.mv(
  yi = log_wMAE,
  V  = log_wMAE_var,
  mods = ~ perc_female,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = filter(perf_epi_sub, cohort != "TAG"),
  method = "REML"
)

###### (5) = (2) excluding TAG: Brain age models: MAEw ~ % female ------
# outcome: MAEw (raw scale) | data: perf_brain_sub without TAG (k = 142) | random: cohort/timepoint + model
brain_wMAE_mv_perc_female_noTAG <- rma.mv(
  yi = wMAE_test,
  V  = wMAE_SE_boot^2,
  mods = ~ perc_female,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = filter(perf_brain_sub, cohort != "TAG"),
  method = "REML"
)

###### (6) = (3) excluding TAG: Epigenetic PAR - brain PAR associations ~ % female ------
# outcome: standardised beta (model 2) | data: mod2_eff2 without TAG (k = 2320) | random: cohort/timepoint + model_combi
mod2_res_mv_perc_female_noTAG <- rma.mv(
  yi = RLM_Estimate_scaled,
  V  = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ perc_female,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model_combi
  ),
  data = filter(mod2_eff2, cohort != "TAG"),
  method = "REML"
)


#### [NEW] results table ----
# b = change in outcome per 1 percentage point more females
# (epi MAEw: log scale; brain MAEw: raw scale; association: standardised beta)
res_epi <- data.frame(
  outcome = "epi MAEw (log)", sample = "all",
  k = epi_wMAE_mv_perc_female$k,
  b = epi_wMAE_mv_perc_female$b[2], ci.lb = epi_wMAE_mv_perc_female$ci.lb[2],
  ci.ub = epi_wMAE_mv_perc_female$ci.ub[2], p = epi_wMAE_mv_perc_female$pval[2]
)
res_brain <- data.frame(
  outcome = "brain MAEw", sample = "all",
  k = brain_wMAE_mv_perc_female$k,
  b = brain_wMAE_mv_perc_female$b[2], ci.lb = brain_wMAE_mv_perc_female$ci.lb[2],
  ci.ub = brain_wMAE_mv_perc_female$ci.ub[2], p = brain_wMAE_mv_perc_female$pval[2]
)
res_assoc <- data.frame(
  outcome = "PAR association", sample = "all",
  k = mod2_res_mv_perc_female$k,
  b = mod2_res_mv_perc_female$b[2], ci.lb = mod2_res_mv_perc_female$ci.lb[2],
  ci.ub = mod2_res_mv_perc_female$ci.ub[2], p = mod2_res_mv_perc_female$pval[2]
)
res_epi_noTAG <- data.frame(
  outcome = "epi MAEw (log)", sample = "excl. TAG",
  k = epi_wMAE_mv_perc_female_noTAG$k,
  b = epi_wMAE_mv_perc_female_noTAG$b[2], ci.lb = epi_wMAE_mv_perc_female_noTAG$ci.lb[2],
  ci.ub = epi_wMAE_mv_perc_female_noTAG$ci.ub[2], p = epi_wMAE_mv_perc_female_noTAG$pval[2]
)
res_brain_noTAG <- data.frame(
  outcome = "brain MAEw", sample = "excl. TAG",
  k = brain_wMAE_mv_perc_female_noTAG$k,
  b = brain_wMAE_mv_perc_female_noTAG$b[2], ci.lb = brain_wMAE_mv_perc_female_noTAG$ci.lb[2],
  ci.ub = brain_wMAE_mv_perc_female_noTAG$ci.ub[2], p = brain_wMAE_mv_perc_female_noTAG$pval[2]
)
res_assoc_noTAG <- data.frame(
  outcome = "PAR association", sample = "excl. TAG",
  k = mod2_res_mv_perc_female_noTAG$k,
  b = mod2_res_mv_perc_female_noTAG$b[2], ci.lb = mod2_res_mv_perc_female_noTAG$ci.lb[2],
  ci.ub = mod2_res_mv_perc_female_noTAG$ci.ub[2], p = mod2_res_mv_perc_female_noTAG$pval[2]
)

sex_moderator <- rbind(res_epi, res_brain, res_assoc, res_epi_noTAG, res_brain_noTAG, res_assoc_noTAG)
sex_moderator

write.csv(sex_moderator, file.path(out_dir, "R3.2_perc_female_moderator.csv"), row.names = FALSE)


#### [NEW] SM figure: associations by % female ----
# Points: mean association per cohort x timepoint (across all model pairs), labelled by cohort and mean age.
# Line:   fitted meta-regression from model (3), all samples, with 95% CI.
library(ggrepel)   # labels that do not overlap

# mean association per cohort x timepoint
sample_means <- mod2_eff2 %>%
  group_by(cohort, timepoint) %>%
  summarise(mean_assoc  = mean(RLM_Estimate_scaled),
            perc_female = first(perc_female),
            mean_age_epi = first(mean_age_epi),
            array       = sub("v[12]$", "", first(array)),   # EPICv1/EPICv2 -> EPIC
            tissue      = first(tissue_test),
            .groups = "drop") %>%
  # label: cohort + mean age at DNAm assessment (per wave); array or tissue added
  # where one wave was split into several samples (e.g. GenR 450k / EPIC, K2H blood / saliva)
  mutate(wave = sub("_(450k|epic|blood|saliva)$", "", timepoint)) %>%
  group_by(cohort, wave) %>%
  mutate(age   = round(mean(mean_age_epi)),
         split = if (n() == 1) "" else if (n_distinct(array) > 1) paste0(" (", array, ")") else paste0(" (", tissue, ")")) %>%
  ungroup() %>%
  mutate(label = paste0(sub(" Childhood$", "", cohort), " ", age, "y", split))   # "K2H Childhood" shortened to "K2H"

# fitted line from model (3)
# 100 evenly spaced % female values from lowest to highest in the data (x-positions of the line)
grid <- seq(min(mod2_eff2$perc_female), max(mod2_eff2$perc_female), length.out = 100)
# predicted association (+ 95% CI) at each of these values from the meta-regression
pred <- predict(mod2_res_mv_perc_female, newmods = grid)
# combine into one table for plotting (geom_line = line, geom_ribbon = CI band)
fit  <- data.frame(perc_female = grid, pred = pred$pred, ci.lb = pred$ci.lb, ci.ub = pred$ci.ub)

# with individual dots
p_sex <- ggplot() +
  # all individual estimates (each model pair in each sample), faint, behind everything
  geom_point(data = mod2_eff2, aes(x = perc_female, y = RLM_Estimate_scaled),
             position = position_jitter(width = 0.4, height = 0, seed = 1),   # seed: same figure on every run
             size = 0.6, alpha = 0.15, colour = "grey60") +
  geom_ribbon(data = fit, aes(x = perc_female, ymin = ci.lb, ymax = ci.ub), fill = "grey85", alpha = 0.7) +
  geom_line(data = fit, aes(x = perc_female, y = pred)) +
  geom_point(data = sample_means, aes(x = perc_female, y = mean_assoc), size = 2.5, colour = "grey20") +
  # all labels shown (max.overlaps = Inf); lines connect labels to their point
  # labels on a white background so they stay readable over the fitted line and dots
  geom_label_repel(data = sample_means, aes(x = perc_female, y = mean_assoc, label = label),
                  size = 2.5, max.overlaps = Inf, min.segment.length = 0,
                  fill = scales::alpha("white", 0.85), label.padding = 0.12,
                  box.padding = 0.5, point.padding = 0.3, force = 3, force_pull = 0.5,
                  segment.colour = "grey50", seed = 1) +   # seed: same label layout on every run
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(x = "Female participants (%)",
       y = "Epigenetic PAR - brain PAR association (standardised beta)") +
  theme_classic()

p_sex

ggsave(file.path(out_dir, "R3.2_assoc_by_perc_female_with_individual_estimates.png"), p_sex, width = 7, height = 5, dpi = 600)
ggsave(file.path(out_dir, "R3.2_assoc_by_perc_female_with_individual_estimates.pdf"), p_sex, width = 7, height = 5)

