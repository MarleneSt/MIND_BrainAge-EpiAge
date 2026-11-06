###### R4.3: offset, slope and MAE after per-cohort linear recalibration (in years) ------

# Per cohort x timepoint x model, recalibration = regress predicted age on
# chronological age (predicted = a + slope * age + residual). This line is fully
# determined by summary statistics each cohort already returned:
#   offset         = mean(predicted - age)                 -> PAD_mean
#   slope          = r * SD(predicted) / SD(age)           -> Pearson, predage_sd, SD age (Table S5)
#   residual SD    = SD of residuals of predicted ~ age    -> PAR_sd
#   recalibrated MAE = sqrt(2/pi) * residual SD (~0.8 x SD; MAE of normal residuals),
#   i.e. the error left after removing offset and slope, in years
# Data: performance data only (no GenR duplicates there).


#### libraries and folders ----
library(tidyverse)
library(metafor)

in_dir  <- "results"   # perf_*.rds # [REPRODUCIBILITY] relative path
out_dir <- "results/revision/R4.3" # [REPRODUCIBILITY] relative path
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE) # [REPRODUCIBILITY] create folder on a fresh run


#### load data ----
perf_epi   <- readRDS(file.path(in_dir, "perf_epi.rds"))
perf_brain <- readRDS(file.path(in_dir, "perf_brain.rds"))

perf <- bind_rows(
  perf_epi   %>% mutate(across(where(is.factor), as.character)),
  perf_brain %>% mutate(across(where(is.factor), as.character))
)

# same exclusions as the original MAEw analyses: these models are not trained
# on chronological age, so offset/slope/MAE are not meaningful for them
perf <- perf %>%
  filter(!model %in% c("DunedinPACE", "DNAmTL", "DunedinPACNI")) %>%
  mutate(model = recode(model, "Bohlin_conv" = "Bohlin", "Knight_conv" = "Knight", "EPIC_conv" = "EPIC"))


#### SD of chronological age per cohort x timepoint (SM Table S5) ----
# Not stored in perf_*.rds. Table S5 gives age as "mean (SD)": age at MRI for
# brain, age at DNAm for epigenetic samples.
# [REPRODUCIBILITY] read from data/cohort_age_sd.csv (extracted once from Table S5
# with archive/make_cohort_age_sd.R) instead of reading SM_Tables.xlsx here
s5 <- read.csv("data/cohort_age_sd.csv")

perf <- perf %>%
  mutate(tp_key = tolower(timepoint)) %>%
  left_join(s5, by = c("cohort", "tp_key", "domain"))

# check: every row should have an SD of age
sum(is.na(perf$sd_age))   # should be 0


#### per-cohort recalibration metrics (years) ----
perf <- perf %>%
  mutate(
    offset      = PAD_mean,
    slope       = Pearson * predage_sd / sd_age,
    resid_sd    = PAR_sd,
    recal_MAE   = sqrt(2 / pi) * resid_sd,
    # Table S5 SDs are rounded to 2 decimals; in very narrow samples (SD < 0.1,
    # e.g. cord blood) this makes the slope imprecise
    sd_age_small = sd_age < 0.1
  )

# check: RMSE^2 = offset^2 + (slope - 1)^2 * SD_age^2 + resid_sd^2 should
# reproduce the reported RMSE (confirms offset/slope/residual SD are consistent)
perf <- perf %>%
  mutate(rmse_check = sqrt(offset^2 + (slope - 1)^2 * sd_age^2 + resid_sd^2))
median(abs(perf$rmse_check - perf$RMSE), na.rm = TRUE)   # should be close to 0

# save per-cohort table (for SM)
per_cohort <- perf %>%
  select(domain, cohort, timepoint, model, N, mean_age, sd_age, sd_age_small,
         offset, slope, recal_MAE, MAE, wMAE_test) %>%
  arrange(domain, model, cohort, timepoint)

write_csv(per_cohort, file.path(out_dir, "R4.3_per_cohort_recalibration.csv"))


#### pooled per model (meta-analysis) ----
# Same model-wise structure as the original MAEw analyses: model as fixed
# effect (no intercept), cohort and timepoint within cohort as random effects.
# Sampling variances:
#   offset:          SD(PAD)^2 / N
#   slope:           resid_sd^2 / ((N - 1) * SD_age^2)
#   log(resid_sd):   1 / (2 * (N - 2))   -> pooled on log scale, back-transformed
perf <- perf %>%
  mutate(
    offset_var       = PAD_sd^2 / N,
    slope_var        = resid_sd^2 / ((N - 1) * sd_age^2),
    log_resid_sd     = log(resid_sd),
    log_resid_sd_var = 1 / (2 * (N - 2))
  )

perf_epi_r   <- perf %>% filter(domain == "epi")
perf_brain_r <- perf %>% filter(domain == "brain")

# --- epigenetic clocks ---
epi_offset_mw <- rma.mv(
  yi = offset, V = offset_var,
  mods = ~ model - 1,
  random = ~ 1 | cohort/timepoint,
  data = perf_epi_r, method = "REML",
  control = list(optimizer = "optim", optmethod = "BFGS")   # default optimizer (nlminb) did not converge for this model only
)

epi_slope_mw <- rma.mv(
  yi = slope, V = slope_var,
  mods = ~ model - 1,
  random = ~ 1 | cohort/timepoint,
  data = perf_epi_r, method = "REML"
)

epi_resid_sd_mw <- rma.mv(
  yi = log_resid_sd, V = log_resid_sd_var,
  mods = ~ model - 1,
  random = ~ 1 | cohort/timepoint,
  data = perf_epi_r, method = "REML"
)

# --- brain age models ---
brain_offset_mw <- rma.mv(
  yi = offset, V = offset_var,
  mods = ~ model - 1,
  random = ~ 1 | cohort/timepoint,
  data = perf_brain_r, method = "REML"
)

brain_slope_mw <- rma.mv(
  yi = slope, V = slope_var,
  mods = ~ model - 1,
  random = ~ 1 | cohort/timepoint,
  data = perf_brain_r, method = "REML"
)

brain_resid_sd_mw <- rma.mv(
  yi = log_resid_sd, V = log_resid_sd_var,
  mods = ~ model - 1,
  random = ~ 1 | cohort/timepoint,
  data = perf_brain_r, method = "REML"
)


#### model-wise results table ----
# recalibrated MAE = sqrt(2/pi) * pooled residual SD (back-transformed from log)
epi_mw <- tibble(
  domain       = "epi",
  model        = sub("^model", "", rownames(epi_offset_mw$b)),
  k            = as.numeric(table(perf_epi_r$model)[sub("^model", "", rownames(epi_offset_mw$b))]),
  offset       = epi_offset_mw$b[, 1],
  offset_lb    = epi_offset_mw$ci.lb,
  offset_ub    = epi_offset_mw$ci.ub,
  slope        = epi_slope_mw$b[, 1],
  slope_lb     = epi_slope_mw$ci.lb,
  slope_ub     = epi_slope_mw$ci.ub,
  recal_MAE    = sqrt(2 / pi) * exp(epi_resid_sd_mw$b[, 1]),
  recal_MAE_lb = sqrt(2 / pi) * exp(epi_resid_sd_mw$ci.lb),
  recal_MAE_ub = sqrt(2 / pi) * exp(epi_resid_sd_mw$ci.ub)
)

brain_mw <- tibble(
  domain       = "brain",
  model        = sub("^model", "", rownames(brain_offset_mw$b)),
  k            = as.numeric(table(perf_brain_r$model)[sub("^model", "", rownames(brain_offset_mw$b))]),
  offset       = brain_offset_mw$b[, 1],
  offset_lb    = brain_offset_mw$ci.lb,
  offset_ub    = brain_offset_mw$ci.ub,
  slope        = brain_slope_mw$b[, 1],
  slope_lb     = brain_slope_mw$ci.lb,
  slope_ub     = brain_slope_mw$ci.ub,
  recal_MAE    = sqrt(2 / pi) * exp(brain_resid_sd_mw$b[, 1]),
  recal_MAE_lb = sqrt(2 / pi) * exp(brain_resid_sd_mw$ci.lb),
  recal_MAE_ub = sqrt(2 / pi) * exp(brain_resid_sd_mw$ci.ub)
)

modelwise <- bind_rows(epi_mw, brain_mw)
print(modelwise %>% mutate(across(where(is.numeric), ~ round(.x, 2))), n = Inf, width = Inf)

write_csv(modelwise, file.path(out_dir, "R4.3_modelwise_offset_slope_recalMAE.csv"))


# meta-analysis

#### overall (global) pooled recalibrated MAE ----
# Same structure as the original global MAEw model: cohort, timepoint within
# cohort and model as random effects.
epi_resid_sd_global <- rma.mv(
  yi = log_resid_sd, V = log_resid_sd_var,
  random = list(~ 1 | cohort/timepoint, ~ 1 | model),
  data = perf_epi_r, method = "REML"
)

brain_resid_sd_global <- rma.mv(
  yi = log_resid_sd, V = log_resid_sd_var,
  random = list(~ 1 | cohort/timepoint, ~ 1 | model),
  data = perf_brain_r, method = "REML"
)

# recalibrated MAE (years) with 95% CI
sqrt(2 / pi) * exp(c(epi_resid_sd_global$b, epi_resid_sd_global$ci.lb, epi_resid_sd_global$ci.ub))
sqrt(2 / pi) * exp(c(brain_resid_sd_global$b, brain_resid_sd_global$ci.lb, brain_resid_sd_global$ci.ub))

# recalibrated MAE (years) with 95% CI
global <- data.frame(
  domain       = c("epi", "brain"),
  k            = c(epi_resid_sd_global$k, brain_resid_sd_global$k),
  recal_MAE    = sqrt(2 / pi) * exp(c(epi_resid_sd_global$b, brain_resid_sd_global$b)),
  recal_MAE_lb = sqrt(2 / pi) * exp(c(epi_resid_sd_global$ci.lb, brain_resid_sd_global$ci.lb)),
  recal_MAE_ub = sqrt(2 / pi) * exp(c(epi_resid_sd_global$ci.ub, brain_resid_sd_global$ci.ub))
)
global

write.csv(global, file.path(out_dir, "R4.3_global_recalMAE.csv"), row.names = FALSE)




#### R4.3: data underlying the SM recalibration tables (epi, brain, per cohort) ----
# [REPRODUCIBILITY] objects and files named by content, not by SM table number
# (was: table_s15/s16/s17 and R4.3_TableS15/S16/S17_data.csv)
library(tidyverse)

out_dir <- "results/revision/R4.3" # [REPRODUCIBILITY] relative path

modelwise  <- read.csv(file.path(out_dir, "R4.3_modelwise_offset_slope_recalMAE.csv"))
global     <- read.csv(file.path(out_dir, "R4.3_global_recalMAE.csv"))
per_cohort <- read.csv(file.path(out_dir, "R4.3_per_cohort_recalibration.csv"))


#### epigenetic clocks and brain age models: pooled ----
# Global row: only recalibrated MAE was pooled globally; offset and slope left empty ("-" in the table)
# log recalibrated MAE = log of the back-transformed estimate (so recalibrated MAE = exp(log), as for MAEw in S9)

# --- epigenetic clocks ---
epi_global <- global %>%
  filter(domain == "epi") %>%
  transmute(Model = "Global", k,
            offset = NA, offset_lb = NA, offset_ub = NA,
            slope = NA, slope_lb = NA, slope_ub = NA,
            log_recal_MAE = log(recal_MAE), recal_MAE, recal_MAE_lb, recal_MAE_ub)

epi_models <- modelwise %>%
  filter(domain == "epi") %>%
  arrange(tolower(model)) %>%
  transmute(Model = model, k,
            offset, offset_lb, offset_ub,
            slope, slope_lb, slope_ub,
            log_recal_MAE = log(recal_MAE), recal_MAE, recal_MAE_lb, recal_MAE_ub)

table_epi <- bind_rows(epi_global, epi_models)
table_epi

# --- brain age models ---
brain_global <- global %>%
  filter(domain == "brain") %>%
  transmute(Model = "Global", k,
            offset = NA, offset_lb = NA, offset_ub = NA,
            slope = NA, slope_lb = NA, slope_ub = NA,
            log_recal_MAE = log(recal_MAE), recal_MAE, recal_MAE_lb, recal_MAE_ub)

brain_models <- modelwise %>%
  filter(domain == "brain") %>%
  arrange(tolower(model)) %>%
  transmute(Model = model, k,
            offset, offset_lb, offset_ub,
            slope, slope_lb, slope_ub,
            log_recal_MAE = log(recal_MAE), recal_MAE, recal_MAE_lb, recal_MAE_ub)

table_brain <- bind_rows(brain_global, brain_models)
table_brain


#### per cohort, timepoint and model ----
# cohort names as in SM Table S5; age as "mean (SD)" as in S5/S6
table_per_cohort <- per_cohort %>%
  mutate(cohort = recode(cohort,
                         "K2H Childhood" = "Kids2Health_childhood",
                         "K2H Infancy"   = "Kids2Health_infancy",
                         "UCI Echo"      = "UCI ECHO"),
         age = sprintf("%.2f (%.2f)", mean_age, sd_age)) %>%
  arrange(tolower(cohort), tolower(timepoint), domain, tolower(model)) %>%
  transmute(cohort, timepoint, modality = domain, model,
            `total n` = N, age,
            offset, slope, `recalibrated MAE` = recal_MAE, MAE, MAEw = wMAE_test)

table_per_cohort


#### save ----
write.csv(table_epi, file.path(out_dir, "R4.3_table_epi_recalibration.csv"), row.names = FALSE)
write.csv(table_brain, file.path(out_dir, "R4.3_table_brain_recalibration.csv"), row.names = FALSE)
write.csv(table_per_cohort, file.path(out_dir, "R4.3_table_per_cohort_recalibration.csv"), row.names = FALSE)
