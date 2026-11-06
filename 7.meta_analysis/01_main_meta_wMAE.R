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

#### perf data loading ----
# renaming performance files to include cohort name
# RUN ONLY ONCE therefore commented out

# # ALSPAC
# folder <- file.path(home, "ALSPAC/analyses/results/performance")
# 
# files <- list.files(folder, pattern = "\\.csv$",full.names = TRUE)
# 
# new_names <- file.path(
#   folder,
#   paste0("ALSPAC_", basename(files))
# )
# 
# file.rename(files, new_names)
# 
# # BHRC
# folder <- file.path(home, "BHRC/analyses/results/performance")
# 
# files <- list.files(folder, pattern = "\\.csv$",full.names = TRUE)
# 
# new_names <- file.path(
#   folder,
#   paste0("BHRC_", basename(files))
# )
# 
# file.rename(files, new_names)
# 
# CannTeen
# folder <- file.path(home, "CannTeen/analyses/results/performance")
# 
# files <- list.files(folder, pattern = "\\.csv$",full.names = TRUE)
# 
# new_names <- file.path(
#   folder,
#   paste0("CannTeen_", basename(files))
# )
# 
# file.rename(files, new_names)
# 
# # DCHS
# folder <- file.path(home, "DCHS/analyses/results/performance")
# 
# files <- list.files(folder, pattern = "\\.csv$",full.names = TRUE)
# 
# new_names <- file.path(
#   folder,
#   paste0("DCHS_", basename(files))
# )
# 
# file.rename(files, new_names)
# 
# # FinnBrain
# folder <- file.path(home, "FinnBrain/analyses/results/performance")
# 
# files <- list.files(folder, pattern = "\\.csv$",full.names = TRUE)
# 
# new_names <- file.path(
#   folder,
#   paste0("FinnBrain_", basename(files))
# )
# 
# file.rename(files, new_names)
# 
# # GenR
# folder <- file.path(home, "GenR/analyses/results/performance")
# 
# files <- list.files(folder, pattern = "\\.csv$",full.names = TRUE)
# 
# new_names <- file.path(
#   folder,
#   paste0("GenR_", basename(files))
# )
# 
# file.rename(files, new_names)
# 
# # GUSTO
# folder <- file.path(home, "/GUSTO/analyses/results/performance")
# 
# files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
# 
# new_names <- file.path(
#   folder,
#   paste0("GUSTO_", basename(files))
# )
# 
# file.rename(files, new_names)
# 
# # Kids2Health_C
# folder <- file.path(home, "Kids2Health_C/analyses/results/performance")
# 
# files <- list.files(folder, pattern = "\\.csv$",full.names = TRUE)
# 
# new_names <- file.path(
#   folder,
#   paste0("Kids2Health_C_", basename(files))
# )
# 
# file.rename(files, new_names)
# 
# # Kids2Health_I
# folder <- file.path(home, "Kids2Health_I/analyses/results/performance")
# 
# files <- list.files(folder, pattern = "\\.csv$",full.names = TRUE)
# 
# new_names <- file.path(
#   folder,
#   paste0("Kids2Health_I_", basename(files))
# )
# 
# file.rename(files, new_names)
# 
# # MTwiNS
# folder <- file.path(home, "MTwiNS/analyses/results/performance")
# 
# files <- list.files(folder, pattern = "\\.csv$",full.names = TRUE)
# 
# new_names <- file.path(
#   folder,
#   paste0("MTwiNS_", basename(files))
# )
# 
# file.rename(files, new_names)
# 
# # NICAP
# folder <- file.path(home, "NICAP/analyses/results/performance")
# 
# files <- list.files(folder, pattern = "\\.csv$",full.names = TRUE)
# 
# new_names <- file.path(
#   folder,
#   paste0("NICAP_", basename(files))
# )
# 
# file.rename(files, new_names)
# 
# # Oregon_ADHD-1000
# folder <- file.path(home, "Oregon_ADHD-1000/analyses/results/performance")
# 
# files <- list.files(folder, pattern = "\\.csv$",full.names = TRUE)
# 
# new_names <- file.path(
#   folder,
#   paste0("Oregon_ADHD-1000_", basename(files))
# )
# 
# file.rename(files, new_names)
# 
# # SAND
# folder <- file.path(home, "SAND/analyses/results/performance")
# 
# files <- list.files(folder, pattern = "\\.csv$",full.names = TRUE)
# 
# new_names <- file.path(
#   folder,
#   paste0("SAND_", basename(files))
# )
# 
# file.rename(files, new_names)
# 
# # TAG
# folder <- file.path(home, "TAG/analyses/results/performance")
# 
# files <- list.files(folder, pattern = "\\.csv$",full.names = TRUE)
# 
# new_names <- file.path(
#   folder,
#   paste0("TAG_", basename(files))
# )
# 
# file.rename(files, new_names)
# 
# # UCI_Echo
# folder <- file.path(home, "UCI_Echo/analyses/results/performance")
# 
# files <- list.files(folder, pattern = "\\.csv$",full.names = TRUE)
# 
# new_names <- file.path(
#   folder,
#   paste0("UCI_Echo_", basename(files))
# )
# 
# file.rename(files, new_names)

# extracting all file paths from data_dir
files_perf <- list.files(
  path = data_dir,
  full.names = TRUE,
  recursive = TRUE
)

# filter to performance files 
files_perf <- files_perf[
  grepl("/performance/.*_performance_(brain|epi)_.*_main\\.csv$", files_perf)
]

# function to parse performance information from file name 
parse_perf <- function(fname) {
  
  base <- basename(fname)
  
  m <- stringr::str_match(
    base,
    "^(.*?)_performance_(brain|epi)_(.+)_main\\.csv$"
  )
  
  tibble::tibble(
    cohort = m[,2],
    domain = m[,3],
    timepoint = m[,4],
    file_path = fname
  )
}
#check the output of this as it might not work if cohorts follow an unusual naming convention 

# apply function 
file_info <- purrr::map_dfr(files_perf, parse_perf)

# read files 
perf_main <- file_info %>%
  mutate(data = purrr::map(file_path, ~ read.csv(.x, check.names = FALSE))) %>%
  select(-file_path) %>%
  tidyr::unnest(data)

# remove predage from model name for cleaner plotting 
perf_main <- perf_main %>%
  mutate(
    model = str_remove(model, "_predage$")
  )

# checkin if there are any timepoints / cohorts with less than 40
small_n <- perf_main %>%
  filter(N < 40) %>%
  select(cohort, timepoint, N, model)

small_n
# only Oregon 38 for pyment --> ok (just because 2 additional ps failed Pyment QC)

# Gen R included some additional models that we're not including so I'm removing them here + CannTeen included Horvath (but only Horvath2013 needed)
perf_main <- perf_main %>%
  filter(!(model %in% c("Mayne", "Lee.RPC", "Lee.CPC", "Lee.refRPC", "Horvath")))

# save this for later reloading without having to rerun the above step 
saveRDS(perf_main, file.path(results_dir, "perf_main.rds")) # [REPRODUCIBILITY] saved to results, not data (was: data_dir), so nothing is written to the data folder

# load it later like 
#perf_main - readRDS(file.path(data_dir, "perf_main.rds"))

#### perf descriptives general ----

model_perf_desc <- perf_main %>%
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

write.csv(model_perf_desc, file.path(results_dir, paste0("model_perf_brain_epi_desc_", "final", ".csv")), row.names = FALSE) # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

#### perf analysis prep ----

# --- Filter to domain ---
perf_brain <- perf_main %>% filter(domain == "brain") %>%
  droplevels()
perf_epi   <- perf_main %>% filter(domain == "epi") %>%
  droplevels()

# --- Fisher z-transform Pearson for secondary analysis ---

perf_brain <- escalc(measure = "ZCOR", ri = Pearson, ni = N, data = perf_brain) %>%
  rename(pearson_z = yi, pearson_var = vi) %>%
  mutate(pearson_se = sqrt(pearson_var))


# for epi, we first invert the correlation with age of DNAmTL so that pooling makes sense as this measure is negatively related to chron age (generally)

perf_epi <- perf_epi %>%
  mutate(
    Pearson = if_else(model == "DNAmTL", -Pearson, Pearson)
  )

perf_epi <- escalc(measure = "ZCOR", ri = Pearson, ni = N, data = perf_epi) %>%
  rename(pearson_z = yi, pearson_var = vi) %>%
  mutate(pearson_se = sqrt(pearson_var))

# --- Check k per model before running ---
perf_brain %>% count(model, name = "k") %>% arrange(k)
perf_epi   %>% count(model, name = "k") %>% arrange(k)

# --- Filter to models with k >= 3 ---
perf_brain <- perf_brain %>% group_by(model) %>% filter(n() >= 3) %>% ungroup() %>% mutate(model = droplevels(factor(model)))
perf_epi   <- perf_epi   %>% group_by(model) %>% filter(n() >= 3) %>% ungroup() %>% mutate(model = droplevels(factor(model)))

# --- make sure model & other vars are factors ---
perf_brain <- perf_brain %>%
  mutate(
    across(
      c(model, cohort, domain, timepoint),
      as.factor
    )
  )

perf_epi <- perf_epi %>%
  mutate(
    across(
      c(model, cohort, domain, timepoint),
      as.factor
    )
  )

# --- check for missing or implausible values ---
perf_brain %>%
  summarise(
    n_missing_wMAE = sum(is.na(wMAE_test)),
    n_missing_se   = sum(is.na(wMAE_SE_boot)),
    n_missing_N    = sum(is.na(N)),
    n_missing_r    = sum(is.na(Pearson)),
    n_nonpos_se    = sum(wMAE_SE_boot <= 0, na.rm = TRUE)
  )

perf_epi %>%
  summarise(
    n_missing_wMAE = sum(is.na(wMAE_test)),
    n_missing_se   = sum(is.na(wMAE_SE_boot)),
    n_missing_N    = sum(is.na(N)),
    n_missing_r    = sum(is.na(Pearson)),
    n_nonpos_se    = sum(wMAE_SE_boot <= 0, na.rm = TRUE)
  )

# MTwiNS tp0_blood highlighted (age = 0 so a lot of outputs cannot be generated), removing it here
perf_epi <- perf_epi %>%
  filter(!(cohort == "MTwiNS" & timepoint == "tp0_blood"))

# --- add in moderators ---

# moderators at the cohort level 
mods_brain <- read.csv(file.path(data_dir,"moderators_cohorts_brain_R.csv"))
mods_epi <- read.csv(file.path(data_dir,"moderators_cohorts_epi_R.csv"))

# check that cohort names correspond & timepoint names 
unique(mods_brain$cohort)
unique(perf_brain$cohort)

# check if joining will work 
perf_brain_keys <- perf_brain %>%
  distinct(cohort, timepoint)

mods_keys <- mods_brain %>%
  distinct(cohort, timepoint)

anti_join(
  perf_brain_keys,
  mods_keys,
  by = c("cohort", "timepoint")
)

anti_join(
  mods_keys,
  perf_brain_keys,
  by = c("cohort", "timepoint")
)

# the latter shows that some rows from the moderator file are not in the perf_brain but that's all fine there arll duplications for later matching with the epi data (e.g., when multiple arrays or tissues are available )

# check if joining will work 
perf_epi_keys <- perf_epi %>%
  distinct(cohort, timepoint)

mods_keys <- mods_epi %>%
  distinct(cohort, timepoint)

anti_join(
  perf_epi_keys,
  mods_keys,
  by = c("cohort", "timepoint")
)

anti_join(
  mods_keys,
  perf_epi_keys,
  by = c("cohort", "timepoint")
)

# Mtwins tp0_blood - excluded so fine to show up here 

perf_brain <- perf_brain %>% left_join(mods_brain, by = c("cohort", "timepoint"))
perf_epi <- perf_epi %>% left_join(mods_epi, by = c("cohort", "timepoint"))

# check if any of the newly added vars are NA 
perf_brain %>% filter(if_any(c(57:60,62,63), is.na)) %>% select(c(1:4,57:63)) # [REPRODUCIBILITY] column positions shifted by -20 after removing unused min/max/whisker columns (was: c(77:80,82,83) / c(1:4,77:83))
perf_epi %>% filter(if_any(c(57:60,62,63), is.na)) %>% select(c(1:4,57:63)) # [REPRODUCIBILITY] column positions shifted by -20 after removing unused min/max/whisker columns (was: c(77:80,82,83) / c(1:4,77:83))


# moderators at the model level 
mods_brain_models <- read.csv(file.path(data_dir,"moderators_brain_models.csv"))
mods_epi_models <- read.csv(file.path(data_dir,"moderators_epi_models.csv"))

# check that cohort names correspond & timepoint names 
unique(mods_brain_models$model)
unique(mods_epi_models$model)

# check if joining will work 
perf_brain_keys <- perf_brain %>%
  distinct(model)

mods_keys <- mods_brain_models %>%
  distinct(model)

anti_join(
  perf_brain_keys,
  mods_keys,
  by = c("model")
)

anti_join(
  mods_keys,
  perf_brain_keys,
  by = c("model")
)

# epi
perf_epi_keys <- perf_epi %>%
  distinct(model)

mods_keys <- mods_epi_models %>%
  distinct(model)

anti_join(
  perf_epi_keys,
  mods_keys,
  by = c("model")
)

anti_join(
  mods_keys,
  perf_epi_keys,
  by = c("model")
)


perf_brain <- perf_brain %>% left_join(mods_brain_models, by = c("model"))
perf_epi <- perf_epi %>% left_join(mods_epi_models, by = c("model"))

# check if any of the newly added vars are NA 
perf_brain %>% filter(if_any(c(64:67), is.na)) %>% select(c(1:4,64:67)) # [REPRODUCIBILITY] column positions shifted by -20 after removing unused min/max/whisker columns (was: c(84:87) / c(1:4,84:87))
perf_epi %>% filter(if_any(c(64:67), is.na)) %>% select(c(1:4,64:67)) # [REPRODUCIBILITY] column positions shifted by -20 after removing unused min/max/whisker columns (was: c(84:87) / c(1:4,84:87))

# --- rename some cohort names for nicer plotting ---

# brain 
perf_brain <- perf_brain %>%
  mutate(
    cohort = recode(
      cohort,
      "SAND" = "FFCWS",
      "Oregon_ADHD-1000" = "Oregon ADHD-1000",
      "Kids2Health_C" = "K2H Childhood",
      "UCI_Echo" = "UCI Echo"
    )
  )

# epi 
perf_epi <- perf_epi %>%
  mutate(
    cohort = recode(
      cohort,
      "SAND" = "FFCWS",
      "Oregon_ADHD-1000" = "Oregon ADHD-1000",
      "Kids2Health_C" = "K2H Childhood",
      "Kids2Health_I" = "K2H Infancy",
      "UCI_Echo" = "UCI Echo"
    )
  )

# save the files 
saveRDS(perf_brain, file.path(results_dir, "perf_brain.rds"))
write.csv(perf_brain, file.path(results_dir, "perf_brain.csv"), row.names = FALSE)

saveRDS(perf_epi, file.path(results_dir, "perf_epi.rds"))
write.csv(perf_epi, file.path(results_dir, "perf_epi.csv"), row.names = FALSE)

#### assoc data loading ----

##### model 2 -----
# extracing model 2 files 
files_mod2 <- list.files(
  path = data_dir,
  pattern = "_model2_PAR\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

parse_fname <- function(fname) {
  base <- basename(fname)
  # e.g. Cohort1_T1_model2_PAR.csv
  m <- str_match(
    base,
    "^([^_]+)_(.+)_model(\\d+)_([^\\.]+)\\.csv$"
  )
  tibble(
    cohort        = m[, 2],
    timepoint     = m[, 3],
    model_num     = m[, 4],
    var_type_short = m[, 5],
    file_path     = fname          # <- avoid name "file"
  )
}
#this might need to be adjusted if people have added in tissue or array 

file_info <- purrr::map_dfr(files_mod2, parse_fname)

mod2_all <- file_info %>% # [REPRODUCIBILITY] renamed (was: mod_all), used as mod2_all below
  mutate(data = purrr::map(file_path, ~ read.csv(.x, check.names = FALSE))) %>%
  select(-file_path) %>%
  unnest(data)

# extract model 2 effect 
mod2_eff <- mod2_all %>%
  filter(
    Model == "Model_2",
    VariableType == "_PAR",
    Predictor == EpiPredictor   # the epi predictor effect
  ) %>%
  mutate(
    brain_model = sub("_PAR$", "", Outcome),
    epi_model   = sub("_PAR$", "", EpiPredictor),
    yi  = RLM_Estimate_scaled, 
    sei = RLM_SE_coeftest_HC3_scaled 
  )


#### assoc data prep ----

# add in model combi 
mod2_eff2 <- mod2_eff %>%
  mutate(
    model_combi = interaction(brain_model, epi_model, drop = TRUE)
  )

# changing sign for DNAmTL 
mod2_eff2 <- mod2_eff2 %>%
  mutate(
    RLM_CI_lower_old = RLM_CI_lower_coeftest_scaled,
    RLM_CI_upper_old = RLM_CI_upper_coeftest_scaled,
    
    tmp_low  = pmin(-RLM_CI_upper_old, -RLM_CI_lower_old),
    tmp_high = pmax(-RLM_CI_upper_old, -RLM_CI_lower_old),
    
    RLM_Estimate_scaled = if_else(
      epi_model == "DNAmTL",
      -RLM_Estimate_scaled,
      RLM_Estimate_scaled
    ),
    
    RLM_CI_lower_coeftest_scaled = if_else(
      epi_model == "DNAmTL",
      tmp_low,
      RLM_CI_lower_coeftest_scaled
    ),
    
    RLM_CI_upper_coeftest_scaled = if_else(
      epi_model == "DNAmTL",
      tmp_high,
      RLM_CI_upper_coeftest_scaled
    )
  ) %>%
  select(-RLM_CI_lower_old, -RLM_CI_upper_old, -tmp_low, -tmp_high)

# add in moderators ()
# cohort level 
mods_association <- read.csv(file.path(data_dir,"moderators_cohorts_overlap_R.csv"))

# moderators at the model level 
mods_brain_models <- read.csv(file.path(data_dir,"moderators_brain_models.csv"))
mods_brain_models$model_brain <- mods_brain_models$model

mods_epi_models <- read.csv(file.path(data_dir,"moderators_epi_models.csv"))
mods_epi_models$model_epi <- mods_epi_models$model

# merge cohort level moderators 

# change names for alignment 
mod2_eff2 <- mod2_eff2 %>%
  mutate(
    cohort = recode(
      cohort,
      "Oregon-ADHD 1000" = "Oregon_ADHD-1000",
      "BHRCS" = "BHRC",
      "K2H" = "Kids2Health_C",
      "GENR" = "GenR",
      "MTwINS" = "MTwiNS"
    )
  )


# check if joining will work (cohort level moderator)
mod2_keys <- mod2_eff2 %>%
  distinct(cohort, timepoint)

mods_keys <- mods_association %>%
  distinct(cohort, timepoint)

anti_join(
  mod2_keys,
  mods_keys,
  by = c("cohort", "timepoint")
)

anti_join(
  mods_keys,
  mod2_keys,
  by = c("cohort", "timepoint")
)

mod2_eff2 <- mod2_eff2 %>% left_join(mods_association, by = c("cohort", "timepoint"))

# we also need to merge in the epi cohort moderators (to get tissue) - mods_epi

mod2_keys <- mod2_eff2 %>%
  distinct(cohort, timepoint)

mods_keys <- mods_epi %>%
  distinct(cohort, timepoint)

anti_join(
  mod2_keys,
  mods_keys,
  by = c("cohort", "timepoint")
)

anti_join(
  mods_keys,
  mod2_keys,
  by = c("cohort", "timepoint")
) # timepoints with no imaging equivalent left over, this is fine 

mod2_eff2 <- mod2_eff2 %>% left_join(mods_epi, by = c("cohort", "timepoint"))

# merge model level moderators 

#adding in variables with the model name for merging 
mod2_eff2$model_brain <- sub("_PAR$", "", mod2_eff2$Outcome)
mod2_eff2$model_epi <- sub("_PAR$", "", mod2_eff2$Predictor)

# check if joining will work (model level moderator - brain)
mod2_keys <- mod2_eff2 %>%
  distinct(model_brain)

mods_keys <- mods_brain_models %>%
  distinct(model_brain)

anti_join(
  mod2_keys,
  mods_keys,
  by = c("model_brain")
)

anti_join(
  mods_keys,
  mod2_keys,
  by = c("model_brain")
)

# check if joining will work (model level moderator - epi)
mod2_keys <- mod2_eff2 %>%
  distinct(model_epi)

mods_keys <- mods_epi_models %>%
  distinct(model_epi)

anti_join(
  mod2_keys,
  mods_keys,
  by = c("model_epi")
)

anti_join(
  mods_keys,
  mod2_keys,
  by = c("model_epi")
) # gestational clocls missing but that's fine

mod2_eff2 <- mod2_eff2 %>% left_join(mods_brain_models, by = c("model_brain"))
mod2_eff2 <- mod2_eff2 %>% left_join(mods_epi_models, by = c("model_epi"))

# rename some cohort names for nicer plotting & to align with performance 
mod2_eff2 <- mod2_eff2 %>%
  mutate(
    cohort = recode(
      cohort,
      "SAND" = "FFCWS",
      "Oregon_ADHD-1000" = "Oregon ADHD-1000",
      "Kids2Health_C" = "K2H Childhood"
    ))

# --- check for missing or implausible values ---
mod2_eff2 %>%
  summarise(
    n_missing_RLM = sum(is.na(RLM_Estimate_scaled)),
    n_missing_SE   = sum(is.na(RLM_SE_coeftest_HC3_scaled)),
    n_missing_N    = sum(is.na(RLM_N)),
    n_missing_age_brain    = sum(is.na(mean_age_brain)),
    n_missing_age_epi    = sum(is.na(mean_age_epi)),
    n_nonpos_se    = sum(RLM_SE_coeftest_HC3_scaled <= 0, na.rm = TRUE)
  ) %>%
  print(width = Inf) # missing is currently canteen 

# check if all is in order
write.csv(mod2_eff2, file.path(results_dir, "mod2_eff2.csv"), row.names = FALSE)
saveRDS(mod2_eff2, file.path(results_dir, "mod2_eff2.rds"))

##### assoc descriptives general -----

assoc_desc <- mod2_eff2 %>%
  group_by(brain_model, epi_model) %>%
  summarise(
    k = n(),
    beta_mean = mean(RLM_Estimate_scaled, na.rm = TRUE),
    beta_sd = sd(RLM_Estimate_scaled, na.rm = TRUE),
    beta_min = min(RLM_Estimate_scaled, na.rm = TRUE),
    beta_max = max(RLM_Estimate_scaled, na.rm = TRUE),
    se_mean = mean(RLM_SE_coeftest_HC3_scaled, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(assoc_desc, file.path(results_dir, paste0("assoc_desc_model2_", "final", ".csv")), row.names = FALSE) # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

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

#### analyses perf brain ----

##### brain wMAE no DunedinPACNI -----

# exclude specific models if needed 
perf_brain_sub <- perf_brain %>% filter(!model %in% c("DunedinPACNI")) %>% droplevels()

###### brain wMAE global ------

# --- set up folder ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_brain <- file.path(results_dir, paste0("perf_brain_wMAE_noDunedin_global_", run_id))
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

cat("Model summary\n")
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

###### wMAE modelwise ------

# --- set-up ---
run_id_modelwise <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_brain_modelwise<- file.path(results_dir, paste0("perf_brain_wMAE_noDunedin_modelwise_", run_id_modelwise))

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

###### wMAE age global ------
# all but DunedinPACNI

# --- set up ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_brain_global_ageint <- file.path(results_dir, paste0("perf_brain_wMAE_noDunedin_global_ageint_", run_id))

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

###### brain wMAE age modelwise ------
# all but DunedinPACNI

# --- set up ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_brain_modelwise_ageint <- file.path(
  results_dir,
  paste0("perf_brain_wMAE_noDunedin_modelwise_ageint_", run_id)
)

dir.create(results_dir_brain_modelwise_ageint, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_brain_modelwise_ageint, "analysis_log_modelwise_ageint.txt")


# --- data prep ---

perf_brain_sub <- perf_brain_sub %>%
  mutate(
    mean_age_c = mean_age - mean(mean_age, na.rm = TRUE),
    model = droplevels(factor(model))
  )


# --- meta ---

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


# optional comparison with global age effect model
age_model_comparison <- anova(
  brain_wMAE_mv_global_ageint,
  brain_wMAE_mv_modelwise_ageint,
  refit = TRUE
)

capture.output(
  age_model_comparison,
  file = file.path(results_dir_brain_modelwise_ageint, "model_comparison_global_vs_modelwise_age.txt")
)

# --- extract modelwise estimates ---

modelwise_age_est <- tibble(
  term = rownames(brain_wMAE_mv_modelwise_ageint$b),
  estimate = as.numeric(brain_wMAE_mv_modelwise_ageint$b),
  se = brain_wMAE_mv_modelwise_ageint$se,
  ci.lb = brain_wMAE_mv_modelwise_ageint$ci.lb,
  ci.ub = brain_wMAE_mv_modelwise_ageint$ci.ub,
  pval = brain_wMAE_mv_modelwise_ageint$pval
) %>%
  mutate(
    type = if_else(str_detect(term, "mean_age_c"), "age_slope", "intercept"),
    model = term %>%
      str_remove("^model") %>%
      str_remove(":mean_age_c") %>%
      str_remove("\\.mean_age_c")
  )

X <- model.matrix(brain_wMAE_mv_modelwise_ageint)

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

# Interpretation:
# intercept = predicted wMAE for that model at the mean sample age
# age_slope = absolute change in wMAE per one-year increase in mean age

write.csv(
  modelwise_age_est,
  file.path(results_dir_brain_modelwise_ageint, "modelwise_age_interaction_results.csv"),
  row.names = FALSE
)

saveRDS(
  brain_wMAE_mv_modelwise_ageint,
  file.path(results_dir_brain_modelwise_ageint, "model_brain_wMAE_mv_modelwise_ageint.rds")
)

capture.output(
  summary(brain_wMAE_mv_modelwise_ageint),
  file = file.path(results_dir_brain_modelwise_ageint, "model_summary_raw_scale.txt")
)

# --- heterogeneity ---

I2 <- calc_I2_mv(
  brain_wMAE_mv_modelwise_ageint,
  component_names = c("cohort", "timepoint_within_cohort")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_brain_modelwise_ageint, "heterogeneity_I2.csv"), row.names = FALSE)


# --- plotting ---

mean_age_mean <- mean(perf_brain_sub$mean_age, na.rm = TRUE)

newdat <- expand.grid(
  mean_age_c = seq(
    min(perf_brain_sub$mean_age_c, na.rm = TRUE),
    max(perf_brain_sub$mean_age_c, na.rm = TRUE),
    length.out = 100
  ),
  model = levels(droplevels(perf_brain_sub$model))
)

X_new <- model.matrix(~ model + model:mean_age_c - 1, newdat)

pred <- predict(
  brain_wMAE_mv_modelwise_ageint,
  newmods = X_new
)

newdat$mean_age <- newdat$mean_age_c + mean_age_mean
newdat$pred  <- as.numeric(pred$pred)
newdat$ci.lb <- as.numeric(pred$ci.lb)
newdat$ci.ub <- as.numeric(pred$ci.ub)

png(
  filename = file.path(
    results_dir_brain_modelwise_ageint,
    "plot_brain_meanage_wMAE_modelwise.png"
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
    y = "Predicted wMAE",
    title = "Brain age model-specific age effects"
  ) +
  theme_minimal(base_size = 14)

dev.off()


# --- LOO ---

cohorts <- unique(perf_brain_sub$cohort)

loo_modelwise_age <- lapply(cohorts, function(cc) {
  
  tmp <- perf_brain_sub %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = wMAE_test,
    V  = wMAE_SE_boot^2,
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
      type = if_else(str_detect(term, "mean_age_c"), "age_slope", "intercept"),
      model = term %>%
        str_remove("^model") %>%
        str_remove(":mean_age_c") %>%
        str_remove("\\.mean_age_c")
    )
}) %>%
  bind_rows()

write.csv(
  loo_modelwise_age,
  file.path(results_dir_brain_modelwise_ageint, "modelwise_leave_one_cohort_out_age_interactions.csv"),
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
  file.path(results_dir_brain_modelwise_ageint, "modelwise_leave_one_cohort_out_age_interactions_summary.csv"),
  row.names = FALSE
)


# --- log ---

sink(log_file)

cat("Brain wMAE multilevel meta-analysis - modelwise age interactions, raw scale\n")
cat("=========================================================================\n\n")
cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_brain_modelwise_ageint, "\n\n")

cat("Model summary, raw wMAE scale\n")
cat("=============================\n")
print(summary(brain_wMAE_mv_modelwise_ageint))

cat("\n\nModelwise age interaction results\n")
cat("=================================\n")
print(modelwise_age_est, n = Inf, width = Inf)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_modelwise_age, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary for age slopes\n")
cat("===========================================\n")
print(loo_modelwise_age_summary, n = Inf, width = Inf)

sink()

#### analyses perf epi ----

##### wMAE no DunedinPACNI & DNAmTL -----

###### epi wMAE global ------

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
    log_wMAE = log(wMAE_test),
    # Delta method: Var(log(X)) ≈ (SE_X / X)^2
    log_wMAE_var = (wMAE_SE_boot / wMAE_test)^2
  )

# --- set up folder ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi <- file.path(results_dir, paste0("perf_epi_wMAE_noDunedinDNAmTL_global_", run_id))

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

cat("Model summary\n")
cat("=============\n")
print(summary(epi_wMAE_mv_global))

cat("Model summary - back transformed \n")
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

###### epi wMAE modelwise ------

# --- set up folder ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi_modelwise <- file.path(results_dir, paste0("perf_epi_wMAE_noDunedinDNAmTL_modelwise_", run_id))

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

###### epi wMAE age global ------
# all but DunedinPACNI & DNAmTL

# --- set up ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi_global_ageint <- file.path(results_dir, paste0("perf_epi_wMAE_noDunedinDNAmTL_global_ageint_", run_id))

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

###### epi wMAE age modelwise ------
# all but DunedinPACNI & DNAmTL

# --- set up ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi_modelwise_ageint <- file.path(results_dir, paste0("perf_epi_wMAE_noDunedinDNAmTL_modelwise_ageint_", run_id))

dir.create(results_dir_epi_modelwise_ageint, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_epi_modelwise_ageint, "analysis_log_modelwise_ageint.txt")

# --- data prep ---

# center age 
perf_epi_sub_nogest <- perf_epi_sub_nogest %>%
  mutate(mean_age_c = mean_age - mean(mean_age, na.rm = TRUE))

# --- meta ---

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

# --- heterogeneity estimates --
I2 <- calc_I2_mv(
  epi_wMAE_mv_modelwise_ageint,
  component_names = c("cohort", "timepoint_within_cohort")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_epi_modelwise_ageint, "heterogeneity_I2.csv"), row.names = FALSE)


# compare model with common slope to model specific splopes
age_model_comparison <- anova(
  epi_wMAE_mv_global_ageint,
  epi_wMAE_mv_modelwise_ageint,
  refit = TRUE
)

age_model_comparison

capture.output(
  age_model_comparison,
  file = file.path(results_dir_epi_modelwise_ageint, "model_comparison_global_vs_modelwise_age.txt")
)

modelwise_age_est <- tibble(
  term = rownames(epi_wMAE_mv_modelwise_ageint$b),
  estimate_log = as.numeric(epi_wMAE_mv_modelwise_ageint$b),
  se_log = epi_wMAE_mv_modelwise_ageint$se,
  ci.lb_log = epi_wMAE_mv_modelwise_ageint$ci.lb,
  ci.ub_log = epi_wMAE_mv_modelwise_ageint$ci.ub,
  pval = epi_wMAE_mv_modelwise_ageint$pval
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

X <- model.matrix(epi_wMAE_mv_modelwise_ageint)

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

# For interpretation:
#   
#   rows with type == "intercept" give geometric mean wMAE at the mean age.
# rows with type == "age_slope" give percent change in geometric mean wMAE per one-year increase in age

write.csv(
  modelwise_age_est,
  file.path(results_dir_epi_modelwise_ageint, "modelwise_age_interaction_backtransformed_results.csv"),
  row.names = FALSE
)

saveRDS(
  epi_wMAE_mv_modelwise_ageint,
  file.path(results_dir_epi_modelwise_ageint, "model_epi_wMAE_mv_modelwise_ageint.rds")
)

capture.output(
  summary(epi_wMAE_mv_modelwise_ageint),
  file = file.path(results_dir_epi_modelwise_ageint, "model_summary_log_scale.txt")
)


# --- plotting (should be adapted) ---
# plotting model-specific age effects, log-transformed model

mean_age_mean <- mean(perf_epi_sub_nogest$mean_age, na.rm = TRUE)

newdat <- expand.grid(
  mean_age_c = seq(
    min(perf_epi_sub_nogest$mean_age_c, na.rm = TRUE),
    max(perf_epi_sub_nogest$mean_age_c, na.rm = TRUE),
    length.out = 100
  ),
  model = unique(perf_epi_sub_nogest$model))

pred <- predict(
  epi_wMAE_mv_modelwise_ageint,
  newmods = model.matrix(~ model + model:mean_age_c - 1, newdat)
)

newdat <- newdat %>%
  mutate(
    # raw age for interpretable x-axis
    mean_age = mean_age_c + mean_age_mean,
    
    # log-scale predictions
    log_pred  = pred$pred,
    log_ci.lb = pred$ci.lb,
    log_ci.ub = pred$ci.ub,
    
    # back-transformed predictions
    pred  = exp(log_pred),
    ci.lb = exp(log_ci.lb),
    ci.ub = exp(log_ci.ub)
  )

png(
  filename = file.path(
    results_dir_epi_modelwise_ageint,
    "plot_epi_meanage_wMAE_modelwise_backtransformed.png"
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
    y = "Predicted geometric mean wMAE",
    title = "Model-specific age effects"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- LOO ---

cohorts <- unique(perf_epi_sub_nogest$cohort)

loo_modelwise_age <- lapply(cohorts, function(cc) {
  
  tmp <- perf_epi_sub_nogest %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = log_wMAE,
    V  = log_wMAE_var,
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

cat("Epi wMAE multilevel meta-analysis - modelwise age interactions, log-transformed\n")
cat("================================================================================\n\n")
cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_epi_modelwise_ageint, "\n\n")

cat("Model summary, log scale\n")
cat("========================\n")
print(summary(epi_wMAE_mv_modelwise_ageint))

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nModel comparison: global age slope vs model-specific age slopes\n")
cat("===============================================================\n")
print(age_model_comparison)

cat("\n\nEffect sizes per model\n")
cat("======================\n")
print(effect_sizes_per_model, n = Inf, width = Inf)

cat("\n\nBack-transformed modelwise age interaction results\n")
cat("==================================================\n")
print(modelwise_age_est, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_modelwise_age, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary for age slopes\n")
cat("===========================================\n")
print(loo_modelwise_age_summary, n = Inf, width = Inf)

sink()

#### analyses associations  ----

#### M2 PAR ----

###### assoc global -----

# --- set up folder ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_assoc <- file.path(
  results_dir,
  paste0("assoc_model2_PAR_global_", run_id)
)

dir.create(results_dir_assoc, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_assoc, "analysis_log.txt")

# --- meta ---
mod2_res_mv_global <- rma.mv(
  yi = RLM_Estimate_scaled,
  V  = RLM_SE_coeftest_HC3_scaled^2,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model_combi
  ),
  data = mod2_eff2,
  method = "REML",
  verbose = TRUE
)

summary(mod2_res_mv_global)

saveRDS(
  mod2_res_mv_global,
  file.path(results_dir_assoc, "model_assoc_mv_global.rds")
)

capture.output(
  summary(mod2_res_mv_global),
  file = file.path(results_dir_assoc, "model_summary.txt")
)

# --- heterogeneity ---
I2 <- calc_I2_mv(
  mod2_res_mv_global,
  component_names = c("cohort", "timepoint_within_cohort", "model_combi")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_assoc, "heterogeneity_I2.csv"), row.names = FALSE)

# --- extracting estimates ---
assoc_global_est <- tibble(
  term = names(coef(mod2_res_mv_global)),
  estimate = as.numeric(coef(mod2_res_mv_global)),
  se = mod2_res_mv_global$se,
  ci.lb = mod2_res_mv_global$ci.lb,
  ci.ub = mod2_res_mv_global$ci.ub,
  pval = mod2_res_mv_global$pval
)

write.csv(
  assoc_global_est,
  file.path(results_dir_assoc, "global_estimates.csv"),
  row.names = FALSE
)

effect_sizes_overall <- tibble(
  k = nrow(mod2_eff2),
  n_cohorts = n_distinct(mod2_eff2$cohort),
  n_timepoints = n_distinct(paste(mod2_eff2$cohort, mod2_eff2$timepoint)),
  n_model_combis = n_distinct(mod2_eff2$model_combi)
)

write.csv(
  effect_sizes_overall,
  file.path(results_dir_assoc, "effect_sizes_overall.csv"),
  row.names = FALSE
)

# --- LOO cohort ---

cohorts <- unique(mod2_eff2$cohort)

loo_cohort <- lapply(cohorts, function(cc) {
  
  tmp <- subset(mod2_eff2, cohort != cc)
  
  fit <- rma.mv(
    yi = RLM_Estimate_scaled,
    V  = RLM_SE_coeftest_HC3_scaled^2,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model_combi
    ),
    data = tmp,
    method = "REML"
  )
  
  data.frame(
    cohort_removed = cc,
    beta = coef(fit)[1],
    se   = sqrt(vcov(fit)[1, 1])
  )
})

loo_cohort <- do.call(rbind, loo_cohort) %>%
  as_tibble()

loo_cohort
range(loo_cohort$beta)
sd(loo_cohort$beta)

loo_summary <- tibble(
  min_beta   = min(loo_cohort$beta, na.rm = TRUE),
  max_beta   = max(loo_cohort$beta, na.rm = TRUE),
  range_beta = max_beta - min_beta,
  sd_beta    = sd(loo_cohort$beta, na.rm = TRUE)
)

write.csv(
  loo_cohort,
  file.path(results_dir_assoc, "leave_one_cohort_out.csv"),
  row.names = FALSE
)

write.csv(
  loo_summary,
  file.path(results_dir_assoc, "leave_one_cohort_out_summary.csv"),
  row.names = FALSE
)


full_beta <- coef(mod2_res_mv_global)[1]

png(
  filename = file.path(results_dir_assoc, "plot_assoc_global_loo.png"),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_cohort,
  aes(
    x = reorder(cohort_removed, beta),
    y = beta
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = beta - 1.96 * se,
      ymax = beta + 1.96 * se
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept = full_beta,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "Pooled standardised beta",
    title = "Leave-one-cohort-out sensitivity of pooled association"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- full analysis log ---

sink(log_file)

cat("Association multilevel meta-analysis\n")
cat("====================================\n\n")
cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_assoc, "\n\n")

cat("Model summary\n")
cat("=============\n")
print(summary(mod2_res_mv_global))

cat("\n\nGlobal association estimate\n")
cat("===========================\n")
print(assoc_global_est, n = Inf, width = Inf)

cat("\n\nEffect-size counts\n")
cat("==================\n")
print(effect_sizes_overall, n = Inf, width = Inf)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_cohort)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_summary)

sink()

###### assoc modelwise ------

# --- set up folder ---

run_id_modelwise <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_assoc_modelwise <- file.path(
  results_dir,
  paste0("assoc_model2_PAR_modelwise_", run_id_modelwise)
)

dir.create(results_dir_assoc_modelwise, recursive = TRUE, showWarnings = FALSE)

log_file_modelwise <- file.path(results_dir_assoc_modelwise, "analysis_log_modelwise.txt")


# --- meta ---

mod2_res_mv_modelwise <- rma.mv(
  yi   = RLM_Estimate_scaled,
  V    = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ model_combi - 1,
  random = ~ 1 | cohort/timepoint,
  data   = mod2_eff2,
  method = "REML"
)

summary(mod2_res_mv_modelwise)

saveRDS(
  mod2_res_mv_modelwise,
  file.path(results_dir_assoc_modelwise, "model_assoc_mv_modelwise.rds")
)

capture.output(
  summary(mod2_res_mv_modelwise),
  file = file.path(results_dir_assoc_modelwise, "modelwise_model_summary.txt")
)


# --- extract modelwise estimates ---

assoc_modelwise_est <- tibble(
  term  = rownames(mod2_res_mv_modelwise$b),
  beta  = as.numeric(mod2_res_mv_modelwise$b),
  se    = mod2_res_mv_modelwise$se,
  ci.lb = mod2_res_mv_modelwise$ci.lb,
  ci.ub = mod2_res_mv_modelwise$ci.ub,
  pval  = mod2_res_mv_modelwise$pval
) %>%
  mutate(model_combi = str_remove(term, "^model_combi")) %>%
  separate(model_combi, into = c("brain_model", "epi_model"), sep = "\\.", remove = FALSE) %>%
  select(-term)

assoc_modelwise_est <- assoc_modelwise_est %>%
  arrange(beta)

X <- model.matrix(mod2_res_mv_modelwise)

model_combi_cols <- colnames(X)[
  str_detect(colnames(X), "^model_combi") &
    !str_detect(colnames(X), ":")
]

effect_sizes_per_model_combi <- tibble(
  model_combi = str_remove(model_combi_cols, "^model_combi"),
  k = colSums(X[, model_combi_cols, drop = FALSE] != 0)
)

assoc_modelwise_est <- assoc_modelwise_est %>%
  left_join(effect_sizes_per_model_combi, by = "model_combi")

write.csv(
  assoc_modelwise_est,
  file.path(results_dir_assoc_modelwise, "modelwise_estimates.csv"),
  row.names = FALSE
)

# --- heterogeneity ---
I2 <- calc_I2_mv(
  mod2_res_mv_modelwise,
  component_names = c("cohort", "timepoint_within_cohort")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_assoc_modelwise, "heterogeneity_I2.csv"), row.names = FALSE)

# --- LOO cohort ---
# slightly different code here as it doesn't run 

# cohorts <- unique(mod2_eff2$cohort)
# 
# loo_cohort_modelwise <- lapply(cohorts, function(cc) {
#   
#   message("Running LOO modelwise: removing ", cc)
#   
#   tmp <- subset(mod2_eff2, cohort != cc)
#   
#   fit <- tryCatch(
#     rma.mv(
#       yi   = RLM_Estimate_scaled,
#       V    = RLM_SE_coeftest_HC3_scaled^2,
#       mods = ~ model_combi - 1,
#       random = ~ 1 | cohort/timepoint,
#       data = tmp,
#       method = "REML"
#     ),
#     error = function(e) {
#       message("LOO failed for cohort: ", cc)
#       message("Error: ", e$message)
#       return(NULL)
#     }
#   )
#   
#   if (is.null(fit)) {
#     return(tibble(
#       cohort_removed = cc,
#       term = NA_character_,
#       brain_model = NA_character_,
#       epi_model = NA_character_,
#       beta = NA_real_,
#       se = NA_real_,
#       ci.lb = NA_real_,
#       ci.ub = NA_real_,
#       pval = NA_real_,
#       converged = FALSE
#     ))
#   }
#   
#   tibble(
#     cohort_removed = cc,
#     term  = rownames(fit$b),
#     beta  = as.numeric(fit$b),
#     se    = fit$se,
#     ci.lb = fit$ci.lb,
#     ci.ub = fit$ci.ub,
#     pval  = fit$pval,
#     converged = TRUE
#   ) %>%
#     mutate(model_combi = str_remove(term, "^model_combi")) %>%
#     separate(model_combi, into = c("brain_model", "epi_model"), sep = "\\.", remove =FALSE) %>%
#     select(cohort_removed, brain_model, epi_model, beta, se, ci.lb, ci.ub, pval, converged)
# }) %>%
#   bind_rows()

# cohorts <- unique(mod2_eff2$cohort)
# 
# loo_cohort_modelwise <- lapply(cohorts, function(cc) {
#   
#   tmp <- subset(mod2_eff2, cohort != cc)
#   
#   fit <- rma.mv(
#     yi   = RLM_Estimate_scaled,
#     V    = RLM_SE_coeftest_HC3_scaled^2,
#     mods = ~ model_combi - 1,
#     random = ~ 1 | cohort/timepoint,
#     data = tmp,
#     method = "REML"
#   )
#   
#   tibble(
#     cohort_removed = cc,
#     term  = rownames(fit$b),
#     beta  = as.numeric(fit$b),
#     se    = fit$se,
#     ci.lb = fit$ci.lb,
#     ci.ub = fit$ci.ub,
#     pval  = fit$pval
#   ) %>%
#     mutate(model_combi = str_remove(term, "^model_combi")) %>%
#     separate(model_combi, into = c("brain_model", "epi_model"), sep = "\\.") %>%
#     select(cohort_removed, brain_model, epi_model, beta, se, ci.lb, ci.ub, pval)
# }) %>%
#   bind_rows()

# had issues with LOCO as models without FFCWS did not converge, with the help of ChatGPT I created a fallback options whereas different/simpler models are tried out 

# function
fit_meta_modelwise <- function(dat) {
  
  dat <- dat %>%
    droplevels()
  
  # remove empty / unsupported model_combi levels
  dat <- dat %>%
    group_by(model_combi) %>%
    filter(n() >= 2) %>%
    ungroup()
  
  # fallback if too little data
  if (nrow(dat) < 3 || n_distinct(dat$model_combi) < 1) {
    return(list(fit = NULL, meta_model = "failed: insufficient data"))
  }
  
  # 1. preferred model
  res <- tryCatch(
    rma.mv(
      yi   = RLM_Estimate_scaled,
      V    = RLM_SE_coeftest_HC3_scaled^2,
      mods = ~ model_combi - 1,
      random = ~ 1 | cohort/timepoint,
      data = dat,
      method = "REML"
    ),
    error = function(e) NULL
  )
  
  if (!is.null(res)) {
    return(list(fit = res, meta_model = "rma.mv REML cohort/timepoint"))
  }
  
  # 2. same model, alternative optimizer
  res <- tryCatch(
    rma.mv(
      yi   = RLM_Estimate_scaled,
      V    = RLM_SE_coeftest_HC3_scaled^2,
      mods = ~ model_combi - 1,
      random = ~ 1 | cohort/timepoint,
      data = dat,
      method = "REML",
      control = list(
        optimizer = "optim",
        optmethod = "BFGS",
        maxit = 10000
      )
    ),
    error = function(e) NULL
  )
  
  if (!is.null(res)) {
    return(list(fit = res, meta_model = "rma.mv REML BFGS cohort/timepoint"))
  }
  
  # 3. simpler multilevel fallback: cohort only
  res <- tryCatch(
    rma.mv(
      yi   = RLM_Estimate_scaled,
      V    = RLM_SE_coeftest_HC3_scaled^2,
      mods = ~ model_combi - 1,
      random = ~ 1 | cohort,
      data = dat,
      method = "REML"
    ),
    error = function(e) NULL
  )
  
  if (!is.null(res)) {
    return(list(fit = res, meta_model = "rma.mv REML cohort only fallback"))
  }
  
  # 4. final fallback: no multilevel structure
  res <- tryCatch(
    rma(
      yi   = RLM_Estimate_scaled,
      vi   = RLM_SE_coeftest_HC3_scaled^2,
      mods = ~ model_combi - 1,
      data = dat,
      method = "REML"
    ),
    error = function(e) NULL
  )
  
  if (!is.null(res)) {
    return(list(fit = res, meta_model = "rma REML fallback"))
  }
  
  return(list(fit = NULL, meta_model = "failed"))
}

# fit model 
cohorts <- unique(mod2_eff2$cohort)

loo_cohort_modelwise <- lapply(cohorts, function(cc) {
  
  message("Running LOO modelwise: removing ", cc)
  
  tmp <- mod2_eff2 %>%
    filter(cohort != cc) %>%
    droplevels()
  
  fit_obj <- fit_meta_modelwise(tmp)
  fit <- fit_obj$fit
  
  if (is.null(fit)) {
    return(tibble(
      cohort_removed = cc,
      model_combi = NA_character_,
      brain_model = NA_character_,
      epi_model = NA_character_,
      beta = NA_real_,
      se = NA_real_,
      ci.lb = NA_real_,
      ci.ub = NA_real_,
      pval = NA_real_,
      converged = FALSE,
      meta_model = fit_obj$meta_model
    ))
  }
  
  tibble(
    cohort_removed = cc,
    term  = rownames(fit$b),
    beta  = as.numeric(fit$b),
    se    = fit$se,
    ci.lb = fit$ci.lb,
    ci.ub = fit$ci.ub,
    pval  = fit$pval,
    converged = TRUE,
    meta_model = fit_obj$meta_model
  ) %>%
    mutate(model_combi = str_remove(term, "^model_combi")) %>%
    separate(
      model_combi,
      into = c("brain_model", "epi_model"),
      sep = "\\.",
      remove = FALSE
    ) %>%
    select(
      cohort_removed,
      model_combi,
      brain_model,
      epi_model,
      beta,
      se,
      ci.lb,
      ci.ub,
      pval,
      converged,
      meta_model
    )
}) %>%
  bind_rows()

loo_cohort_modelwise %>%
  count(cohort_removed, meta_model)

write.csv(
  loo_cohort_modelwise,
  file.path(results_dir_assoc_modelwise, "modelwise_leave_one_cohort_out.csv"),
  row.names = FALSE
)

modelwise_loo_summary <- loo_cohort_modelwise %>%
  filter(converged) %>%
  group_by(brain_model, epi_model) %>%
  summarise(
    n_loo      = n(),
    meta_models_used = paste(unique(meta_model), collapse = "; "),
    min_beta   = min(beta, na.rm = TRUE),
    max_beta   = max(beta, na.rm = TRUE),
    range_beta = max_beta - min_beta,
    sd_beta    = sd(beta, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  modelwise_loo_summary,
  file.path(results_dir_assoc_modelwise, "modelwise_leave_one_cohort_out_summary.csv"),
  row.names = FALSE
)


# --- plot heatmap ---

png(
  filename = file.path(results_dir_assoc_modelwise, "heatmap_assoc_modelwise.png"),
  width = 1800,
  height = 1800,
  res = 300
)

assoc_modelwise_est_plot <- assoc_modelwise_est %>%
  mutate(
    epi_model = factor(epi_model, levels = rev(sort(unique(epi_model))))
  )

ggplot(assoc_modelwise_est_plot, aes(x = brain_model, y = epi_model, fill = beta)) +
  geom_tile() +
  scale_fill_gradient2(
    mid = "white",
    midpoint = 0,
    limits = c(-0.1, 0.1),
    oob = scales::squish
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Brain age model",
    y = "Epigenetic age model",
    fill = "Beta",
    title = "Modelwise brain age - epigenetic age associations"
  )

dev.off()

# --- full analysis log ---

sink(log_file_modelwise)

cat("Association modelwise multilevel meta-analysis\n")
cat("==============================================\n\n")
cat("Run ID:", run_id_modelwise, "\n")
cat("Output directory:", results_dir_assoc_modelwise, "\n\n")

cat("Modelwise model summary\n")
cat("=======================\n")
print(summary(mod2_res_mv_modelwise))

cat("\n\nModelwise estimates\n")
cat("===================\n")
print(assoc_modelwise_est, n = Inf, width = Inf)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nModelwise leave-one-cohort-out results\n")
cat("=====================================\n")
print(loo_cohort_modelwise, n = Inf, width = Inf)

cat("\n\nModelwise leave-one-cohort-out summary\n")
cat("=====================================\n")
print(modelwise_loo_summary, n = Inf, width = Inf)

sink()

###### assoc age global (brain) ------
# --- set up ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_assoc_global_ageint_brain <- file.path(results_dir, paste0("assoc_model2_PAR_global_ageint_brain_", run_id))

dir.create(results_dir_assoc_global_ageint_brain, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_assoc_global_ageint_brain, "analysis_log_global_ageint.txt")

# --- data prep ---

# center age 
mod2_eff2 <- mod2_eff2 %>%
  mutate(
    mean_age_brain_c = mean_age_brain - mean(mean_age_brain, na.rm = TRUE)
  )

# --- meta ---

mod2_res_mv_global_ageint_brain <- rma.mv(
  yi = RLM_Estimate_scaled,
  V  = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ mean_age_brain_c,
  random = list(
    ~ 1 | cohort/timepoint, 
    ~ 1 | model_combi  
  ),
  data   = mod2_eff2,
  method = "REML",
  verbose = TRUE
)

summary(mod2_res_mv_global_ageint_brain)

saveRDS(
  mod2_res_mv_global_ageint_brain,
  file.path(results_dir_assoc_global_ageint_brain, "model_assoc_mv_global_ageint_brain.rds")
)

capture.output(
  summary(mod2_res_mv_global_ageint_brain),
  file = file.path(results_dir_assoc_global_ageint_brain, "global_ageint_brain_model_summary.txt")
)

# --- heterogeneity ---
I2 <- calc_I2_mv(
  mod2_res_mv_global_ageint_brain,
  component_names = c("cohort", "timepoint_within_cohort", "model_combi")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_assoc_global_ageint_brain, "heterogeneity_I2.csv"), row.names = FALSE)

# --- extract estimates ---
assoc_global_ageint_brain_est <- tibble(
  term = names(coef(mod2_res_mv_global_ageint_brain)),
  estimate = as.numeric(coef(mod2_res_mv_global_ageint_brain)),
  se = mod2_res_mv_global_ageint_brain$se,
  ci.lb = mod2_res_mv_global_ageint_brain$ci.lb,
  ci.ub = mod2_res_mv_global_ageint_brain$ci.ub,
  pval = mod2_res_mv_global_ageint_brain$pval
)

write.csv(
  assoc_global_ageint_brain_est,
  file.path(results_dir_assoc_global_ageint_brain, "age_moderator_results.csv"),
  row.names = FALSE
)

# --- plotting (should be adapted) ---

# plotting 

newdat <- data.frame(
  mean_age_brain_c = seq(
    min(mod2_eff2$mean_age_brain_c, na.rm = TRUE),
    max(mod2_eff2$mean_age_brain_c, na.rm = TRUE),
    length.out = 100
  )
)

pred <- predict(mod2_res_mv_global_ageint_brain, newmods = newdat$mean_age_brain_c)

newdat$pred  <- pred$pred
newdat$ci.lb <- pred$ci.lb
newdat$ci.ub <- pred$ci.ub

newdat$mean_age_brain <- newdat$mean_age_brain_c + mean(mod2_eff2$mean_age_brain, na.rm = TRUE)

png(
  filename = file.path(results_dir_assoc_global_ageint_brain, "association_global_meanage_meta_regression.png"),
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
    mods = ~ mean_age_brain_c,
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


loo_assoc_ageint_brain_summary <- loo_assoc_ageint %>%
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
    results_dir_assoc_global_ageint_brain,
    "leave_one_cohort_out_age_interaction.csv"
  ),
  row.names = FALSE
)

png(
  filename = file.path(
    results_dir_assoc_global_ageint_brain,
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
    yintercept = coef(mod2_res_mv_global_ageint_brain)[2],
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
cat("Output directory:", results_dir_assoc_global_ageint_brain, "\n\n")

cat("Model summary\n")
cat("=============\n")
print(summary(mod2_res_mv_global_ageint_brain))

cat("\n\nGlobal age interaction estimates\n")
cat("================================\n")
print(assoc_global_ageint_brain_est)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_assoc_ageint, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_assoc_ageint_brain_summary, n = Inf, width = Inf)

cat("\n\nSession info\n")
cat("============\n")
print(sessionInfo())

sink()

###### assoc age global (epi) ------

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_assoc_global_ageint_epi <- file.path(results_dir, paste0("assoc_model2_PAR_global_ageint_epi_", run_id))

dir.create(results_dir_assoc_global_ageint_epi, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_assoc_global_ageint_epi, "analysis_log_global_ageint.txt")

# --- data prep ---

# center age 
mod2_eff2 <- mod2_eff2 %>%
  mutate(
    mean_age_epi_c = mean_age_epi - mean(mean_age_epi, na.rm = TRUE)
  )

# --- meta ---

mod2_res_mv_global_ageint_epi <- rma.mv(
  yi = RLM_Estimate_scaled,
  V  = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ mean_age_epi_c,
  random = list(
    ~ 1 | cohort/timepoint, 
    ~ 1 | model_combi  
  ),
  data   = mod2_eff2,
  method = "REML",
  verbose = TRUE
)

summary(mod2_res_mv_global_ageint_epi)

saveRDS(
  mod2_res_mv_global_ageint_epi,
  file.path(results_dir_assoc_global_ageint_epi, "model_assoc_mv_global_ageint_epi.rds")
)

capture.output(
  summary(mod2_res_mv_global_ageint_epi),
  file = file.path(results_dir_assoc_global_ageint_epi, "global_ageint_epi_model_summary.txt")
)

# --- heterogeneity ---
I2 <- calc_I2_mv(
  mod2_res_mv_global_ageint_epi,
  component_names = c("cohort", "timepoint_within_cohort", "model_combi")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_assoc_global_ageint_epi, "heterogeneity_I2.csv"), row.names = FALSE)

# --- extract estimates ---
assoc_global_ageint_epi_est <- tibble(
  term = names(coef(mod2_res_mv_global_ageint_epi)),
  estimate = as.numeric(coef(mod2_res_mv_global_ageint_epi)),
  se = mod2_res_mv_global_ageint_epi$se,
  ci.lb = mod2_res_mv_global_ageint_epi$ci.lb,
  ci.ub = mod2_res_mv_global_ageint_epi$ci.ub,
  pval = mod2_res_mv_global_ageint_epi$pval
)

write.csv(
  assoc_global_ageint_epi_est,
  file.path(results_dir_assoc_global_ageint_epi, "age_moderator_results.csv"),
  row.names = FALSE
)

# --- plotting (should be adapted) ---

# plotting 

newdat <- data.frame(
  mean_age_epi_c = seq(
    min(mod2_eff2$mean_age_epi_c, na.rm = TRUE),
    max(mod2_eff2$mean_age_epi_c, na.rm = TRUE),
    length.out = 100
  )
)

pred <- predict(mod2_res_mv_global_ageint_epi, newmods = newdat$mean_age_epi_c)

newdat$pred  <- pred$pred
newdat$ci.lb <- pred$ci.lb
newdat$ci.ub <- pred$ci.ub

newdat$mean_age_epi <- newdat$mean_age_epi_c + mean(mod2_eff2$mean_age_epi, na.rm = TRUE)

png(
  filename = file.path(results_dir_assoc_global_ageint_epi, "association_global_meanage_meta_regression.png"),
  width = 1800, height = 1200, res = 300
)

ggplot() +
  geom_point(
    data = mod2_eff2,
    aes(x = mean_age_epi, y = RLM_Estimate_scaled, colour = cohort),
    alpha = .5
  ) +
  geom_line(
    data = newdat,
    aes(x = mean_age_epi, y = pred),
    colour = "black",
    linewidth = 1
  ) +
  geom_ribbon(
    data = newdat,
    aes(x = mean_age_epi, ymin = ci.lb, ymax = ci.ub),
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
    mods = ~ mean_age_epi_c,
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

loo_assoc_ageint_epi_summary <- loo_assoc_ageint %>%
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
    results_dir_assoc_global_ageint_epi,
    "leave_one_cohort_out_age_interaction.csv"
  ),
  row.names = FALSE
)

png(
  filename = file.path(
    results_dir_assoc_global_ageint_epi,
    "plot_assoc_meanage_epi_loo.png"
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
    yintercept = coef(mod2_res_mv_global_ageint_epi)[2],
    linetype = "dashed"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "Estimated epi age slope",
    title = "Leave-one-cohort-out sensitivity of epi age effect on associations"
  ) +
  theme_minimal()

dev.off()


# --- full analysis log ---

sink(log_file)

cat("Association global age interaction meta-analysis: epi mean age\n")
cat("===============================================================\n\n")

cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_assoc_global_ageint_epi, "\n\n")

cat("Model summary\n")
cat("=============\n")
print(summary(mod2_res_mv_global_ageint_epi))

cat("\n\nGlobal age interaction estimates\n")
cat("================================\n")
print(assoc_global_ageint_epi_est)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_assoc_ageint, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_assoc_ageint_epi_summary, n = Inf, width = Inf)

sink()

##### assoc age modelwise - brain mean age -----

# --- set up ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_assoc_age_brain_modelwise <- file.path(
  results_dir,
  paste0("assoc_model2_PAR_modelwise_ageint_brain_", run_id)
)

dir.create(results_dir_assoc_age_brain_modelwise, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  results_dir_assoc_age_brain_modelwise,
  "analysis_log_modelwise_ageint_brain.txt"
)


# --- data prep ---

mod2_eff2_age_brain <- mod2_eff2 %>%
  mutate(
    model_combi = droplevels(factor(model_combi)),
    mean_age_brain_c = mean_age_brain - mean(mean_age_brain, na.rm = TRUE)
  )


# --- meta ---

mod2_res_mv_pairs_ageint_brain <- rma.mv(
  yi = RLM_Estimate_scaled,
  V  = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ model_combi + model_combi:mean_age_brain_c - 1,
  random = ~ 1 | cohort/timepoint,
  data = mod2_eff2_age_brain,
  method = "REML"
)

summary(mod2_res_mv_pairs_ageint_brain)

saveRDS(
  mod2_res_mv_pairs_ageint_brain,
  file.path(
    results_dir_assoc_age_brain_modelwise,
    "model_assoc_mv_modelwise_ageint_brain.rds"
  )
)

capture.output(
  summary(mod2_res_mv_pairs_ageint_brain),
  file = file.path(
    results_dir_assoc_age_brain_modelwise,
    "model_summary_modelwise_ageint_brain.txt"
  )
)

# --- heterogeneity ---
I2 <- calc_I2_mv(
  mod2_res_mv_pairs_ageint_brain,
  component_names = c("cohort", "timepoint_within_cohort")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_assoc_age_brain_modelwise, "heterogeneity_I2.csv"), row.names = FALSE)



# --- extract estimates ---

pair_age_est_brain <- tibble(
  term  = rownames(mod2_res_mv_pairs_ageint_brain$b),
  b     = as.numeric(mod2_res_mv_pairs_ageint_brain$b),
  se    = mod2_res_mv_pairs_ageint_brain$se,
  ci.lb = mod2_res_mv_pairs_ageint_brain$ci.lb,
  ci.ub = mod2_res_mv_pairs_ageint_brain$ci.ub,
  pval  = mod2_res_mv_pairs_ageint_brain$pval
)

pair_age_est_brain_parsed <- pair_age_est_brain %>%
  mutate(
    term = str_remove(term, "^model_combi"),
    
    param_type = if_else(
      str_detect(term, ":mean_age_brain_c"),
      "age_slope",
      "intercept"
    ),
    
    pair_term = str_remove(term, ":mean_age_brain_c")
  ) %>%
  separate(
    pair_term,
    into = c("brain_model", "epi_model"),
    sep = "\\.",
    remove = FALSE
  )

# --- effect sizes ---
X <- model.matrix(mod2_res_mv_pairs_ageint_brain)

model_combi_cols <- colnames(X)[
  str_detect(colnames(X), "^model_combi") &
    !str_detect(colnames(X), ":")
]

effect_sizes_per_pair_brain <- tibble(
  pair_term = str_remove(model_combi_cols, "^model_combi"),
  k = colSums(X[, model_combi_cols, drop = FALSE] != 0)
) %>%
  separate(
    pair_term,
    into = c("brain_model", "epi_model"),
    sep = "\\.",
    remove = FALSE
  )

pair_age_est_brain_parsed <- pair_age_est_brain_parsed %>%
  left_join(effect_sizes_per_pair_brain, by = c("brain_model", "epi_model"))

pair_intercepts_brain_age <- pair_age_est_brain_parsed %>%
  filter(param_type == "intercept")

pair_slopes_brain_age <- pair_age_est_brain_parsed %>%
  filter(param_type == "age_slope")

write.csv(
  pair_age_est_brain_parsed,
  file.path(
    results_dir_assoc_age_brain_modelwise,
    "associations_pairwise_modelwise_ageint_brain_all.csv"
  ),
  row.names = FALSE
)

write.csv(
  pair_intercepts_brain_age,
  file.path(
    results_dir_assoc_age_brain_modelwise,
    "associations_pairwise_intercepts_mean_age_brain.csv"
  ),
  row.names = FALSE
)

write.csv(
  pair_slopes_brain_age,
  file.path(
    results_dir_assoc_age_brain_modelwise,
    "associations_pairwise_age_slopes_mean_age_brain.csv"
  ),
  row.names = FALSE
)



# --- heatmap ---

gen1_order <- sort(c(
  "ZhangEN", "ZhangBLUP", "Wu", "skinHorvath", "PedBE",
  "PCBrainAge", "Horvath2013", "Hannum", "CorticalClock",
  "cAge", "AltumAge"
))

not_gen1_order <- sort(c(
  "PCGrimAge", "PhenoAge", "DamAge", "AdaptAge",
  "DNAmTL", "DunedinPACE"
))

epi_order_gen <- c(gen1_order, not_gen1_order)

pair_slopes_brain_age <- pair_slopes_brain_age %>%
  mutate(
    epi_model = factor(epi_model, levels = rev(epi_order_gen))
  )

png(
  file.path(
    results_dir_assoc_age_brain_modelwise,
    "heatmap_pairwise_age_slopes_mean_age_brain.png"
  ),
  width = 1800,
  height = 1800,
  res = 300
)

ggplot(
  pair_slopes_brain_age,
  aes(
    x = brain_model,
    y = epi_model,
    fill = b
  )
) +
  geom_tile() +
  geom_hline(
    yintercept = length(epi_order_gen) - length(gen1_order) + 0.5,
    linetype = "dashed"
  ) +
  scale_fill_gradient2(
    mid = "white",
    midpoint = 0,
    limits = c(-0.016, 0.016),
    oob = scales::squish
  ) +
  theme_minimal() +
  labs(
    title = "Age moderation of brain-epi associations: brain mean age",
    fill = "Slope"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

dev.off()


# --- leave-one-cohort-out ---

# slightly different version here as previously failed for modelwise analyses but with the age interaction it seems to be fine 

# --- leave-one-cohort-out ---

cohorts <- unique(mod2_eff2_age_brain$cohort)

loo_pairs_age_brain <- lapply(cohorts, function(cc) {
  
  message("Running LOO modelwise age interaction: removing ", cc)
  
  tmp <- mod2_eff2_age_brain %>%
    filter(cohort != cc)
  
  fit <- tryCatch(
    rma.mv(
      yi = RLM_Estimate_scaled,
      V  = RLM_SE_coeftest_HC3_scaled^2,
      mods = ~ model_combi + model_combi:mean_age_brain_c - 1,
      random = ~ 1 | cohort/timepoint,
      data = tmp,
      method = "REML"
    ),
    error = function(e) {
      message("LOO failed for cohort: ", cc)
      message("Error: ", e$message)
      return(e)
    }
  )
  
  if (inherits(fit, "error")) {
    return(
      tibble(
        cohort_removed = cc,
        term = NA_character_,
        b = NA_real_,
        se = NA_real_,
        ci.lb = NA_real_,
        ci.ub = NA_real_,
        pval = NA_real_,
        converged = FALSE,
        error_message = fit$message
      )
    )
  }
  
  tibble(
    cohort_removed = cc,
    term = rownames(fit$b),
    b = as.numeric(fit$b),
    se = fit$se,
    ci.lb = fit$ci.lb,
    ci.ub = fit$ci.ub,
    pval = fit$pval,
    converged = TRUE,
    error_message = NA_character_
  )
}) %>%
  bind_rows()

# check failures
loo_pairs_age_brain %>%
  filter(!converged) %>%
  distinct(cohort_removed, error_message)

# cohorts <- unique(mod2_eff2_age_brain$cohort)
# 
# loo_pairs_age_brain <- lapply(cohorts, function(cc) {
#   
#   tmp <- mod2_eff2_age_brain %>%
#     filter(cohort != cc)
#   
#   fit <- rma.mv(
#     yi = RLM_Estimate_scaled,
#     V  = RLM_SE_coeftest_HC3_scaled^2,
#     mods = ~ model_combi + model_combi:mean_age_brain_c - 1,
#     random = ~ 1 | cohort/timepoint,
#     data = tmp,
#     method = "REML"
#   )
#   
#   tibble(
#     cohort_removed = cc,
#     term = rownames(fit$b),
#     b = as.numeric(fit$b),
#     se = fit$se,
#     ci.lb = fit$ci.lb,
#     ci.ub = fit$ci.ub,
#     pval = fit$pval
#   )
# }) %>%
#   bind_rows()

full_pairs_age_brain <- tibble(
  term = rownames(mod2_res_mv_pairs_ageint_brain$b),
  b_full = as.numeric(mod2_res_mv_pairs_ageint_brain$b)
)

loo_compare_brain_age <- loo_pairs_age_brain %>%
  left_join(full_pairs_age_brain, by = "term") %>%
  mutate(diff = b - b_full)

loo_summary_brain_age <- loo_pairs_age_brain %>%
  filter(converged) %>%
  mutate(
    term_clean = str_remove(term, "^model_combi"),
    param_type = if_else(
      str_detect(term_clean, ":mean_age_brain_c"),
      "age_slope",
      "intercept"
    ),
    pair_term = str_remove(term_clean, ":mean_age_brain_c")
  ) %>%
  separate(
    pair_term,
    into = c("brain_model", "epi_model"),
    sep = "\\.",
    remove = FALSE
  ) %>%
  group_by(brain_model, epi_model, param_type) %>%
  summarise(
    min_b = min(b, na.rm = TRUE),
    max_b = max(b, na.rm = TRUE),
    range_b = max_b - min_b,
    sd_b = sd(b, na.rm = TRUE),
    min_pval = min(pval, na.rm = TRUE),
    max_pval = max(pval, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  loo_pairs_age_brain,
  file.path(
    results_dir_assoc_age_brain_modelwise,
    "modelwise_leave_one_cohort_out_ageint_brain.csv"
  ),
  row.names = FALSE
)

write.csv(
  loo_summary_brain_age,
  file.path(
    results_dir_assoc_age_brain_modelwise,
    "modelwise_leave_one_cohort_out_ageint_brain_summary.csv"
  ),
  row.names = FALSE
)


# --- log ---

sink(log_file)

cat("Association modelwise age interaction meta-analysis: brain mean age\n")
cat("===================================================================\n\n")
cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_assoc_age_brain_modelwise, "\n\n")

cat("Model summary\n")
cat("=============\n")
print(summary(mod2_res_mv_pairs_ageint_brain))

cat("\n\nPairwise intercepts\n")
cat("===================\n")
print(pair_intercepts_brain_age, n = Inf, width = Inf)

cat("\n\nPairwise age slopes\n")
cat("===================\n")
print(pair_slopes_brain_age, n = Inf, width = Inf)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_summary_brain_age, n = Inf, width = Inf)

cat("\n\nSession info\n")
cat("============\n")
print(sessionInfo())

sink()

##### assoc age modelwise - epi mean age -----

# --- set up ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_assoc_age_epi_modelwise <- file.path(
  results_dir,
  paste0("assoc_model2_PAR_modelwise_ageint_epi_", run_id)
)

dir.create(results_dir_assoc_age_epi_modelwise, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  results_dir_assoc_age_epi_modelwise,
  "analysis_log_modelwise_ageint_epi.txt"
)


# --- data prep ---

mod2_eff2_age_epi <- mod2_eff2 %>%
  mutate(
    model_combi = droplevels(factor(model_combi)),
    mean_age_epi_c = mean_age_epi - mean(mean_age_epi, na.rm = TRUE)
  )


# --- meta ---

mod2_res_mv_pairs_ageint_epi <- rma.mv(
  yi = RLM_Estimate_scaled,
  V  = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ model_combi + model_combi:mean_age_epi_c - 1,
  random = ~ 1 | cohort/timepoint,
  data = mod2_eff2_age_epi,
  method = "REML"
)

summary(mod2_res_mv_pairs_ageint_epi)

saveRDS(
  mod2_res_mv_pairs_ageint_epi,
  file.path(
    results_dir_assoc_age_epi_modelwise,
    "model_assoc_mv_modelwise_ageint_epi.rds"
  )
)

capture.output(
  summary(mod2_res_mv_pairs_ageint_epi),
  file = file.path(
    results_dir_assoc_age_epi_modelwise,
    "model_summary_modelwise_ageint_epi.txt"
  )
)

# --- heterogeneity ---
I2 <- calc_I2_mv(
  mod2_res_mv_pairs_ageint_epi,
  component_names = c("cohort", "timepoint_within_cohort")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_assoc_age_epi_modelwise, "heterogeneity_I2.csv"), row.names = FALSE)

# --- extract estimates ---

pair_age_est_epi <- tibble(
  term  = rownames(mod2_res_mv_pairs_ageint_epi$b),
  b     = as.numeric(mod2_res_mv_pairs_ageint_epi$b),
  se    = mod2_res_mv_pairs_ageint_epi$se,
  ci.lb = mod2_res_mv_pairs_ageint_epi$ci.lb,
  ci.ub = mod2_res_mv_pairs_ageint_epi$ci.ub,
  pval  = mod2_res_mv_pairs_ageint_epi$pval
)

pair_age_est_epi_parsed <- pair_age_est_epi %>%
  mutate(
    term = str_remove(term, "^model_combi"),
    
    param_type = if_else(
      str_detect(term, ":mean_age_epi_c"),
      "age_slope",
      "intercept"
    ),
    
    pair_term = str_remove(term, ":mean_age_epi_c")
  ) %>%
  separate(
    pair_term,
    into = c("brain_model", "epi_model"),
    sep = "\\.",
    remove = FALSE
  )

# --- effect sizes ---
X <- model.matrix(mod2_res_mv_pairs_ageint_epi)

model_combi_cols <- colnames(X)[
  str_detect(colnames(X), "^model_combi") &
    !str_detect(colnames(X), ":")
]

effect_sizes_per_pair_epi<- tibble(
  pair_term = str_remove(model_combi_cols, "^model_combi"),
  k = colSums(X[, model_combi_cols, drop = FALSE] != 0)
) %>%
  separate(
    pair_term,
    into = c("brain_model", "epi_model"),
    sep = "\\.",
    remove = FALSE
  )

pair_age_est_epi_parsed <- pair_age_est_epi_parsed %>%
  left_join(effect_sizes_per_pair_epi, by = c("brain_model", "epi_model"))

pair_intercepts_epi_age <- pair_age_est_epi_parsed %>%
  filter(param_type == "intercept")

pair_slopes_epi_age <- pair_age_est_epi_parsed %>%
  filter(param_type == "age_slope")

write.csv(
  pair_age_est_epi_parsed,
  file.path(
    results_dir_assoc_age_epi_modelwise,
    "associations_pairwise_modelwise_ageint_epi_all.csv"
  ),
  row.names = FALSE
)

write.csv(
  pair_intercepts_epi_age,
  file.path(
    results_dir_assoc_age_epi_modelwise,
    "associations_pairwise_intercepts_mean_age_epi.csv"
  ),
  row.names = FALSE
)

write.csv(
  pair_slopes_epi_age,
  file.path(
    results_dir_assoc_age_epi_modelwise,
    "associations_pairwise_age_slopes_mean_age_epi.csv"
  ),
  row.names = FALSE
)



# --- heatmap ---

gen1_order <- sort(c(
  "ZhangEN", "ZhangBLUP", "Wu", "skinHorvath", "PedBE",
  "PCBrainAge", "Horvath2013", "Hannum", "CorticalClock",
  "cAge", "AltumAge"
))

not_gen1_order <- sort(c(
  "PCGrimAge", "PhenoAge", "DamAge", "AdaptAge",
  "DNAmTL", "DunedinPACE"
))

epi_order_gen <- c(gen1_order, not_gen1_order)

pair_slopes_epi_age <- pair_slopes_epi_age %>%
  mutate(
    epi_model = factor(epi_model, levels = rev(epi_order_gen))
  )

png(
  file.path(
    results_dir_assoc_age_epi_modelwise,
    "heatmap_pairwise_age_slopes_mean_age_epi.png"
  ),
  width = 1800,
  height = 1800,
  res = 300
)

ggplot(
  pair_slopes_epi_age,
  aes(
    x = brain_model,
    y = epi_model,
    fill = b
  )
) +
  geom_tile() +
  geom_hline(
    yintercept = length(epi_order_gen) - length(gen1_order) + 0.5,
    linetype = "dashed"
  ) +
  scale_fill_gradient2(
    mid = "white",
    midpoint = 0,
    limits = c(-0.016, 0.016),
    oob = scales::squish
  ) +
  theme_minimal() +
  labs(
    title = "Age moderation of brain-epi associations: epi mean age",
    fill = "Slope"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

dev.off()


# --- leave-one-cohort-out ---

# slightly different version here as previously failed for modelwise analyses but with the age interaction it seems to be fine 

# --- leave-one-cohort-out ---

cohorts <- unique(mod2_eff2_age_epi$cohort)

loo_pairs_age_epi <- lapply(cohorts, function(cc) {
  
  message("Running LOO modelwise age interaction: removing ", cc)
  
  tmp <- mod2_eff2_age_epi %>%
    filter(cohort != cc)
  
  fit <- tryCatch(
    rma.mv(
      yi = RLM_Estimate_scaled,
      V  = RLM_SE_coeftest_HC3_scaled^2,
      mods = ~ model_combi + model_combi:mean_age_epi_c - 1,
      random = ~ 1 | cohort/timepoint,
      data = tmp,
      method = "REML"
    ),
    error = function(e) {
      message("LOO failed for cohort: ", cc)
      message("Error: ", e$message)
      return(e)
    }
  )
  
  if (inherits(fit, "error")) {
    return(
      tibble(
        cohort_removed = cc,
        term = NA_character_,
        b = NA_real_,
        se = NA_real_,
        ci.lb = NA_real_,
        ci.ub = NA_real_,
        pval = NA_real_,
        converged = FALSE,
        error_message = fit$message
      )
    )
  }
  
  tibble(
    cohort_removed = cc,
    term = rownames(fit$b),
    b = as.numeric(fit$b),
    se = fit$se,
    ci.lb = fit$ci.lb,
    ci.ub = fit$ci.ub,
    pval = fit$pval,
    converged = TRUE,
    error_message = NA_character_
  )
}) %>%
  bind_rows()

# check failures
loo_pairs_age_epi %>%
  filter(!converged) %>%
  distinct(cohort_removed, error_message)

# cohorts <- unique(mod2_eff2_age_epi$cohort)
# 
# loo_pairs_age_epi <- lapply(cohorts, function(cc) {
#   
#   tmp <- mod2_eff2_age_epi %>%
#     filter(cohort != cc)
#   
#   fit <- rma.mv(
#     yi = RLM_Estimate_scaled,
#     V  = RLM_SE_coeftest_HC3_scaled^2,
#     mods = ~ model_combi + model_combi:mean_age_epi_c - 1,
#     random = ~ 1 | cohort/timepoint,
#     data = tmp,
#     method = "REML"
#   )
#   
#   tibble(
#     cohort_removed = cc,
#     term = rownames(fit$b),
#     b = as.numeric(fit$b),
#     se = fit$se,
#     ci.lb = fit$ci.lb,
#     ci.ub = fit$ci.ub,
#     pval = fit$pval
#   )
# }) %>%
#   bind_rows()

full_pairs_age_epi <- tibble(
  term = rownames(mod2_res_mv_pairs_ageint_epi$b),
  b_full = as.numeric(mod2_res_mv_pairs_ageint_epi$b)
)

loo_compare_epi_age <- loo_pairs_age_epi %>%
  left_join(full_pairs_age_epi, by = "term") %>%
  mutate(diff = b - b_full)

loo_summary_epi_age <- loo_pairs_age_epi %>%
  filter(converged) %>%
  mutate(
    term_clean = str_remove(term, "^model_combi"),
    param_type = if_else(
      str_detect(term_clean, ":mean_age_epi_c"),
      "age_slope",
      "intercept"
    ),
    pair_term = str_remove(term_clean, ":mean_age_epi_c")
  ) %>%
  separate(
    pair_term,
    into = c("brain_model", "epi_model"),
    sep = "\\.",
    remove = FALSE
  ) %>%
  group_by(brain_model, epi_model, param_type) %>%
  summarise(
    min_b = min(b, na.rm = TRUE),
    max_b = max(b, na.rm = TRUE),
    range_b = max_b - min_b,
    sd_b = sd(b, na.rm = TRUE),
    min_pval = min(pval, na.rm = TRUE),
    max_pval = max(pval, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  loo_pairs_age_epi,
  file.path(
    results_dir_assoc_age_epi_modelwise,
    "modelwise_leave_one_cohort_out_ageint_epi.csv"
  ),
  row.names = FALSE
)

write.csv(
  loo_summary_epi_age,
  file.path(
    results_dir_assoc_age_epi_modelwise,
    "modelwise_leave_one_cohort_out_ageint_epi_summary.csv"
  ),
  row.names = FALSE
)


# --- log ---

sink(log_file)

cat("Association modelwise age interaction meta-analysis: epi mean age\n")
cat("===================================================================\n\n")
cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_assoc_age_epi_modelwise, "\n\n")

cat("Model summary\n")
cat("=============\n")
print(summary(mod2_res_mv_pairs_ageint_epi))

cat("\n\nPairwise intercepts\n")
cat("===================\n")
print(pair_intercepts_epi_age, n = Inf, width = Inf)

cat("\n\nPairwise age slopes\n")
cat("===================\n")
print(pair_slopes_epi_age, n = Inf, width = Inf)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_summary_epi_age, n = Inf, width = Inf)

sink()

#### age overlap ----

##### perf brain ----

# the below does the following
# assesses overlap as no / yes including boundaries i.e., if training is 0-16 and sample is 16-18 that would count as overlap 
# assessed % overlap, here an overlap as in the above example is counted as 0%, % is defined as overlap / study age range --> “How much of the study population is covered by the training population?”

perf_brain_sub <- perf_brain_sub %>%
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

###### brain wMAE overlap_brain global ------

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_brain_ageoverlap_binary <- file.path(
  results_dir,
  paste0("perf_brain_wMAE_noDunedin_global_ageoverlap_binary_", run_id)
)

dir.create(results_dir_brain_ageoverlap_binary, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_brain_ageoverlap_binary, "analysis_log_ageoverlap_binary.txt")


# --- meta ---

brain_wMAE_mv_ageoverlap_binary <- rma.mv(
  yi = wMAE_test,
  V  = wMAE_SE_boot^2,
  mods = ~ overlap_brain,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = perf_brain_sub,
  method = "REML"
)

summary(brain_wMAE_mv_ageoverlap_binary)

# --- heterogeneity ---
I2 <- calc_I2_mv(
  brain_wMAE_mv_ageoverlap_binary,
  component_names = c("cohort", "timepoint_within_cohort", "model")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_brain_ageoverlap_binary, "heterogeneity_I2.csv"), row.names = FALSE)

# --- predicted values ---

newdat_ageoverlap_binary <- data.frame(
  overlap_brain = levels(factor(perf_brain_sub$overlap_brain))
)

newdat_ageoverlap_binary <- newdat_ageoverlap_binary %>%
  mutate(
    overlap_brain = factor(
      overlap_brain,
      levels = levels(factor(perf_brain_sub$overlap_brain))
    )
  )

X_new_ageoverlap_binary <- model.matrix(
  ~ overlap_brain,
  data = newdat_ageoverlap_binary
)

pred_ageoverlap_binary <- predict(
  brain_wMAE_mv_ageoverlap_binary,
  newmods = X_new_ageoverlap_binary[, -1, drop = FALSE]
)

newdat_ageoverlap_binary$pred  <- pred_ageoverlap_binary$pred
newdat_ageoverlap_binary$ci.lb <- pred_ageoverlap_binary$ci.lb
newdat_ageoverlap_binary$ci.ub <- pred_ageoverlap_binary$ci.ub


# --- plot ---

png(
  filename = file.path(
    results_dir_brain_ageoverlap_binary,
    "plot_ageoverlap_binary_wMAE_meta.png"
  ),
  width = 1600,
  height = 1200,
  res = 300
)

ggplot() +
  geom_jitter(
    data = perf_brain_sub,
    aes(
      x = overlap_brain,
      y = wMAE_test
    ),
    width = 0.12,
    alpha = 0.25
  ) +
  geom_point(
    data = newdat_ageoverlap_binary,
    aes(
      x = overlap_brain,
      y = pred
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
  labs(
    x = "Brain age training overlap",
    y = "Predicted wMAE",
    title = "Meta-regression of brain age performance by training age overlap"
  ) +
  theme_minimal(base_size = 14)

dev.off()


# --- save outputs ---

saveRDS(
  brain_wMAE_mv_ageoverlap_binary,
  file.path(
    results_dir_brain_ageoverlap_binary,
    "model_brain_wMAE_mv_ageoverlap_binary.rds"
  )
)

write.csv(
  newdat_ageoverlap_binary,
  file.path(
    results_dir_brain_ageoverlap_binary,
    "ageoverlap_binary_predicted_values.csv"
  ),
  row.names = FALSE
)

capture.output(
  summary(brain_wMAE_mv_ageoverlap_binary),
  file = file.path(
    results_dir_brain_ageoverlap_binary,
    "model_summary.txt"
  )
)

# --- LOO cohort ---

cohorts <- unique(perf_brain_sub$cohort)

loo_ageoverlap_binary <- lapply(cohorts, function(cc) {
  
  tmp <- perf_brain_sub %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = wMAE_test,
    V  = wMAE_SE_boot^2,
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
    intercept = coef(fit)[1],
    intercept_se = sqrt(vcov(fit)[1, 1]),
    overlap_effect = coef(fit)[2],
    overlap_effect_se = sqrt(vcov(fit)[2, 2])
  )
}) %>%
  bind_rows()

loo_ageoverlap_binary_summary <- tibble(
  min_overlap_effect = min(loo_ageoverlap_binary$overlap_effect, na.rm = TRUE),
  max_overlap_effect = max(loo_ageoverlap_binary$overlap_effect, na.rm = TRUE),
  range_overlap_effect = max_overlap_effect - min_overlap_effect,
  sd_overlap_effect = sd(loo_ageoverlap_binary$overlap_effect, na.rm = TRUE)
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
    x = reorder(cohort_removed, overlap_effect),
    y = overlap_effect
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = overlap_effect - 1.96 * overlap_effect_se,
      ymax = overlap_effect + 1.96 * overlap_effect_se
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept = coef(brain_wMAE_mv_ageoverlap_binary)[2],
    linetype = "dashed"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "Estimated overlap effect",
    title = "Leave-one-cohort-out sensitivity of age-overlap effect"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- full analysis log ---

sink(log_file)

cat("Brain wMAE multilevel meta-analysis - age overlap binary moderator\n")
cat("==================================================================\n\n")

cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_brain_ageoverlap_binary, "\n\n")

cat("Model summary\n")
cat("=============\n")
print(summary(brain_wMAE_mv_ageoverlap_binary))

cat("\n\nPredicted values\n")
cat("================\n")
print(newdat_ageoverlap_binary)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

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

###### brain wMAE overlap_perc_brain global ------

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_brain_ageoverlap_perc <- file.path(
  results_dir,
  paste0("perf_brain_wMAE_noDunedin_global_ageoverlap_perc_", run_id)
)

dir.create(results_dir_brain_ageoverlap_perc, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  results_dir_brain_ageoverlap_perc,
  "analysis_log_ageoverlap_perc.txt"
)

# --- meta ---

brain_wMAE_mv_ageoverlap_perc <- rma.mv(
  yi = wMAE_test,
  V  = wMAE_SE_boot^2,
  mods = ~ overlap_perc_brain,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = perf_brain_sub,
  method = "REML"
)

summary(brain_wMAE_mv_ageoverlap_perc)

# --- heterogeneity ---
I2 <- calc_I2_mv(
  brain_wMAE_mv_ageoverlap_perc,
  component_names = c("cohort", "timepoint_within_cohort", "model")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_brain_ageoverlap_perc, "heterogeneity_I2.csv"), row.names = FALSE)

# --- extracting estimates ---
global_est <- tibble(
  term = names(coef(brain_wMAE_mv_ageoverlap_perc)),
  estimate = as.numeric(coef(brain_wMAE_mv_ageoverlap_perc)),
  se = brain_wMAE_mv_ageoverlap_perc$se,
  ci.lb = brain_wMAE_mv_ageoverlap_perc$ci.lb,
  ci.ub = brain_wMAE_mv_ageoverlap_perc$ci.ub,
  pval = brain_wMAE_mv_ageoverlap_perc$pval
)

write.csv(
  global_est,
  file.path(results_dir_brain_ageoverlap_perc, "global_estimates.csv"),
  row.names = FALSE
)

# --- predicted values ---
newdat_ageoverlap_perc <- data.frame(
  overlap_perc_brain = seq(
    min(perf_brain_sub$overlap_perc_brain, na.rm = TRUE),
    max(perf_brain_sub$overlap_perc_brain, na.rm = TRUE),
    length.out = 100
  )
) 

pred_ageoverlap_perc <- predict(
  brain_wMAE_mv_ageoverlap_perc,
  newmods = newdat_ageoverlap_perc$overlap_perc_brain
)

newdat_ageoverlap_perc$pred  <- pred_ageoverlap_perc$pred
newdat_ageoverlap_perc$ci.lb <- pred_ageoverlap_perc$ci.lb
newdat_ageoverlap_perc$ci.ub <- pred_ageoverlap_perc$ci.ub


# --- cohort-timepoint means ---

cohort_time_means_ageoverlap_perc <- perf_brain_sub %>%
  group_by(cohort, timepoint) %>%
  summarise(
    overlap_perc_brain =
      mean(overlap_perc_brain, na.rm = TRUE),
    
    mean_wMAE =
      mean(wMAE_test, na.rm = TRUE),
    
    .groups = "drop"
  )


# --- plot ---

png(
  filename = file.path(
    results_dir_brain_ageoverlap_perc,
    "plot_ageoverlap_perc_wMAE_meta.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot() +
  
  geom_point(
    data = perf_brain_sub,
    aes(
      x = overlap_perc_brain,
      y = wMAE_test
    ),
    alpha = 0.25
  ) +
  
  geom_point(
    data = cohort_time_means_ageoverlap_perc,
    aes(
      x = overlap_perc_brain,
      y = mean_wMAE
    ),
    size = 3
  ) +
  
  geom_line(
    data = newdat_ageoverlap_perc,
    aes(
      x = overlap_perc_brain,
      y = pred
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
  
  labs(
    x = "Brain age training overlap (%)",
    y = "Predicted wMAE",
    title = "Meta-regression of brain age performance by proportional age overlap"
  ) +
  
  theme_minimal(base_size = 14)

dev.off()


# --- save outputs ---

saveRDS(
  brain_wMAE_mv_ageoverlap_perc,
  file.path(
    results_dir_brain_ageoverlap_perc,
    "model_brain_wMAE_mv_ageoverlap_perc.rds"
  )
)

write.csv(
  newdat_ageoverlap_perc,
  file.path(
    results_dir_brain_ageoverlap_perc,
    "ageoverlap_perc_predicted_values.csv"
  ),
  row.names = FALSE
)

capture.output(
  summary(brain_wMAE_mv_ageoverlap_perc),
  file = file.path(
    results_dir_brain_ageoverlap_perc,
    "model_summary.txt"
  )
)


# --- LOO cohort ---

cohorts <- unique(perf_brain_sub$cohort)

loo_ageoverlap_perc <- lapply(cohorts, function(cc) {
  
  tmp <- perf_brain_sub %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = wMAE_test,
    V  = wMAE_SE_boot^2,
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
    
    intercept = coef(fit)[1],
    intercept_se = sqrt(vcov(fit)[1, 1]),
    
    overlap_perc_effect = coef(fit)[2],
    overlap_perc_effect_se = sqrt(vcov(fit)[2, 2])
  )
}) %>%
  bind_rows()

loo_ageoverlap_perc_summary <- tibble(
  min_overlap_perc_effect =
    min(loo_ageoverlap_perc$overlap_perc_effect, na.rm = TRUE),
  
  max_overlap_perc_effect =
    max(loo_ageoverlap_perc$overlap_perc_effect, na.rm = TRUE),
  
  range_overlap_perc_effect =
    max_overlap_perc_effect - min_overlap_perc_effect,
  
  sd_overlap_perc_effect =
    sd(loo_ageoverlap_perc$overlap_perc_effect, na.rm = TRUE)
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
    x = reorder(cohort_removed, overlap_perc_effect),
    y = overlap_perc_effect
  )
) +
  
  geom_point(size = 3) +
  
  geom_errorbar(
    aes(
      ymin =
        overlap_perc_effect -
        1.96 * overlap_perc_effect_se,
      
      ymax =
        overlap_perc_effect +
        1.96 * overlap_perc_effect_se
    ),
    width = 0.15
  ) +
  
  geom_hline(
    yintercept = coef(brain_wMAE_mv_ageoverlap_perc)[2],
    linetype = "dashed"
  ) +
  
  coord_flip() +
  
  labs(
    x = "Cohort removed",
    y = "Estimated overlap percentage effect",
    title = "Leave-one-cohort-out sensitivity of proportional age-overlap effect"
  ) +
  
  theme_minimal(base_size = 14)

dev.off()


# --- full analysis log ---

sink(log_file)

cat("Brain wMAE multilevel meta-analysis - age overlap percentage moderator\n")
cat("======================================================================\n\n")

cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_brain_ageoverlap_perc, "\n\n")

cat("Model summary\n")
cat("=============\n")
print(summary(brain_wMAE_mv_ageoverlap_perc))

cat("\n\nPredicted values\n")
cat("================\n")
print(newdat_ageoverlap_perc)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_ageoverlap_perc, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_ageoverlap_perc_summary, n = Inf, width = Inf)

sink()

###### brain wMAE processing global ------

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_brain_processing <- file.path(
  results_dir,
  paste0("perf_brain_wMAE_noDunedin_global_processing_", run_id)
)

dir.create(results_dir_brain_processing, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  results_dir_brain_processing,
  "analysis_log_processing.txt"
)


# --- data prep ---

perf_brain_sub <- perf_brain_sub %>%
  mutate(
    processing = factor(
      processing,
      levels = c("minimal", "extracted")
    )
  )


# --- meta ---

brain_wMAE_mv_processing <- rma.mv(
  yi = wMAE_test,
  V  = wMAE_SE_boot^2,
  mods = ~ processing,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = perf_brain_sub,
  method = "REML"
)

summary(brain_wMAE_mv_processing)

# --- heterogeneity ---
I2 <- calc_I2_mv(
  brain_wMAE_mv_processing,
  component_names = c("cohort", "timepoint_within_cohort", "model")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_brain_processing, "heterogeneity_I2.csv"), row.names = FALSE)

# --- extracting estimates ---
global_est <- tibble(
  term = names(coef(brain_wMAE_mv_processing)),
  estimate = as.numeric(coef(brain_wMAE_mv_processing)),
  se = brain_wMAE_mv_processing$se,
  ci.lb = brain_wMAE_mv_processing$ci.lb,
  ci.ub = brain_wMAE_mv_processing$ci.ub,
  pval = brain_wMAE_mv_processing$pval
)

write.csv(
  global_est,
  file.path(results_dir_brain_processing, "global_estimates.csv"),
  row.names = FALSE
)


# --- predicted values ---

newdat_processing <- data.frame(
  processing = levels(perf_brain_sub$processing)
)

newdat_processing <- newdat_processing %>%
  mutate(
    processing = factor(
      processing,
      levels = levels(perf_brain_sub$processing)
    )
  )

X_new_processing <- model.matrix(
  ~ processing,
  data = newdat_processing
)

pred_processing <- predict(
  brain_wMAE_mv_processing,
  newmods = X_new_processing[, -1, drop = FALSE]
)

newdat_processing$pred  <- pred_processing$pred
newdat_processing$ci.lb <- pred_processing$ci.lb
newdat_processing$ci.ub <- pred_processing$ci.ub


# --- plot ---

png(
  filename = file.path(
    results_dir_brain_processing,
    "plot_processing_wMAE_meta.png"
  ),
  width = 1600,
  height = 1200,
  res = 300
)

ggplot() +
  geom_jitter(
    data = perf_brain_sub,
    aes(
      x = processing,
      y = wMAE_test
    ),
    width = 0.12,
    alpha = 0.25
  ) +
  geom_point(
    data = newdat_processing,
    aes(
      x = processing,
      y = pred
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
  labs(
    x = "Processing level",
    y = "Predicted wMAE",
    title = "Meta-regression of brain age performance by processing level"
  ) +
  theme_minimal(base_size = 14)

dev.off()


# --- save outputs ---

saveRDS(
  brain_wMAE_mv_processing,
  file.path(
    results_dir_brain_processing,
    "model_brain_wMAE_mv_processing.rds"
  )
)

write.csv(
  newdat_processing,
  file.path(
    results_dir_brain_processing,
    "processing_predicted_values.csv"
  ),
  row.names = FALSE
)

capture.output(
  summary(brain_wMAE_mv_processing),
  file = file.path(
    results_dir_brain_processing,
    "model_summary.txt"
  )
)


# --- LOO cohort ---

cohorts <- unique(perf_brain_sub$cohort)

loo_processing <- lapply(cohorts, function(cc) {
  
  tmp <- perf_brain_sub %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = wMAE_test,
    V  = wMAE_SE_boot^2,
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
    intercept = coef(fit)[1],
    intercept_se = sqrt(vcov(fit)[1, 1]),
    processing_effect = coef(fit)[2],
    processing_effect_se = sqrt(vcov(fit)[2, 2])
  )
}) %>%
  bind_rows()

loo_processing_summary <- tibble(
  min_processing_effect = min(loo_processing$processing_effect, na.rm = TRUE),
  max_processing_effect = max(loo_processing$processing_effect, na.rm = TRUE),
  range_processing_effect = max_processing_effect - min_processing_effect,
  sd_processing_effect = sd(loo_processing$processing_effect, na.rm = TRUE)
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


# --- LOO plot ---

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
    x = reorder(cohort_removed, processing_effect),
    y = processing_effect
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = processing_effect - 1.96 * processing_effect_se,
      ymax = processing_effect + 1.96 * processing_effect_se
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept = coef(brain_wMAE_mv_processing)[2],
    linetype = "dashed"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "Estimated processing effect",
    title = "Leave-one-cohort-out sensitivity of processing effect"
  ) +
  theme_minimal(base_size = 14)

dev.off()


# --- full analysis log ---

sink(log_file)

cat("Brain wMAE multilevel meta-analysis - processing moderator\n")
cat("==========================================================\n\n")

cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_brain_processing, "\n\n")

cat("Model summary\n")
cat("=============\n")
print(summary(brain_wMAE_mv_processing))

cat("\n\nPredicted values\n")
cat("================\n")
print(newdat_processing)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_processing, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_processing_summary, n = Inf, width = Inf)

sink()

##### perf epi ----

# the below does the following
# assesses overlap as no / yes including boundaries i.e., if training is 0-16 and sample is 16-18 that would count as overlap 
# assessed % overlap, here an overlap as in the above example is counted as 0%, % is defined as overlap / study age range --> “How much of the study population is covered by the training population?”

# need to change age for pheno Age as in the current loaded sample it's not accurate 
perf_epi_sub <- perf_epi_sub %>%
  mutate(
    training_age_min_epi = if_else(
      model == "PhenoAge",
      20,
      training_age_min_epi
    ),
    
    training_age_max_epi = if_else(
      model == "PhenoAge",
      85,
      training_age_max_epi
    )
  )

perf_epi_sub <- perf_epi_sub %>%
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

###### epi wMAE overlap_epi global ------

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi_ageoverlap_binary <- file.path(
  results_dir,
  paste0("perf_epi_wMAE_noDunedinDNAmTL_global_ageoverlap_binary_", run_id)
)

dir.create(results_dir_epi_ageoverlap_binary, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_epi_ageoverlap_binary, "analysis_log_ageoverlap_binary.txt")


# --- data prep ---

perf_epi_sub <- perf_epi_sub %>%
  mutate(
    overlap_epi = factor(overlap_epi, levels = c("no", "yes"))
  )

perf_epi_sub_nogest <- perf_epi_sub %>%
  filter(!model %in% c("Bohlin", "Knight", "EPIC")) %>%
  droplevels()

# --- meta ---

epi_wMAE_mv_ageoverlap_binary <- rma.mv(
  yi = log_wMAE,
  V  = log_wMAE_var,
  mods = ~ overlap_epi,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = perf_epi_sub_nogest,
  method = "REML"
)

summary(epi_wMAE_mv_ageoverlap_binary)

# --- heterogeneity ---
I2 <- calc_I2_mv(
  epi_wMAE_mv_ageoverlap_binary,
  component_names = c("cohort", "timepoint_within_cohort", "model")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_epi_ageoverlap_binary, "heterogeneity_I2.csv"), row.names = FALSE)

# --- back-transform model results ---

ageoverlap_binary_est <- tibble(
  term = names(coef(epi_wMAE_mv_ageoverlap_binary)),
  estimate_log = as.numeric(coef(epi_wMAE_mv_ageoverlap_binary)),
  se_log = epi_wMAE_mv_ageoverlap_binary$se,
  ci.lb_log = epi_wMAE_mv_ageoverlap_binary$ci.lb,
  ci.ub_log = epi_wMAE_mv_ageoverlap_binary$ci.ub,
  pval = epi_wMAE_mv_ageoverlap_binary$pval
) %>%
  mutate(
    estimate_exp = exp(estimate_log),
    ci.lb_exp = exp(ci.lb_log),
    ci.ub_exp = exp(ci.ub_log),
    percent_change = (exp(estimate_log) - 1) * 100,
    percent_ci.lb = (exp(ci.lb_log) - 1) * 100,
    percent_ci.ub = (exp(ci.ub_log) - 1) * 100
  )

write.csv(
  ageoverlap_binary_est,
  file.path(results_dir_epi_ageoverlap_binary, "ageoverlap_binary_backtransformed_results.csv"),
  row.names = FALSE
)


# --- predicted values ---

newdat_ageoverlap_binary <- data.frame(
  overlap_epi = levels(perf_epi_sub_nogest$overlap_epi)
) %>%
  mutate(
    overlap_epi = factor(
      overlap_epi,
      levels = levels(perf_epi_sub_nogest$overlap_epi)
    )
  )

X_new_ageoverlap_binary <- model.matrix(
  ~ overlap_epi,
  data = newdat_ageoverlap_binary
)

pred_ageoverlap_binary <- predict(
  epi_wMAE_mv_ageoverlap_binary,
  newmods = X_new_ageoverlap_binary[, -1, drop = FALSE]
)

newdat_ageoverlap_binary <- newdat_ageoverlap_binary %>%
  mutate(
    log_pred = pred_ageoverlap_binary$pred,
    log_ci.lb = pred_ageoverlap_binary$ci.lb,
    log_ci.ub = pred_ageoverlap_binary$ci.ub,
    
    pred = exp(log_pred),
    ci.lb = exp(log_ci.lb),
    ci.ub = exp(log_ci.ub)
  )

write.csv(
  newdat_ageoverlap_binary,
  file.path(results_dir_epi_ageoverlap_binary, "ageoverlap_binary_predicted_values.csv"),
  row.names = FALSE
)


# --- plot ---

png(
  filename = file.path(
    results_dir_epi_ageoverlap_binary,
    "plot_ageoverlap_binary_wMAE_meta_backtransformed.png"
  ),
  width = 1600,
  height = 1200,
  res = 300
)

ggplot() +
  geom_jitter(
    data = perf_epi_sub_nogest,
    aes(x = overlap_epi, y = wMAE_test),
    width = 0.12,
    alpha = 0.25
  ) +
  geom_point(
    data = newdat_ageoverlap_binary,
    aes(x = overlap_epi, y = pred),
    size = 4
  ) +
  geom_errorbar(
    data = newdat_ageoverlap_binary,
    aes(x = overlap_epi, ymin = ci.lb, ymax = ci.ub),
    width = 0.15
  ) +
  labs(
    x = "Epi clock training age overlap",
    y = "Geometric mean weighted MAE",
    title = "Meta-regression of epi clock performance by training age overlap"
  ) +
  theme_minimal(base_size = 14)

dev.off()


# --- save outputs ---

saveRDS(
  epi_wMAE_mv_ageoverlap_binary,
  file.path(
    results_dir_epi_ageoverlap_binary,
    "model_epi_wMAE_mv_ageoverlap_binary_logscale.rds"
  )
)

capture.output(
  summary(epi_wMAE_mv_ageoverlap_binary),
  file = file.path(
    results_dir_epi_ageoverlap_binary,
    "model_summary_log_scale.txt"
  )
)


# --- LOO cohort ---

cohorts <- unique(perf_epi_sub_nogest$cohort)

loo_ageoverlap_binary <- lapply(cohorts, function(cc) {
  
  tmp <- perf_epi_sub_nogest %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = log_wMAE,
    V  = log_wMAE_var,
    mods = ~ overlap_epi,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model
    ),
    data = tmp,
    method = "REML"
  )
  
  log_intercept <- coef(fit)[1]
  log_overlap_effect <- coef(fit)[2]
  
  log_intercept_se <- sqrt(vcov(fit)[1, 1])
  log_overlap_effect_se <- sqrt(vcov(fit)[2, 2])
  
  tibble(
    cohort_removed = cc,
    
    log_intercept = log_intercept,
    log_intercept_se = log_intercept_se,
    intercept_wMAE = exp(log_intercept),
    intercept_ci.lb = exp(log_intercept - 1.96 * log_intercept_se),
    intercept_ci.ub = exp(log_intercept + 1.96 * log_intercept_se),
    
    log_overlap_effect = log_overlap_effect,
    log_overlap_effect_se = log_overlap_effect_se,
    overlap_ratio = exp(log_overlap_effect),
    overlap_percent_change = (exp(log_overlap_effect) - 1) * 100,
    overlap_percent_ci.lb = (exp(log_overlap_effect - 1.96 * log_overlap_effect_se) - 1) * 100,
    overlap_percent_ci.ub = (exp(log_overlap_effect + 1.96 * log_overlap_effect_se) - 1) * 100
  )
}) %>%
  bind_rows()

loo_ageoverlap_binary_summary <- tibble(
  min_overlap_percent_change = min(loo_ageoverlap_binary$overlap_percent_change, na.rm = TRUE),
  max_overlap_percent_change = max(loo_ageoverlap_binary$overlap_percent_change, na.rm = TRUE),
  range_overlap_percent_change = max_overlap_percent_change - min_overlap_percent_change,
  sd_overlap_percent_change = sd(loo_ageoverlap_binary$overlap_percent_change, na.rm = TRUE),
  
  min_log_overlap_effect = min(loo_ageoverlap_binary$log_overlap_effect, na.rm = TRUE),
  max_log_overlap_effect = max(loo_ageoverlap_binary$log_overlap_effect, na.rm = TRUE),
  range_log_overlap_effect = max_log_overlap_effect - min_log_overlap_effect,
  sd_log_overlap_effect = sd(loo_ageoverlap_binary$log_overlap_effect, na.rm = TRUE)
)

write.csv(
  loo_ageoverlap_binary,
  file.path(results_dir_epi_ageoverlap_binary, "leave_one_cohort_out_ageoverlap_binary.csv"),
  row.names = FALSE
)

write.csv(
  loo_ageoverlap_binary_summary,
  file.path(results_dir_epi_ageoverlap_binary, "leave_one_cohort_out_ageoverlap_binary_summary.csv"),
  row.names = FALSE
)


# --- LOO plot ---

full_overlap_percent_change <- (
  exp(coef(epi_wMAE_mv_ageoverlap_binary)[2]) - 1
) * 100

png(
  filename = file.path(
    results_dir_epi_ageoverlap_binary,
    "plot_ageoverlap_binary_loo_percent_change.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_ageoverlap_binary,
  aes(
    x = reorder(cohort_removed, overlap_percent_change),
    y = overlap_percent_change
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = overlap_percent_ci.lb,
      ymax = overlap_percent_ci.ub
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept = full_overlap_percent_change,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 0,
    colour = "grey60"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "% change in geometric mean wMAE",
    title = "Leave-one-cohort-out sensitivity of age-overlap effect"
  ) +
  theme_minimal(base_size = 14)

dev.off()


# --- full analysis log ---

sink(log_file)

cat("Epi wMAE multilevel meta-analysis - age overlap binary moderator, log-transformed\n")
cat("================================================================================\n\n")

cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_epi_ageoverlap_binary, "\n\n")

cat("Model summary, log scale\n")
cat("========================\n")
print(summary(epi_wMAE_mv_ageoverlap_binary))

cat("\n\nBack-transformed model results\n")
cat("==============================\n")
print(ageoverlap_binary_est, n = Inf, width = Inf)

cat("\n\nPredicted values\n")
cat("================\n")
print(newdat_ageoverlap_binary)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_ageoverlap_binary, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_ageoverlap_binary_summary, n = Inf, width = Inf)

sink()

###### epi wMAE overlap_perc_epi global ------

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi_ageoverlap_perc <- file.path(
  results_dir,
  paste0("perf_epi_wMAE_noDunedinDNAmTL_global_ageoverlap_perc_", run_id)
)

dir.create(results_dir_epi_ageoverlap_perc, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  results_dir_epi_ageoverlap_perc,
  "analysis_log_ageoverlap_perc.txt"
)


# --- data prep ---

perf_epi_sub_nogest <- perf_epi_sub %>%
  filter(!model %in% c("Bohlin", "Knight", "EPIC")) %>%
  droplevels()


# --- meta ---

epi_wMAE_mv_ageoverlap_perc <- rma.mv(
  yi = log_wMAE,
  V  = log_wMAE_var,
  mods = ~ overlap_perc_epi,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = perf_epi_sub_nogest,
  method = "REML"
)

summary(epi_wMAE_mv_ageoverlap_perc)

# --- heterogeneity ---
I2 <- calc_I2_mv(
  epi_wMAE_mv_ageoverlap_perc,
  component_names = c("cohort", "timepoint_within_cohort", "model")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_epi_ageoverlap_perc, "heterogeneity_I2.csv"), row.names = FALSE)

# --- back-transform model results ---

ageoverlap_perc_est <- tibble(
  term = names(coef(epi_wMAE_mv_ageoverlap_perc)),
  estimate_log = as.numeric(coef(epi_wMAE_mv_ageoverlap_perc)),
  se_log = epi_wMAE_mv_ageoverlap_perc$se,
  ci.lb_log = epi_wMAE_mv_ageoverlap_perc$ci.lb,
  ci.ub_log = epi_wMAE_mv_ageoverlap_perc$ci.ub,
  pval = epi_wMAE_mv_ageoverlap_perc$pval
) %>%
  mutate(
    estimate_exp = exp(estimate_log),
    ci.lb_exp = exp(ci.lb_log),
    ci.ub_exp = exp(ci.ub_log),
    percent_change = (exp(estimate_log) - 1) * 100,
    percent_ci.lb = (exp(ci.lb_log) - 1) * 100,
    percent_ci.ub = (exp(ci.ub_log) - 1) * 100
  )

write.csv(
  ageoverlap_perc_est,
  file.path(results_dir_epi_ageoverlap_perc, "ageoverlap_perc_backtransformed_results.csv"),
  row.names = FALSE
)


# --- predicted values ---

newdat_ageoverlap_perc <- data.frame(
  overlap_perc_epi = seq(
    min(perf_epi_sub_nogest$overlap_perc_epi, na.rm = TRUE),
    max(perf_epi_sub_nogest$overlap_perc_epi, na.rm = TRUE),
    length.out = 100
  )
)

pred_ageoverlap_perc <- predict(
  epi_wMAE_mv_ageoverlap_perc,
  newmods = newdat_ageoverlap_perc$overlap_perc_epi
)

newdat_ageoverlap_perc <- newdat_ageoverlap_perc %>%
  mutate(
    log_pred = pred_ageoverlap_perc$pred,
    log_ci.lb = pred_ageoverlap_perc$ci.lb,
    log_ci.ub = pred_ageoverlap_perc$ci.ub,
    pred = exp(log_pred),
    ci.lb = exp(log_ci.lb),
    ci.ub = exp(log_ci.ub)
  )

write.csv(
  newdat_ageoverlap_perc,
  file.path(results_dir_epi_ageoverlap_perc, "ageoverlap_perc_predicted_values.csv"),
  row.names = FALSE
)


# --- plot ---

png(
  filename = file.path(
    results_dir_epi_ageoverlap_perc,
    "plot_ageoverlap_perc_wMAE_meta_backtransformed.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot() +
  geom_point(
    data = perf_epi_sub_nogest,
    aes(x = overlap_perc_epi, y = wMAE_test),
    alpha = 0.25
  ) +
  geom_line(
    data = newdat_ageoverlap_perc,
    aes(x = overlap_perc_epi, y = pred),
    linewidth = 1
  ) +
  geom_ribbon(
    data = newdat_ageoverlap_perc,
    aes(x = overlap_perc_epi, ymin = ci.lb, ymax = ci.ub),
    alpha = 0.15
  ) +
  labs(
    x = "Epi clock training overlap (%)",
    y = "Geometric mean weighted MAE",
    title = "Meta-regression of epi clock performance by proportional age overlap"
  ) +
  theme_minimal(base_size = 14)

dev.off()


# --- save outputs ---

saveRDS(
  epi_wMAE_mv_ageoverlap_perc,
  file.path(
    results_dir_epi_ageoverlap_perc,
    "model_epi_wMAE_mv_ageoverlap_perc_logscale.rds"
  )
)

capture.output(
  summary(epi_wMAE_mv_ageoverlap_perc),
  file = file.path(
    results_dir_epi_ageoverlap_perc,
    "model_summary_log_scale.txt"
  )
)


# --- LOO cohort ---

cohorts <- unique(perf_epi_sub_nogest$cohort)

loo_ageoverlap_perc <- lapply(cohorts, function(cc) {
  
  tmp <- perf_epi_sub_nogest %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = log_wMAE,
    V  = log_wMAE_var,
    mods = ~ overlap_perc_epi,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model
    ),
    data = tmp,
    method = "REML"
  )
  
  log_intercept <- coef(fit)[1]
  log_overlap_perc_effect <- coef(fit)[2]
  
  log_intercept_se <- sqrt(vcov(fit)[1, 1])
  log_overlap_perc_effect_se <- sqrt(vcov(fit)[2, 2])
  
  tibble(
    cohort_removed = cc,
    
    log_intercept = log_intercept,
    log_intercept_se = log_intercept_se,
    intercept_wMAE = exp(log_intercept),
    intercept_ci.lb = exp(log_intercept - 1.96 * log_intercept_se),
    intercept_ci.ub = exp(log_intercept + 1.96 * log_intercept_se),
    
    log_overlap_perc_effect = log_overlap_perc_effect,
    log_overlap_perc_effect_se = log_overlap_perc_effect_se,
    overlap_perc_ratio_per_1pct = exp(log_overlap_perc_effect),
    overlap_perc_percent_change_per_1pct =
      (exp(log_overlap_perc_effect) - 1) * 100,
    overlap_perc_percent_ci.lb =
      (exp(log_overlap_perc_effect - 1.96 * log_overlap_perc_effect_se) - 1) * 100,
    overlap_perc_percent_ci.ub =
      (exp(log_overlap_perc_effect + 1.96 * log_overlap_perc_effect_se) - 1) * 100
  )
}) %>%
  bind_rows()

loo_ageoverlap_perc_summary <- tibble(
  min_overlap_perc_percent_change_per_1pct =
    min(loo_ageoverlap_perc$overlap_perc_percent_change_per_1pct, na.rm = TRUE),
  max_overlap_perc_percent_change_per_1pct =
    max(loo_ageoverlap_perc$overlap_perc_percent_change_per_1pct, na.rm = TRUE),
  range_overlap_perc_percent_change_per_1pct =
    max_overlap_perc_percent_change_per_1pct -
    min_overlap_perc_percent_change_per_1pct,
  sd_overlap_perc_percent_change_per_1pct =
    sd(loo_ageoverlap_perc$overlap_perc_percent_change_per_1pct, na.rm = TRUE),
  
  min_log_overlap_perc_effect =
    min(loo_ageoverlap_perc$log_overlap_perc_effect, na.rm = TRUE),
  max_log_overlap_perc_effect =
    max(loo_ageoverlap_perc$log_overlap_perc_effect, na.rm = TRUE),
  range_log_overlap_perc_effect =
    max_log_overlap_perc_effect - min_log_overlap_perc_effect,
  sd_log_overlap_perc_effect =
    sd(loo_ageoverlap_perc$log_overlap_perc_effect, na.rm = TRUE)
)

write.csv(
  loo_ageoverlap_perc,
  file.path(results_dir_epi_ageoverlap_perc, "leave_one_cohort_out_ageoverlap_perc.csv"),
  row.names = FALSE
)

write.csv(
  loo_ageoverlap_perc_summary,
  file.path(results_dir_epi_ageoverlap_perc, "leave_one_cohort_out_ageoverlap_perc_summary.csv"),
  row.names = FALSE
)


# --- LOO plot ---

full_overlap_perc_percent_change <- (
  exp(coef(epi_wMAE_mv_ageoverlap_perc)[2]) - 1
) * 100

png(
  filename = file.path(
    results_dir_epi_ageoverlap_perc,
    "plot_ageoverlap_perc_loo_percent_change.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_ageoverlap_perc,
  aes(
    x = reorder(cohort_removed, overlap_perc_percent_change_per_1pct),
    y = overlap_perc_percent_change_per_1pct
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = overlap_perc_percent_ci.lb,
      ymax = overlap_perc_percent_ci.ub
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept = full_overlap_perc_percent_change,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 0,
    colour = "grey60"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "% change in geometric mean wMAE per 1% overlap",
    title = "Leave-one-cohort-out sensitivity of proportional age-overlap effect"
  ) +
  theme_minimal(base_size = 14)

dev.off()


# --- full analysis log ---

sink(log_file)

cat("Epi wMAE multilevel meta-analysis - age overlap percentage moderator, log-transformed\n")
cat("====================================================================================\n\n")

cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_epi_ageoverlap_perc, "\n\n")

cat("Model summary, log scale\n")
cat("========================\n")
print(summary(epi_wMAE_mv_ageoverlap_perc))

cat("\n\nBack-transformed model results\n")
cat("==============================\n")
print(ageoverlap_perc_est, width = Inf)

cat("\n\nPredicted values\n")
cat("================\n")
print(newdat_ageoverlap_perc)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_ageoverlap_perc, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_ageoverlap_perc_summary,width = Inf)

sink()

###### epi wMAE tissue  ------

# create tissue overlap variable 
perf_epi_sub <- perf_epi_sub %>%
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
  perf_epi_sub$tissue_training,
  perf_epi_sub$tissue_test,
  perf_epi_sub$tissue_overlap
)

###### epi wMAE tissue_overlap global ------

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi_tissue_overlap <- file.path(
  results_dir,
  paste0("perf_epi_wMAE_noDunedinDNAmTL_global_tissue_overlap_", run_id)
)

dir.create(results_dir_epi_tissue_overlap, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_epi_tissue_overlap, "analysis_log_tissue_overlap.txt")


# --- data prep ---

perf_epi_sub <- perf_epi_sub %>%
  mutate(
    tissue_overlap = factor(tissue_overlap, levels = c("no", "yes"))
  )

# --- meta ---

epi_wMAE_mv_tissue_overlap <- rma.mv(
  yi = log_wMAE,
  V  = log_wMAE_var,
  mods = ~ tissue_overlap,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = perf_epi_sub,
  method = "REML"
)

summary(epi_wMAE_mv_tissue_overlap)

# --- heterogeneity ---
I2 <- calc_I2_mv(
  epi_wMAE_mv_tissue_overlap,
  component_names = c("cohort", "timepoint_within_cohort", "model")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_epi_tissue_overlap, "heterogeneity_I2.csv"), row.names = FALSE)


# --- back-transform model results ---

tissue_overlap_est <- tibble(
  term = names(coef(epi_wMAE_mv_tissue_overlap)),
  estimate_log = as.numeric(coef(epi_wMAE_mv_tissue_overlap)),
  se_log = epi_wMAE_mv_tissue_overlap$se,
  ci.lb_log = epi_wMAE_mv_tissue_overlap$ci.lb,
  ci.ub_log = epi_wMAE_mv_tissue_overlap$ci.ub,
  pval = epi_wMAE_mv_tissue_overlap$pval
) %>%
  mutate(
    estimate_exp = exp(estimate_log),
    ci.lb_exp = exp(ci.lb_log),
    ci.ub_exp = exp(ci.ub_log),
    percent_change = (exp(estimate_log) - 1) * 100,
    percent_ci.lb = (exp(ci.lb_log) - 1) * 100,
    percent_ci.ub = (exp(ci.ub_log) - 1) * 100
  )

write.csv(
  tissue_overlap_est,
  file.path(results_dir_epi_tissue_overlap, "tissue_overlap_backtransformed_results.csv"),
  row.names = FALSE
)


# --- predicted values ---

newdat_tissue_overlap <- data.frame(
  tissue_overlap = levels(perf_epi_sub$tissue_overlap)
) %>%
  mutate(
    tissue_overlap = factor(
      tissue_overlap,
      levels = levels(perf_epi_sub$tissue_overlap)
    )
  )

X_new_tissue_overlap <- model.matrix(
  ~ tissue_overlap,
  data = newdat_tissue_overlap
)

pred_tissue_overlap <- predict(
  epi_wMAE_mv_tissue_overlap,
  newmods = X_new_tissue_overlap[, -1, drop = FALSE]
)

newdat_tissue_overlap <- newdat_tissue_overlap %>%
  mutate(
    log_pred = pred_tissue_overlap$pred,
    log_ci.lb = pred_tissue_overlap$ci.lb,
    log_ci.ub = pred_tissue_overlap$ci.ub,
    
    pred = exp(log_pred),
    ci.lb = exp(log_ci.lb),
    ci.ub = exp(log_ci.ub)
  )

write.csv(
  newdat_tissue_overlap,
  file.path(results_dir_epi_tissue_overlap, "tissue_overlap_predicted_values.csv"),
  row.names = FALSE
)


# --- plot ---

png(
  filename = file.path(
    results_dir_epi_tissue_overlap,
    "plot_tissue_overlap_wMAE_meta_backtransformed.png"
  ),
  width = 1600,
  height = 1200,
  res = 300
)

ggplot() +
  geom_jitter(
    data = perf_epi_sub,
    aes(x = tissue_overlap, y = wMAE_test),
    width = 0.12,
    alpha = 0.25
  ) +
  geom_point(
    data = newdat_tissue_overlap,
    aes(x = tissue_overlap, y = pred),
    size = 4
  ) +
  geom_errorbar(
    data = newdat_tissue_overlap,
    aes(x = tissue_overlap, ymin = ci.lb, ymax = ci.ub),
    width = 0.15
  ) +
  labs(
    x = "Epi clock tissue overlap",
    y = "Geometric mean weighted MAE",
    title = "Meta-regression of epi clock performance by tissue overlap"
  ) +
  theme_minimal(base_size = 14)

dev.off()


# --- save outputs ---

saveRDS(
  epi_wMAE_mv_tissue_overlap,
  file.path(
    results_dir_epi_tissue_overlap,
    "model_epi_wMAE_mv_tissue_overlap_logscale.rds"
  )
)

capture.output(
  summary(epi_wMAE_mv_tissue_overlap),
  file = file.path(
    results_dir_epi_tissue_overlap,
    "model_summary_log_scale.txt"
  )
)


# --- LOO cohort ---

cohorts <- unique(perf_epi_sub$cohort)

loo_tissue_overlap <- lapply(cohorts, function(cc) {
  
  tmp <- perf_epi_sub %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = log_wMAE,
    V  = log_wMAE_var,
    mods = ~ tissue_overlap,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model
    ),
    data = tmp,
    method = "REML"
  )
  
  log_intercept <- coef(fit)[1]
  log_overlap_effect <- coef(fit)[2]
  
  log_intercept_se <- sqrt(vcov(fit)[1, 1])
  log_overlap_effect_se <- sqrt(vcov(fit)[2, 2])
  
  tibble(
    cohort_removed = cc,
    
    log_intercept = log_intercept,
    log_intercept_se = log_intercept_se,
    intercept_wMAE = exp(log_intercept),
    intercept_ci.lb = exp(log_intercept - 1.96 * log_intercept_se),
    intercept_ci.ub = exp(log_intercept + 1.96 * log_intercept_se),
    
    log_overlap_effect = log_overlap_effect,
    log_overlap_effect_se = log_overlap_effect_se,
    overlap_ratio = exp(log_overlap_effect),
    overlap_percent_change = (exp(log_overlap_effect) - 1) * 100,
    overlap_percent_ci.lb = (exp(log_overlap_effect - 1.96 * log_overlap_effect_se) - 1) * 100,
    overlap_percent_ci.ub = (exp(log_overlap_effect + 1.96 * log_overlap_effect_se) - 1) * 100
  )
}) %>%
  bind_rows()

loo_tissue_overlap_summary <- tibble(
  min_overlap_percent_change = min(loo_tissue_overlap$overlap_percent_change, na.rm = TRUE),
  max_overlap_percent_change = max(loo_tissue_overlap$overlap_percent_change, na.rm = TRUE),
  range_overlap_percent_change = max_overlap_percent_change - min_overlap_percent_change,
  sd_overlap_percent_change = sd(loo_tissue_overlap$overlap_percent_change, na.rm = TRUE),
  
  min_log_overlap_effect = min(loo_tissue_overlap$log_overlap_effect, na.rm = TRUE),
  max_log_overlap_effect = max(loo_tissue_overlap$log_overlap_effect, na.rm = TRUE),
  range_log_overlap_effect = max_log_overlap_effect - min_log_overlap_effect,
  sd_log_overlap_effect = sd(loo_tissue_overlap$log_overlap_effect, na.rm = TRUE)
)

write.csv(
  loo_tissue_overlap,
  file.path(results_dir_epi_tissue_overlap, "leave_one_cohort_out_tissue_overlap.csv"),
  row.names = FALSE
)

write.csv(
  loo_tissue_overlap_summary,
  file.path(results_dir_epi_tissue_overlap, "leave_one_cohort_out_tissue_overlap_summary.csv"),
  row.names = FALSE
)


# --- LOO plot ---

full_overlap_percent_change <- (
  exp(coef(epi_wMAE_mv_tissue_overlap)[2]) - 1
) * 100

png(
  filename = file.path(
    results_dir_epi_tissue_overlap,
    "plot_tissue_overlap_loo_percent_change.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_tissue_overlap,
  aes(
    x = reorder(cohort_removed, overlap_percent_change),
    y = overlap_percent_change
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = overlap_percent_ci.lb,
      ymax = overlap_percent_ci.ub
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept = full_overlap_percent_change,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 0,
    colour = "grey60"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "% change in geometric mean wMAE",
    title = "Leave-one-cohort-out sensitivity of tissue_overlap effect"
  ) +
  theme_minimal(base_size = 14)

dev.off()


# --- full analysis log ---

sink(log_file)

cat("Epi wMAE multilevel meta-analysis - tissue overlap, log-transformed\n")
cat("================================================================================\n\n")

cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_epi_tissue_overlap, "\n\n")

cat("Model summary, log scale\n")
cat("========================\n")
print(summary(epi_wMAE_mv_tissue_overlap))

cat("\n\nBack-transformed model results\n")
cat("==============================\n")
print(tissue_overlap_est, n = Inf, width = Inf)

cat("\n\nPredicted values\n")
cat("================\n")
print(newdat_tissue_overlap)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_tissue_overlap, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_tissue_overlap_summary, n = Inf, width = Inf)

sink()

###### epi wMAE tissue_training global ------

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi_tissue_training <- file.path(
  results_dir,
  paste0("perf_epi_wMAE_noDunedinDNAmTL_global_tissue_training_", run_id)
)

dir.create(results_dir_epi_tissue_training, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_epi_tissue_training, "analysis_log_tissue_training.txt")


# --- data prep ---

perf_epi_sub <- perf_epi_sub %>%
  mutate(
    tissue_training = factor(tissue_training)
  )

# --- meta ---

epi_wMAE_mv_tissue_training <- rma.mv(
  yi = log_wMAE,
  V  = log_wMAE_var,
  mods = ~ tissue_training -1 ,#taking out intercept so that there is no reference category 
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = perf_epi_sub,
  method = "REML"
)

summary(epi_wMAE_mv_tissue_training)

# --- heterogeneity ---
I2 <- calc_I2_mv(
  epi_wMAE_mv_tissue_training,
  component_names = c("cohort", "timepoint_within_cohort", "model")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_epi_tissue_training, "heterogeneity_I2.csv"), row.names = FALSE)

# --- back-transform model results ---

tissue_training_est <- tibble(
  term = names(coef(epi_wMAE_mv_tissue_training)),
  estimate_log = as.numeric(coef(epi_wMAE_mv_tissue_training)),
  se_log = epi_wMAE_mv_tissue_training$se,
  ci.lb_log = epi_wMAE_mv_tissue_training$ci.lb,
  ci.ub_log = epi_wMAE_mv_tissue_training$ci.ub,
  pval = epi_wMAE_mv_tissue_training$pval
) %>%
  mutate(
    estimate_exp = exp(estimate_log),
    ci.lb_exp = exp(ci.lb_log),
    ci.ub_exp = exp(ci.ub_log),
    percent_change = (exp(estimate_log) - 1) * 100,
    percent_ci.lb = (exp(ci.lb_log) - 1) * 100,
    percent_ci.ub = (exp(ci.ub_log) - 1) * 100
  )

write.csv(
  tissue_training_est,
  file.path(results_dir_epi_tissue_training, "tissue_training_backtransformed_results.csv"),
  row.names = FALSE
)


# --- predicted values ---

newdat_tissue_training <- data.frame(
  tissue_training = levels(perf_epi_sub$tissue_training)
) %>%
  mutate(
    tissue_training = factor(
      tissue_training,
      levels = levels(perf_epi_sub$tissue_training)
    )
  )

X_new_tissue_training <- model.matrix(
  ~ tissue_training - 1,
  data = newdat_tissue_training
)

pred_tissue_training <- predict(
  epi_wMAE_mv_tissue_training,
  newmods = X_new_tissue_training
)

newdat_tissue_training <- newdat_tissue_training %>%
  mutate(
    log_pred = pred_tissue_training$pred,
    log_ci.lb = pred_tissue_training$ci.lb,
    log_ci.ub = pred_tissue_training$ci.ub,
    
    pred = exp(log_pred),
    ci.lb = exp(log_ci.lb),
    ci.ub = exp(log_ci.ub)
  )

write.csv(
  newdat_tissue_training,
  file.path(results_dir_epi_tissue_training, "tissue_training_predicted_values.csv"),
  row.names = FALSE
)


# --- plot ---

png(
  filename = file.path(
    results_dir_epi_tissue_training,
    "plot_tissue_training_wMAE_meta_backtransformed.png"
  ),
  width = 1600,
  height = 1200,
  res = 300
)

ggplot() +
  geom_jitter(
    data = perf_epi_sub,
    aes(x = tissue_training, y = wMAE_test),
    width = 0.12,
    alpha = 0.25
  ) +
  geom_point(
    data = newdat_tissue_training,
    aes(x = tissue_training, y = pred),
    size = 4
  ) +
  geom_errorbar(
    data = newdat_tissue_training,
    aes(x = tissue_training, ymin = ci.lb, ymax = ci.ub),
    width = 0.15
  ) +
  labs(
    x = "Epi clock tissue overlap",
    y = "Geometric mean weighted MAE",
    title = "Meta-regression of epi clock performance by training tissue"
  ) +
  theme_minimal(base_size = 14)

dev.off()


# --- save outputs ---

saveRDS(
  epi_wMAE_mv_tissue_training,
  file.path(
    results_dir_epi_tissue_training,
    "model_epi_wMAE_mv_tissue_training_logscale.rds"
  )
)

capture.output(
  summary(epi_wMAE_mv_tissue_training),
  file = file.path(
    results_dir_epi_tissue_training,
    "model_summary_log_scale.txt"
  )
)


# --- LOO cohort ---

cohorts <- unique(perf_epi_sub$cohort)

loo_tissue_training <- lapply(cohorts, function(cc) {
  
  tmp <- perf_epi_sub %>%
    filter(cohort != cc) %>%
    droplevels()
  
  fit <- rma.mv(
    yi = log_wMAE,
    V  = log_wMAE_var,
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
    log_est = as.numeric(coef(fit)),
    log_se = fit$se,
    log_ci.lb = fit$ci.lb,
    log_ci.ub = fit$ci.ub
  ) %>%
    mutate(
      tissue_training = str_remove(term, "^tissue_training"),
      pred_wMAE = exp(log_est),
      ci.lb = exp(log_ci.lb),
      ci.ub = exp(log_ci.ub)
    )
}) %>%
  bind_rows()


# --- LOO summary by tissue ---

loo_tissue_training_summary <- loo_tissue_training %>%
  group_by(tissue_training) %>%
  summarise(
    min_wMAE = min(pred_wMAE, na.rm = TRUE),
    max_wMAE = max(pred_wMAE, na.rm = TRUE),
    range_wMAE = max_wMAE - min_wMAE,
    sd_wMAE = sd(pred_wMAE, na.rm = TRUE),
    
    min_log_est = min(log_est, na.rm = TRUE),
    max_log_est = max(log_est, na.rm = TRUE),
    range_log_est = max_log_est - min_log_est,
    sd_log_est = sd(log_est, na.rm = TRUE),
    
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
  filename = file.path(
    results_dir_epi_tissue_training,
    "plot_tissue_training_loo.png"
  ),
  width = 2000,
  height = 1400,
  res = 300
)

ggplot(
  loo_tissue_training,
  aes(
    x = reorder(cohort_removed, pred_wMAE),
    y = pred_wMAE
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
    y = "Predicted geometric mean wMAE",
    title = "Leave-one-cohort-out sensitivity by training tissue"
  ) +
  theme_minimal(base_size = 14)

dev.off()


# --- full analysis log ---

sink(log_file)

cat("Epi wMAE multilevel meta-analysis - tissue training, log-transformed\n")
cat("================================================================================\n\n")

cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_epi_tissue_training, "\n\n")

cat("Model summary, log scale\n")
cat("========================\n")
print(summary(epi_wMAE_mv_tissue_training))

cat("\n\nBack-transformed model results\n")
cat("==============================\n")
print(tissue_training_est, n = Inf, width = Inf)

cat("\n\nPredicted values\n")
cat("================\n")
print(newdat_tissue_training)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_tissue_training, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_tissue_training_summary, n = Inf, width = Inf)

sink()

###### epi wMAE tissue_test global ------

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi_tissue_test <- file.path(
  results_dir,
  paste0("perf_epi_wMAE_noDunedinDNAmTL_global_tissue_test_", run_id)
)

dir.create(results_dir_epi_tissue_test, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_epi_tissue_test, "analysis_log_tissue_test.txt")


# --- data prep ---

perf_epi_sub <- perf_epi_sub %>%
  mutate(
    tissue_test = factor(tissue_test)
  )

# --- meta ---

epi_wMAE_mv_tissue_test <- rma.mv(
  yi = log_wMAE,
  V  = log_wMAE_var,
  mods = ~ tissue_test -1 ,#taking out intercept so that there is no reference category 
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = perf_epi_sub,
  method = "REML"
)

summary(epi_wMAE_mv_tissue_test)

# --- heterogeneity ---
I2 <- calc_I2_mv(
  epi_wMAE_mv_tissue_test,
  component_names = c("cohort", "timepoint_within_cohort", "model")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_epi_tissue_test, "heterogeneity_I2.csv"), row.names = FALSE)

# --- back-transform model results ---

tissue_test_est <- tibble(
  term = names(coef(epi_wMAE_mv_tissue_test)),
  estimate_log = as.numeric(coef(epi_wMAE_mv_tissue_test)),
  se_log = epi_wMAE_mv_tissue_test$se,
  ci.lb_log = epi_wMAE_mv_tissue_test$ci.lb,
  ci.ub_log = epi_wMAE_mv_tissue_test$ci.ub,
  pval = epi_wMAE_mv_tissue_test$pval
) %>%
  mutate(
    estimate_exp = exp(estimate_log),
    ci.lb_exp = exp(ci.lb_log),
    ci.ub_exp = exp(ci.ub_log),
    percent_change = (exp(estimate_log) - 1) * 100,
    percent_ci.lb = (exp(ci.lb_log) - 1) * 100,
    percent_ci.ub = (exp(ci.ub_log) - 1) * 100
  )

write.csv(
  tissue_test_est,
  file.path(results_dir_epi_tissue_test, "tissue_test_backtransformed_results.csv"),
  row.names = FALSE
)


# --- predicted values ---

newdat_tissue_test <- data.frame(
  tissue_test = levels(perf_epi_sub$tissue_test)
) %>%
  mutate(
    tissue_test = factor(
      tissue_test,
      levels = levels(perf_epi_sub$tissue_test)
    )
  )

X_new_tissue_test <- model.matrix(
  ~ tissue_test - 1,
  data = newdat_tissue_test
)

pred_tissue_test <- predict(
  epi_wMAE_mv_tissue_test,
  newmods = X_new_tissue_test
)

newdat_tissue_test <- newdat_tissue_test %>%
  mutate(
    log_pred = pred_tissue_test$pred,
    log_ci.lb = pred_tissue_test$ci.lb,
    log_ci.ub = pred_tissue_test$ci.ub,
    
    pred = exp(log_pred),
    ci.lb = exp(log_ci.lb),
    ci.ub = exp(log_ci.ub)
  )

write.csv(
  newdat_tissue_test,
  file.path(results_dir_epi_tissue_test, "tissue_test_predicted_values.csv"),
  row.names = FALSE
)


# --- plot ---

png(
  filename = file.path(
    results_dir_epi_tissue_test,
    "plot_tissue_test_wMAE_meta_backtransformed.png"
  ),
  width = 1600,
  height = 1200,
  res = 300
)

ggplot() +
  geom_jitter(
    data = perf_epi_sub,
    aes(x = tissue_test, y = wMAE_test),
    width = 0.12,
    alpha = 0.25
  ) +
  geom_point(
    data = newdat_tissue_test,
    aes(x = tissue_test, y = pred),
    size = 4
  ) +
  geom_errorbar(
    data = newdat_tissue_test,
    aes(x = tissue_test, ymin = ci.lb, ymax = ci.ub),
    width = 0.15
  ) +
  labs(
    x = "Epi clock tissue overlap",
    y = "Geometric mean weighted MAE",
    title = "Meta-regression of epi clock performance by test tissue"
  ) +
  theme_minimal(base_size = 14)

dev.off()


# --- save outputs ---

saveRDS(
  epi_wMAE_mv_tissue_test,
  file.path(
    results_dir_epi_tissue_test,
    "model_epi_wMAE_mv_tissue_test_logscale.rds"
  )
)

capture.output(
  summary(epi_wMAE_mv_tissue_test),
  file = file.path(
    results_dir_epi_tissue_test,
    "model_summary_log_scale.txt"
  )
)


# --- LOO cohort ---

cohorts <- unique(perf_epi_sub$cohort)

loo_tissue_test <- lapply(cohorts, function(cc) {
  
  tmp <- perf_epi_sub %>%
    filter(cohort != cc) %>%
    droplevels()
  
  fit <- rma.mv(
    yi = log_wMAE,
    V  = log_wMAE_var,
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
    log_est = as.numeric(coef(fit)),
    log_se = fit$se,
    log_ci.lb = fit$ci.lb,
    log_ci.ub = fit$ci.ub
  ) %>%
    mutate(
      tissue_test = str_remove(term, "^tissue_test"),
      pred_wMAE = exp(log_est),
      ci.lb = exp(log_ci.lb),
      ci.ub = exp(log_ci.ub)
    )
}) %>%
  bind_rows()


# --- LOO summary by tissue ---

loo_tissue_test_summary <- loo_tissue_test %>%
  group_by(tissue_test) %>%
  summarise(
    min_wMAE = min(pred_wMAE, na.rm = TRUE),
    max_wMAE = max(pred_wMAE, na.rm = TRUE),
    range_wMAE = max_wMAE - min_wMAE,
    sd_wMAE = sd(pred_wMAE, na.rm = TRUE),
    
    min_log_est = min(log_est, na.rm = TRUE),
    max_log_est = max(log_est, na.rm = TRUE),
    range_log_est = max_log_est - min_log_est,
    sd_log_est = sd(log_est, na.rm = TRUE),
    
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
  filename = file.path(
    results_dir_epi_tissue_test,
    "plot_tissue_test_loo.png"
  ),
  width = 2000,
  height = 1400,
  res = 300
)

ggplot(
  loo_tissue_test,
  aes(
    x = reorder(cohort_removed, pred_wMAE),
    y = pred_wMAE
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
    y = "Predicted geometric mean wMAE",
    title = "Leave-one-cohort-out sensitivity by test tissue"
  ) +
  theme_minimal(base_size = 14)

dev.off()


# --- full analysis log ---

sink(log_file)

cat("Epi wMAE multilevel meta-analysis - tissue test, log-transformed\n")
cat("================================================================================\n\n")

cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_epi_tissue_test, "\n\n")

cat("Model summary, log scale\n")
cat("========================\n")
print(summary(epi_wMAE_mv_tissue_test))

cat("\n\nBack-transformed model results\n")
cat("==============================\n")
print(tissue_test_est, n = Inf, width = Inf)

cat("\n\nPredicted values\n")
cat("================\n")
print(newdat_tissue_test)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

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

###### epi wMAE generation_epi global ------

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi_generation_epi <- file.path(
  results_dir,
  paste0("perf_epi_wMAE_noDunedinDNAmTL_global_generation_epi_", run_id)
)

dir.create(results_dir_epi_generation_epi, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_epi_generation_epi, "analysis_log_generation_epi.txt")


# --- data prep ---

perf_epi_sub <- perf_epi_sub %>%
  mutate(
    generation_epi = factor(generation_epi)
  )

# --- meta ---

epi_wMAE_mv_generation_epi <- rma.mv(
  yi = log_wMAE,
  V  = log_wMAE_var,
  mods = ~ generation_epi -1 ,#taking out intercept so that there is no reference category 
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = perf_epi_sub,
  method = "REML"
)

summary(epi_wMAE_mv_generation_epi)


# --- back-transform model results ---

generation_epi_est <- tibble(
  term = names(coef(epi_wMAE_mv_generation_epi)),
  estimate_log = as.numeric(coef(epi_wMAE_mv_generation_epi)),
  se_log = epi_wMAE_mv_generation_epi$se,
  ci.lb_log = epi_wMAE_mv_generation_epi$ci.lb,
  ci.ub_log = epi_wMAE_mv_generation_epi$ci.ub,
  pval = epi_wMAE_mv_generation_epi$pval
) %>%
  mutate(
    estimate_exp = exp(estimate_log),
    ci.lb_exp = exp(ci.lb_log),
    ci.ub_exp = exp(ci.ub_log),
    percent_change = (exp(estimate_log) - 1) * 100,
    percent_ci.lb = (exp(ci.lb_log) - 1) * 100,
    percent_ci.ub = (exp(ci.ub_log) - 1) * 100
  )

write.csv(
  generation_epi_est,
  file.path(results_dir_epi_generation_epi, "generation_epi_backtransformed_results.csv"),
  row.names = FALSE
)


# --- predicted values ---

newdat_generation_epi <- data.frame(
  generation_epi = levels(perf_epi_sub$generation_epi)
) %>%
  mutate(
    generation_epi = factor(
      generation_epi,
      levels = levels(perf_epi_sub$generation_epi)
    )
  )

X_new_generation_epi <- model.matrix(
  ~ generation_epi - 1,
  data = newdat_generation_epi
)

pred_generation_epi <- predict(
  epi_wMAE_mv_generation_epi,
  newmods = X_new_generation_epi
)

newdat_generation_epi <- newdat_generation_epi %>%
  mutate(
    log_pred = pred_generation_epi$pred,
    log_ci.lb = pred_generation_epi$ci.lb,
    log_ci.ub = pred_generation_epi$ci.ub,
    
    pred = exp(log_pred),
    ci.lb = exp(log_ci.lb),
    ci.ub = exp(log_ci.ub)
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
    "plot_generation_epi_wMAE_meta_backtransformed.png"
  ),
  width = 1600,
  height = 1200,
  res = 300
)

ggplot() +
  geom_jitter(
    data = perf_epi_sub,
    aes(x = generation_epi, y = wMAE_test),
    width = 0.12,
    alpha = 0.25
  ) +
  geom_point(
    data = newdat_generation_epi,
    aes(x = generation_epi, y = pred),
    size = 4
  ) +
  geom_errorbar(
    data = newdat_generation_epi,
    aes(x = generation_epi, ymin = ci.lb, ymax = ci.ub),
    width = 0.15
  ) +
  labs(
    x = "Epi clock tissue overlap",
    y = "Geometric mean weighted MAE",
    title = "Meta-regression of epi clock performance by training tissue"
  ) +
  theme_minimal(base_size = 14)

dev.off()


# --- save outputs ---

saveRDS(
  epi_wMAE_mv_generation_epi,
  file.path(
    results_dir_epi_generation_epi,
    "model_epi_wMAE_mv_generation_epi_logscale.rds"
  )
)

capture.output(
  summary(epi_wMAE_mv_generation_epi),
  file = file.path(
    results_dir_epi_generation_epi,
    "model_summary_log_scale.txt"
  )
)


# --- LOO cohort ---

cohorts <- unique(perf_epi_sub$cohort)

loo_generation_epi <- lapply(cohorts, function(cc) {
  
  tmp <- perf_epi_sub %>%
    filter(cohort != cc) %>%
    droplevels()
  
  fit <- rma.mv(
    yi = log_wMAE,
    V  = log_wMAE_var,
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
    log_est = as.numeric(coef(fit)),
    log_se = fit$se,
    log_ci.lb = fit$ci.lb,
    log_ci.ub = fit$ci.ub
  ) %>%
    mutate(
      generation_epi = str_remove(term, "^generation_epi"),
      pred_wMAE = exp(log_est),
      ci.lb = exp(log_ci.lb),
      ci.ub = exp(log_ci.ub)
    )
}) %>%
  bind_rows()


# --- LOO summary by tissue ---

loo_generation_epi_summary <- loo_generation_epi %>%
  group_by(generation_epi) %>%
  summarise(
    min_wMAE = min(pred_wMAE, na.rm = TRUE),
    max_wMAE = max(pred_wMAE, na.rm = TRUE),
    range_wMAE = max_wMAE - min_wMAE,
    sd_wMAE = sd(pred_wMAE, na.rm = TRUE),
    
    min_log_est = min(log_est, na.rm = TRUE),
    max_log_est = max(log_est, na.rm = TRUE),
    range_log_est = max_log_est - min_log_est,
    sd_log_est = sd(log_est, na.rm = TRUE),
    
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
  filename = file.path(
    results_dir_epi_generation_epi,
    "plot_generation_epi_loo.png"
  ),
  width = 2000,
  height = 1400,
  res = 300
)

ggplot(
  loo_generation_epi,
  aes(
    x = reorder(cohort_removed, pred_wMAE),
    y = pred_wMAE
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
    y = "Predicted geometric mean wMAE",
    title = "Leave-one-cohort-out sensitivity by training tissue"
  ) +
  theme_minimal(base_size = 14)

dev.off()


# --- full analysis log ---

sink(log_file)

cat("Epi wMAE multilevel meta-analysis - generation, log-transformed\n")
cat("================================================================================\n\n")

cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_epi_generation_epi, "\n\n")

cat("Model summary, log scale\n")
cat("========================\n")
print(summary(epi_wMAE_mv_generation_epi))

cat("\n\nBack-transformed model results\n")
cat("==============================\n")
print(generation_epi_est, n = Inf, width = Inf)

cat("\n\nPredicted values\n")
cat("================\n")
print(newdat_generation_epi)

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

###### epi wMAE age generation_epi ------
# all but DunedinPACNI & DNAmTL

# --- set up ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_epi_generation_epi_ageint <- file.path(results_dir, paste0("perf_epi_wMAE_noDunedin_generation_epi_ageint_", run_id))

dir.create(results_dir_epi_generation_epi_ageint, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_epi_generation_epi_ageint, "analysis_log_generation_epi_ageint.txt")

# --- data prep ---

# center age 
perf_epi_sub_nogest <- perf_epi_sub %>%
  filter(generation_epi != "1_gest") %>%
  droplevels() %>%
  mutate(
    generation_epi = factor(generation_epi),
    mean_age_c = mean_age - mean(mean_age, na.rm = TRUE)
  )

# --- meta ---

epi_wMAE_mv_generation_epi_ageint <- rma.mv(
  yi = log_wMAE,
  V  = log_wMAE_var,
  mods = ~ generation_epi + generation_epi:mean_age_c - 1,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = perf_epi_sub_nogest,
  method = "REML"
)

summary(epi_wMAE_mv_generation_epi_ageint)


generation_epi_age_est <- tibble(
  term = rownames(epi_wMAE_mv_generation_epi_ageint$b),
  estimate_log = as.numeric(epi_wMAE_mv_generation_epi_ageint$b),
  se_log = epi_wMAE_mv_generation_epi_ageint$se,
  ci.lb_log = epi_wMAE_mv_generation_epi_ageint$ci.lb,
  ci.ub_log = epi_wMAE_mv_generation_epi_ageint$ci.ub,
  pval = epi_wMAE_mv_generation_epi_ageint$pval
) %>%
  mutate(
    type = if_else(str_detect(term, "mean_age_c"), "age_slope", "intercept"),
    generation_epi = term %>%
      str_remove("^generation_epi") %>%
      str_remove(":mean_age_c")
  ) %>%
  mutate(
    estimate_exp = exp(estimate_log),
    ci.lb_exp = exp(ci.lb_log),
    ci.ub_exp = exp(ci.ub_log),
    percent_change = (exp(estimate_log) - 1) * 100,
    percent_ci.lb = (exp(ci.lb_log) - 1) * 100,
    percent_ci.ub = (exp(ci.ub_log) - 1) * 100
  )

# For interpretation:
#   
#   rows with type == "intercept" give geometric mean wMAE at the mean age.
# rows with type == "age_slope" give percent change in geometric mean wMAE per one-year increase in age

write.csv(
  generation_epi_age_est,
  file.path(results_dir_epi_generation_epi_ageint, "generation_epi_age_interaction_backtransformed_results.csv"),
  row.names = FALSE
)

saveRDS(
  epi_wMAE_mv_generation_epi_ageint,
  file.path(results_dir_epi_generation_epi_ageint, "model_epi_wMAE_mv_generation_epi_ageint.rds")
)

capture.output(
  summary(epi_wMAE_mv_generation_epi_ageint),
  file = file.path(results_dir_epi_generation_epi_ageint, "model_summary_log_scale.txt")
)


###
# also runthe model once with intercept 
epi_wMAE_mv_generation_epi_ageint_withintercept <- rma.mv(
  yi = log_wMAE,
  V  = log_wMAE_var,
  mods = ~ generation_epi*mean_age_c,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model
  ),
  data = perf_epi_sub_nogest,
  method = "REML"
)

summary(epi_wMAE_mv_generation_epi_ageint_withintercept)

generation_epi_age_est_withintercept <- tibble(
  term = rownames(epi_wMAE_mv_generation_epi_ageint_withintercept$b),
  estimate_log = as.numeric(epi_wMAE_mv_generation_epi_ageint_withintercept$b),
  se_log = epi_wMAE_mv_generation_epi_ageint_withintercept$se,
  ci.lb_log = epi_wMAE_mv_generation_epi_ageint_withintercept$ci.lb,
  ci.ub_log = epi_wMAE_mv_generation_epi_ageint_withintercept$ci.ub,
  pval = epi_wMAE_mv_generation_epi_ageint_withintercept$pval
) %>%
  mutate(
    type = if_else(str_detect(term, "mean_age_c"), "age_slope", "intercept"),
    generation_epi = term %>%
      str_remove("^generation_epi") %>%
      str_remove(":mean_age_c")
  ) %>%
  mutate(
    estimate_exp = exp(estimate_log),
    ci.lb_exp = exp(ci.lb_log),
    ci.ub_exp = exp(ci.ub_log),
    percent_change = (exp(estimate_log) - 1) * 100,
    percent_ci.lb = (exp(ci.lb_log) - 1) * 100,
    percent_ci.ub = (exp(ci.ub_log) - 1) * 100
  )

# For interpretation:
#   
#   rows with type == "intercept" give geometric mean wMAE at the mean age.
# rows with type == "age_slope" give percent change in geometric mean wMAE per one-year increase in age

write.csv(
  generation_epi_age_est,
  file.path(results_dir_epi_generation_epi_ageint, "generation_epi_age_interaction_backtransformed_results_withintercept.csv"),
  row.names = FALSE
)

saveRDS(
  epi_wMAE_mv_generation_epi_ageint_withintercept,
  file.path(results_dir_epi_generation_epi_ageint, "model_epi_wMAE_mv_generation_epi_ageint_withintercept.rds")
)

capture.output(
  summary(epi_wMAE_mv_generation_epi_ageint_withintercept),
  file = file.path(results_dir_epi_generation_epi_ageint, "model_summary_log_scale_withintercept.txt")
)

# --- plotting (should be adapted) ---
# plotting model-specific age effects, log-transformed model

mean_age_mean <- mean(perf_epi_sub_nogest$mean_age, na.rm = TRUE)

newdat <- expand.grid(
  mean_age_c = seq(
    min(perf_epi_sub_nogest$mean_age_c, na.rm = TRUE),
    max(perf_epi_sub_nogest$mean_age_c, na.rm = TRUE),
    length.out = 100
  ),
  generation_epi = levels(perf_epi_sub_nogest$generation_epi)
)

X_new <- model.matrix(
  ~ generation_epi + generation_epi:mean_age_c - 1,
  data = newdat
)

pred <- predict(
  epi_wMAE_mv_generation_epi_ageint,
  newmods = X_new
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
    results_dir_epi_generation_epi_ageint,
    "plot_epi_meanage_wMAE_generation_epi_backtransformed.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(newdat, aes(mean_age, pred, colour = generation_epi)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = ci.lb, ymax = ci.ub, fill = generation_epi),
    alpha = .12,
    colour = NA
  ) +
  labs(
    x = "Mean sample age",
    y = "Predicted geometric mean wMAE",
    title = "Age effects on epi clock performance by generation"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- LOO ---

cohorts <- unique(perf_epi_sub_nogest$cohort)

loo_generation_epi_age <- lapply(cohorts, function(cc) {
  
  tmp <- perf_epi_sub_nogest %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = log_wMAE,
    V  = log_wMAE_var,
    mods = ~ generation_epi + generation_epi:mean_age_c - 1,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model
    ),
    data = tmp,
    method = "REML"
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
      
      generation_epi = term %>%
        str_remove("^generation_epi") %>%
        str_remove(":mean_age_c"),
      
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
  loo_generation_epi_age,
  file.path(results_dir_epi_generation_epi_ageint, "generation_epi_leave_one_cohort_out_age_interactions.csv"),
  row.names = FALSE
)

loo_generation_epi_age_summary <- loo_generation_epi_age %>%
  filter(type == "age_slope") %>%
  group_by(generation_epi) %>%
  summarise(
    min_percent_per_year = min(percent_change, na.rm = TRUE),
    max_percent_per_year = max(percent_change, na.rm = TRUE),
    range_percent_per_year = max_percent_per_year - min_percent_per_year,
    sd_percent_per_year = sd(percent_change, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  loo_generation_epi_age_summary,
  file.path(results_dir_epi_generation_epi_ageint, "generation_epi_leave_one_cohort_out_age_interactions_summary.csv"),
  row.names = FALSE
)
# --- log ---

sink(log_file)

cat("Epi wMAE multilevel meta-analysis - generation_epi age interactions, log-transformed\n")
cat("================================================================================\n\n")
cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_epi_generation_epi_ageint, "\n\n")

cat("Model summary, log scale\n")
cat("========================\n")
print(summary(epi_wMAE_mv_generation_epi_ageint))

cat("\n\nBack-transformed generation_epi age interaction results\n")
cat("==================================================\n")
print(generation_epi_age_est, n = Inf, width = Inf)

cat("Model summary with intercept, log scale\n")
cat("========================\n")
print(summary(epi_wMAE_mv_generation_epi_ageint_withintercept))

cat("\n\nBack-transformed generation_epi age interaction results with intercept\n")
cat("==================================================\n")
print(generation_epi_age_est_withintercept, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_generation_epi_age, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary for age slopes\n")
cat("===========================================\n")
print(loo_generation_epi_age_summary, n = Inf, width = Inf)

sink()

##### association ----

###### assoc generation_epi global (4-level) ------

# --- set up ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_assoc_generation_epi <- file.path(
  results_dir,
  paste0("assoc_model2_PAR_global_generation_epi_", run_id)
)

dir.create(
  results_dir_assoc_generation_epi,
  recursive = TRUE,
  showWarnings = FALSE
)

log_file <- file.path(
  results_dir_assoc_generation_epi,
  "analysis_log_generation_epi.txt"
)


# --- data prep ---

mod2_eff2 <- mod2_eff2 %>%
  mutate(
    generation_epi = factor(generation_epi)
  )


# --- meta ---

mod2_res_mv_generation_epi <- rma.mv(
  yi = RLM_Estimate_scaled,
  V  = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ generation_epi - 1,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model_combi
  ),
  data = mod2_eff2,
  method = "REML"
)

summary(mod2_res_mv_generation_epi)


# --- save model ---

saveRDS(
  mod2_res_mv_generation_epi,
  file.path(
    results_dir_assoc_generation_epi,
    "model_assoc_mv_generation_epi.rds"
  )
)

capture.output(
  summary(mod2_res_mv_generation_epi),
  file = file.path(
    results_dir_assoc_generation_epi,
    "generation_epi_model_summary.txt"
  )
)

# --- heterogeneity ---
I2 <- calc_I2_mv(
  mod2_res_mv_generation_epi,
  component_names = c("cohort", "timepoint_within_cohort", "model_combi")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_assoc_generation_epi, "heterogeneity_I2.csv"), row.names = FALSE)

# --- estimates ---

assoc_generation_epi_est <- tibble(
  term = names(coef(mod2_res_mv_generation_epi)),
  estimate = as.numeric(coef(mod2_res_mv_generation_epi)),
  se = mod2_res_mv_generation_epi$se,
  ci.lb = mod2_res_mv_generation_epi$ci.lb,
  ci.ub = mod2_res_mv_generation_epi$ci.ub,
  pval = mod2_res_mv_generation_epi$pval
)

write.csv(
  assoc_generation_epi_est,
  file.path(
    results_dir_assoc_generation_epi,
    "generation_epi_estimates.csv"
  ),
  row.names = FALSE
)


# --- predicted values ---

newdat_generation_epi <- data.frame(
  generation_epi = levels(mod2_eff2$generation_epi)
) %>%
  mutate(
    generation_epi = factor(
      generation_epi,
      levels = levels(mod2_eff2$generation_epi)
    )
  )

X_new_generation_epi <- model.matrix(
  ~ generation_epi - 1,
  data = newdat_generation_epi
)

pred_generation_epi <- predict(
  mod2_res_mv_generation_epi,
  newmods = X_new_generation_epi
)

newdat_generation_epi <- newdat_generation_epi %>%
  mutate(
    pred = pred_generation_epi$pred,
    ci.lb = pred_generation_epi$ci.lb,
    ci.ub = pred_generation_epi$ci.ub
  )

write.csv(
  newdat_generation_epi,
  file.path(
    results_dir_assoc_generation_epi,
    "generation_epi_predicted_values.csv"
  ),
  row.names = FALSE
)


# --- plot ---

png(
  filename = file.path(
    results_dir_assoc_generation_epi,
    "plot_generation_epi_meta.png"
  ),
  width = 1600,
  height = 1200,
  res = 300
)

ggplot() +
  geom_jitter(
    data = mod2_eff2,
    aes(
      x = generation_epi,
      y = RLM_Estimate_scaled
    ),
    width = 0.12,
    alpha = 0.20
  ) +
  geom_point(
    data = newdat_generation_epi,
    aes(
      x = generation_epi,
      y = pred
    ),
    size = 4
  ) +
  geom_errorbar(
    data = newdat_generation_epi,
    aes(
      x = generation_epi,
      ymin = ci.lb,
      ymax = ci.ub
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    x = "Epigenetic clock generation",
    y = "Standardised association beta",
    title = "Brain-PAR – Epi-PAR associations by epigenetic clock generation"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- LOO cohort ---

cohorts <- unique(mod2_eff2$cohort)

loo_generation_epi <- lapply(cohorts, function(cc) {
  
  tmp <- mod2_eff2 %>%
    filter(cohort != cc) %>%
    droplevels()
  
  fit <- rma.mv(
    yi = RLM_Estimate_scaled,
    V  = RLM_SE_coeftest_HC3_scaled^2,
    mods = ~ generation_epi - 1,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model_combi
    ),
    data = tmp,
    method = "REML"
  )
  
  tibble(
    cohort_removed = cc,
    term = names(coef(fit)),
    estimate = as.numeric(coef(fit)),
    se = fit$se,
    ci.lb = fit$ci.lb,
    ci.ub = fit$ci.ub,
    pval = fit$pval
  ) %>%
    mutate(
      generation_epi = str_remove(term, "^generation_epi")
    )
}) %>%
  bind_rows()

loo_generation_epi_summary <- loo_generation_epi %>%
  group_by(generation_epi) %>%
  summarise(
    min_estimate = min(estimate, na.rm = TRUE),
    max_estimate = max(estimate, na.rm = TRUE),
    range_estimate = max_estimate - min_estimate,
    sd_estimate = sd(estimate, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  loo_generation_epi,
  file.path(results_dir_assoc_generation_epi,
            "leave_one_cohort_out_generation_epi.csv"),
  row.names = FALSE
)

write.csv(
  loo_generation_epi_summary,
  file.path(results_dir_assoc_generation_epi,
            "leave_one_cohort_out_generation_epi_summary.csv"),
  row.names = FALSE
)

# --- LOO plot ---

png(
  filename = file.path(
    results_dir_assoc_generation_epi,
    "plot_generation_epi_loo.png"
  ),
  width = 2000,
  height = 1400,
  res = 300
)

ggplot(
  loo_generation_epi,
  aes(
    x = reorder(cohort_removed, estimate),
    y = estimate
  )
) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = ci.lb, ymax = ci.ub),
    width = 0.15
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  facet_wrap(~ generation_epi, scales = "free_x") +
  labs(
    x = "Cohort removed",
    y = "Estimated standardised association beta",
    title = "Leave-one-cohort-out sensitivity by epigenetic clock generation"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- full analysis log ---

sink(log_file)

cat("Association meta-analysis: epigenetic clock generation\n")
cat("======================================================\n\n")

cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_assoc_generation_epi, "\n\n")

cat("Model summary\n")
cat("=============\n")
print(summary(mod2_res_mv_generation_epi))

cat("\n\nGeneration-specific estimates\n")
cat("=============================\n")
print(assoc_generation_epi_est, width = Inf)

cat("\n\nPredicted values\n")
cat("================\n")
print(newdat_generation_epi)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_generation_epi, width = Inf, n = Inf )

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_generation_epi_summary, width = Inf, n = Inf )

sink()

###### assoc generation_epi global (Gen1 vs others) ------

# --- set up ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_assoc_generation_epi_binary <- file.path(
  results_dir,
  paste0("assoc_model2_PAR_global_generation_epi_binary_", run_id)
)

dir.create(
  results_dir_assoc_generation_epi_binary,
  recursive = TRUE,
  showWarnings = FALSE
)

log_file <- file.path(
  results_dir_assoc_generation_epi_binary,
  "analysis_log_generation_epi_binary.txt"
)


# --- data prep ---

mod2_eff2 <- mod2_eff2 %>%
  mutate(
    generation_epi_binary = if_else(
      generation_epi == 1,
      "Gen1",
      "GenNext"
    ),
    
    generation_epi_binary = factor(
      generation_epi_binary,
      levels = c("Gen1", "GenNext")
    )
  )


# --- meta ---

mod2_res_mv_generation_epi_binary <- rma.mv(
  yi = RLM_Estimate_scaled,
  V  = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ generation_epi_binary,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model_combi
  ),
  data = mod2_eff2,
  method = "REML"
)

summary(mod2_res_mv_generation_epi_binary)


# --- save model ---

saveRDS(
  mod2_res_mv_generation_epi_binary,
  file.path(
    results_dir_assoc_generation_epi_binary,
    "model_assoc_mv_generation_epi_binary.rds"
  )
)

capture.output(
  summary(mod2_res_mv_generation_epi_binary),
  file = file.path(
    results_dir_assoc_generation_epi_binary,
    "generation_epi_binary_model_summary.txt"
  )
)

# --- heterogeneity ---
I2 <- calc_I2_mv(
  mod2_res_mv_generation_epi_binary,
  component_names = c("cohort", "timepoint_within_cohort", "model_combi")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_assoc_generation_epi_binary, "heterogeneity_I2.csv"), row.names = FALSE)



# --- estimates ---

assoc_generation_epi_binary_est <- tibble(
  term = names(coef(mod2_res_mv_generation_epi_binary)),
  estimate = as.numeric(coef(mod2_res_mv_generation_epi_binary)),
  se = mod2_res_mv_generation_epi_binary$se,
  ci.lb = mod2_res_mv_generation_epi_binary$ci.lb,
  ci.ub = mod2_res_mv_generation_epi_binary$ci.ub,
  pval = mod2_res_mv_generation_epi_binary$pval
)

write.csv(
  assoc_generation_epi_binary_est,
  file.path(
    results_dir_assoc_generation_epi_binary,
    "generation_epi_binary_estimates.csv"
  ),
  row.names = FALSE
)


# --- predicted values ---

newdat_generation_epi_binary <- data.frame(
  generation_epi_binary = levels(mod2_eff2$generation_epi_binary)
)

X_new_generation_epi_binary <- model.matrix(
  ~ generation_epi_binary,
  data = newdat_generation_epi_binary
)

pred_generation_epi_binary <- predict(
  mod2_res_mv_generation_epi_binary,
  newmods = X_new_generation_epi_binary[, -1, drop = FALSE]
)

newdat_generation_epi_binary <- newdat_generation_epi_binary %>%
  mutate(
    pred = pred_generation_epi_binary$pred,
    ci.lb = pred_generation_epi_binary$ci.lb,
    ci.ub = pred_generation_epi_binary$ci.ub
  )

write.csv(
  newdat_generation_epi_binary,
  file.path(
    results_dir_assoc_generation_epi_binary,
    "generation_epi_binary_predicted_values.csv"
  ),
  row.names = FALSE
)


# --- plot ---

png(
  filename = file.path(
    results_dir_assoc_generation_epi_binary,
    "plot_generation_epi_binary_meta.png"
  ),
  width = 1600,
  height = 1200,
  res = 300
)

ggplot() +
  geom_jitter(
    data = mod2_eff2,
    aes(
      x = generation_epi_binary,
      y = RLM_Estimate_scaled
    ),
    width = 0.12,
    alpha = 0.20
  ) +
  geom_point(
    data = newdat_generation_epi_binary,
    aes(
      x = generation_epi_binary,
      y = pred
    ),
    size = 4
  ) +
  geom_errorbar(
    data = newdat_generation_epi_binary,
    aes(
      x = generation_epi_binary,
      ymin = ci.lb,
      ymax = ci.ub
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    x = "Epigenetic clock generation",
    y = "Standardised association beta",
    title = "Brain-PAR – Epi-PAR associations: Gen1 vs later-generation clocks"
  ) +
  theme_minimal(base_size = 14)

dev.off()


# --- LOO cohort ---

cohorts <- unique(mod2_eff2$cohort)

loo_generation_epi_binary <- lapply(cohorts, function(cc) {
  
  tmp <- mod2_eff2 %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = RLM_Estimate_scaled,
    V  = RLM_SE_coeftest_HC3_scaled^2,
    mods = ~ generation_epi_binary,
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
    
    generation_difference = coef(fit)[2],
    generation_difference_se = sqrt(vcov(fit)[2,2]),
    
    generation_difference_ci.lb =
      coef(fit)[2] - 1.96 * sqrt(vcov(fit)[2,2]),
    
    generation_difference_ci.ub =
      coef(fit)[2] + 1.96 * sqrt(vcov(fit)[2,2])
  )
}) %>%
  bind_rows()

loo_generation_epi_binary_summary <- tibble(
  min_generation_difference =
    min(loo_generation_epi_binary$generation_difference, na.rm = TRUE),
  
  max_generation_difference =
    max(loo_generation_epi_binary$generation_difference, na.rm = TRUE),
  
  range_generation_difference =
    max_generation_difference - min_generation_difference,
  
  sd_generation_difference =
    sd(loo_generation_epi_binary$generation_difference, na.rm = TRUE)
)

write.csv(
  loo_generation_epi_binary,
  file.path(
    results_dir_assoc_generation_epi_binary,
    "leave_one_cohort_out_generation_epi_binary.csv"
  ),
  row.names = FALSE
)

write.csv(
  loo_generation_epi_binary_summary,
  file.path(
    results_dir_assoc_generation_epi_binary,
    "leave_one_cohort_out_generation_epi_binary_summary.csv"
  ),
  row.names = FALSE
)


# --- LOO plot ---

png(
  filename = file.path(
    results_dir_assoc_generation_epi_binary,
    "plot_generation_epi_binary_loo.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_generation_epi_binary,
  aes(
    x = reorder(cohort_removed, generation_difference),
    y = generation_difference
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = generation_difference_ci.lb,
      ymax = generation_difference_ci.ub
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept = coef(mod2_res_mv_generation_epi_binary)[2],
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 0,
    colour = "grey60"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "Difference vs Gen1",
    title = "Leave-one-cohort-out sensitivity: GenNext vs Gen1"
  ) +
  theme_minimal(base_size = 14)

dev.off()


# --- full analysis log ---

sink(log_file)

cat("Association meta-analysis: Gen1 vs later-generation epigenetic clocks\n")
cat("====================================================================\n\n")

cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_assoc_generation_epi_binary, "\n\n")

cat("Model summary\n")
cat("=============\n")
print(summary(mod2_res_mv_generation_epi_binary))

cat("\n\nGeneration binary estimates\n")
cat("===========================\n")
print(assoc_generation_epi_binary_est, width = Inf, n = Inf )

cat("\n\nPredicted values\n")
cat("================\n")
print(newdat_generation_epi_binary)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_generation_epi_binary, width = Inf, n = Inf )

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_generation_epi_binary_summary, width = Inf, n = Inf )

cat("\n\nSession info\n")
cat("============\n")
print(sessionInfo())

sink()

###### assoc tissue saliva training or test  ------
# --- set up ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_assoc_tissue_saliva_testortrain <- file.path(
  results_dir,
  paste0("assoc_model2_PAR_global_tissue_saliva_testortrain_", run_id)
)

dir.create(
  results_dir_assoc_tissue_saliva_testortrain,
  recursive = TRUE,
  showWarnings = FALSE
)

log_file <- file.path(
  results_dir_assoc_tissue_saliva_testortrain,
  "analysis_log_tissue_saliva_testortrain.txt"
)


# --- data prep ---

mod2_eff2 <- mod2_eff2 %>%
  mutate(
    tissue_saliva_testortrain = case_when(
      
      tissue_test == "saliva" ~ "yes",
      
      tissue_saliva == "yes" ~ "yes",
      
      TRUE ~ "no"
    ),
    
    tissue_saliva_testortrain = factor(
      tissue_saliva_testortrain,
      levels = c("no", "yes")
    )
  )

# check if it worked 
table(mod2_eff2$tissue_saliva_testortrain)

table(
  mod2_eff2$tissue_test,
  mod2_eff2$tissue_saliva,
  mod2_eff2$tissue_saliva_testortrain
)



# --- meta ---

mod2_res_mv_tissue_saliva_testortrain <- rma.mv(
  yi = RLM_Estimate_scaled,
  V  = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ tissue_saliva_testortrain,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model_combi
  ),
  data = mod2_eff2,
  method = "REML"
)

summary(mod2_res_mv_tissue_saliva_testortrain)


# --- save model ---

saveRDS(
  mod2_res_mv_tissue_saliva_testortrain,
  file.path(
    results_dir_assoc_tissue_saliva_testortrain,
    "model_assoc_mv_tissue_saliva_testortrain.rds"
  )
)

capture.output(
  summary(mod2_res_mv_tissue_saliva_testortrain),
  file = file.path(
    results_dir_assoc_tissue_saliva_testortrain,
    "tissue_saliva_testortrain_model_summary.txt"
  )
)

# --- heterogeneity ---
I2 <- calc_I2_mv(
  mod2_res_mv_tissue_saliva_testortrain,
  component_names = c("cohort", "timepoint_within_cohort", "model_combi")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_assoc_tissue_saliva_testortrain, "heterogeneity_I2.csv"), row.names = FALSE)



# --- estimates ---

assoc_tissue_saliva_testortrain_est <- tibble(
  term = names(coef(mod2_res_mv_tissue_saliva_testortrain)),
  estimate = as.numeric(coef(mod2_res_mv_tissue_saliva_testortrain)),
  se = mod2_res_mv_tissue_saliva_testortrain$se,
  ci.lb = mod2_res_mv_tissue_saliva_testortrain$ci.lb,
  ci.ub = mod2_res_mv_tissue_saliva_testortrain$ci.ub,
  pval = mod2_res_mv_tissue_saliva_testortrain$pval
)

write.csv(
  assoc_tissue_saliva_testortrain_est,
  file.path(
    results_dir_assoc_tissue_saliva_testortrain,
    "tissue_saliva_testortrain_estimates.csv"
  ),
  row.names = FALSE
)


# --- predicted values ---

newdat_tissue_saliva_testortrain <- data.frame(
  tissue_saliva_testortrain =
    levels(mod2_eff2$tissue_saliva_testortrain)
)

X_new_tissue_saliva_testortrain <- model.matrix(
  ~ tissue_saliva_testortrain,
  data = newdat_tissue_saliva_testortrain
)

pred_tissue_saliva_testortrain <- predict(
  mod2_res_mv_tissue_saliva_testortrain,
  newmods =
    X_new_tissue_saliva_testortrain[, -1, drop = FALSE]
)

newdat_tissue_saliva_testortrain <- 
  newdat_tissue_saliva_testortrain %>%
  mutate(
    pred = pred_tissue_saliva_testortrain$pred,
    ci.lb = pred_tissue_saliva_testortrain$ci.lb,
    ci.ub = pred_tissue_saliva_testortrain$ci.ub
  )

write.csv(
  newdat_tissue_saliva_testortrain,
  file.path(
    results_dir_assoc_tissue_saliva_testortrain,
    "tissue_saliva_testortrain_predicted_values.csv"
  ),
  row.names = FALSE
)


# --- plot ---

png(
  filename = file.path(
    results_dir_assoc_tissue_saliva_testortrain,
    "plot_tissue_saliva_testortrain_meta.png"
  ),
  width = 1600,
  height = 1200,
  res = 300
)

ggplot() +
  geom_jitter(
    data = mod2_eff2,
    aes(
      x = tissue_saliva_testortrain,
      y = RLM_Estimate_scaled
    ),
    width = 0.12,
    alpha = 0.20
  ) +
  geom_point(
    data = newdat_tissue_saliva_testortrain,
    aes(
      x = tissue_saliva_testortrain,
      y = pred
    ),
    size = 4
  ) +
  geom_errorbar(
    data = newdat_tissue_saliva_testortrain,
    aes(
      x = tissue_saliva_testortrain,
      ymin = ci.lb,
      ymax = ci.ub
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    x = "Saliva involved in training or testing",
    y = "Standardised association beta",
    title = "Brain-PAR – Epi-PAR associations by saliva involvement"
  ) +
  theme_minimal(base_size = 14)

dev.off()


# --- LOO cohort ---

cohorts <- unique(mod2_eff2$cohort)

loo_tissue_saliva_testortrain <- lapply(cohorts, function(cc) {
  
  tmp <- mod2_eff2 %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = RLM_Estimate_scaled,
    V  = RLM_SE_coeftest_HC3_scaled^2,
    mods = ~ tissue_saliva_testortrain,
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
    
    saliva_difference = coef(fit)[2],
    saliva_difference_se = sqrt(vcov(fit)[2,2]),
    
    saliva_difference_ci.lb =
      coef(fit)[2] - 1.96 * sqrt(vcov(fit)[2,2]),
    
    saliva_difference_ci.ub =
      coef(fit)[2] + 1.96 * sqrt(vcov(fit)[2,2])
  )
}) %>%
  bind_rows()

loo_tissue_saliva_testortrain_summary <- tibble(
  min_saliva_difference =
    min(
      loo_tissue_saliva_testortrain$saliva_difference,
      na.rm = TRUE
    ),
  
  max_saliva_difference =
    max(
      loo_tissue_saliva_testortrain$saliva_difference,
      na.rm = TRUE
    ),
  
  range_saliva_difference =
    max_saliva_difference - min_saliva_difference,
  
  sd_saliva_difference =
    sd(
      loo_tissue_saliva_testortrain$saliva_difference,
      na.rm = TRUE
    )
)

write.csv(
  loo_tissue_saliva_testortrain,
  file.path(
    results_dir_assoc_tissue_saliva_testortrain,
    "leave_one_cohort_out_tissue_saliva_testortrain.csv"
  ),
  row.names = FALSE
)

write.csv(
  loo_tissue_saliva_testortrain_summary,
  file.path(
    results_dir_assoc_tissue_saliva_testortrain,
    "leave_one_cohort_out_tissue_saliva_testortrain_summary.csv"
  ),
  row.names = FALSE
)


# --- LOO plot ---

png(
  filename = file.path(
    results_dir_assoc_tissue_saliva_testortrain,
    "plot_tissue_saliva_testortrain_loo.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_tissue_saliva_testortrain,
  aes(
    x = reorder(cohort_removed, saliva_difference),
    y = saliva_difference
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = saliva_difference_ci.lb,
      ymax = saliva_difference_ci.ub
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept =
      coef(mod2_res_mv_tissue_saliva_testortrain)[2],
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 0,
    colour = "grey60"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "Difference vs non-saliva",
    title = "Leave-one-cohort-out sensitivity: saliva involvement"
  ) +
  theme_minimal(base_size = 14)

dev.off()


# --- full analysis log ---

sink(log_file)

cat("Association meta-analysis: saliva involvement\n")
cat("=============================================\n\n")

cat("Run ID:", run_id, "\n")
cat(
  "Output directory:",
  results_dir_assoc_tissue_saliva_testortrain,
  "\n\n"
)

cat("Model summary\n")
cat("=============\n")
print(summary(mod2_res_mv_tissue_saliva_testortrain))

cat("\n\nSaliva involvement estimates\n")
cat("================================\n")
print(
  assoc_tissue_saliva_testortrain_est,
  width = Inf,
  n = Inf
)

cat("\n\nPredicted values\n")
cat("================\n")
print(newdat_tissue_saliva_testortrain)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(
  loo_tissue_saliva_testortrain,
  width = Inf,
  n = Inf
)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(
  loo_tissue_saliva_testortrain_summary,
  width = Inf,
  n = Inf
)

sink()

###### assoc tissue saliva training AND test  ------
# --- set up ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_assoc_tissue_saliva_testANDtrain <- file.path(
  results_dir,
  paste0("assoc_model2_PAR_global_tissue_saliva_testANDtrain_", run_id)
)

dir.create(
  results_dir_assoc_tissue_saliva_testANDtrain,
  recursive = TRUE,
  showWarnings = FALSE
)

log_file <- file.path(
  results_dir_assoc_tissue_saliva_testANDtrain,
  "analysis_log_tissue_saliva_testANDtrain.txt"
)


# --- data prep ---

mod2_eff2 <- mod2_eff2 %>%
  mutate(
    tissue_saliva_testANDtrain = case_when(
      tissue_test == "saliva" & tissue_saliva == "yes" ~ "yes",
      TRUE ~ "no"
    ),
    
    tissue_saliva_testANDtrain = factor(
      tissue_saliva_testANDtrain,
      levels = c("no", "yes")
    )
  )

# check if it worked 
table(mod2_eff2$tissue_saliva_testANDtrain)

table(
  mod2_eff2$tissue_test,
  mod2_eff2$tissue_saliva,
  mod2_eff2$tissue_saliva_testANDtrain
)



# --- meta ---

mod2_res_mv_tissue_saliva_testANDtrain <- rma.mv(
  yi = RLM_Estimate_scaled,
  V  = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ tissue_saliva_testANDtrain,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model_combi
  ),
  data = mod2_eff2,
  method = "REML"
)

summary(mod2_res_mv_tissue_saliva_testANDtrain)


# --- save model ---

saveRDS(
  mod2_res_mv_tissue_saliva_testANDtrain,
  file.path(
    results_dir_assoc_tissue_saliva_testANDtrain,
    "model_assoc_mv_tissue_saliva_testANDtrain.rds"
  )
)

capture.output(
  summary(mod2_res_mv_tissue_saliva_testANDtrain),
  file = file.path(
    results_dir_assoc_tissue_saliva_testANDtrain,
    "tissue_saliva_testANDtrain_model_summary.txt"
  )
)

# --- heterogeneity ---
I2 <- calc_I2_mv(
  mod2_res_mv_tissue_saliva_testANDtrain,
  component_names = c("cohort", "timepoint_within_cohort", "model_combi")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_assoc_tissue_saliva_testANDtrain, "heterogeneity_I2.csv"), row.names = FALSE)



# --- estimates ---

assoc_tissue_saliva_testANDtrain_est <- tibble(
  term = names(coef(mod2_res_mv_tissue_saliva_testANDtrain)),
  estimate = as.numeric(coef(mod2_res_mv_tissue_saliva_testANDtrain)),
  se = mod2_res_mv_tissue_saliva_testANDtrain$se,
  ci.lb = mod2_res_mv_tissue_saliva_testANDtrain$ci.lb,
  ci.ub = mod2_res_mv_tissue_saliva_testANDtrain$ci.ub,
  pval = mod2_res_mv_tissue_saliva_testANDtrain$pval
)

write.csv(
  assoc_tissue_saliva_testANDtrain_est,
  file.path(
    results_dir_assoc_tissue_saliva_testANDtrain,
    "tissue_saliva_testANDtrain_estimates.csv"
  ),
  row.names = FALSE
)


# --- predicted values ---

newdat_tissue_saliva_testANDtrain <- data.frame(
  tissue_saliva_testANDtrain =
    levels(mod2_eff2$tissue_saliva_testANDtrain)
)

X_new_tissue_saliva_testANDtrain <- model.matrix(
  ~ tissue_saliva_testANDtrain,
  data = newdat_tissue_saliva_testANDtrain
)

pred_tissue_saliva_testANDtrain <- predict(
  mod2_res_mv_tissue_saliva_testANDtrain,
  newmods =
    X_new_tissue_saliva_testANDtrain[, -1, drop = FALSE]
)

newdat_tissue_saliva_testANDtrain <- 
  newdat_tissue_saliva_testANDtrain %>%
  mutate(
    pred = pred_tissue_saliva_testANDtrain$pred,
    ci.lb = pred_tissue_saliva_testANDtrain$ci.lb,
    ci.ub = pred_tissue_saliva_testANDtrain$ci.ub
  )

write.csv(
  newdat_tissue_saliva_testANDtrain,
  file.path(
    results_dir_assoc_tissue_saliva_testANDtrain,
    "tissue_saliva_testANDtrain_predicted_values.csv"
  ),
  row.names = FALSE
)


# --- plot ---

png(
  filename = file.path(
    results_dir_assoc_tissue_saliva_testANDtrain,
    "plot_tissue_saliva_testANDtrain_meta.png"
  ),
  width = 1600,
  height = 1200,
  res = 300
)

ggplot() +
  geom_jitter(
    data = mod2_eff2,
    aes(
      x = tissue_saliva_testANDtrain,
      y = RLM_Estimate_scaled
    ),
    width = 0.12,
    alpha = 0.20
  ) +
  geom_point(
    data = newdat_tissue_saliva_testANDtrain,
    aes(
      x = tissue_saliva_testANDtrain,
      y = pred
    ),
    size = 4
  ) +
  geom_errorbar(
    data = newdat_tissue_saliva_testANDtrain,
    aes(
      x = tissue_saliva_testANDtrain,
      ymin = ci.lb,
      ymax = ci.ub
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    x = "Saliva involved in training and testing",
    y = "Standardised association beta",
    title = "Brain-PAR – Epi-PAR associations by saliva involvement"
  ) +
  theme_minimal(base_size = 14)

dev.off()


# --- LOO cohort ---

cohorts <- unique(mod2_eff2$cohort)

loo_tissue_saliva_testANDtrain <- lapply(cohorts, function(cc) {
  
  tmp <- mod2_eff2 %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = RLM_Estimate_scaled,
    V  = RLM_SE_coeftest_HC3_scaled^2,
    mods = ~ tissue_saliva_testANDtrain,
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
    
    saliva_difference = coef(fit)[2],
    saliva_difference_se = sqrt(vcov(fit)[2,2]),
    
    saliva_difference_ci.lb =
      coef(fit)[2] - 1.96 * sqrt(vcov(fit)[2,2]),
    
    saliva_difference_ci.ub =
      coef(fit)[2] + 1.96 * sqrt(vcov(fit)[2,2])
  )
}) %>%
  bind_rows()

loo_tissue_saliva_testANDtrain_summary <- tibble(
  min_saliva_difference =
    min(
      loo_tissue_saliva_testANDtrain$saliva_difference,
      na.rm = TRUE
    ),
  
  max_saliva_difference =
    max(
      loo_tissue_saliva_testANDtrain$saliva_difference,
      na.rm = TRUE
    ),
  
  range_saliva_difference =
    max_saliva_difference - min_saliva_difference,
  
  sd_saliva_difference =
    sd(
      loo_tissue_saliva_testANDtrain$saliva_difference,
      na.rm = TRUE
    )
)

write.csv(
  loo_tissue_saliva_testANDtrain,
  file.path(
    results_dir_assoc_tissue_saliva_testANDtrain,
    "leave_one_cohort_out_tissue_saliva_testANDtrain.csv"
  ),
  row.names = FALSE
)

write.csv(
  loo_tissue_saliva_testANDtrain_summary,
  file.path(
    results_dir_assoc_tissue_saliva_testANDtrain,
    "leave_one_cohort_out_tissue_saliva_testANDtrain_summary.csv"
  ),
  row.names = FALSE
)


# --- LOO plot ---

png(
  filename = file.path(
    results_dir_assoc_tissue_saliva_testANDtrain,
    "plot_tissue_saliva_testANDtrain_loo.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_tissue_saliva_testANDtrain,
  aes(
    x = reorder(cohort_removed, saliva_difference),
    y = saliva_difference
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = saliva_difference_ci.lb,
      ymax = saliva_difference_ci.ub
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept =
      coef(mod2_res_mv_tissue_saliva_testANDtrain)[2],
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 0,
    colour = "grey60"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "Difference vs non-saliva",
    title = "Leave-one-cohort-out sensitivity: saliva involvement"
  ) +
  theme_minimal(base_size = 14)

dev.off()


# --- full analysis log ---

sink(log_file)

cat("Association meta-analysis: saliva involvement\n")
cat("=============================================\n\n")

cat("Run ID:", run_id, "\n")
cat(
  "Output directory:",
  results_dir_assoc_tissue_saliva_testANDtrain,
  "\n\n"
)

cat("Model summary\n")
cat("=============\n")
print(summary(mod2_res_mv_tissue_saliva_testANDtrain))

cat("\n\nSaliva involvement estimates\n")
cat("================================\n")
print(
  assoc_tissue_saliva_testANDtrain_est,
  width = Inf,
  n = Inf
)

cat("\n\nPredicted values\n")
cat("================\n")
print(newdat_tissue_saliva_testANDtrain)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(
  loo_tissue_saliva_testANDtrain,
  width = Inf,
  n = Inf
)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(
  loo_tissue_saliva_testANDtrain_summary,
  width = Inf,
  n = Inf
)

sink()

###### assoc tissue buccal training or test  ------
# --- set up ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_assoc_tissue_buccal_testortrain <- file.path(
  results_dir,
  paste0("assoc_model2_PAR_global_tissue_buccal_testortrain_", run_id)
)

dir.create(
  results_dir_assoc_tissue_buccal_testortrain,
  recursive = TRUE,
  showWarnings = FALSE
)

log_file <- file.path(
  results_dir_assoc_tissue_buccal_testortrain,
  "analysis_log_tissue_buccal_testortrain.txt"
)


# --- data prep ---

mod2_eff2 <- mod2_eff2 %>%
  mutate(
    tissue_buccal_testortrain = case_when(
      
      tissue_test == "buccal" ~ "yes",
      
      tissue_buccal == "yes" ~ "yes",
      
      TRUE ~ "no"
    ),
    
    tissue_buccal_testortrain = factor(
      tissue_buccal_testortrain,
      levels = c("no", "yes")
    )
  )

# check if it worked 
table(mod2_eff2$tissue_buccal_testortrain)

table(
  mod2_eff2$tissue_test,
  mod2_eff2$tissue_buccal,
  mod2_eff2$tissue_buccal_testortrain
)



# --- meta ---

mod2_res_mv_tissue_buccal_testortrain <- rma.mv(
  yi = RLM_Estimate_scaled,
  V  = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ tissue_buccal_testortrain,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model_combi
  ),
  data = mod2_eff2,
  method = "REML"
)

summary(mod2_res_mv_tissue_buccal_testortrain)


# --- save model ---

saveRDS(
  mod2_res_mv_tissue_buccal_testortrain,
  file.path(
    results_dir_assoc_tissue_buccal_testortrain,
    "model_assoc_mv_tissue_buccal_testortrain.rds"
  )
)

capture.output(
  summary(mod2_res_mv_tissue_buccal_testortrain),
  file = file.path(
    results_dir_assoc_tissue_buccal_testortrain,
    "tissue_buccal_testortrain_model_summary.txt"
  )
)

# --- heterogeneity ---
I2 <- calc_I2_mv(
  mod2_res_mv_tissue_buccal_testortrain,
  component_names = c("cohort", "timepoint_within_cohort", "model_combi")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_assoc_tissue_buccal_testortrain, "heterogeneity_I2.csv"), row.names = FALSE)

# --- estimates ---

assoc_tissue_buccal_testortrain_est <- tibble(
  term = names(coef(mod2_res_mv_tissue_buccal_testortrain)),
  estimate = as.numeric(coef(mod2_res_mv_tissue_buccal_testortrain)),
  se = mod2_res_mv_tissue_buccal_testortrain$se,
  ci.lb = mod2_res_mv_tissue_buccal_testortrain$ci.lb,
  ci.ub = mod2_res_mv_tissue_buccal_testortrain$ci.ub,
  pval = mod2_res_mv_tissue_buccal_testortrain$pval
)

write.csv(
  assoc_tissue_buccal_testortrain_est,
  file.path(
    results_dir_assoc_tissue_buccal_testortrain,
    "tissue_buccal_testortrain_estimates.csv"
  ),
  row.names = FALSE
)


# --- predicted values ---

newdat_tissue_buccal_testortrain <- data.frame(
  tissue_buccal_testortrain =
    levels(mod2_eff2$tissue_buccal_testortrain)
)

X_new_tissue_buccal_testortrain <- model.matrix(
  ~ tissue_buccal_testortrain,
  data = newdat_tissue_buccal_testortrain
)

pred_tissue_buccal_testortrain <- predict(
  mod2_res_mv_tissue_buccal_testortrain,
  newmods =
    X_new_tissue_buccal_testortrain[, -1, drop = FALSE]
)

newdat_tissue_buccal_testortrain <- 
  newdat_tissue_buccal_testortrain %>%
  mutate(
    pred = pred_tissue_buccal_testortrain$pred,
    ci.lb = pred_tissue_buccal_testortrain$ci.lb,
    ci.ub = pred_tissue_buccal_testortrain$ci.ub
  )

write.csv(
  newdat_tissue_buccal_testortrain,
  file.path(
    results_dir_assoc_tissue_buccal_testortrain,
    "tissue_buccal_testortrain_predicted_values.csv"
  ),
  row.names = FALSE
)


# --- plot ---

png(
  filename = file.path(
    results_dir_assoc_tissue_buccal_testortrain,
    "plot_tissue_buccal_testortrain_meta.png"
  ),
  width = 1600,
  height = 1200,
  res = 300
)

ggplot() +
  geom_jitter(
    data = mod2_eff2,
    aes(
      x = tissue_buccal_testortrain,
      y = RLM_Estimate_scaled
    ),
    width = 0.12,
    alpha = 0.20
  ) +
  geom_point(
    data = newdat_tissue_buccal_testortrain,
    aes(
      x = tissue_buccal_testortrain,
      y = pred
    ),
    size = 4
  ) +
  geom_errorbar(
    data = newdat_tissue_buccal_testortrain,
    aes(
      x = tissue_buccal_testortrain,
      ymin = ci.lb,
      ymax = ci.ub
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    x = "buccal involved in training or testing",
    y = "Standardised association beta",
    title = "Brain-PAR – Epi-PAR associations by buccal involvement"
  ) +
  theme_minimal(base_size = 14)

dev.off()


# --- LOO cohort ---

cohorts <- unique(mod2_eff2$cohort)

loo_tissue_buccal_testortrain <- lapply(cohorts, function(cc) {
  
  tmp <- mod2_eff2 %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = RLM_Estimate_scaled,
    V  = RLM_SE_coeftest_HC3_scaled^2,
    mods = ~ tissue_buccal_testortrain,
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
    
    buccal_difference = coef(fit)[2],
    buccal_difference_se = sqrt(vcov(fit)[2,2]),
    
    buccal_difference_ci.lb =
      coef(fit)[2] - 1.96 * sqrt(vcov(fit)[2,2]),
    
    buccal_difference_ci.ub =
      coef(fit)[2] + 1.96 * sqrt(vcov(fit)[2,2])
  )
}) %>%
  bind_rows()

loo_tissue_buccal_testortrain_summary <- tibble(
  min_buccal_difference =
    min(
      loo_tissue_buccal_testortrain$buccal_difference,
      na.rm = TRUE
    ),
  
  max_buccal_difference =
    max(
      loo_tissue_buccal_testortrain$buccal_difference,
      na.rm = TRUE
    ),
  
  range_buccal_difference =
    max_buccal_difference - min_buccal_difference,
  
  sd_buccal_difference =
    sd(
      loo_tissue_buccal_testortrain$buccal_difference,
      na.rm = TRUE
    )
)

write.csv(
  loo_tissue_buccal_testortrain,
  file.path(
    results_dir_assoc_tissue_buccal_testortrain,
    "leave_one_cohort_out_tissue_buccal_testortrain.csv"
  ),
  row.names = FALSE
)

write.csv(
  loo_tissue_buccal_testortrain_summary,
  file.path(
    results_dir_assoc_tissue_buccal_testortrain,
    "leave_one_cohort_out_tissue_buccal_testortrain_summary.csv"
  ),
  row.names = FALSE
)


# --- LOO plot ---

png(
  filename = file.path(
    results_dir_assoc_tissue_buccal_testortrain,
    "plot_tissue_buccal_testortrain_loo.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_tissue_buccal_testortrain,
  aes(
    x = reorder(cohort_removed, buccal_difference),
    y = buccal_difference
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = buccal_difference_ci.lb,
      ymax = buccal_difference_ci.ub
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept =
      coef(mod2_res_mv_tissue_buccal_testortrain)[2],
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 0,
    colour = "grey60"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "Difference vs non-buccal",
    title = "Leave-one-cohort-out sensitivity: buccal involvement"
  ) +
  theme_minimal(base_size = 14)

dev.off()


# --- full analysis log ---

sink(log_file)

cat("Association meta-analysis: buccal involvement\n")
cat("=============================================\n\n")

cat("Run ID:", run_id, "\n")
cat(
  "Output directory:",
  results_dir_assoc_tissue_buccal_testortrain,
  "\n\n"
)

cat("Model summary\n")
cat("=============\n")
print(summary(mod2_res_mv_tissue_buccal_testortrain))

cat("\n\nbuccal involvement estimates\n")
cat("================================\n")
print(
  assoc_tissue_buccal_testortrain_est,
  width = Inf,
  n = Inf
)

cat("\n\nPredicted values\n")
cat("================\n")
print(newdat_tissue_buccal_testortrain)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(
  loo_tissue_buccal_testortrain,
  width = Inf,
  n = Inf
)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(
  loo_tissue_buccal_testortrain_summary,
  width = Inf,
  n = Inf
)

cat("\n\nSession info\n")
cat("============\n")
print(sessionInfo())

sink()

###### assoc tissue buccal training AND test  ------
# --- set up ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_assoc_tissue_buccal_testANDtrain <- file.path(
  results_dir,
  paste0("assoc_model2_PAR_global_tissue_buccal_testANDtrain_", run_id)
)

dir.create(
  results_dir_assoc_tissue_buccal_testANDtrain,
  recursive = TRUE,
  showWarnings = FALSE
)

log_file <- file.path(
  results_dir_assoc_tissue_buccal_testANDtrain,
  "analysis_log_tissue_buccal_testANDtrain.txt"
)


# --- data prep ---

mod2_eff2 <- mod2_eff2 %>%
  mutate(
    tissue_buccal_testANDtrain = case_when(
      tissue_test == "buccal" & tissue_buccal == "yes" ~ "yes",
      TRUE ~ "no"
    ),
    
    tissue_buccal_testANDtrain = factor(
      tissue_buccal_testANDtrain,
      levels = c("no", "yes")
    )
  )

# check if it worked 
table(mod2_eff2$tissue_buccal_testANDtrain)

table(
  mod2_eff2$tissue_test,
  mod2_eff2$tissue_buccal,
  mod2_eff2$tissue_buccal_testANDtrain
)



# --- meta ---

mod2_res_mv_tissue_buccal_testANDtrain <- rma.mv(
  yi = RLM_Estimate_scaled,
  V  = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ tissue_buccal_testANDtrain,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model_combi
  ),
  data = mod2_eff2,
  method = "REML"
)

summary(mod2_res_mv_tissue_buccal_testANDtrain)


# --- save model ---

saveRDS(
  mod2_res_mv_tissue_buccal_testANDtrain,
  file.path(
    results_dir_assoc_tissue_buccal_testANDtrain,
    "model_assoc_mv_tissue_buccal_testANDtrain.rds"
  )
)

capture.output(
  summary(mod2_res_mv_tissue_buccal_testANDtrain),
  file = file.path(
    results_dir_assoc_tissue_buccal_testANDtrain,
    "tissue_buccal_testANDtrain_model_summary.txt"
  )
)

# --- heterogeneity ---
I2 <- calc_I2_mv(
  mod2_res_mv_tissue_buccal_testANDtrain,
  component_names = c("cohort", "timepoint_within_cohort", "model_combi")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_assoc_tissue_buccal_testANDtrain, "heterogeneity_I2.csv"), row.names = FALSE)

# --- estimates ---

assoc_tissue_buccal_testANDtrain_est <- tibble(
  term = names(coef(mod2_res_mv_tissue_buccal_testANDtrain)),
  estimate = as.numeric(coef(mod2_res_mv_tissue_buccal_testANDtrain)),
  se = mod2_res_mv_tissue_buccal_testANDtrain$se,
  ci.lb = mod2_res_mv_tissue_buccal_testANDtrain$ci.lb,
  ci.ub = mod2_res_mv_tissue_buccal_testANDtrain$ci.ub,
  pval = mod2_res_mv_tissue_buccal_testANDtrain$pval
)

write.csv(
  assoc_tissue_buccal_testANDtrain_est,
  file.path(
    results_dir_assoc_tissue_buccal_testANDtrain,
    "tissue_buccal_testANDtrain_estimates.csv"
  ),
  row.names = FALSE
)


# --- predicted values ---

newdat_tissue_buccal_testANDtrain <- data.frame(
  tissue_buccal_testANDtrain =
    levels(mod2_eff2$tissue_buccal_testANDtrain)
)

X_new_tissue_buccal_testANDtrain <- model.matrix(
  ~ tissue_buccal_testANDtrain,
  data = newdat_tissue_buccal_testANDtrain
)

pred_tissue_buccal_testANDtrain <- predict(
  mod2_res_mv_tissue_buccal_testANDtrain,
  newmods =
    X_new_tissue_buccal_testANDtrain[, -1, drop = FALSE]
)

newdat_tissue_buccal_testANDtrain <- 
  newdat_tissue_buccal_testANDtrain %>%
  mutate(
    pred = pred_tissue_buccal_testANDtrain$pred,
    ci.lb = pred_tissue_buccal_testANDtrain$ci.lb,
    ci.ub = pred_tissue_buccal_testANDtrain$ci.ub
  )

write.csv(
  newdat_tissue_buccal_testANDtrain,
  file.path(
    results_dir_assoc_tissue_buccal_testANDtrain,
    "tissue_buccal_testANDtrain_predicted_values.csv"
  ),
  row.names = FALSE
)


# --- plot ---

png(
  filename = file.path(
    results_dir_assoc_tissue_buccal_testANDtrain,
    "plot_tissue_buccal_testANDtrain_meta.png"
  ),
  width = 1600,
  height = 1200,
  res = 300
)

ggplot() +
  geom_jitter(
    data = mod2_eff2,
    aes(
      x = tissue_buccal_testANDtrain,
      y = RLM_Estimate_scaled
    ),
    width = 0.12,
    alpha = 0.20
  ) +
  geom_point(
    data = newdat_tissue_buccal_testANDtrain,
    aes(
      x = tissue_buccal_testANDtrain,
      y = pred
    ),
    size = 4
  ) +
  geom_errorbar(
    data = newdat_tissue_buccal_testANDtrain,
    aes(
      x = tissue_buccal_testANDtrain,
      ymin = ci.lb,
      ymax = ci.ub
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    x = "buccal involved in training AND testing",
    y = "Standardised association beta",
    title = "Brain-PAR – Epi-PAR associations by buccal involvement"
  ) +
  theme_minimal(base_size = 14)

dev.off()


# --- LOO cohort ---
# removed all LOO analyses as only one cohort has buccal data 

# --- full analysis log ---

sink(log_file)

cat("Association meta-analysis: buccal involvement\n")
cat("=============================================\n\n")

cat("Run ID:", run_id, "\n")
cat(
  "Output directory:",
  results_dir_assoc_tissue_buccal_testANDtrain,
  "\n\n"
)

cat("Model summary\n")
cat("=============\n")
print(summary(mod2_res_mv_tissue_buccal_testANDtrain))

cat("\n\nbuccal involvement estimates\n")
cat("================================\n")
print(
  assoc_tissue_buccal_testANDtrain_est,
  width = Inf,
  n = Inf
)

cat("\n\nPredicted values\n")
cat("================\n")
print(newdat_tissue_buccal_testANDtrain)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
cat("left out as only one cohort contributes effect sizes to this")

cat("\n\nSession info\n")
cat("============\n")
print(sessionInfo())

sink()

###### assoc tissue cortex training or test  ------
# --- set up ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_assoc_tissue_cortical <- file.path(
  results_dir,
  paste0("assoc_model2_PAR_global_tissue_cortical_", run_id)
)

dir.create(
  results_dir_assoc_tissue_cortical,
  recursive = TRUE,
  showWarnings = FALSE
)

log_file <- file.path(
  results_dir_assoc_tissue_cortical,
  "analysis_log_tissue_cortical.txt"
)


# --- data prep ---

mod2_eff2 <- mod2_eff2 %>%
  mutate(
    tissue_cortical = factor(tissue_cortical)
  )

# check if it worked 
table(mod2_eff2$tissue_cortical)

table(
  mod2_eff2$tissue_test,
  mod2_eff2$tissue_training,
  mod2_eff2$tissue_cortical
)

# --- meta ---

mod2_res_mv_tissue_cortical <- rma.mv(
  yi = RLM_Estimate_scaled,
  V  = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ tissue_cortical,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model_combi
  ),
  data = mod2_eff2,
  method = "REML"
)

summary(mod2_res_mv_tissue_cortical)


# --- save model ---

saveRDS(
  mod2_res_mv_tissue_cortical,
  file.path(
    results_dir_assoc_tissue_cortical,
    "model_assoc_mv_tissue_cortical.rds"
  )
)

capture.output(
  summary(mod2_res_mv_tissue_cortical),
  file = file.path(
    results_dir_assoc_tissue_cortical,
    "tissue_cortical_model_summary.txt"
  )
)

# --- heterogeneity ---
I2 <- calc_I2_mv(
  mod2_res_mv_tissue_cortical,
  component_names = c("cohort", "timepoint_within_cohort", "model_combi")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_assoc_tissue_cortical, "heterogeneity_I2.csv"), row.names = FALSE)

# --- estimates ---

assoc_tissue_cortical_est <- tibble(
  term = names(coef(mod2_res_mv_tissue_cortical)),
  estimate = as.numeric(coef(mod2_res_mv_tissue_cortical)),
  se = mod2_res_mv_tissue_cortical$se,
  ci.lb = mod2_res_mv_tissue_cortical$ci.lb,
  ci.ub = mod2_res_mv_tissue_cortical$ci.ub,
  pval = mod2_res_mv_tissue_cortical$pval
)

write.csv(
  assoc_tissue_cortical_est,
  file.path(
    results_dir_assoc_tissue_cortical,
    "tissue_cortical_estimates.csv"
  ),
  row.names = FALSE
)


# --- predicted values ---

newdat_tissue_cortical <- data.frame(
  tissue_cortical =
    levels(mod2_eff2$tissue_cortical)
)

X_new_tissue_cortical <- model.matrix(
  ~ tissue_cortical,
  data = newdat_tissue_cortical
)

pred_tissue_cortical <- predict(
  mod2_res_mv_tissue_cortical,
  newmods =
    X_new_tissue_cortical[, -1, drop = FALSE]
)

newdat_tissue_cortical <- 
  newdat_tissue_cortical %>%
  mutate(
    pred = pred_tissue_cortical$pred,
    ci.lb = pred_tissue_cortical$ci.lb,
    ci.ub = pred_tissue_cortical$ci.ub
  )

write.csv(
  newdat_tissue_cortical,
  file.path(
    results_dir_assoc_tissue_cortical,
    "tissue_cortical_predicted_values.csv"
  ),
  row.names = FALSE
)


# --- plot ---

png(
  filename = file.path(
    results_dir_assoc_tissue_cortical,
    "plot_tissue_cortical_meta.png"
  ),
  width = 1600,
  height = 1200,
  res = 300
)

ggplot() +
  geom_jitter(
    data = mod2_eff2,
    aes(
      x = tissue_cortical,
      y = RLM_Estimate_scaled
    ),
    width = 0.12,
    alpha = 0.20
  ) +
  geom_point(
    data = newdat_tissue_cortical,
    aes(
      x = tissue_cortical,
      y = pred
    ),
    size = 4
  ) +
  geom_errorbar(
    data = newdat_tissue_cortical,
    aes(
      x = tissue_cortical,
      ymin = ci.lb,
      ymax = ci.ub
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    x = "cortex involved in training or testing",
    y = "Standardised association beta",
    title = "Brain-PAR – Epi-PAR associations by cortex involvement"
  ) +
  theme_minimal(base_size = 14)

dev.off()


# --- LOO cohort ---

cohorts <- unique(mod2_eff2$cohort)

loo_tissue_cortical <- lapply(cohorts, function(cc) {
  
  tmp <- mod2_eff2 %>%
    filter(cohort != cc)
  
  fit <- rma.mv(
    yi = RLM_Estimate_scaled,
    V  = RLM_SE_coeftest_HC3_scaled^2,
    mods = ~ tissue_cortical,
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
    
    cortex_difference = coef(fit)[2],
    cortex_difference_se = sqrt(vcov(fit)[2,2]),
    
    cortex_difference_ci.lb =
      coef(fit)[2] - 1.96 * sqrt(vcov(fit)[2,2]),
    
    cortex_difference_ci.ub =
      coef(fit)[2] + 1.96 * sqrt(vcov(fit)[2,2])
  )
}) %>%
  bind_rows()

loo_tissue_cortical_summary <- tibble(
  min_cortex_difference =
    min(
      loo_tissue_cortical$cortex_difference,
      na.rm = TRUE
    ),
  
  max_cortex_difference =
    max(
      loo_tissue_cortical$cortex_difference,
      na.rm = TRUE
    ),
  
  range_cortex_difference =
    max_cortex_difference - min_cortex_difference,
  
  sd_cortex_difference =
    sd(
      loo_tissue_cortical$cortex_difference,
      na.rm = TRUE
    )
)

write.csv(
  loo_tissue_cortical,
  file.path(
    results_dir_assoc_tissue_cortical,
    "leave_one_cohort_out_tissue_cortical.csv"
  ),
  row.names = FALSE
)

write.csv(
  loo_tissue_cortical_summary,
  file.path(
    results_dir_assoc_tissue_cortical,
    "leave_one_cohort_out_tissue_cortical_summary.csv"
  ),
  row.names = FALSE
)


# --- LOO plot ---

png(
  filename = file.path(
    results_dir_assoc_tissue_cortical,
    "plot_tissue_cortical_loo.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_tissue_cortical,
  aes(
    x = reorder(cohort_removed, cortex_difference),
    y = cortex_difference
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = cortex_difference_ci.lb,
      ymax = cortex_difference_ci.ub
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept =
      coef(mod2_res_mv_tissue_cortical)[2],
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 0,
    colour = "grey60"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "Difference vs non-cortex",
    title = "Leave-one-cohort-out sensitivity: cortex involvement"
  ) +
  theme_minimal(base_size = 14)

dev.off()


# --- full analysis log ---

sink(log_file)

cat("Association meta-analysis: cortex involvement\n")
cat("=============================================\n\n")

cat("Run ID:", run_id, "\n")
cat(
  "Output directory:",
  results_dir_assoc_tissue_cortical,
  "\n\n"
)

cat("Model summary\n")
cat("=============\n")
print(summary(mod2_res_mv_tissue_cortical))

cat("\n\ncortex involvement estimates\n")
cat("================================\n")
print(
  assoc_tissue_cortical_est,
  width = Inf,
  n = Inf
)

cat("\n\nPredicted values\n")
cat("================\n")
print(newdat_tissue_cortical)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(
  loo_tissue_cortical,
  width = Inf,
  n = Inf
)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(
  loo_tissue_cortical_summary,
  width = Inf,
  n = Inf
)

cat("\n\nSession info\n")
cat("============\n")
print(sessionInfo())

sink()

##### association - purer comparisons ----

mod2_eff2 <- mod2_eff2 %>%
  mutate(
    # epi overlap
    overlap_lower_epi = pmax(min_age_epi, training_age_min_epi),
    overlap_upper_epi = pmin(max_age_epi, training_age_max_epi),
    
    overlap_epi = if_else(
      overlap_upper_epi >= overlap_lower_epi,
      "yes",
      "no"
    ),
    
    overlap_width_epi = if_else(
      overlap_epi == "yes",
      overlap_upper_epi - overlap_lower_epi,
      0
    ),
    
    study_width_epi = max_age_epi - min_age_epi,
    
    overlap_perc_epi = 100 * overlap_width_epi / study_width_epi,
    
    
    # brain overlap
    overlap_lower_brain = pmax(min_age_brain, training_age_min_brain),
    overlap_upper_brain = pmin(max_age_brain, training_age_max_brain),
    
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
    
    study_width_brain = max_age_brain - min_age_brain,
    
    overlap_perc_brain = 100 * overlap_width_brain / study_width_brain
  ) %>%
  select(
    -overlap_lower_epi, -overlap_upper_epi,
    -overlap_width_epi, -study_width_epi,
    -overlap_lower_brain, -overlap_upper_brain,
    -overlap_width_brain, -study_width_brain
  )

# create tissue overlap variable 
mod2_eff2 <- mod2_eff2 %>%
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
  mod2_eff2$tissue_training,
  mod2_eff2$tissue_test,
  mod2_eff2$tissue_overlap
)

mod2_eff2_full_overlap <- mod2_eff2 %>%
  filter(
    overlap_perc_brain == 100,
    overlap_perc_epi == 100,
    tissue_overlap == "yes"
  )


# descriptives: full-overlap subset

desc_overlap <- tibble(
  
  n_associations = nrow(mod2_eff2_full_overlap),
  
  n_brain_models =
    n_distinct(mod2_eff2_full_overlap$brain_model),
  
  n_epi_models =
    n_distinct(mod2_eff2_full_overlap$epi_model),
  
  n_model_combinations =
    n_distinct(mod2_eff2_full_overlap$model_combi),
  
  n_cohorts =
    n_distinct(mod2_eff2_full_overlap$cohort),
  
  n_timepoints =
    n_distinct(
      paste(mod2_eff2_full_overlap$cohort,
            mod2_eff2_full_overlap$timepoint)
    ),
  
  mean_n =
    mean(mod2_eff2_full_overlap$RLM_N, na.rm = TRUE),
  
  median_n =
    median(mod2_eff2_full_overlap$RLM_N, na.rm = TRUE),
  
  min_n =
    min(mod2_eff2_full_overlap$RLM_N, na.rm = TRUE),
  
  max_n =
    max(mod2_eff2_full_overlap$RLM_N, na.rm = TRUE)
)

desc_overlap

tibble(
  metric = c(
    "Associations retained",
    "Brain models retained",
    "Epi models retained",
    "Model combinations retained"
  ),
  
  value = c(
    nrow(mod2_eff2_full_overlap) / nrow(mod2_eff2),
    
    n_distinct(mod2_eff2_full_overlap$brain_model) /
      n_distinct(mod2_eff2$brain_model),
    
    n_distinct(mod2_eff2_full_overlap$epi_model) /
      n_distinct(mod2_eff2$epi_model),
    
    n_distinct(mod2_eff2_full_overlap$model_combi) /
      n_distinct(mod2_eff2$model_combi)
  ) * 100
)

# which models disappear? 

setdiff(
  unique(mod2_eff2$brain_model),
  unique(mod2_eff2_full_overlap$brain_model)
)

setdiff(
  unique(mod2_eff2$epi_model),
  unique(mod2_eff2_full_overlap$epi_model)
)

###### assoc global -----

# --- set up folder ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_assoc <- file.path(
  results_dir,
  paste0("assoc_model2_PAR_global_pure", run_id)
)

dir.create(results_dir_assoc, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_assoc, "analysis_log.txt")

# --- descriptives for 'pure' dataset ---

desc_overlap <- tibble(
  n_associations = nrow(mod2_eff2_full_overlap),
  n_brain_models = n_distinct(mod2_eff2_full_overlap$brain_model),
  n_epi_models = n_distinct(mod2_eff2_full_overlap$epi_model),
  n_model_combinations = n_distinct(mod2_eff2_full_overlap$model_combi),
  n_cohorts = n_distinct(mod2_eff2_full_overlap$cohort),
  n_timepoints = n_distinct(
    paste(mod2_eff2_full_overlap$cohort, mod2_eff2_full_overlap$timepoint)
  ),
  mean_n = mean(mod2_eff2_full_overlap$RLM_N, na.rm = TRUE),
  median_n = median(mod2_eff2_full_overlap$RLM_N, na.rm = TRUE),
  min_n = min(mod2_eff2_full_overlap$RLM_N, na.rm = TRUE),
  max_n = max(mod2_eff2_full_overlap$RLM_N, na.rm = TRUE)
)

retention_summary <- tibble(
  metric = c(
    "Associations retained",
    "Brain models retained",
    "Epi models retained",
    "Model combinations retained"
  ),
  value_percent = c(
    nrow(mod2_eff2_full_overlap) / nrow(mod2_eff2),
    n_distinct(mod2_eff2_full_overlap$brain_model) / n_distinct(mod2_eff2$brain_model),
    n_distinct(mod2_eff2_full_overlap$epi_model) / n_distinct(mod2_eff2$epi_model),
    n_distinct(mod2_eff2_full_overlap$model_combi) / n_distinct(mod2_eff2$model_combi)
  ) * 100
)

missing_brain_models <- tibble(
  brain_model = setdiff(
    unique(mod2_eff2$brain_model),
    unique(mod2_eff2_full_overlap$brain_model)
  )
)

missing_epi_models <- tibble(
  epi_model = setdiff(
    unique(mod2_eff2$epi_model),
    unique(mod2_eff2_full_overlap$epi_model)
  )
)

write.csv(
  desc_overlap,
  file.path(results_dir_assoc, "descriptives_full_overlap.csv"),
  row.names = FALSE
)

write.csv(
  retention_summary,
  file.path(results_dir_assoc, "retention_summary_full_overlap.csv"),
  row.names = FALSE
)

write.csv(
  missing_brain_models,
  file.path(results_dir_assoc, "excluded_brain_models_full_overlap.csv"),
  row.names = FALSE
)

write.csv(
  missing_epi_models,
  file.path(results_dir_assoc, "excluded_epi_models_full_overlap.csv"),
  row.names = FALSE
)


# --- meta ---
mod2_res_mv_global <- rma.mv(
  yi = RLM_Estimate_scaled,
  V  = RLM_SE_coeftest_HC3_scaled^2,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model_combi
  ),
  data = mod2_eff2_full_overlap,
  method = "REML",
  verbose = TRUE
)

summary(mod2_res_mv_global)

saveRDS(
  mod2_res_mv_global,
  file.path(results_dir_assoc, "model_assoc_mv_global.rds")
)

capture.output(
  summary(mod2_res_mv_global),
  file = file.path(results_dir_assoc, "model_summary.txt")
)

# --- heterogeneity ---
I2 <- calc_I2_mv(
  mod2_res_mv_global,
  component_names = c("cohort", "timepoint_within_cohort", "model_combi")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_assoc, "heterogeneity_I2.csv"), row.names = FALSE)


# --- LOO cohort ---

cohorts <- unique(mod2_eff2_full_overlap$cohort)

loo_cohort <- lapply(cohorts, function(cc) {
  
  tmp <- subset(mod2_eff2_full_overlap, cohort != cc)
  
  fit <- rma.mv(
    yi = RLM_Estimate_scaled,
    V  = RLM_SE_coeftest_HC3_scaled^2,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | model_combi
    ),
    data = tmp,
    method = "REML"
  )
  
  data.frame(
    cohort_removed = cc,
    beta = coef(fit)[1],
    se   = sqrt(vcov(fit)[1, 1])
  )
})

loo_cohort <- do.call(rbind, loo_cohort) %>%
  as_tibble()

loo_cohort
range(loo_cohort$beta)
sd(loo_cohort$beta)

loo_summary <- tibble(
  min_beta   = min(loo_cohort$beta, na.rm = TRUE),
  max_beta   = max(loo_cohort$beta, na.rm = TRUE),
  range_beta = max_beta - min_beta,
  sd_beta    = sd(loo_cohort$beta, na.rm = TRUE)
)

write.csv(
  loo_cohort,
  file.path(results_dir_assoc, "leave_one_cohort_out.csv"),
  row.names = FALSE
)

write.csv(
  loo_summary,
  file.path(results_dir_assoc, "leave_one_cohort_out_summary.csv"),
  row.names = FALSE
)


full_beta <- coef(mod2_res_mv_global)[1]

png(
  filename = file.path(results_dir_assoc, "plot_assoc_global_loo.png"),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_cohort,
  aes(
    x = reorder(cohort_removed, beta),
    y = beta
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = beta - 1.96 * se,
      ymax = beta + 1.96 * se
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept = full_beta,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "Pooled standardised beta",
    title = "Leave-one-cohort-out sensitivity of pooled association"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- full analysis log ---

sink(log_file)

cat("Association multilevel meta-analysis - tissue & age overlap subsample \n")
cat("====================================\n\n")
cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_assoc, "\n\n")

cat("\n\nDescriptives for full-overlap dataset\n")
cat("=====================================\n")
print(desc_overlap, width = Inf)

cat("\n\nRetention summary\n")
cat("=================\n")
print(retention_summary)

cat("\n\nBrain models excluded by full-overlap filter\n")
cat("============================================\n")
print(missing_brain_models)

cat("\n\nEpi models excluded by full-overlap filter\n")
cat("==========================================\n")
print(missing_epi_models)

cat("Model summary\n")
cat("=============\n")
print(summary(mod2_res_mv_global))

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_cohort)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_summary)

sink()

###### assoc modelwise ------

# --- set up folder ---

run_id_modelwise <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_assoc_modelwise <- file.path(
  results_dir,
  paste0("assoc_model2_PAR_modelwise_pure_", run_id_modelwise)
)

dir.create(results_dir_assoc_modelwise, recursive = TRUE, showWarnings = FALSE)

log_file_modelwise <- file.path(results_dir_assoc_modelwise, "analysis_log_modelwise.txt")


# --- meta ---

# remove all empty factor levels
mod2_eff2_full_overlap <- mod2_eff2_full_overlap %>%
  mutate(
    model_combi = droplevels(factor(model_combi)),
    cohort = droplevels(factor(cohort)),
    timepoint = droplevels(factor(timepoint))
  )

mod2_res_mv_modelwise <- rma.mv(
  yi   = RLM_Estimate_scaled,
  V    = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ model_combi - 1,
  random = ~ 1 | cohort/timepoint,
  data   = mod2_eff2_full_overlap,
  method = "REML"
)

summary(mod2_res_mv_modelwise)

saveRDS(
  mod2_res_mv_modelwise,
  file.path(results_dir_assoc_modelwise, "model_assoc_mv_modelwise.rds")
)

capture.output(
  summary(mod2_res_mv_modelwise),
  file = file.path(results_dir_assoc_modelwise, "modelwise_model_summary.txt")
)

# --- heterogeneity ---

I2 <- calc_I2_mv(
  mod2_res_mv_modelwise,
  component_names = c("cohort", "timepoint_within_cohort")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_assoc_modelwise, "heterogeneity_I2.csv"), row.names = FALSE)


# --- extract modelwise estimates ---

assoc_modelwise_est <- tibble(
  term  = rownames(mod2_res_mv_modelwise$b),
  beta  = as.numeric(mod2_res_mv_modelwise$b),
  se    = mod2_res_mv_modelwise$se,
  ci.lb = mod2_res_mv_modelwise$ci.lb,
  ci.ub = mod2_res_mv_modelwise$ci.ub,
  pval  = mod2_res_mv_modelwise$pval
) %>%
  mutate(model_combi = str_remove(term, "^model_combi")) %>%
  separate(model_combi, into = c("brain_model", "epi_model"), sep = "\\.", remove = FALSE) %>%
  select(-term)

assoc_modelwise_est <- assoc_modelwise_est  %>%
  arrange(beta)

X <- model.matrix(mod2_res_mv_modelwise)

model_combi_cols <- colnames(X)[
  str_detect(colnames(X), "^model_combi") &
    !str_detect(colnames(X), ":")
]

effect_sizes_per_model_combi <- tibble(
  model_combi = str_remove(model_combi_cols, "^model_combi"),
  k = colSums(X[, model_combi_cols, drop = FALSE] != 0))

assoc_modelwise_est <- assoc_modelwise_est %>%
  left_join(effect_sizes_per_model_combi, by = "model_combi")

write.csv(
  assoc_modelwise_est,
  file.path(results_dir_assoc_modelwise, "modelwise_estimates.csv"),
  row.names = FALSE
)


# --- LOO cohort ---
# slightly different code here as it doesn't run 

cohorts <- unique(mod2_eff2_full_overlap$cohort)

loo_cohort_modelwise <- lapply(cohorts, function(cc) {
  
  message("Running LOO modelwise: removing ", cc)
  
  tmp <- subset(mod2_eff2_full_overlap, cohort != cc)
  
  fit <- tryCatch(
    rma.mv(
      yi   = RLM_Estimate_scaled,
      V    = RLM_SE_coeftest_HC3_scaled^2,
      mods = ~ model_combi - 1,
      random = ~ 1 | cohort/timepoint,
      data = tmp,
      method = "REML"
    ),
    error = function(e) {
      message("LOO failed for cohort: ", cc)
      message("Error: ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(fit)) {
    return(tibble(
      cohort_removed = cc,
      term = NA_character_,
      brain_model = NA_character_,
      epi_model = NA_character_,
      beta = NA_real_,
      se = NA_real_,
      ci.lb = NA_real_,
      ci.ub = NA_real_,
      pval = NA_real_,
      converged = FALSE
    ))
  }
  
  tibble(
    cohort_removed = cc,
    term  = rownames(fit$b),
    beta  = as.numeric(fit$b),
    se    = fit$se,
    ci.lb = fit$ci.lb,
    ci.ub = fit$ci.ub,
    pval  = fit$pval,
    converged = TRUE
  ) %>%
    mutate(model_combi = str_remove(term, "^model_combi")) %>%
    separate(model_combi, into = c("brain_model", "epi_model"), sep = "\\.", remove = FALSE) %>%
    select(cohort_removed, brain_model, epi_model, beta, se, ci.lb, ci.ub, pval, converged)
}) %>%
  bind_rows()

# cohorts <- unique(mod2_eff2_full_overlap$cohort)
# 
# loo_cohort_modelwise <- lapply(cohorts, function(cc) {
#   
#   tmp <- subset(mod2_eff2_full_overlap, cohort != cc)
#   
#   fit <- rma.mv(
#     yi   = RLM_Estimate_scaled,
#     V    = RLM_SE_coeftest_HC3_scaled^2,
#     mods = ~ model_combi - 1,
#     random = ~ 1 | cohort/timepoint,
#     data = tmp,
#     method = "REML"
#   )
#   
#   tibble(
#     cohort_removed = cc,
#     term  = rownames(fit$b),
#     beta  = as.numeric(fit$b),
#     se    = fit$se,
#     ci.lb = fit$ci.lb,
#     ci.ub = fit$ci.ub,
#     pval  = fit$pval
#   ) %>%
#     mutate(model_combi = str_remove(term, "^model_combi")) %>%
#     separate(model_combi, into = c("brain_model", "epi_model"), sep = "\\.") %>%
#     select(cohort_removed, brain_model, epi_model, beta, se, ci.lb, ci.ub, pval)
# }) %>%
#   bind_rows()

loo_cohort_modelwise %>%
  filter(!converged) %>%
  distinct(cohort_removed)

write.csv(
  loo_cohort_modelwise,
  file.path(results_dir_assoc_modelwise, "modelwise_leave_one_cohort_out.csv"),
  row.names = FALSE
)

modelwise_loo_summary <- loo_cohort_modelwise %>%
  filter(converged) %>%
  group_by(brain_model, epi_model) %>%
  summarise(
    min_beta   = min(beta, na.rm = TRUE),
    max_beta   = max(beta, na.rm = TRUE),
    range_beta = max_beta - min_beta,
    sd_beta    = sd(beta, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  modelwise_loo_summary,
  file.path(results_dir_assoc_modelwise, "modelwise_leave_one_cohort_out_summary.csv"),
  row.names = FALSE
)


# --- plot heatmap ---

png(
  filename = file.path(results_dir_assoc_modelwise, "heatmap_assoc_modelwise.png"),
  width = 1800,
  height = 1800,
  res = 300
)

assoc_modelwise_est_plot <- assoc_modelwise_est %>%
  mutate(
    epi_model = factor(epi_model, levels = rev(sort(unique(epi_model))))
  )

ggplot(assoc_modelwise_est_plot, aes(x = brain_model, y = epi_model, fill = beta)) +
  geom_tile() +
  scale_fill_gradient2(
    mid = "white",
    midpoint = 0,
    limits = c(-0.1, 0.1),
    oob = scales::squish
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Brain age model",
    y = "Epigenetic age model",
    fill = "Beta",
    title = "Modelwise brain age - epigenetic age associations"
  )

dev.off()

# --- full analysis log ---

sink(log_file_modelwise)

cat("Association modelwise multilevel meta-analysis\n")
cat("==============================================\n\n")
cat("Run ID:", run_id_modelwise, "\n")
cat("Output directory:", results_dir_assoc_modelwise, "\n\n")

cat("Modelwise model summary\n")
cat("=======================\n")
print(summary(mod2_res_mv_modelwise))

cat("\n\nModelwise estimates\n")
cat("===================\n")
print(assoc_modelwise_est, n = Inf, width = Inf)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nModelwise leave-one-cohort-out results\n")
cat("=====================================\n")
print(loo_cohort_modelwise, n = Inf, width = Inf)

cat("\n\nModelwise leave-one-cohort-out summary\n")
cat("=====================================\n")
print(modelwise_loo_summary, n = Inf, width = Inf)

sink()

##### assoc age global - pure sample/brain age -----

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_assoc_global_ageint_brain <- file.path(
  results_dir,
  paste0("assoc_model2_PAR_global_ageint_brain_pure_", run_id)
)

dir.create(results_dir_assoc_global_ageint_brain, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  results_dir_assoc_global_ageint_brain,
  "analysis_log_global_ageint_brain.txt"
)

# --- data prep ---

mod2_eff2_full_overlap_age_brain <- mod2_eff2_full_overlap %>%
  mutate(
    cohort = droplevels(factor(cohort)),
    timepoint = droplevels(factor(timepoint)),
    model_combi = droplevels(factor(model_combi)),
    mean_age_brain_c = mean_age_brain - mean(mean_age_brain, na.rm = TRUE)
  )

# --- meta ---

mod2_res_mv_global_ageint_brain <- rma.mv(
  yi = RLM_Estimate_scaled,
  V  = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ mean_age_brain_c,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model_combi
  ),
  data = mod2_eff2_full_overlap_age_brain,
  method = "REML",
  verbose = TRUE
)

summary(mod2_res_mv_global_ageint_brain)

saveRDS(
  mod2_res_mv_global_ageint_brain,
  file.path(
    results_dir_assoc_global_ageint_brain,
    "model_assoc_mv_global_ageint_brain.rds"
  )
)

capture.output(
  summary(mod2_res_mv_global_ageint_brain),
  file = file.path(
    results_dir_assoc_global_ageint_brain,
    "global_ageint_brain_model_summary.txt"
  )
)

# --- heterogeneity ---
I2 <- calc_I2_mv(
  mod2_res_mv_global_ageint_brain,
  component_names = c("cohort", "timepoint_within_cohort", "model_combi")
)

I2_tbl <- tibble(
  component = c("total", names(I2$I2_components)),
  I2_percent = c(I2$I2_total, as.numeric(I2$I2_components))
)

write.csv(I2_tbl, file.path(results_dir_assoc_global_ageint_brain, "heterogeneity_I2.csv"), row.names = FALSE)


# --- estimates ---

assoc_global_ageint_brain_est <- tibble(
  term = names(coef(mod2_res_mv_global_ageint_brain)),
  estimate = as.numeric(coef(mod2_res_mv_global_ageint_brain)),
  se = mod2_res_mv_global_ageint_brain$se,
  ci.lb = mod2_res_mv_global_ageint_brain$ci.lb,
  ci.ub = mod2_res_mv_global_ageint_brain$ci.ub,
  pval = mod2_res_mv_global_ageint_brain$pval
)

write.csv(
  assoc_global_ageint_brain_est,
  file.path(results_dir_assoc_global_ageint_brain, "global_ageint_brain_estimates.csv"),
  row.names = FALSE
)


# --- prediction plot ---

newdat <- data.frame(
  mean_age_brain_c = seq(
    min(mod2_eff2_full_overlap_age_brain$mean_age_brain_c, na.rm = TRUE),
    max(mod2_eff2_full_overlap_age_brain$mean_age_brain_c, na.rm = TRUE),
    length.out = 100
  )
)

pred <- predict(
  mod2_res_mv_global_ageint_brain,
  newmods = newdat$mean_age_brain_c
)

newdat$pred <- pred$pred
newdat$ci.lb <- pred$ci.lb
newdat$ci.ub <- pred$ci.ub

newdat$mean_age_brain <- newdat$mean_age_brain_c +
  mean(mod2_eff2_full_overlap_age_brain$mean_age_brain, na.rm = TRUE)

png(
  filename = file.path(
    results_dir_assoc_global_ageint_brain,
    "association_global_meanage_brain_meta_regression.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot() +
  geom_point(
    data = mod2_eff2_full_overlap_age_brain,
    aes(
      x = mean_age_brain,
      y = RLM_Estimate_scaled,
      colour = cohort
    ),
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
  theme_minimal(base_size = 14) +
  labs(
    title = "Brain-PAR - Epi-PAR association by brain mean age",
    x = "Mean brain age",
    y = "Standardised beta (RLM)",
    colour = "Cohort"
  )

dev.off()

# --- LOO cohort ---

cohorts <- unique(mod2_eff2_full_overlap_age_brain$cohort)

loo_assoc_ageint_brain <- lapply(cohorts, function(cc) {
  
  tmp <- mod2_eff2_full_overlap_age_brain %>%
    filter(cohort != cc) %>%
    mutate(
      cohort = droplevels(factor(cohort)),
      timepoint = droplevels(factor(timepoint)),
      model_combi = droplevels(factor(model_combi))
    )
  
  fit <- rma.mv(
    yi = RLM_Estimate_scaled,
    V  = RLM_SE_coeftest_HC3_scaled^2,
    mods = ~ mean_age_brain_c,
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
    intercept_se = sqrt(vcov(fit)[1, 1]),
    age_slope = coef(fit)[2],
    age_slope_se = sqrt(vcov(fit)[2, 2])
  )
}) %>%
  bind_rows()

loo_assoc_ageint_brain_summary <- loo_assoc_ageint_brain %>%
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
  loo_assoc_ageint_brain,
  file.path(results_dir_assoc_global_ageint_brain, "leave_one_cohort_out_age_interaction.csv"),
  row.names = FALSE
)

write.csv(
  loo_assoc_ageint_brain_summary,
  file.path(results_dir_assoc_global_ageint_brain, "leave_one_cohort_out_age_interaction_summary.csv"),
  row.names = FALSE
)

png(
  filename = file.path(results_dir_assoc_global_ageint_brain, "plot_assoc_meanage_brain_loo.png"),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_assoc_ageint_brain,
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
    yintercept = coef(mod2_res_mv_global_ageint_brain)[2],
    linetype = "dashed"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "Estimated brain mean-age slope",
    title = "Leave-one-cohort-out sensitivity of brain mean-age effect"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- log ---

sink(log_file)

cat("Association global age interaction meta-analysis: brain mean age, filtered subset\n")
cat("================================================================================\n\n")

cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_assoc_global_ageint_brain, "\n\n")

cat("Model summary\n")
cat("=============\n")
print(summary(mod2_res_mv_global_ageint_brain))

cat("\n\nGlobal age interaction estimates\n")
cat("================================\n")
print(assoc_global_ageint_brain_est)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_tbl, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_assoc_ageint_brain, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_assoc_ageint_brain_summary, n = Inf, width = Inf)

cat("\n\nSession info\n")
cat("============\n")
print(sessionInfo())

sink()


##### assoc age global - pure sample/epi age -----


run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_assoc_global_ageint_epi <- file.path(
  results_dir,
  paste0("assoc_model2_PAR_global_ageint_epi_pure_", run_id)
)

dir.create(results_dir_assoc_global_ageint_epi, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  results_dir_assoc_global_ageint_epi,
  "analysis_log_global_ageint_epi.txt"
)

# --- data prep ---

mod2_eff2_full_overlap_age_epi <- mod2_eff2_full_overlap %>%
  mutate(
    cohort = droplevels(factor(cohort)),
    timepoint = droplevels(factor(timepoint)),
    model_combi = droplevels(factor(model_combi)),
    mean_age_epi_c = mean_age_epi - mean(mean_age_epi, na.rm = TRUE)
  )

# --- meta ---

mod2_res_mv_global_ageint_epi <- rma.mv(
  yi = RLM_Estimate_scaled,
  V  = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ mean_age_epi_c,
  random = list(
    ~ 1 | cohort/timepoint,
    ~ 1 | model_combi
  ),
  data = mod2_eff2_full_overlap_age_epi,
  method = "REML",
  verbose = TRUE
)

summary(mod2_res_mv_global_ageint_epi)

saveRDS(
  mod2_res_mv_global_ageint_epi,
  file.path(
    results_dir_assoc_global_ageint_epi,
    "model_assoc_mv_global_ageint_epi.rds"
  )
)

capture.output(
  summary(mod2_res_mv_global_ageint_epi),
  file = file.path(
    results_dir_assoc_global_ageint_epi,
    "global_ageint_epi_model_summary.txt"
  )
)

# --- estimates ---

assoc_global_ageint_epi_est <- tibble(
  term = names(coef(mod2_res_mv_global_ageint_epi)),
  estimate = as.numeric(coef(mod2_res_mv_global_ageint_epi)),
  se = mod2_res_mv_global_ageint_epi$se,
  ci.lb = mod2_res_mv_global_ageint_epi$ci.lb,
  ci.ub = mod2_res_mv_global_ageint_epi$ci.ub,
  pval = mod2_res_mv_global_ageint_epi$pval
)

write.csv(
  assoc_global_ageint_epi_est,
  file.path(results_dir_assoc_global_ageint_epi, "global_ageint_epi_estimates.csv"),
  row.names = FALSE
)

# --- I2 ---

I2_epi_age <- calc_I2_mv(
  mod2_res_mv_global_ageint_epi,
  component_names = c(
    "cohort",
    "timepoint_within_cohort",
    "model_combi"
  )
)

I2_epi_age_tbl <- tibble(
  component = c("total", names(I2_epi_age$I2_components)),
  I2_percent = c(
    I2_epi_age$I2_total,
    as.numeric(I2_epi_age$I2_components)
  )
)

write.csv(
  I2_epi_age_tbl,
  file.path(results_dir_assoc_global_ageint_epi, "I2_heterogeneity.csv"),
  row.names = FALSE
)

# --- prediction plot ---

newdat <- data.frame(
  mean_age_epi_c = seq(
    min(mod2_eff2_full_overlap_age_epi$mean_age_epi_c, na.rm = TRUE),
    max(mod2_eff2_full_overlap_age_epi$mean_age_epi_c, na.rm = TRUE),
    length.out = 100
  )
)

pred <- predict(
  mod2_res_mv_global_ageint_epi,
  newmods = newdat$mean_age_epi_c
)

newdat$pred <- pred$pred
newdat$ci.lb <- pred$ci.lb
newdat$ci.ub <- pred$ci.ub

newdat$mean_age_epi <- newdat$mean_age_epi_c +
  mean(mod2_eff2_full_overlap_age_epi$mean_age_epi, na.rm = TRUE)

png(
  filename = file.path(
    results_dir_assoc_global_ageint_epi,
    "association_global_meanage_epi_meta_regression.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot() +
  geom_point(
    data = mod2_eff2_full_overlap_age_epi,
    aes(
      x = mean_age_epi,
      y = RLM_Estimate_scaled,
      colour = cohort
    ),
    alpha = .5
  ) +
  geom_line(
    data = newdat,
    aes(x = mean_age_epi, y = pred),
    colour = "black",
    linewidth = 1
  ) +
  geom_ribbon(
    data = newdat,
    aes(x = mean_age_epi, ymin = ci.lb, ymax = ci.ub),
    alpha = .15
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Brain-PAR - Epi-PAR association by epi mean age",
    x = "Mean epigenetic age",
    y = "Standardised beta (RLM)",
    colour = "Cohort"
  )

dev.off()

# --- LOO cohort ---

cohorts <- unique(mod2_eff2_full_overlap_age_epi$cohort)

loo_assoc_ageint_epi <- lapply(cohorts, function(cc) {
  
  tmp <- mod2_eff2_full_overlap_age_epi %>%
    filter(cohort != cc) %>%
    mutate(
      cohort = droplevels(factor(cohort)),
      timepoint = droplevels(factor(timepoint)),
      model_combi = droplevels(factor(model_combi))
    )
  
  fit <- rma.mv(
    yi = RLM_Estimate_scaled,
    V  = RLM_SE_coeftest_HC3_scaled^2,
    mods = ~ mean_age_epi_c,
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
    intercept_se = sqrt(vcov(fit)[1, 1]),
    age_slope = coef(fit)[2],
    age_slope_se = sqrt(vcov(fit)[2, 2])
  )
}) %>%
  bind_rows()

loo_assoc_ageint_epi_summary <- loo_assoc_ageint_epi %>%
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
  loo_assoc_ageint_epi,
  file.path(results_dir_assoc_global_ageint_epi, "leave_one_cohort_out_age_interaction.csv"),
  row.names = FALSE
)

write.csv(
  loo_assoc_ageint_epi_summary,
  file.path(results_dir_assoc_global_ageint_epi, "leave_one_cohort_out_age_interaction_summary.csv"),
  row.names = FALSE
)

png(
  filename = file.path(results_dir_assoc_global_ageint_epi, "plot_assoc_meanage_epi_loo.png"),
  width = 1800,
  height = 1200,
  res = 300
)

ggplot(
  loo_assoc_ageint_epi,
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
    yintercept = coef(mod2_res_mv_global_ageint_epi)[2],
    linetype = "dashed"
  ) +
  coord_flip() +
  labs(
    x = "Cohort removed",
    y = "Estimated epi mean-age slope",
    title = "Leave-one-cohort-out sensitivity of epi mean-age effect"
  ) +
  theme_minimal(base_size = 14)

dev.off()

# --- log ---

sink(log_file)

cat("Association global age interaction meta-analysis: epi mean age, filtered subset\n")
cat("==============================================================================\n\n")

cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_assoc_global_ageint_epi, "\n\n")

cat("Model summary\n")
cat("=============\n")
print(summary(mod2_res_mv_global_ageint_epi))

cat("\n\nGlobal age interaction estimates\n")
cat("================================\n")
print(assoc_global_ageint_epi_est)

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_epi_age_tbl)

cat("\n\nLeave-one-cohort-out results\n")
cat("============================\n")
print(loo_assoc_ageint_epi, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_assoc_ageint_epi_summary, n = Inf, width = Inf)

cat("\n\nSession info\n")
cat("============\n")
print(sessionInfo())

sink()

##### assoc age modelwise - pure sample/brain age -----

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_assoc_age_brain_modelwise <- file.path(
  results_dir,
  paste0("assoc_model2_PAR_modelwise_ageint_brain_pure_", run_id)
)

dir.create(results_dir_assoc_age_brain_modelwise, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  results_dir_assoc_age_brain_modelwise,
  "analysis_log_modelwise_ageint_brain.txt"
)

# --- data prep ---

mod2_eff2_full_overlap_age_brain <- mod2_eff2_full_overlap %>%
  mutate(
    cohort = droplevels(factor(cohort)),
    timepoint = droplevels(factor(timepoint)),
    model_combi = droplevels(factor(model_combi)),
    mean_age_brain_c = mean_age_brain - mean(mean_age_brain, na.rm = TRUE)
  )

# --- meta ---

mod2_res_mv_pairs_ageint_brain <- rma.mv(
  yi = RLM_Estimate_scaled,
  V  = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ model_combi + model_combi:mean_age_brain_c - 1,
  random = ~ 1 | cohort/timepoint,
  data = mod2_eff2_full_overlap_age_brain,
  method = "REML"
)

summary(mod2_res_mv_pairs_ageint_brain)

saveRDS(
  mod2_res_mv_pairs_ageint_brain,
  file.path(
    results_dir_assoc_age_brain_modelwise,
    "model_assoc_mv_modelwise_ageint_brain.rds"
  )
)

capture.output(
  summary(mod2_res_mv_pairs_ageint_brain),
  file = file.path(
    results_dir_assoc_age_brain_modelwise,
    "model_summary_modelwise_ageint_brain.txt"
  )
)

# --- I2 ---

I2_pairs_age_brain <- calc_I2_mv(
  mod2_res_mv_pairs_ageint_brain,
  component_names = c(
    "cohort",
    "timepoint_within_cohort"
  )
)

I2_pairs_age_brain_tbl <- tibble(
  component = c("total", names(I2_pairs_age_brain$I2_components)),
  I2_percent = c(
    I2_pairs_age_brain$I2_total,
    as.numeric(I2_pairs_age_brain$I2_components)
  )
)

write.csv(
  I2_pairs_age_brain_tbl,
  file.path(results_dir_assoc_age_brain_modelwise, "I2_heterogeneity.csv"),
  row.names = FALSE
)

# --- extract estimates ---

pair_age_est_brain <- tibble(
  term = rownames(mod2_res_mv_pairs_ageint_brain$b),
  b = as.numeric(mod2_res_mv_pairs_ageint_brain$b),
  se = mod2_res_mv_pairs_ageint_brain$se,
  ci.lb = mod2_res_mv_pairs_ageint_brain$ci.lb,
  ci.ub = mod2_res_mv_pairs_ageint_brain$ci.ub,
  pval = mod2_res_mv_pairs_ageint_brain$pval
)

pair_age_est_brain_parsed <- pair_age_est_brain %>%
  mutate(
    term = str_remove(term, "^model_combi"),
    param_type = if_else(
      str_detect(term, ":mean_age_brain_c"),
      "age_slope",
      "intercept"
    ),
    pair_term = str_remove(term, ":mean_age_brain_c")
  ) %>%
  separate(
    pair_term,
    into = c("brain_model", "epi_model"),
    sep = "\\.",
    remove = FALSE
  )

X <- model.matrix(mod2_res_mv_pairs_ageint_brain)

model_combi_cols <- colnames(X)[
  str_detect(colnames(X), "^model_combi") &
    !str_detect(colnames(X), ":")
]

effect_sizes_per_pair_brain <- tibble(
  pair_term = str_remove(model_combi_cols, "^model_combi"),
  k = colSums(X[, model_combi_cols, drop = FALSE] != 0)
) %>%
  separate(
    pair_term,
    into = c("brain_model", "epi_model"),
    sep = "\\.",
    remove = FALSE
  )

pair_age_est_brain_parsed <- pair_age_est_brain_parsed %>%
  left_join(
    effect_sizes_per_pair_brain,
    by = c("brain_model", "epi_model")
  )

pair_intercepts_brain_age <- pair_age_est_brain_parsed %>%
  filter(param_type == "intercept")

pair_slopes_brain_age <- pair_age_est_brain_parsed %>%
  filter(param_type == "age_slope")

write.csv(
  pair_age_est_brain_parsed,
  file.path(
    results_dir_assoc_age_brain_modelwise,
    "associations_pairwise_modelwise_ageint_brain_all.csv"
  ),
  row.names = FALSE
)

write.csv(
  pair_intercepts_brain_age,
  file.path(
    results_dir_assoc_age_brain_modelwise,
    "associations_pairwise_intercepts_mean_age_brain.csv"
  ),
  row.names = FALSE
)

write.csv(
  pair_slopes_brain_age,
  file.path(
    results_dir_assoc_age_brain_modelwise,
    "associations_pairwise_age_slopes_mean_age_brain.csv"
  ),
  row.names = FALSE
)


# --- heatmap ---

epi_order <- pair_slopes_brain_age %>%
  distinct(epi_model) %>%
  arrange(epi_model) %>%
  pull(epi_model)

brain_order <- pair_slopes_brain_age %>%
  distinct(brain_model) %>%
  arrange(brain_model) %>%
  pull(brain_model)

pair_slopes_brain_age_plot <- pair_slopes_brain_age %>%
  mutate(
    epi_model = factor(epi_model, levels = rev(epi_order)),
    brain_model = factor(brain_model, levels = brain_order)
  )

png(
  file.path(
    results_dir_assoc_age_brain_modelwise,
    "heatmap_pairwise_age_slopes_mean_age_brain.png"
  ),
  width = 1800,
  height = 1800,
  res = 300
)

ggplot(
  pair_slopes_brain_age_plot,
  aes(
    x = brain_model,
    y = epi_model,
    fill = b
  )
) +
  geom_tile() +
  scale_fill_gradient2(
    mid = "white",
    midpoint = 0,
    limits = c(-0.016, 0.016),
    oob = scales::squish
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Age moderation of brain-epi associations: brain mean age",
    x = "Brain age model",
    y = "Epigenetic age model",
    fill = "Slope"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

dev.off()

# --- leave-one-cohort-out ---

cohorts <- unique(mod2_eff2_full_overlap_age_brain$cohort)

loo_pairs_age_brain <- lapply(cohorts, function(cc) {
  
  message("Running LOO modelwise age interaction: removing ", cc)
  
  tmp <- mod2_eff2_full_overlap_age_brain %>%
    filter(cohort != cc) %>%
    mutate(
      cohort = droplevels(factor(cohort)),
      timepoint = droplevels(factor(timepoint)),
      model_combi = droplevels(factor(model_combi))
    )
  
  fit <- tryCatch(
    rma.mv(
      yi = RLM_Estimate_scaled,
      V  = RLM_SE_coeftest_HC3_scaled^2,
      mods = ~ model_combi + model_combi:mean_age_brain_c - 1,
      random = ~ 1 | cohort/timepoint,
      data = tmp,
      method = "REML"
    ),
    error = function(e) {
      message("LOO failed for cohort: ", cc)
      message("Error: ", e$message)
      return(e)
    }
  )
  
  if (inherits(fit, "error")) {
    return(
      tibble(
        cohort_removed = cc,
        term = NA_character_,
        b = NA_real_,
        se = NA_real_,
        ci.lb = NA_real_,
        ci.ub = NA_real_,
        pval = NA_real_,
        converged = FALSE,
        error_message = fit$message
      )
    )
  }
  
  tibble(
    cohort_removed = cc,
    term = rownames(fit$b),
    b = as.numeric(fit$b),
    se = fit$se,
    ci.lb = fit$ci.lb,
    ci.ub = fit$ci.ub,
    pval = fit$pval,
    converged = TRUE,
    error_message = NA_character_
  )
}) %>%
  bind_rows()

full_pairs_age_brain <- tibble(
  term = rownames(mod2_res_mv_pairs_ageint_brain$b),
  b_full = as.numeric(mod2_res_mv_pairs_ageint_brain$b)
)

loo_compare_brain_age <- loo_pairs_age_brain %>%
  left_join(full_pairs_age_brain, by = "term") %>%
  mutate(diff = b - b_full)

loo_summary_brain_age <- loo_compare_brain_age %>%
  filter(converged) %>%
  group_by(cohort_removed) %>%
  summarise(
    mean_abs_diff = mean(abs(diff), na.rm = TRUE),
    max_abs_diff = max(abs(diff), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abs_diff))

write.csv(
  loo_pairs_age_brain,
  file.path(
    results_dir_assoc_age_brain_modelwise,
    "modelwise_leave_one_cohort_out_ageint_brain.csv"
  ),
  row.names = FALSE
)

write.csv(
  loo_summary_brain_age,
  file.path(
    results_dir_assoc_age_brain_modelwise,
    "modelwise_leave_one_cohort_out_ageint_brain_summary.csv"
  ),
  row.names = FALSE
)

# --- log ---

sink(log_file)

cat("Association modelwise age interaction meta-analysis: brain mean age, filtered subset\n")
cat("===================================================================================\n\n")

cat("Run ID:", run_id, "\n")
cat("Output directory:", results_dir_assoc_age_brain_modelwise, "\n\n")

cat("Model summary\n")
cat("=============\n")
print(summary(mod2_res_mv_pairs_ageint_brain))

cat("\n\nHeterogeneity statistics\n")
cat("========================\n")
print(I2_pairs_age_brain_tbl)

cat("\n\nPairwise intercepts\n")
cat("===================\n")
print(pair_intercepts_brain_age, n = Inf, width = Inf)

cat("\n\nPairwise age slopes\n")
cat("===================\n")
print(pair_slopes_brain_age, n = Inf, width = Inf)

cat("\n\nLeave-one-cohort-out failures\n")
cat("=============================\n")
print(
  loo_pairs_age_brain %>%
    filter(!converged) %>%
    distinct(cohort_removed, error_message),
  n = Inf,
  width = Inf
)

cat("\n\nLeave-one-cohort-out summary\n")
cat("============================\n")
print(loo_summary_brain_age, n = Inf, width = Inf)

cat("\n\nSession info\n")
cat("============\n")
print(sessionInfo())

sink()

