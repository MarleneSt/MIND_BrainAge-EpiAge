#### extracting association files for model 1 and model 3 sensitivity analyses 

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
results_dir <- "results" #[REPRODUCIBILITY] relative path

##### model 1 -----
# extracting model 1 files 
files_mod1 <- list.files(
  path = data_dir,
  pattern = "_model1_PAR\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

parse_fname <- function(fname) {
  base <- basename(fname)
  # e.g. Cohort1_T1_model1_PAR.csv
  m <- str_match(
    base,
    "^([^_]+)_(.+)_model(\\d+)_([^\\.]+)\\.csv$"
  )
  tibble(
    cohort        = m[, 2],
    timepoint     = m[, 3],
    model_num     = m[, 4],
    var_type_short = m[, 5],
    file_path     = fname          
  )
}
#this might need to be adjusted if people have added in tissue or array 

file_info <- purrr::map_dfr(files_mod1, parse_fname)

mod1_all <- file_info %>%
  mutate(data = purrr::map(file_path, ~ read.csv(.x, check.names = FALSE))) %>%
  select(-file_path) %>%
  unnest(data)

# extract model effect 
mod1_eff <- mod1_all %>%
  filter(
    Model == "Model_1",
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
mod1_eff2 <- mod1_eff %>%
  mutate(
    model_combi = interaction(brain_model, epi_model, drop = TRUE)
  )

# changing sign for DNAmTL 
mod1_eff2 <- mod1_eff2 %>%
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
mods_brain <- read.csv(file.path(data_dir,"moderators_cohorts_brain_R.csv"))
mods_epi <- read.csv(file.path(data_dir,"moderators_cohorts_epi_R.csv"))


# change names for alignment 
mod1_eff2 <- mod1_eff2 %>%
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
mod1_keys <- mod1_eff2 %>%
  distinct(cohort, timepoint)

mods_keys <- mods_association %>%
  distinct(cohort, timepoint)

anti_join(
  mod1_keys,
  mods_keys,
  by = c("cohort", "timepoint")
)

anti_join(
  mods_keys,
  mod1_keys,
  by = c("cohort", "timepoint")
)

mod1_eff2 <- mod1_eff2 %>% left_join(mods_association, by = c("cohort", "timepoint"))

# we also need to merge in the epi cohort moderators (to get tissue) - mods_epi

mod1_keys <- mod1_eff2 %>%
  distinct(cohort, timepoint)

mods_keys <- mods_epi %>%
  distinct(cohort, timepoint)

anti_join(
  mod1_keys,
  mods_keys,
  by = c("cohort", "timepoint")
)

anti_join(
  mods_keys,
  mod1_keys,
  by = c("cohort", "timepoint")
) # timepoints with no imaging equivalent left over, this is fine 

mod1_eff2 <- mod1_eff2 %>% left_join(mods_epi, by = c("cohort", "timepoint"))

# merge model level moderators 

#adding in variables with the model name for merging 
mod1_eff2$model_brain <- sub("_PAR$", "", mod1_eff2$Outcome)
mod1_eff2$model_epi <- sub("_PAR$", "", mod1_eff2$Predictor)

# check if joining will work (model level moderator - brain)
mod1_keys <- mod1_eff2 %>%
  distinct(model_brain)

mods_keys <- mods_brain_models %>%
  distinct(model_brain)

anti_join(
  mod1_keys,
  mods_keys,
  by = c("model_brain")
)

anti_join(
  mods_keys,
  mod1_keys,
  by = c("model_brain")
)

# check if joining will work (model level moderator - epi)
mod1_keys <- mod1_eff2 %>%
  distinct(model_epi)

mods_keys <- mods_epi_models %>%
  distinct(model_epi)

anti_join(
  mod1_keys,
  mods_keys,
  by = c("model_epi")
)

anti_join(
  mods_keys,
  mod1_keys,
  by = c("model_epi")
) # gestational clocls missing but that's fine

mod1_eff2 <- mod1_eff2 %>% left_join(mods_brain_models, by = c("model_brain"))
mod1_eff2 <- mod1_eff2 %>% left_join(mods_epi_models, by = c("model_epi"))

# rename some cohort names for nicer plotting & to align with performance 
mod1_eff2 <- mod1_eff2 %>%
  mutate(
    cohort = recode(
      cohort,
      "SAND" = "FFCWS",
      "Oregon_ADHD-1000" = "Oregon ADHD-1000",
      "Kids2Health_C" = "K2H Childhood"
    ))

# --- check for missing or implausible values ---
mod1_eff2 %>%
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
write.csv(mod1_eff2, file.path(results_dir, "mod1_eff2.csv"), row.names = FALSE)
saveRDS(mod1_eff2, file.path(results_dir, "mod1_eff2.rds"))

##### assoc descriptives general -----

assoc_desc <- mod1_eff2 %>%
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

write.csv(assoc_desc, file.path(results_dir, paste0("assoc_desc_model1_", "final", ".csv")), row.names = FALSE) # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))



##### model 3 -----

# extracing model 3 files 
files_mod3 <- list.files(
  path = data_dir,
  pattern = "_model3_PAR\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

parse_fname <- function(fname) {
  base <- basename(fname)
  # e.g. Cohort1_T1_model1_PAR.csv
  m <- str_match(
    base,
    "^([^_]+)_(.+)_model(\\d+)_([^\\.]+)\\.csv$"
  )
  tibble(
    cohort        = m[, 2],
    timepoint     = m[, 3],
    model_num     = m[, 4],
    var_type_short = m[, 5],
    file_path     = fname          
  )
}
#this might need to be adjusted if people have added in tissue or array 

file_info <- purrr::map_dfr(files_mod3, parse_fname)

mod3_all <- file_info %>%
  mutate(data = purrr::map(file_path, ~ read.csv(.x, check.names = FALSE))) %>%
  select(-file_path) %>%
  unnest(data)

# extract model 3 effect 
mod3_eff <- mod3_all %>%
  filter(
    Model == "Model_3",
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
mod3_eff2 <- mod3_eff %>%
  mutate(
    model_combi = interaction(brain_model, epi_model, drop = TRUE)
  )

# changing sign for DNAmTL 
mod3_eff2 <- mod3_eff2 %>%
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
mods_brain <- read.csv(file.path(data_dir,"moderators_cohorts_brain_R.csv"))
mods_epi <- read.csv(file.path(data_dir,"moderators_cohorts_epi_R.csv"))


# merge cohort level moderators 

# change names for alignment 
mod3_eff2 <- mod3_eff2 %>%
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
mod3_keys <- mod3_eff2 %>%
  distinct(cohort, timepoint)

mods_keys <- mods_association %>%
  distinct(cohort, timepoint)

anti_join(
  mod3_keys,
  mods_keys,
  by = c("cohort", "timepoint")
)

anti_join(
  mods_keys,
  mod3_keys,
  by = c("cohort", "timepoint")
)

mod3_eff2 <- mod3_eff2 %>% left_join(mods_association, by = c("cohort", "timepoint"))

# we also need to merge in the epi cohort moderators (to get tissue) - mods_epi

mod3_keys <- mod3_eff2 %>%
  distinct(cohort, timepoint)

mods_keys <- mods_epi %>%
  distinct(cohort, timepoint)

anti_join(
  mod3_keys,
  mods_keys,
  by = c("cohort", "timepoint")
)

anti_join(
  mods_keys,
  mod3_keys,
  by = c("cohort", "timepoint")
) # timepoints with no imaging equivalent left over, this is fine 

mod3_eff2 <- mod3_eff2 %>% left_join(mods_epi, by = c("cohort", "timepoint"))

# merge model level moderators 

#adding in variables with the model name for merging 
mod3_eff2$model_brain <- sub("_PAR$", "", mod3_eff2$Outcome)
mod3_eff2$model_epi <- sub("_PAR$", "", mod3_eff2$Predictor)

# check if joining will work (model level moderator - brain)
mod3_keys <- mod3_eff2 %>%
  distinct(model_brain)

mods_keys <- mods_brain_models %>%
  distinct(model_brain)

anti_join(
  mod3_keys,
  mods_keys,
  by = c("model_brain")
)

anti_join(
  mods_keys,
  mod3_keys,
  by = c("model_brain")
)

# check if joining will work (model level moderator - epi)
mod3_keys <- mod3_eff2 %>%
  distinct(model_epi)

mods_keys <- mods_epi_models %>%
  distinct(model_epi)

anti_join(
  mod3_keys,
  mods_keys,
  by = c("model_epi")
)

anti_join(
  mods_keys,
  mod3_keys,
  by = c("model_epi")
) # gestational clocls missing but that's fine

mod3_eff2 <- mod3_eff2 %>% left_join(mods_brain_models, by = c("model_brain"))
mod3_eff2 <- mod3_eff2 %>% left_join(mods_epi_models, by = c("model_epi"))

# rename some cohort names for nicer plotting & to align with performance 
mod3_eff2 <- mod3_eff2 %>%
  mutate(
    cohort = recode(
      cohort,
      "SAND" = "FFCWS",
      "Oregon_ADHD-1000" = "Oregon ADHD-1000",
      "Kids2Health_C" = "K2H Childhood"
    ))

# --- check for missing or implausible values ---
mod3_eff2 %>%
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
write.csv(mod3_eff2, file.path(results_dir, "mod3_eff2.csv"), row.names = FALSE)
saveRDS(mod3_eff2, file.path(results_dir, "mod3_eff2.rds"))

##### assoc descriptives general -----

assoc_desc <- mod3_eff2 %>%
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

write.csv(assoc_desc, file.path(results_dir, paste0("assoc_desc_model3_", "final", ".csv")), row.names = FALSE) # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

