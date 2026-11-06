###### make_cohort_age_sd.R: one-off extraction of the SD of chronological age ------
# Creates data/cohort_age_sd.csv (input for 14_R4.3_recalibrated_metrics.R).
# The SD of age per cohort x timepoint is not stored in perf_*.rds; it comes from
# SM Table S5 (age as "mean (SD)": age at MRI for brain, age at DNAm for epi).
# Code copied from 14_R4.3_recalibrated_metrics.R (Table S5 section).
# Run once from the Reproducibilty folder, with SM_Tables.xlsx in that folder.

library(tidyverse)
library(readxl)

s5 <- read_excel("SM_Tables.xlsx", sheet = "TableS5", skip = 3) %>%
  filter(modality %in% c("brain", "epi")) %>%
  mutate(
    age_str = if_else(modality == "brain", `age at MRI`, `age at DNAm`),
    sd_age  = as.numeric(str_match(age_str, "\\(\\s*([0-9.]+)\\s*\\)")[, 2]),
    # cohort names that are spelled differently in Table S5
    cohort  = recode(cohort,
                     "Kids2Health_childhood" = "K2H Childhood",
                     "Kids2Health_infancy"   = "K2H Infancy",
                     "UCI ECHO"              = "UCI Echo"),
    tp_key  = tolower(timepoint)
  ) %>%
  select(cohort, tp_key, domain = modality, sd_age)

write.csv(s5, "data/cohort_age_sd.csv", row.names = FALSE)
