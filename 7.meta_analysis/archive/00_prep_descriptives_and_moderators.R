library(tidyverse)
library(readr)
library(openxlsx)

base_dir <- "data" # [REPRODUCIBILITY] relative path
results_path <- "results/" # [REPRODUCIBILITY] relative path

files <- list.files(
  path = base_dir,
  pattern = "descriptives_.*_wins\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

# function below written by ChatGPT will require some double checking / spot checks 

extract_info <- function(file) {
  fname <- basename(file)
  
  relative_path <- str_remove(file, paste0("^", base_dir, "/"))
  cohort <- str_split(relative_path, "/", simplify = TRUE)[,1]
  
  dat <- read_csv(file, show_col_types = FALSE)
  
  out <- tibble(
    cohort = cohort,
    full_file_name = fname,
    full_path = file
  )
  
  # categorical brain/epi sex files
  if (str_detect(fname, "^descriptives_categorical_(brain|epi)_.*_wins\\.csv$")) {
    group <- str_match(fname, "^descriptives_categorical_(brain|epi)_(.*)_wins\\.csv$")[,2]
    timepoint <- str_match(fname, "^descriptives_categorical_(brain|epi)_(.*)_wins\\.csv$")[,3]
    
    sex_var <- paste0("SEX_", group)
    
    sex_counts <- dat %>%
      filter(variable == sex_var) %>%
      mutate(
        sex = case_when(
          level == "0" ~ "male",
          level == "1" ~ "female",
          TRUE ~ paste0("level_", level)
        )
      ) %>%
      select(sex, count) %>%
      pivot_wider(names_from = sex, values_from = count, names_prefix = "n_")
    
    return(bind_cols(out, tibble(file_type = "sex", group = group, timepoint = timepoint), sex_counts))
  }
  
  # categorical overlap sex files
  if (str_detect(fname, "^descriptives_categorical_overlap_brain_epi_.*_wins\\.csv$")) {
    
    timepoint <- str_match(
      fname,
      "^descriptives_categorical_overlap_brain_epi_(.*)_wins\\.csv$"
    )[,2]
    
    sex_counts <- dat %>%
      filter(variable %in% c("SEX_brain", "SEX_epi")) %>%
      mutate(
        sex = case_when(
          level == "0" ~ "male",
          level == "1" ~ "female",
          TRUE ~ paste0("level_", level)
        )
      ) %>%
      mutate(variable = paste0(variable, "_", sex)) %>%
      select(variable, count) %>%
      pivot_wider(
        names_from = variable,
        values_from = count
      )
    
    return(
      bind_cols(
        out,
        tibble(
          file_type = "overlap_sex",
          group = "brain_epi",
          timepoint = timepoint
        ),
        sex_counts
      )
    )
  }
  
  # numeric brain/epi age files
  if (str_detect(fname, "^descriptives_numeric_(brain|epi)_.*_wins\\.csv$")) {
    group <- str_match(fname, "^descriptives_numeric_(brain|epi)_(.*)_wins\\.csv$")[,2]
    timepoint <- str_match(fname, "^descriptives_numeric_(brain|epi)_(.*)_wins\\.csv$")[,3]
    
    age_var <- paste0("AGE_", group)
    
    age_stats <- dat %>%
      filter(variable == age_var) %>%
      transmute(
        age_mean = mean,
        age_min = min,
        age_max = max
      )
    
    return(bind_cols(out, tibble(file_type = "age", group = group, timepoint = timepoint), age_stats))
  }
  
  # overlap numeric files
  if (str_detect(fname, "^descriptives_numeric_overlap_brain_epi_.*_wins\\.csv$")) {
    timepoint <- str_match(fname, "^descriptives_numeric_overlap_brain_epi_(.*)_wins\\.csv$")[,2]
    
    overlap_stats <- dat %>%
      filter(variable %in% c(
        "AGE_brain",
        "AGE_epi",
        "AGE_diff_brain_minus_epi",
        "AGE_diff_abs"
      )) %>%
      select(variable, mean, sd, min, max) %>%
      pivot_wider(
        names_from = variable,
        values_from = c(mean, sd, min, max),
        names_glue = "{variable}_{.value}"
      )
    
    return(bind_cols(out, tibble(file_type = "overlap_age", group = "brain_epi", timepoint = timepoint), overlap_stats))
  }
  
  NULL
}

results <- map_dfr(files, extract_info)

write.xlsx(results, file.path(results_path, "collated_descriptives.xlsx"))


# 6. helper function
add_missing_cols <- function(df, cols) {
  
  missing <- setdiff(cols, names(df))
  
  df[missing] <- NA
  
  df %>%
    select(all_of(cols))
}


# 7. build brain table

brain_cols <- c(
  "cohort",
  "timepoint",
  "mean_age",
  "min_age",
  "max_age",
  "age_range",
  "overlap_age",
  "n_male",
  "n_female"
)

moderators_cohorts_brain <- results %>%
  filter(file_type %in% c("age", "sex"),
         group == "brain") %>%
  group_by(cohort, timepoint) %>%
  summarise(
    mean_age = first(na.omit(age_mean)),
    min_age = first(na.omit(age_min)),
    max_age = first(na.omit(age_max)),
    n_male = first(na.omit(n_male)),
    n_female = first(na.omit(n_female)),
    .groups = "drop"
  ) %>%
  mutate(
    age_range = max_age - min_age,
    overlap_age = NA
  ) %>%
  add_missing_cols(brain_cols)


# 8. build epi table

epi_cols <- c(
  "cohort",
  "timepoint",
  "mean_age",
  "min_age",
  "max_age",
  "age_range",
  "overlap_age",
  "tissue_test",
  "n_male",
  "n_female"
)

moderators_cohorts_epi <- results %>%
  filter(file_type %in% c("age", "sex"),
         group == "epi") %>%
  group_by(cohort, timepoint) %>%
  summarise(
    mean_age = first(na.omit(age_mean)),
    min_age = first(na.omit(age_min)),
    max_age = first(na.omit(age_max)),
    n_male = first(na.omit(n_male)),
    n_female = first(na.omit(n_female)),
    .groups = "drop"
  ) %>%
  mutate(
    age_range = max_age - min_age,
    overlap_age = NA,
    tissue_test = NA
  ) %>%
  add_missing_cols(epi_cols)


# 9. build overlap table

overlap_cols <- c(
  "cohort",
  "timepoint",
  "mean_age_brain",
  "min_age_brain",
  "max_age_brain",
  "age_range_brain",
  "mean_age_epi",
  "min_age_epi",
  "max_age_epi",
  "age_range_epi",
  "age_range",
  "age_diff_brain_minus_epi_mean",
  "age_diff_abs_mean",
  "n_male",
  "n_female",
  "cortical",
  "saliva",
  "buccal"
)

moderators_cohorts_overlap <- results %>%
  filter(file_type %in% c("overlap_age", "overlap_sex")) %>%
  group_by(cohort, timepoint) %>%
  summarise(
    mean_age_brain = first(na.omit(AGE_brain_mean)),
    min_age_brain = first(na.omit(AGE_brain_min)),
    max_age_brain = first(na.omit(AGE_brain_max)),
    
    mean_age_epi = first(na.omit(AGE_epi_mean)),
    min_age_epi = first(na.omit(AGE_epi_min)),
    max_age_epi = first(na.omit(AGE_epi_max)),
    
    age_diff_brain_minus_epi_mean =
      first(na.omit(AGE_diff_brain_minus_epi_mean)),
    
    age_diff_abs_mean =
      first(na.omit(AGE_diff_abs_mean)),
    
    n_male =
      first(na.omit(SEX_brain_male)),
    
    n_female =
      first(na.omit(SEX_brain_female)),
    
    .groups = "drop"
  ) %>%
  mutate(
    age_range_brain = max_age_brain - min_age_brain,
    age_range_epi = max_age_epi - min_age_epi,
    age_range = NA,
    cortical = NA,
    saliva = NA,
    buccal = NA
  ) %>%
  add_missing_cols(overlap_cols)


# 10. export csv files

write_csv(
  moderators_cohorts_brain, file.path(results_path,
  "moderators_cohorts_brain.csv")
)

write_csv(
  moderators_cohorts_epi, file.path(results_path,
  "moderators_cohorts_epi.csv")
)

write_csv(
  moderators_cohorts_overlap, file.path(results_path,
  "moderators_cohorts_overlap.csv")
)

# note: currently does not extract sex/gender info when they are not coded 0/1 - could add something to also extract when coded F/M or Female/Male
