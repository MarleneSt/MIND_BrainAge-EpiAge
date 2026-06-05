library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# Load data
df_surf <- read_csv("CorticalMeasuresENIGMA_SurfAvg_NAtmp.csv")
df_rh <- read_table("aparc_stats_vol_rh.txt")
df_lh <- read_table("aparc_stats_vol_lh.txt")

# align subjID 
names(df_rh)[names(df_rh) == "rh.aparc.volume"] <- "SubjID"
names(df_lh)[names(df_lh) == "lh.aparc.volume"] <- "SubjID"

# Identify NA cells
na_df <- df_surf %>%
  mutate(SubjID = as.character(SubjID)) %>%
  pivot_longer(-SubjID, names_to = "region", values_to = "value") %>%
  filter(is.na(value))

# Convert region names to match right/left hemisphere volume data
na_df <- na_df %>%
  mutate(
    hemisphere = case_when(
      str_starts(region, "L_") ~ "lh",
      str_starts(region, "R_") ~ "rh",
      TRUE ~ NA_character_
    ),
    region_clean = region %>%
      str_remove("^[LR]_") %>%
      str_remove("_surfavg$"),  # <- THIS removes the '_surfavg' suffix
    region_volume = paste0(hemisphere, "_", region_clean, "_volume")
  )


# Apply NA to right hemisphere volume data
for (i in 1:nrow(na_df)) {
  subj <- na_df$SubjID[i]
  region <- na_df$region_volume[i]
  
  if (!is.na(region) && region %in% names(df_rh)) {
    df_rh[df_rh$rh.aparc.volume == subj, region] <- NA
  }
}

# Apply NA to left hemisphere volume data
for (i in 1:nrow(na_df)) {
  subj <- na_df$SubjID[i]
  region <- na_df$region_volume[i]
  
  if (!is.na(region) && region %in% names(df_lh)) {
    df_lh[df_lh$lh.aparc.volume == subj, region] <- NA
  }
}

# Optionally write updated files
write.table(df_rh, "aparc_stats_vol_rh_withNAs.txt", row.names = FALSE)
write.table(df_lh, "TMP_aparc_stats_vol_lh_withNAs.txt", row.names = FALSE)
