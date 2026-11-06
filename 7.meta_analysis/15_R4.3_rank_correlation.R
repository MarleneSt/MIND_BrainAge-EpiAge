###### R4.3: do model rankings by MAEw and recalibrated MAE align? ------
# Spearman rank correlation between model-wise MAEw and recalibrated MAE, per
# domain (reported: rho = 0.90 for epigenetic clocks, 1 for brain age models).
# MAEw: model-specific estimates from 01_main_meta_wMAE.R (same values as SM Tables S9 and S11).
# Recalibrated MAE: R4.3_modelwise_offset_slope_recalMAE.csv (from 14_R4.3_recalibrated_metrics.R).

library(tidyverse)

in_dir  <- "results" # [REPRODUCIBILITY] relative path
out_dir <- "results/revision/R4.3" # [REPRODUCIBILITY] relative path
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE) # [REPRODUCIBILITY] create folder on a fresh run


#### data ----
# MAEw per model (as reported in SM Tables S9 and S11)
# [REPRODUCIBILITY] read from the model-wise output of 01_main_meta_wMAE.R
# (was: read_excel("SM_Tables.xlsx", sheet = "TableS9"/"TableS11", fixed cell ranges))
maew_epi <- read.csv(file.path(in_dir, "perf_epi_wMAE_noDunedinDNAmTL_modelwise_final/modelwise_wMAE_estimates.csv")) %>%
  transmute(domain = "epi", model, MAEw = wMAE)

maew_brain <- read.csv(file.path(in_dir, "perf_brain_wMAE_noDunedin_modelwise_final/modelwise_wMAE_estimates.csv")) %>%
  transmute(domain = "brain", model, MAEw = wMAE)

recal <- read.csv(file.path(out_dir, "R4.3_modelwise_offset_slope_recalMAE.csv")) %>%
  select(domain, model, recal_MAE)

dat <- bind_rows(maew_epi, maew_brain) %>%
  inner_join(recal, by = c("domain", "model"))

# check: all models matched (18 epi, 7 brain)
count(dat, domain)

# Spearman rank correlation per domain
rho <- dat %>%
  group_by(domain) %>%
  summarise(rho = cor(MAEw, recal_MAE, method = "spearman"), n = n())
rho

write.csv(rho, file.path(out_dir, "R4.3_rank_correlation_MAEw_vs_recalMAE.csv"), row.names = FALSE) # [REPRODUCIBILITY] saved (previously only shown in the plot title)
