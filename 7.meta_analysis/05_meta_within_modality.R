#### META-ANALYSES SCRIPT WITHIN MODALITY
#### MIND BrainAge~EpiAge
#### last updated: 06/07/2026

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
# [REPRODUCIBILITY] relative paths; run from the Reproducibilty folder
# (was: /Volumes/.../data/analysis_210526 and /Volumes/.../results/results_020726/within_modality)
home <- "data"
data_dir <- file.path(home)
results_dir <- file.path("results", "within_modality")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE) # [REPRODUCIBILITY] create folder on a fresh run (needed for "../perf_epi.rds" below)

#### some helpers ----

# Helper: I2 for multilevel metafor models
# code from Viechtbauer's metafor website: https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate

# note: I2 is currently not used throughout this script as I started to implement it from the correlation meta-analyses onwars 

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

# load performance files as needed for N 
perf_epi <- readRDS(file.path(results_dir, "../perf_epi.rds"))
perf_brain <- readRDS(file.path(results_dir, "../perf_brain.rds"))

perf_epi <- perf_epi %>%
  mutate(
    model = recode(
      model,
      "Bohlin_conv" = "Bohlin",
      "Knight_conv" = "Knight",
      "EPIC_conv" = "EPIC"
    )
  )

#### within-modality correlation data loading ----

# --- rename files to include cohort name ---
# RUN ONLY ONCE

# cohort_dirs <- list.dirs(data_dir, recursive = FALSE, full.names = TRUE)
# 
# rename_corr_files <- function(cohort_dir) {
#   
#   cohort <- basename(cohort_dir)
#   
#   corr_dir <- file.path(
#     cohort_dir,
#     "analyses/results/performance/corr_cov"
#   )
#   
#   if (!dir.exists(corr_dir)) return(NULL)
#   
#   files <- list.files(
#     corr_dir,
#     pattern = "\\.csv$",
#     full.names = TRUE
#   )
#   
#   # only rename files that do not already start with cohort name
#   files_to_rename <- files[
#     !str_starts(basename(files), paste0(cohort, "_"))
#   ]
#   
#   new_names <- file.path(
#     dirname(files_to_rename),
#     paste0(cohort, "_", basename(files_to_rename))
#   )
#   
#   file.rename(files_to_rename, new_names)
# }
# 
# purrr::walk(cohort_dirs, rename_corr_files)

# --- extract all corr/cov files ---

files_corr <- list.files(
  path = data_dir,
  full.names = TRUE,
  recursive = TRUE
)

# keep only standard performance/corr_cov files, not overlap/corr_cov
files_corr <- files_corr[
  grepl("/performance/corr_cov/", files_corr) &
    !grepl("/performance/overlap/corr_cov/", files_corr)
]

# keep only correlation csvs of interest
files_corr <- files_corr[
  grepl("_correlations_(PAR|predage)_pearson\\.csv$", files_corr)
]

# exclude already-processed long files
files_corr <- files_corr[
  !grepl("_long\\.csv$", files_corr)
]

# quick check
basename(files_corr)[1:20]

# parser 

parse_corr_file <- function(fname) {
  
  base <- basename(fname)
  
  m <- stringr::str_match(
    base,
    "^(.*?)_(brain|epi)_(.+?)_correlations_(PAR|predage)_pearson\\.csv$"
  )
  
  tibble::tibble(
    cohort = m[, 2],
    domain = m[, 3],
    timepoint = m[, 4],
    metric = m[, 5],
    file_path = fname
  )
}

# matrix to long helper 

read_corr_matrix_long <- function(fname) {
  
  raw <- read.csv(
    fname,
    check.names = FALSE
  )
  
  # first column contains row/model names, even if its header is blank
  names(raw)[1] <- "model_1"
  
  if (anyDuplicated(raw$model_1)) {
    warning("Duplicate row labels in: ", basename(fname))
    
    raw <- raw %>%
      mutate(model_1 = make.unique(as.character(model_1)))
  }
  
  raw %>%
    pivot_longer(
      cols = -model_1,
      names_to = "model_2",
      values_to = "r"
    ) %>%
    mutate(r = as.numeric(r)) %>%
    filter(model_1 != model_2) %>%
    mutate(
      pair_id = map2_chr(
        model_1,
        model_2,
        ~ paste(sort(c(.x, .y)), collapse = "___")
      )
    ) %>%
    distinct(pair_id, .keep_all = TRUE) %>%
    select(-pair_id)
}
# read 
corr_file_info <- purrr::map_dfr(
  files_corr,
  parse_corr_file
)

corr_file_info %>%
  filter(
    is.na(cohort) |
      is.na(domain) |
      is.na(timepoint) |
      is.na(metric)
  )

corr_all <- corr_file_info %>%
  mutate(
    data = purrr::map(file_path, read_corr_matrix_long)
  ) %>%
  select(-file_path) %>%
  tidyr::unnest(data)

corr_all_clean <- corr_all %>%
  mutate(
    model_1 = str_remove(model_1, "_PAR$"),
    model_2 = str_remove(model_2, "_PAR$"),
    model_1 = str_remove(model_1, "_predage$"),
    model_2 = str_remove(model_2, "_predage$")
  ) %>%
  filter(
    !model_1 %in% c("AGE_brain", "AGE_epi"),
    !model_2 %in% c("AGE_brain", "AGE_epi")
  )

# rename for nicer plotting / consistent model matching
corr_all_clean <- corr_all_clean %>%
  mutate(
    cohort = recode(
      cohort,
      "SAND" = "FFCWS",
      "Oregon_ADHD-1000" = "Oregon ADHD-1000",
      "Kids2Health_C" = "K2H Childhood",
      "Kids2Health_I" = "K2H Infancy",
      "UCI_Echo" = "UCI Echo"
    ),
    across(
      c(model_1, model_2),
      ~ recode(
        .x,
        "Bohlin_conv" = "Bohlin",
        "Knight_conv" = "Knight",
        "EPIC_conv" = "EPIC"
      )
    )
  )


# remove additional gestational models that are not needed 
corr_all_clean <- corr_all_clean %>%
  filter(!(model_1 %in% c("Mayne", "Lee.RPC", "Lee.CPC", "Lee.refRPC")))

corr_all_clean <- corr_all_clean %>%
  filter(!(model_2 %in% c("Mayne", "Lee.RPC", "Lee.CPC", "Lee.refRPC")))

# ensure that each model pair is only in there once per cohort & timepoint 
corr_all_clean <- corr_all_clean %>%
  mutate(
    model_a = map2_chr(
      model_1,
      model_2,
      ~ sort(c(.x, .y))[1]
    ),
    
    model_b = map2_chr(
      model_1,
      model_2,
      ~ sort(c(.x, .y))[2]
    )
  )

# now check if each model pair appears only ince 
corr_all_clean %>%
  count(cohort, domain, timepoint, metric, model_a, model_b) %>%
  filter(n > 1)

# check for duplicates 
dup_pairs <- corr_all_clean %>%
  count(
    cohort,
    domain,
    timepoint,
    metric,
    model_a,
    model_b
  ) %>%
  filter(n > 1)

dup_pairs


write.csv(
  corr_all_clean,
  file.path(results_dir, "correlations_all.csv"),
  row.names = FALSE
)

saveRDS(
  corr_all_clean,
  file.path(results_dir, "correlations_all.rds")
)

# some checks
corr_all_clean%>%
  count(cohort, timepoint, domain) %>%
  print(with = Inf, n = Inf)

# check for duplicates 
dup_check <- corr_all_clean %>%
  count(cohort, domain, timepoint, metric, model_a, model_b) %>%
  filter(n > 1)

dup_check

tibble(file_path = files_corr, file = basename(files_corr)) %>%
  count(file) %>%
  filter(n > 1)

corr_all_clean %>%
  count(domain, metric, timepoint) %>%
  arrange(domain, metric, timepoint)

# check impossible correlations 

# correlations outside -1 / 1 
corr_all_clean %>%
  filter(
    !is.na(r),
    (r < -1 | r > 1)
  )

# correlations = 1 
corr_all_clean %>%
  filter(abs(r) == 1)

# correlations = NA
corr_all_clean %>%
  summarise(
    n_missing_r = sum(is.na(r))
  )

# duplicated rows 
dup_corr <- corr_all_clean %>%
  count(
    cohort,
    domain,
    timepoint,
    metric,
    model_a,
    model_b
  ) %>%
  filter(n > 1)

dup_corr


##### reverse DNAmTL -----

corr_all_clean <- corr_all_clean %>%
  mutate(
    r = case_when(
      str_detect(model_a, "DNAmTL") |
        str_detect(model_b, "DNAmTL") ~ -r,
      TRUE ~ r
    )
  )


##### add in N ----

# I'm getting this from the performance dfs 

# --- prepare model-specific N lookup ---

# load perf_epi and perf_brain if not already in here

n_lookup <- bind_rows(
  perf_brain %>%
    select(cohort, domain, timepoint, model, N, mean_age),
  
  perf_epi %>%
    select(cohort, domain, timepoint, model, N, mean_age)
) %>%
  mutate(
    domain = as.character(domain),
    model = as.character(model),
    timepoint = as.character(timepoint),
    cohort = as.character(cohort)
  ) %>%
  distinct(cohort, domain, timepoint, model, .keep_all = TRUE)

# check correspondenc ebefore merging 
setdiff(unique(corr_all_clean$cohort), unique(n_lookup$cohort))
setdiff(unique(n_lookup$cohort), unique(corr_all_clean$cohort))

anti_join(
  corr_all_clean %>% distinct(cohort, domain, timepoint),
  n_lookup %>% distinct(cohort, domain, timepoint),
  by = c("cohort", "domain", "timepoint")
)

# check exact cohort/domain/timepoint/model matching
anti_join(
  corr_all_clean %>%
    select(cohort, domain, timepoint, model = model_a) %>%
    bind_rows(
      corr_all_clean %>%
        select(cohort, domain, timepoint, model = model_b)
    ) %>%
    distinct(),
  n_lookup %>%
    distinct(cohort, domain, timepoint, model),
  by = c("cohort", "domain", "timepoint", "model")
)


#  merge N into correlation dataframe 
# for now just N but could also have mean age? 
corr_all_meta <- corr_all_clean %>%
  left_join(
    n_lookup %>%
      rename(model_a = model, N_a = N),
    by = c("cohort", "domain", "timepoint", "model_a")
  ) %>%
  left_join(
    n_lookup %>%
      rename(model_b = model, N_b = N),
    by = c("cohort", "domain", "timepoint", "model_b")
  ) %>%
  mutate(
    N_pair = pmin(N_a, N_b, na.rm = TRUE)
  )

# any failed merges? 
corr_all_meta %>%
  filter(is.na(N_a) | is.na(N_b)) %>%
  count(domain, metric, model_a, model_b, sort = TRUE) %>% 
  print(width = Inf, n = Inf)

corr_all_meta %>%
  count(cohort, domain, timepoint, metric, model_a, model_b) %>%
  filter(n > 1)

# pair counts
pair_counts <- corr_all_meta %>%
  group_by(
    domain,
    metric,
    model_a,
    model_b
  ) %>%
  summarise(
    n_correlations = n(),
    n_cohorts = n_distinct(cohort),
    n_timepoints = n_distinct(paste(cohort, timepoint)),
    mean_N = mean(N_pair, na.rm = TRUE),
    min_N = min(N_pair, na.rm = TRUE),
    max_N = max(N_pair, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(domain, metric, model_a, model_b)

write.csv(
  pair_counts,
  file.path(results_dir, "within_modality_pair_counts.csv"),
  row.names = FALSE
)

##### Fisher's z -----
corr_all_meta <- corr_all_meta %>%
  filter(!is.na(N_pair), N_pair > 3, is.finite(r), abs(r) < 1)

corr_all_meta <- metafor::escalc(
  measure = "ZCOR",
  ri = r,
  ni = N_pair,
  data = corr_all_meta
) %>%
  rename(
    fisher_z = yi,
    fisher_var = vi
  ) %>%
  mutate(
    fisher_se = sqrt(fisher_var),
    pair = paste(model_a, model_b, sep = "___")
  )

##### create pair ID -----
corr_all_meta <- corr_all_meta %>%
  mutate(
    pair = paste(model_a, model_b, sep = "___")
  )

#### split datasets ----
corr_brain_PAR <- corr_all_meta %>%
  filter(domain == "brain", metric == "PAR")

corr_brain_predage <- corr_all_meta %>%
  filter(domain == "brain", metric == "predage")

corr_epi_PAR <- corr_all_meta %>%
  filter(domain == "epi", metric == "PAR")

corr_epi_predage <- corr_all_meta %>%
  filter(domain == "epi", metric == "predage")

#### ANALYSIS ----

##### function -----
run_within_modality_corr_meta <- function(dat, analysis_name, results_dir) {
  
  message("starting function")
  
  run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), "%Y%m%d_%H%M%S"))
  
  out_dir <- file.path(
    results_dir,
    paste0("within_modality_corr_", analysis_name, "_", run_id)
  )
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  log_file <- file.path(out_dir, "analysis_log.txt")
  
  message("global model")
  
  # global model
  global_fit <- rma.mv(
    yi = fisher_z,
    V  = fisher_var,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | pair
    ),
    data = dat,
    method = "REML"
  )
  
  message("global estimate")
  
  global_est <- tibble(
    analysis = analysis_name,
    fisher_z = as.numeric(coef(global_fit)[1]),
    fisher_se = global_fit$se[1],
    z_ci.lb = global_fit$ci.lb[1],
    z_ci.ub = global_fit$ci.ub[1],
    pooled_r = transf.ztor(fisher_z),
    ci.lb = transf.ztor(z_ci.lb),
    ci.ub = transf.ztor(z_ci.ub),
    pval = global_fit$pval[1]
  )
  
  message("pairwise model")
  
  # pair/model-combination model
  pair_fit <- rma.mv(
    yi = fisher_z,
    V  = fisher_var,
    mods = ~ pair - 1,
    random = list(
      ~ 1 | cohort/timepoint
    ),
    data = dat,
    method = "REML"
  )
  
  message("pair estimate")
  
  pair_est <- tibble(
    pair = str_remove(rownames(pair_fit$b), "^pair"),
    fisher_z = as.numeric(pair_fit$b),
    fisher_se = pair_fit$se,
    z_ci.lb = pair_fit$ci.lb,
    z_ci.ub = pair_fit$ci.ub,
    pval = pair_fit$pval
  ) %>%
    separate(pair, into = c("model_a", "model_b"), sep = "___", remove = FALSE) %>%
    mutate(
      pooled_r = transf.ztor(fisher_z),
      ci.lb = transf.ztor(z_ci.lb),
      ci.ub = transf.ztor(z_ci.ub)
    )
  
  message("pair counts")
  
  pair_counts <- dat %>%
    group_by(domain, metric, model_a, model_b, pair) %>%
    summarise(
      k = n(),
      n_cohorts = n_distinct(cohort),
      n_timepoints = n_distinct(paste(cohort, timepoint)),
      mean_N = mean(N_pair, na.rm = TRUE),
      min_N = min(N_pair, na.rm = TRUE),
      max_N = max(N_pair, na.rm = TRUE),
      .groups = "drop"
    )
  
  pair_est <- pair_est %>%
    left_join(pair_counts, by = c("pair", "model_a", "model_b"))
  
  message("writing csvs")
  
  write.csv(dat, file.path(out_dir, "input_data.csv"), row.names = FALSE)
  write.csv(global_est, file.path(out_dir, "global_pooled_correlation.csv"), row.names = FALSE)
  write.csv(pair_est, file.path(out_dir, "pairwise_pooled_correlations.csv"), row.names = FALSE)
  write.csv(pair_counts, file.path(out_dir, "pair_counts.csv"), row.names = FALSE)
  
  message("saving rds")
  
  saveRDS(global_fit, file.path(out_dir, "global_model.rds"))
  saveRDS(pair_fit, file.path(out_dir, "pairwise_model.rds"))
  
  message("writing log")
  
  sink(log_file)
  cat("Within-modality correlation meta-analysis\n")
  cat("========================================\n\n")
  cat("Analysis:", analysis_name, "\n")
  cat("Run ID:", run_id, "\n\n")
  
  cat("Global model\n")
  cat("============\n")
  print(summary(global_fit))
  print(global_est, width = Inf, n = Inf)
  
  cat("\n\nPairwise model\n")
  cat("==============\n")
  print(summary(pair_fit))
  
  cat("\n\nPairwise pooled correlations\n")
  cat("============================\n")
  print(pair_est, n = Inf, width = Inf)
  
  sink()
  
  list(
    out_dir = out_dir,
    global_fit = global_fit,
    pair_fit = pair_fit,
    global_est = global_est,
    pair_est = pair_est
  )
  message("done")
}

##### apply function ----
res_predage_brain <- run_within_modality_corr_meta(dat = corr_brain_predage, analysis_name = "brain_predage", results_dir)

res_predage_epi <- run_within_modality_corr_meta(dat = corr_epi_predage, analysis_name = "epi_predage", results_dir)

res_PAR_brain <- run_within_modality_corr_meta(dat = corr_brain_PAR, analysis_name = "brain_PAR", results_dir)

res_PAR_epi <- run_within_modality_corr_meta(dat = corr_epi_PAR, analysis_name = "epi_PAR", results_dir)

##### longer version function -----

run_within_modality_corr_meta_detailed <- function(dat, analysis_name, results_dir) {
  
  run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), "%Y%m%d_%H%M%S"))
  
  out_dir <- file.path(
    results_dir,
    paste0("longer_within_modality_corr_", analysis_name, "_", run_id)
  )
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(out_dir, "analysis_log.txt")
  
  # --- global model ---
  
  global_fit <- rma.mv(
    yi = fisher_z,
    V  = fisher_var,
    random = list(
      ~ 1 | cohort/timepoint,
      ~ 1 | pair
    ),
    data = dat,
    method = "REML"
  )
  
  I2_global <- calc_I2_mv(
    global_fit,
    component_names = c("cohort_timepoint", "pair")
  )
  
  I2_global_df <- tibble(
    component = c("total", names(I2_global$I2_components)),
    I2 = c(
      I2_global$I2_total,
      as.numeric(I2_global$I2_components)
    )
  )
  
  global_est <- tibble(
    analysis = analysis_name,
    fisher_z = as.numeric(coef(global_fit)[1]),
    fisher_se = global_fit$se[1],
    z_ci.lb = global_fit$ci.lb[1],
    z_ci.ub = global_fit$ci.ub[1],
    pooled_r = transf.ztor(fisher_z),
    ci.lb = transf.ztor(z_ci.lb),
    ci.ub = transf.ztor(z_ci.ub),
    pval = global_fit$pval[1]
  )
  
  # --- pair/model-combination model ---
  
  pair_fit <- rma.mv(
    yi = fisher_z,
    V  = fisher_var,
    mods = ~ pair - 1,
    random = list(
      ~ 1 | cohort/timepoint
    ),
    data = dat,
    method = "REML"
  )
  
  I2_pair <- calc_I2_mv(
    pair_fit,
    component_names = c("cohort_timepoint")
  )
  
  I2_pair_df <- tibble(
    component = c("total", names(I2_pair$I2_components)),
    I2 = c(
      I2_pair$I2_total,
      as.numeric(I2_pair$I2_components)
    )
  )
  
  pair_est <- tibble(
    pair = str_remove(rownames(pair_fit$b), "^pair"),
    fisher_z = as.numeric(pair_fit$b),
    fisher_se = pair_fit$se,
    z_ci.lb = pair_fit$ci.lb,
    z_ci.ub = pair_fit$ci.ub,
    pval = pair_fit$pval
  ) %>%
    separate(pair, into = c("model_a", "model_b"), sep = "___", remove = FALSE) %>%
    mutate(
      pooled_r = transf.ztor(fisher_z),
      ci.lb = transf.ztor(z_ci.lb),
      ci.ub = transf.ztor(z_ci.ub)
    )
  
  pair_counts <- dat %>%
    group_by(domain, metric, model_a, model_b, pair) %>%
    summarise(
      k = n(),
      n_cohorts = n_distinct(cohort),
      n_timepoints = n_distinct(paste(cohort, timepoint)),
      mean_N = mean(N_pair, na.rm = TRUE),
      min_N = min(N_pair, na.rm = TRUE),
      max_N = max(N_pair, na.rm = TRUE),
      .groups = "drop"
    )
  
  pair_est <- pair_est %>%
    left_join(pair_counts, by = c("pair", "model_a", "model_b"))
  
  # --- leave-one-cohort-out: global ---
  
  cohorts <- unique(dat$cohort)
  
  loo_global <- lapply(cohorts, function(cc) {
    
    tmp <- dat %>%
      filter(cohort != cc) %>%
      droplevels()
    
    fit <- rma.mv(
      yi = fisher_z,
      V  = fisher_var,
      random = list(
        ~ 1 | cohort/timepoint,
        ~ 1 | pair
      ),
      data = tmp,
      method = "REML"
    )
    
    z_est <- coef(fit)[1]
    z_se <- sqrt(vcov(fit)[1, 1])
    
    tibble(
      cohort_removed = cc,
      fisher_z = z_est,
      fisher_se = z_se,
      z_ci.lb = z_est - 1.96 * z_se,
      z_ci.ub = z_est + 1.96 * z_se,
      pooled_r = transf.ztor(z_est),
      ci.lb = transf.ztor(z_est - 1.96 * z_se),
      ci.ub = transf.ztor(z_est + 1.96 * z_se)
    )
  }) %>%
    bind_rows()
  
  loo_global_summary <- tibble(
    min_r = min(loo_global$pooled_r, na.rm = TRUE),
    max_r = max(loo_global$pooled_r, na.rm = TRUE),
    range_r = max_r - min_r,
    sd_r = sd(loo_global$pooled_r, na.rm = TRUE),
    min_z = min(loo_global$fisher_z, na.rm = TRUE),
    max_z = max(loo_global$fisher_z, na.rm = TRUE),
    range_z = max_z - min_z,
    sd_z = sd(loo_global$fisher_z, na.rm = TRUE)
  )
  
  # --- leave-one-cohort-out: pairwise ---
  
  loo_pair <- lapply(cohorts, function(cc) {
    
    tmp <- dat %>%
      filter(cohort != cc) %>%
      droplevels()
    
    fit <- rma.mv(
      yi = fisher_z,
      V  = fisher_var,
      mods = ~ pair - 1,
      random = list(
        ~ 1 | cohort/timepoint
      ),
      data = tmp,
      method = "REML"
    )
    
    tibble(
      cohort_removed = cc,
      pair = str_remove(rownames(fit$b), "^pair"),
      fisher_z = as.numeric(fit$b),
      fisher_se = fit$se,
      z_ci.lb = fit$ci.lb,
      z_ci.ub = fit$ci.ub,
      pval = fit$pval
    ) %>%
      separate(pair, into = c("model_a", "model_b"), sep = "___", remove = FALSE) %>%
      mutate(
        pooled_r = transf.ztor(fisher_z),
        ci.lb = transf.ztor(z_ci.lb),
        ci.ub = transf.ztor(z_ci.ub)
      )
  }) %>%
    bind_rows()
  
  loo_pair_summary <- loo_pair %>%
    group_by(pair, model_a, model_b) %>%
    summarise(
      min_r = min(pooled_r, na.rm = TRUE),
      max_r = max(pooled_r, na.rm = TRUE),
      range_r = max_r - min_r,
      sd_r = sd(pooled_r, na.rm = TRUE),
      min_z = min(fisher_z, na.rm = TRUE),
      max_z = max(fisher_z, na.rm = TRUE),
      range_z = max_z - min_z,
      sd_z = sd(fisher_z, na.rm = TRUE),
      .groups = "drop"
    )
  
  # --- plots ---
  
  full_r <- transf.ztor(coef(global_fit)[1])
  
  png(
    filename = file.path(out_dir, "plot_global_loo.png"),
    width = 1800,
    height = 1200,
    res = 300
  )
  
  print(
    ggplot(
      loo_global,
      aes(
        x = reorder(cohort_removed, pooled_r),
        y = pooled_r
      )
    ) +
      geom_point(size = 3) +
      geom_errorbar(
        aes(ymin = ci.lb, ymax = ci.ub),
        width = 0.15
      ) +
      geom_hline(
        yintercept = full_r,
        linetype = "dashed"
      ) +
      coord_flip() +
      labs(
        x = "Cohort removed",
        y = "Pooled Pearson r",
        title = paste0("Leave-one-cohort-out sensitivity: ", analysis_name)
      ) +
      theme_minimal(base_size = 14)
  )
  
  dev.off()
  
  png(
    filename = file.path(out_dir, "plot_pairwise_loo.png"),
    width = 2600,
    height = 1800,
    res = 300
  )
  
  print(
    ggplot(
      loo_pair,
      aes(
        x = reorder(cohort_removed, pooled_r),
        y = pooled_r
      )
    ) +
      geom_point(size = 1.8) +
      geom_errorbar(
        aes(ymin = ci.lb, ymax = ci.ub),
        width = 0.12
      ) +
      facet_wrap(~ pair, scales = "free_y") +
      coord_flip() +
      labs(
        x = "Cohort removed",
        y = "Pooled Pearson r",
        title = paste0("Pairwise leave-one-cohort-out sensitivity: ", analysis_name)
      ) +
      theme_minimal(base_size = 9)
  )
  
  dev.off()
  
  # --- save outputs ---
  
  write.csv(dat, file.path(out_dir, "input_data.csv"), row.names = FALSE)
  write.csv(global_est, file.path(out_dir, "global_pooled_correlation.csv"), row.names = FALSE)
  write.csv(pair_est, file.path(out_dir, "pairwise_pooled_correlations.csv"), row.names = FALSE)
  write.csv(pair_counts, file.path(out_dir, "pair_counts.csv"), row.names = FALSE)
  write.csv(I2_global_df, file.path(out_dir, "I2_global.csv"), row.names = FALSE)
  write.csv(I2_pair_df, file.path(out_dir, "I2_pairwise.csv"), row.names = FALSE)
  write.csv(loo_global, file.path(out_dir, "leave_one_cohort_out_global.csv"), row.names = FALSE)
  write.csv(loo_global_summary, file.path(out_dir, "leave_one_cohort_out_global_summary.csv"), row.names = FALSE)
  write.csv(loo_pair, file.path(out_dir, "leave_one_cohort_out_pairwise.csv"), row.names = FALSE)
  write.csv(loo_pair_summary, file.path(out_dir, "leave_one_cohort_out_pairwise_summary.csv"), row.names = FALSE)
  
  saveRDS(global_fit, file.path(out_dir, "global_model.rds"))
  saveRDS(pair_fit, file.path(out_dir, "pairwise_model.rds"))
  
  # --- log ---
  
  sink(log_file)
  
  cat("Within-modality correlation meta-analysis\n")
  cat("========================================\n\n")
  cat("Analysis:", analysis_name, "\n")
  cat("Run ID:", run_id, "\n")
  cat("Output directory:", out_dir, "\n\n")
  
  cat("Global model\n")
  cat("============\n")
  print(summary(global_fit))
  
  cat("\n\nGlobal estimate, back-transformed to Pearson r\n")
  cat("==============================================\n")
  print(global_est)
  
  cat("\n\nGlobal I2\n")
  cat("=========\n")
  print(I2_global_df)
  
  cat("\n\nPairwise model\n")
  cat("==============\n")
  print(summary(pair_fit))
  
  cat("\n\nPairwise I2\n")
  cat("===========\n")
  print(I2_pair_df)
  
  cat("\n\nPairwise pooled correlations\n")
  cat("============================\n")
  print(pair_est, n = Inf, width = Inf)
  
  cat("\n\nPair counts\n")
  cat("===========\n")
  print(pair_counts, n = Inf, width = Inf)
  
  cat("\n\nLeave-one-cohort-out global\n")
  cat("===========================\n")
  print(loo_global, n = Inf, width = Inf)
  
  cat("\n\nLeave-one-cohort-out global summary\n")
  cat("===================================\n")
  print(loo_global_summary)
  
  cat("\n\nLeave-one-cohort-out pairwise\n")
  cat("=============================\n")
  print(loo_pair, n = Inf, width = Inf)
  
  cat("\n\nLeave-one-cohort-out pairwise summary\n")
  cat("=====================================\n")
  print(loo_pair_summary, n = Inf, width = Inf)
  
  cat("\n\nSession info\n")
  cat("============\n")
  print(sessionInfo())
  
  sink()
  
  list(
    out_dir = out_dir,
    global_fit = global_fit,
    pair_fit = pair_fit,
    global_est = global_est,
    pair_est = pair_est,
    pair_counts = pair_counts,
    I2_global = I2_global_df,
    I2_pair = I2_pair_df,
    loo_global = loo_global,
    loo_pair = loo_pair
  )
}

##### run detailed function -----
res_predage_brain <- run_within_modality_corr_meta_detailed(
  corr_all_meta %>% filter(domain == "brain", metric == "predage"),
  "predage_brain",
  results_dir
)

res_predage_epi <- run_within_modality_corr_meta_detailed(
  corr_all_meta %>% filter(domain == "epi", metric == "predage"),
  "predage_epi",
  results_dir
)

res_PAR_brain <- run_within_modality_corr_meta_detailed(
  corr_all_meta %>% filter(domain == "brain", metric == "PAR"),
  "PAR_brain",
  results_dir
)

res_PAR_epi <- run_within_modality_corr_meta_detailed(
  corr_all_meta %>% filter(domain == "epi", metric == "PAR"),
  "PAR_epi",
  results_dir
)
