#### META-ANALYSES SCRIPT Model 1 + 3
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
home <- "data" # [REPRODUCIBILITY] relative path
data_dir <- file.path(home)
results_dir <- "results/" # [REPRODUCIBILITY] relative path

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

#### load data ----


# load association files based on model 1 and 3 
mod1_eff2 <- readRDS(file.path(results_dir, "mod1_eff2.rds"))
mod3_eff2 <- readRDS(file.path(results_dir, "mod3_eff2.rds"))

#### M3 PAR ----

results_dir <- "results/assoc_model3/" # [REPRODUCIBILITY] relative path

###### assoc global -----

# --- set up folder ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_assoc <- file.path(
  results_dir,
  paste0("assoc_model3_PAR_global_", run_id)
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
  data = mod3_eff2,
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
  k = nrow(mod3_eff2),
  n_cohorts = n_distinct(mod3_eff2$cohort),
  n_timepoints = n_distinct(paste(mod3_eff2$cohort, mod3_eff2$timepoint)),
  n_model_combis = n_distinct(mod3_eff2$model_combi)
)

write.csv(
  effect_sizes_overall,
  file.path(results_dir_assoc, "effect_sizes_overall.csv"),
  row.names = FALSE
)

# --- LOO cohort ---

cohorts <- unique(mod3_eff2$cohort)

loo_cohort <- lapply(cohorts, function(cc) {
  
  tmp <- subset(mod3_eff2, cohort != cc)
  
  # [REPRODUCIBILITY] if the default optimizer (nlminb) does not converge for a
  # refit, retry with BFGS (otherwise one failed refit stops the whole script)
  fit <- tryCatch(
    rma.mv(
      yi = RLM_Estimate_scaled,
      V  = RLM_SE_coeftest_HC3_scaled^2,
      random = list(
        ~ 1 | cohort/timepoint,
        ~ 1 | model_combi
      ),
      data = tmp,
      method = "REML"
    ),
    error = function(e) rma.mv(
      yi = RLM_Estimate_scaled,
      V  = RLM_SE_coeftest_HC3_scaled^2,
      random = list(
        ~ 1 | cohort/timepoint,
        ~ 1 | model_combi
      ),
      data = tmp,
      method = "REML",
      control = list(optimizer = "optim", optmethod = "BFGS")
    )
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
  paste0("assoc_model3_PAR_modelwise_", run_id_modelwise)
)

dir.create(results_dir_assoc_modelwise, recursive = TRUE, showWarnings = FALSE)

log_file_modelwise <- file.path(results_dir_assoc_modelwise, "analysis_log_modelwise.txt")


# --- meta ---

mod2_res_mv_modelwise <- rma.mv(
  yi   = RLM_Estimate_scaled,
  V    = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ model_combi - 1,
  random = ~ 1 | cohort/timepoint,
  data   = mod3_eff2,
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

# cohorts <- unique(mod3_eff2$cohort)
# 
# loo_cohort_modelwise <- lapply(cohorts, function(cc) {
#   
#   message("Running LOO modelwise: removing ", cc)
#   
#   tmp <- subset(mod3_eff2, cohort != cc)
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

# cohorts <- unique(mod3_eff2$cohort)
# 
# loo_cohort_modelwise <- lapply(cohorts, function(cc) {
#   
#   tmp <- subset(mod3_eff2, cohort != cc)
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
cohorts <- unique(mod3_eff2$cohort)

loo_cohort_modelwise <- lapply(cohorts, function(cc) {
  
  message("Running LOO modelwise: removing ", cc)
  
  tmp <- mod3_eff2 %>%
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

results_dir_assoc_global_ageint_brain <- file.path(results_dir, paste0("assoc_model3_PAR_global_ageint_brain_", run_id))

dir.create(results_dir_assoc_global_ageint_brain, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_assoc_global_ageint_brain, "analysis_log_global_ageint.txt")

# --- data prep ---

# center age 
mod3_eff2 <- mod3_eff2 %>%
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
  data   = mod3_eff2,
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
    min(mod3_eff2$mean_age_brain_c, na.rm = TRUE),
    max(mod3_eff2$mean_age_brain_c, na.rm = TRUE),
    length.out = 100
  )
)

pred <- predict(mod2_res_mv_global_ageint_brain, newmods = newdat$mean_age_brain_c)

newdat$pred  <- pred$pred
newdat$ci.lb <- pred$ci.lb
newdat$ci.ub <- pred$ci.ub

newdat$mean_age_brain <- newdat$mean_age_brain_c + mean(mod3_eff2$mean_age_brain, na.rm = TRUE)

png(
  filename = file.path(results_dir_assoc_global_ageint_brain, "association_global_meanage_meta_regression.png"),
  width = 1800, height = 1200, res = 300
)

ggplot() +
  geom_point(
    data = mod3_eff2,
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

cohorts <- unique(mod3_eff2$cohort)

loo_assoc_ageint <- lapply(cohorts, function(cc) {
  
  tmp <- subset(mod3_eff2, cohort != cc)
  
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

results_dir_assoc_global_ageint_epi <- file.path(results_dir, paste0("assoc_model3_PAR_global_ageint_epi_", run_id))

dir.create(results_dir_assoc_global_ageint_epi, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_assoc_global_ageint_epi, "analysis_log_global_ageint.txt")

# --- data prep ---

# center age 
mod3_eff2 <- mod3_eff2 %>%
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
  data   = mod3_eff2,
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
    min(mod3_eff2$mean_age_epi_c, na.rm = TRUE),
    max(mod3_eff2$mean_age_epi_c, na.rm = TRUE),
    length.out = 100
  )
)

pred <- predict(mod2_res_mv_global_ageint_epi, newmods = newdat$mean_age_epi_c)

newdat$pred  <- pred$pred
newdat$ci.lb <- pred$ci.lb
newdat$ci.ub <- pred$ci.ub

newdat$mean_age_epi <- newdat$mean_age_epi_c + mean(mod3_eff2$mean_age_epi, na.rm = TRUE)

png(
  filename = file.path(results_dir_assoc_global_ageint_epi, "association_global_meanage_meta_regression.png"),
  width = 1800, height = 1200, res = 300
)

ggplot() +
  geom_point(
    data = mod3_eff2,
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

cohorts <- unique(mod3_eff2$cohort)

loo_assoc_ageint <- lapply(cohorts, function(cc) {
  
  tmp <- subset(mod3_eff2, cohort != cc)
  
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
  paste0("assoc_model3_PAR_modelwise_ageint_brain_", run_id)
)

dir.create(results_dir_assoc_age_brain_modelwise, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  results_dir_assoc_age_brain_modelwise,
  "analysis_log_modelwise_ageint_brain.txt"
)


# --- data prep ---

mod3_eff2_age_brain <- mod3_eff2 %>%
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
  data = mod3_eff2_age_brain,
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

cohorts <- unique(mod3_eff2_age_brain$cohort)

loo_pairs_age_brain <- lapply(cohorts, function(cc) {
  
  message("Running LOO modelwise age interaction: removing ", cc)
  
  tmp <- mod3_eff2_age_brain %>%
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

# cohorts <- unique(mod3_eff2_age_brain$cohort)
# 
# loo_pairs_age_brain <- lapply(cohorts, function(cc) {
#   
#   tmp <- mod3_eff2_age_brain %>%
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
  paste0("assoc_model3_PAR_modelwise_ageint_epi_", run_id)
)

dir.create(results_dir_assoc_age_epi_modelwise, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  results_dir_assoc_age_epi_modelwise,
  "analysis_log_modelwise_ageint_epi.txt"
)


# --- data prep ---

mod3_eff2_age_epi <- mod3_eff2 %>%
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
  data = mod3_eff2_age_epi,
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

cohorts <- unique(mod3_eff2_age_epi$cohort)

loo_pairs_age_epi <- lapply(cohorts, function(cc) {
  
  message("Running LOO modelwise age interaction: removing ", cc)
  
  tmp <- mod3_eff2_age_epi %>%
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

# cohorts <- unique(mod3_eff2_age_epi$cohort)
# 
# loo_pairs_age_epi <- lapply(cohorts, function(cc) {
#   
#   tmp <- mod3_eff2_age_epi %>%
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

#### M1 PAR ----

results_dir <- "results/assoc_model1/" # [REPRODUCIBILITY] relative path

###### assoc global -----

# --- set up folder ---

run_id <- "final" # [REPRODUCIBILITY] fixed name (was: format(Sys.time(), ...))

results_dir_assoc <- file.path(
  results_dir,
  paste0("assoc_model1_PAR_global_", run_id)
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
  data = mod1_eff2,
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
  k = nrow(mod1_eff2),
  n_cohorts = n_distinct(mod1_eff2$cohort),
  n_timepoints = n_distinct(paste(mod1_eff2$cohort, mod1_eff2$timepoint)),
  n_model_combis = n_distinct(mod1_eff2$model_combi)
)

write.csv(
  effect_sizes_overall,
  file.path(results_dir_assoc, "effect_sizes_overall.csv"),
  row.names = FALSE
)

# --- LOO cohort ---

cohorts <- unique(mod1_eff2$cohort)

loo_cohort <- lapply(cohorts, function(cc) {
  
  tmp <- subset(mod1_eff2, cohort != cc)
  
  # [REPRODUCIBILITY] if the default optimizer (nlminb) does not converge for a
  # refit, retry with BFGS (otherwise one failed refit stops the whole script)
  fit <- tryCatch(
    rma.mv(
      yi = RLM_Estimate_scaled,
      V  = RLM_SE_coeftest_HC3_scaled^2,
      random = list(
        ~ 1 | cohort/timepoint,
        ~ 1 | model_combi
      ),
      data = tmp,
      method = "REML"
    ),
    error = function(e) rma.mv(
      yi = RLM_Estimate_scaled,
      V  = RLM_SE_coeftest_HC3_scaled^2,
      random = list(
        ~ 1 | cohort/timepoint,
        ~ 1 | model_combi
      ),
      data = tmp,
      method = "REML",
      control = list(optimizer = "optim", optmethod = "BFGS")
    )
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
  paste0("assoc_model1_PAR_modelwise_", run_id_modelwise)
)

dir.create(results_dir_assoc_modelwise, recursive = TRUE, showWarnings = FALSE)

log_file_modelwise <- file.path(results_dir_assoc_modelwise, "analysis_log_modelwise.txt")


# --- meta ---

mod2_res_mv_modelwise <- rma.mv(
  yi   = RLM_Estimate_scaled,
  V    = RLM_SE_coeftest_HC3_scaled^2,
  mods = ~ model_combi - 1,
  random = ~ 1 | cohort/timepoint,
  data   = mod1_eff2,
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

# cohorts <- unique(mod1_eff2$cohort)
# 
# loo_cohort_modelwise <- lapply(cohorts, function(cc) {
#   
#   message("Running LOO modelwise: removing ", cc)
#   
#   tmp <- subset(mod1_eff2, cohort != cc)
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

# cohorts <- unique(mod1_eff2$cohort)
# 
# loo_cohort_modelwise <- lapply(cohorts, function(cc) {
#   
#   tmp <- subset(mod1_eff2, cohort != cc)
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
cohorts <- unique(mod1_eff2$cohort)

loo_cohort_modelwise <- lapply(cohorts, function(cc) {
  
  message("Running LOO modelwise: removing ", cc)
  
  tmp <- mod1_eff2 %>%
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

results_dir_assoc_global_ageint_brain <- file.path(results_dir, paste0("assoc_model1_PAR_global_ageint_brain_", run_id))

dir.create(results_dir_assoc_global_ageint_brain, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_assoc_global_ageint_brain, "analysis_log_global_ageint.txt")

# --- data prep ---

# center age 
mod1_eff2 <- mod1_eff2 %>%
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
  data   = mod1_eff2,
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
    min(mod1_eff2$mean_age_brain_c, na.rm = TRUE),
    max(mod1_eff2$mean_age_brain_c, na.rm = TRUE),
    length.out = 100
  )
)

pred <- predict(mod2_res_mv_global_ageint_brain, newmods = newdat$mean_age_brain_c)

newdat$pred  <- pred$pred
newdat$ci.lb <- pred$ci.lb
newdat$ci.ub <- pred$ci.ub

newdat$mean_age_brain <- newdat$mean_age_brain_c + mean(mod1_eff2$mean_age_brain, na.rm = TRUE)

png(
  filename = file.path(results_dir_assoc_global_ageint_brain, "association_global_meanage_meta_regression.png"),
  width = 1800, height = 1200, res = 300
)

ggplot() +
  geom_point(
    data = mod1_eff2,
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

cohorts <- unique(mod1_eff2$cohort)

loo_assoc_ageint <- lapply(cohorts, function(cc) {
  
  tmp <- subset(mod1_eff2, cohort != cc)
  
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

results_dir_assoc_global_ageint_epi <- file.path(results_dir, paste0("assoc_model1_PAR_global_ageint_epi_", run_id))

dir.create(results_dir_assoc_global_ageint_epi, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_dir_assoc_global_ageint_epi, "analysis_log_global_ageint.txt")

# --- data prep ---

# center age 
mod1_eff2 <- mod1_eff2 %>%
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
  data   = mod1_eff2,
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
    min(mod1_eff2$mean_age_epi_c, na.rm = TRUE),
    max(mod1_eff2$mean_age_epi_c, na.rm = TRUE),
    length.out = 100
  )
)

pred <- predict(mod2_res_mv_global_ageint_epi, newmods = newdat$mean_age_epi_c)

newdat$pred  <- pred$pred
newdat$ci.lb <- pred$ci.lb
newdat$ci.ub <- pred$ci.ub

newdat$mean_age_epi <- newdat$mean_age_epi_c + mean(mod1_eff2$mean_age_epi, na.rm = TRUE)

png(
  filename = file.path(results_dir_assoc_global_ageint_epi, "association_global_meanage_meta_regression.png"),
  width = 1800, height = 1200, res = 300
)

ggplot() +
  geom_point(
    data = mod1_eff2,
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

cohorts <- unique(mod1_eff2$cohort)

loo_assoc_ageint <- lapply(cohorts, function(cc) {
  
  tmp <- subset(mod1_eff2, cohort != cc)
  
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
  paste0("assoc_model1_PAR_modelwise_ageint_brain_", run_id)
)

dir.create(results_dir_assoc_age_brain_modelwise, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  results_dir_assoc_age_brain_modelwise,
  "analysis_log_modelwise_ageint_brain.txt"
)


# --- data prep ---

mod1_eff2_age_brain <- mod1_eff2 %>%
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
  data = mod1_eff2_age_brain,
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

cohorts <- unique(mod1_eff2_age_brain$cohort)

loo_pairs_age_brain <- lapply(cohorts, function(cc) {
  
  message("Running LOO modelwise age interaction: removing ", cc)
  
  tmp <- mod1_eff2_age_brain %>%
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

# cohorts <- unique(mod1_eff2_age_brain$cohort)
# 
# loo_pairs_age_brain <- lapply(cohorts, function(cc) {
#   
#   tmp <- mod1_eff2_age_brain %>%
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
  paste0("assoc_model1_PAR_modelwise_ageint_epi_", run_id)
)

dir.create(results_dir_assoc_age_epi_modelwise, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  results_dir_assoc_age_epi_modelwise,
  "analysis_log_modelwise_ageint_epi.txt"
)


# --- data prep ---

mod1_eff2_age_epi <- mod1_eff2 %>%
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
  data = mod1_eff2_age_epi,
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

cohorts <- unique(mod1_eff2_age_epi$cohort)

loo_pairs_age_epi <- lapply(cohorts, function(cc) {
  
  message("Running LOO modelwise age interaction: removing ", cc)
  
  tmp <- mod1_eff2_age_epi %>%
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

# cohorts <- unique(mod1_eff2_age_epi$cohort)
# 
# loo_pairs_age_epi <- lapply(cohorts, function(cc) {
#   
#   tmp <- mod1_eff2_age_epi %>%
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





