###### 8.fdr: FDR correction for the revision (association analyses) ------

# FDR:
#   (A) moderator analyses: corrected across moderators per outcome ("p_assoc")
#   (B) model-wise results: corrected across the number of model pairs
#

#### folders ----
assoc_dir <- "results" # [REPRODUCIBILITY] relative path
perf_dir  <- "results" # [REPRODUCIBILITY] relative path
rev_dir   <- "results/revision" # [REPRODUCIBILITY] relative path
out_dir   <- file.path(rev_dir, "fdr")

library(metafor)   # needed to read the saved rma.mv objects


#### (A) moderator family: associations (model 2) ----

# --- moderators read from the saved models (deduplicated data, 22/10/2026 and 03/11/2026) ---
m_age_brain <- readRDS(file.path(assoc_dir, "assoc_model2_PAR_global_ageint_brain_final/model_assoc_mv_global_ageint_brain.rds"))
m_gen_2     <- readRDS(file.path(assoc_dir, "assoc_model2_PAR_global_generation_epi_binary_final/model_assoc_mv_generation_epi_binary.rds"))
m_saliva    <- readRDS(file.path(assoc_dir, "assoc_model2_PAR_global_tissue_saliva_testortrain_final/model_assoc_mv_tissue_saliva_testortrain.rds"))
m_buccal    <- readRDS(file.path(assoc_dir, "assoc_model2_PAR_global_tissue_buccal_testortrain_final/model_assoc_mv_tissue_buccal_testortrain.rds"))
m_cortex    <- readRDS(file.path(assoc_dir, "assoc_model2_PAR_global_tissue_cortical_final/model_assoc_mv_tissue_cortical.rds"))
m_agediff   <- readRDS(file.path(assoc_dir, "assoc_model2_PAR_global_agediff_final/model_assoc_mv_global_agediff.rds"))

# --- model pair and generation (4 level): omnibus test WITH intercept (R2.5) ---
omnibus <- read.csv(file.path(rev_dir, "R2.5/assoc_model2_PAR_modelwise_final/omnibus_differences_main_text_R2.5.csv"))
omni_pairs <- omnibus[omnibus$test == "assoc_pairs", ]
omni_gen4  <- omnibus[omnibus$test == "assoc_generation", ]

# --- moderator table (one row per moderator test) ---
moderators <- data.frame(
  moderator = c("age (brain)",
                "model pair",
                "generation (2 level)",
                "generation (4 level)",
                "saliva test or training",
                "buccal test or training",
                "cortex training",
                "time interval DNAm-MRI"),
  Q  = c(m_age_brain$QM, omni_pairs$QM, m_gen_2$QM, omni_gen4$QM,
         m_saliva$QM, m_buccal$QM, m_cortex$QM, m_agediff$QM),
  df = c(m_age_brain$QMdf[1], omni_pairs$df, m_gen_2$QMdf[1], omni_gen4$df,
         m_saliva$QMdf[1], m_buccal$QMdf[1], m_cortex$QMdf[1], m_agediff$QMdf[1]),
  p_assoc = c(m_age_brain$QMp, omni_pairs$p, m_gen_2$QMp, omni_gen4$p,
              m_saliva$QMp, m_buccal$QMp, m_cortex$QMp, m_agediff$QMp)
)

# --- new moderator: proportion of female participants (R3.2) ---
female <- read.csv(file.path(rev_dir, "R3.2/R3.2_perc_female_moderator.csv"))
female <- female[female$outcome == "PAR association" & female$sample == "all", ]
# single moderator with 1 df: QM = (b / se)^2, se from the 95% CI
se_female <- (female$ci.ub - female$ci.lb) / (2 * qnorm(0.975))
moderators <- rbind(moderators,
                    data.frame(moderator = "proportion female", Q = (female$b / se_female)^2, df = 1, p_assoc = female$p))

# --- FDR correction across moderators ---
moderators$p_assoc_fdr <- p.adjust(moderators$p_assoc, method = "fdr")
moderators


#### (B) model-wise results: associations per model pair (136 pairs) ----

# --- main effect per pair (model without intercept, deduplicated data) ---
modelwise <- read.csv(file.path(assoc_dir, "assoc_model2_PAR_modelwise_final/modelwise_estimates.csv"))

# --- age interaction per pair (brain mean age) ---
m_pairs_age <- readRDS(file.path(assoc_dir, "assoc_model2_PAR_modelwise_ageint_brain_final/model_assoc_mv_modelwise_ageint_brain.rds"))
i_age <- grep(":mean_age_brain_c", rownames(m_pairs_age$b))
age_terms <- data.frame(
  model_combi       = sub(":mean_age_brain_c", "", sub("^model_combi", "", rownames(m_pairs_age$b)[i_age]), fixed = TRUE),
  b_age             = m_pairs_age$b[i_age, 1],
  p_assoc_model_age = m_pairs_age$pval[i_age]
)

# --- combine by pair (check: all 136 pairs matched) ---
pairs <- merge(modelwise[, c("model_combi", "brain_model", "epi_model", "beta", "pval")],
               age_terms, by = "model_combi")
names(pairs)[names(pairs) == "pval"] <- "p_assoc_model"
stopifnot(nrow(pairs) == 136)

# --- FDR correction across the 136 pairs (separately for main effect and age interaction) ---
pairs$p_assoc_model_fdr     <- p.adjust(pairs$p_assoc_model, method = "fdr")
pairs$p_assoc_model_age_fdr <- p.adjust(pairs$p_assoc_model_age, method = "fdr")

pairs <- pairs[order(pairs$brain_model, pairs$epi_model),
               c("brain_model", "epi_model", "beta", "p_assoc_model", "p_assoc_model_fdr",
                 "b_age", "p_assoc_model_age", "p_assoc_model_age_fdr")]
head(pairs)


#### (C) moderator families: performance (MAEw) ----

# --- proportion of female participants (R3.2), MAEw rows ---
female_all   <- read.csv(file.path(rev_dir, "R3.2/R3.2_perc_female_moderator.csv"))
female_epi   <- female_all[female_all$outcome == "epi MAEw (log)" & female_all$sample == "all", ]
female_brain <- female_all[female_all$outcome == "brain MAEw"     & female_all$sample == "all", ]
se_female_epi   <- (female_epi$ci.ub   - female_epi$ci.lb)   / (2 * qnorm(0.975))
se_female_brain <- (female_brain$ci.ub - female_brain$ci.lb) / (2 * qnorm(0.975))

# --- model: omnibus test WITH intercept (R2.5) ---
omni_epi   <- omnibus[omnibus$test == "epi_wMAE_models", ]
omni_brain <- omnibus[omnibus$test == "brain_wMAE_models", ]

###### (C1) epi MAEw ------
# moderators read from the saved models (results_020726, as in the paper)
m_epi_age     <- readRDS(file.path(perf_dir, "perf_epi_wMAE_noDunedinDNAmTL_global_ageint_final/model_epi_wMAE_mv_global_ageint_logscale.rds"))
m_epi_overlap <- readRDS(file.path(perf_dir, "perf_epi_wMAE_noDunedinDNAmTL_global_ageoverlap_perc_final/model_epi_wMAE_mv_ageoverlap_perc_logscale.rds"))
m_epi_tis_ov  <- readRDS(file.path(perf_dir, "perf_epi_wMAE_noDunedinDNAmTL_global_tissue_overlap_final/model_epi_wMAE_mv_tissue_overlap_logscale.rds"))
m_epi_tis_te  <- readRDS(file.path(perf_dir, "perf_epi_wMAE_noDunedinDNAmTL_global_tissue_test_final/model_epi_wMAE_mv_tissue_test_logscale.rds"))
m_epi_tis_tr  <- readRDS(file.path(perf_dir, "perf_epi_wMAE_noDunedinDNAmTL_global_tissue_training_final/model_epi_wMAE_mv_tissue_training_logscale.rds"))
m_epi_gen     <- readRDS(file.path(perf_dir, "perf_epi_wMAE_noDunedinDNAmTL_global_generation_epi_final/model_epi_wMAE_mv_generation_epi_logscale.rds"))

moderators_epi <- data.frame(
  moderator = c("age",
                "model",
                "age range overlap",
                "tissue overlap",
                "tissue test",
                "tissue training",
                "generation (4 level)",
                "proportion female"),
  Q  = c(m_epi_age$QM, omni_epi$QM, m_epi_overlap$QM, m_epi_tis_ov$QM,
         m_epi_tis_te$QM, m_epi_tis_tr$QM, m_epi_gen$QM, (female_epi$b / se_female_epi)^2),
  df = c(m_epi_age$QMdf[1], omni_epi$df, m_epi_overlap$QMdf[1], m_epi_tis_ov$QMdf[1],
         m_epi_tis_te$QMdf[1], m_epi_tis_tr$QMdf[1], m_epi_gen$QMdf[1], 1),
  p_epi_maew = c(m_epi_age$QMp, omni_epi$p, m_epi_overlap$QMp, m_epi_tis_ov$QMp,
                 m_epi_tis_te$QMp, m_epi_tis_tr$QMp, m_epi_gen$QMp, female_epi$p)
)
moderators_epi$p_epi_maew_fdr <- p.adjust(moderators_epi$p_epi_maew, method = "fdr")
moderators_epi

###### (C2) brain MAEw ------
m_brain_age     <- readRDS(file.path(perf_dir, "perf_brain_wMAE_noDunedin_global_ageint_final/model_brain_wMAE_mv_global_ageint.rds"))
m_brain_overlap <- readRDS(file.path(perf_dir, "perf_brain_wMAE_noDunedin_global_ageoverlap_perc_final/model_brain_wMAE_mv_ageoverlap_perc.rds"))
m_brain_proc    <- readRDS(file.path(perf_dir, "perf_brain_wMAE_noDunedin_global_processing_final/model_brain_wMAE_mv_processing.rds"))

moderators_brain <- data.frame(
  moderator = c("age",
                "model",
                "age range overlap",
                "processing level",
                "proportion female"),
  Q  = c(m_brain_age$QM, omni_brain$QM, m_brain_overlap$QM, m_brain_proc$QM,
         (female_brain$b / se_female_brain)^2),
  df = c(m_brain_age$QMdf[1], omni_brain$df, m_brain_overlap$QMdf[1], m_brain_proc$QMdf[1], 1),
  p_brain_maew = c(m_brain_age$QMp, omni_brain$p, m_brain_overlap$QMp, m_brain_proc$QMp, female_brain$p)
)
moderators_brain$p_brain_maew_fdr <- p.adjust(moderators_brain$p_brain_maew, method = "fdr")
moderators_brain


#### (D) moderator families: performance (Pearson r) ----
# Unchanged from the paper (proportion female was only tested on MAEw); added here
# so that all moderator pFDRs can be traced to this script.
# All moderators, including model, are read from the saved models (results_020726,
# as in the paper; model omnibus test without intercept, as in the paper).

###### (D1) epi Pearson r ------
m_epi_r_age     <- readRDS(file.path(perf_dir, "perf_epi_pearson_global_ageint_final/model_epi_pearson_mv_global_ageint.rds"))
m_epi_r_model   <- readRDS(file.path(perf_dir, "perf_epi_pearson_modelwise_final/model_epi_pearson_mv_modelwise.rds"))
m_epi_r_overlap <- readRDS(file.path(perf_dir, "perf_epi_pearson_noGest_global_ageoverlap_perc_final/model_epi_pearson_mv_ageoverlap_perc.rds"))
m_epi_r_tis_ov  <- readRDS(file.path(perf_dir, "perf_epi_pearson_global_tissue_overlap_final/model_epi_pearson_mv_tissue_overlap.rds"))
m_epi_r_tis_te  <- readRDS(file.path(perf_dir, "perf_epi_pearson_global_tissue_test_final/model_epi_pearson_mv_tissue_test.rds"))
m_epi_r_tis_tr  <- readRDS(file.path(perf_dir, "perf_epi_pearson_global_tissue_training_final/model_epi_pearson_mv_tissue_training.rds"))
m_epi_r_gen     <- readRDS(file.path(perf_dir, "perf_epi_pearson_global_generation_epi_final/model_epi_pearson_mv_generation_epi.rds"))

moderators_epi_r <- data.frame(
  moderator = c("age",
                "model",
                "age range overlap",
                "tissue overlap",
                "tissue test",
                "tissue training",
                "generation (4 level)"),
  Q  = c(m_epi_r_age$QM, m_epi_r_model$QM, m_epi_r_overlap$QM, m_epi_r_tis_ov$QM,
         m_epi_r_tis_te$QM, m_epi_r_tis_tr$QM, m_epi_r_gen$QM),
  df = c(m_epi_r_age$QMdf[1], m_epi_r_model$QMdf[1], m_epi_r_overlap$QMdf[1], m_epi_r_tis_ov$QMdf[1],
         m_epi_r_tis_te$QMdf[1], m_epi_r_tis_tr$QMdf[1], m_epi_r_gen$QMdf[1]),
  p_epi_r = c(m_epi_r_age$QMp, m_epi_r_model$QMp, m_epi_r_overlap$QMp, m_epi_r_tis_ov$QMp,
              m_epi_r_tis_te$QMp, m_epi_r_tis_tr$QMp, m_epi_r_gen$QMp)
)
moderators_epi_r$p_epi_r_fdr <- p.adjust(moderators_epi_r$p_epi_r, method = "fdr")
moderators_epi_r

###### (D2) brain Pearson r ------
m_brain_r_age     <- readRDS(file.path(perf_dir, "perf_brain_pearson_global_ageint_final/model_brain_pearson_mv_global_ageint.rds"))
m_brain_r_model   <- readRDS(file.path(perf_dir, "perf_brain_pearson_modelwise_final/model_brain_pearson_mv_modelwise.rds"))
m_brain_r_overlap <- readRDS(file.path(perf_dir, "perf_brain_pearson_global_ageoverlap_perc_final/model_brain_pearson_mv_ageoverlap_perc.rds"))
m_brain_r_proc    <- readRDS(file.path(perf_dir, "perf_brain_pearson_global_processing_final/model_brain_pearson_mv_processing.rds"))

moderators_brain_r <- data.frame(
  moderator = c("age",
                "model",
                "age range overlap",
                "processing level"),
  Q  = c(m_brain_r_age$QM, m_brain_r_model$QM, m_brain_r_overlap$QM, m_brain_r_proc$QM),
  df = c(m_brain_r_age$QMdf[1], m_brain_r_model$QMdf[1], m_brain_r_overlap$QMdf[1], m_brain_r_proc$QMdf[1]),
  p_brain_r = c(m_brain_r_age$QMp, m_brain_r_model$QMp, m_brain_r_overlap$QMp, m_brain_r_proc$QMp)
)
moderators_brain_r$p_brain_r_fdr <- p.adjust(moderators_brain_r$p_brain_r, method = "fdr")
moderators_brain_r


#### save ----
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(moderators, file.path(out_dir, "moderator_fdr_assoc_revision.csv"), row.names = FALSE)
write.csv(pairs,      file.path(out_dir, "modelwise_res_FDR_assoc_revision.csv"), row.names = FALSE)
write.csv(moderators_epi,   file.path(out_dir, "moderator_fdr_epi_maew_revision.csv"), row.names = FALSE)
write.csv(moderators_brain, file.path(out_dir, "moderator_fdr_brain_maew_revision.csv"), row.names = FALSE)
write.csv(moderators_epi_r,   file.path(out_dir, "moderator_fdr_epi_r_revision.csv"), row.names = FALSE)
write.csv(moderators_brain_r, file.path(out_dir, "moderator_fdr_brain_r_revision.csv"), row.names = FALSE)
