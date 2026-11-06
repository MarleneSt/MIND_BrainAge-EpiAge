# MIND BrainAge–EpiAge: meta-analysis scripts

Scripts to reproduce the meta-analyses, tables and figures of *"Epigenetic and brain age across development: Performance and associations in the MIND consortium"* from the cohort-level results contributed by each cohort.

The cohort-level analyses (data preprocessing, epigenetic clock and brain age estimation, per-cohort performance and association models) were run locally by each cohort with the scripts in `brain_epi/` of this repository. The scripts in this folder pool those cohort-level results.

`data/` contains the cohort-level summary results contributed by each cohort (no individual-level data); the same summary statistics can also be explored in the interactive dashboard (https://epi-brain-age-dashboard.streamlit.app/).

## 1. System requirements

**Operating system.** Tested on macOS 26.4 (Apple M4 Pro, 24 GB RAM). The scripts use only standard R and Python packages and should run on any recent Windows, macOS or Linux system.

**Hardware.** No non-standard hardware required. About 8 GB of free memory is recommended.

**Software** (versions used for the reported results):

| Software | Version |
|---|---|
| R | 4.5.2 |
| metafor | 4.8-0 |
| tidyverse | 2.0.0 (dplyr 1.2.1, tidyr 1.3.2, ggplot2 4.0.2, purrr 1.2.1, tibble 3.3.1, readr 2.2.0, stringr 1.6.0, forcats 1.0.1) |
| pheatmap | 1.0.13 |
| Polychrome | 1.5.4 |
| patchwork | 1.3.2 |
| scales | 1.4.0 |
| ggrepel | 0.9.8 |
| corrplot | 0.95 |
| cowplot | 1.2.0 |
| Python (Figures 1 and 4B only) | 3.13.9, with plotly 6.3.0, kaleido 1.5.0, pandas 2.3.3, numpy 2.3.5, scipy 1.16.3, matplotlib 3.10.6, seaborn 0.13.2 |

## 2. Installation

1. Install R (≥ 4.5; https://cran.r-project.org) and Python 3 (for Figures 1 and 4B only).
2. Install the R packages:
   ```
   Rscript install_packages.R
   ```
3. Install the Python packages:
   ```
   pip install plotly kaleido pandas numpy scipy matplotlib seaborn
   ```

Typical install time on a standard desktop computer: about 10 minutes.

## 3. Demo / running the pipeline

Run all scripts in order from this folder:

```
Rscript run_all.R
```

`run_all.R` runs scripts `01`–`21` in order, each in a fresh R (or Python) session, and stops at the first script that fails. To resume from a given script after a failure, e.g. script 05: `Rscript run_all.R 05`. Scripts can also be run one by one in the same order (working directory = this folder).

**Expected output.** All results are written to `results/` (see section 5 for which file contains which result). At the end, `run_all.R` prints a short table comparing key results with the values reported in the paper. All rows should show `match = yes` (also saved in `results/logs/key_checks.csv`).

Logs for each script, run times (`run_times.csv`) and R session information (`sessionInfo.txt`) are saved in `results/logs/`.

**Run time.** About 3 hours on a standard desktop computer (Apple M4 Pro, 24 GB RAM). The longest scripts are `05` (≈ 85 min), `01` (≈ 50 min) and `07` (≈ 40 min); all other scripts take less than 5 minutes each. Most of the run time comes from the leave-one-cohort-out sensitivity analyses, which refit each model once per cohort (12–15 times); the main results do not depend on them.

## 4. Instructions for use

**Folder structure**

```
data/        cohort-level input files (one folder per cohort, as returned by the cohort-level scripts)
             + moderator tables (moderators_*_R.csv, moderators_*_models.csv)
             + cohort_age_sd.csv (SD of age per cohort x timepoint; used by 14)
             + descriptives_table_20260625.csv (cohort descriptives; used by 21)
results/     all outputs (created by the scripts)
archive/     superseded or one-off scripts (not part of the pipeline)
```

The moderator tables were first generated from the cohort descriptives (`archive/00_prep_descriptives_and_moderators.R`) and then completed by hand (tissue moderators; files ending in `_R.csv`). `cohort_age_sd.csv` was extracted once from SM Table S5 with `archive/make_cohort_age_sd.R`. `descriptives_table_20260625.csv` contains the cohort descriptives per cohort and timepoint (sample size, % female, mean age and mean predicted ages); the same descriptives are reported in SM Table S5.

**Scripts**

| Script | Content |
|---|---|
| `01_main_meta_wMAE.R` | Reads and cleans all cohort-level files; performance (MAEw) meta-analyses for epigenetic clocks and brain age models, including moderators (see below); main epigenetic PAR – brain PAR association meta-analyses and moderators; restricted subsample of associations |
| `02_meta_pure_performance.R` | MAEw in the restricted subsample (full age range and tissue overlap) |
| `03_meta_correlation.R` | Performance based on Pearson's r (correlation of predicted with chronological age), including moderators |
| `04_meta_unweightedMAE.R` | Performance based on unweighted MAE |
| `05_meta_within_modality.R` | Within-modality correlations (between epigenetic clocks; between brain age models) |
| `06_extracting_assoc_mod1_mod3.R` | Reads cohort-level association results without covariates (model 1) and with cell type adjustment (model 3) |
| `07_meta_assoc_model1_and_3.R` | Meta-analyses of models 1 and 3 (sensitivity analyses for covariate adjustment) |
| `08_R1_time_interval.R` | Time interval between DNAm and MRI assessment as moderator of the associations (revision) |
| `09_R2.1_age_effect_unweighted_MAE.R` | Age effect on unweighted MAE, adjusted for test-sample age range (revision) |
| `10_R2.1_per_clock_age_effects.R` | Per-clock age effects: MAEw vs unweighted MAE adjusted for age range (revision) |
| `11_R2.3_clock_fixed_effect.R` | Age range and tissue overlap moderators with epigenetic clock as fixed effect (revision) |
| `12_R2.5_meta_with_intercept.R` | Omnibus tests of differences between models, model pairs and clock generations, from models with intercept (revision) |
| `13_R3.2_sex_moderator.R` | Proportion of female participants as moderator (revision) |
| `14_R4.3_recalibrated_metrics.R` | Offset, slope and MAE after per-cohort linear recalibration (revision) |
| `15_R4.3_rank_correlation.R` | Rank correlation between MAEw and recalibrated MAE (revision) |
| `16_fdr.R` | FDR correction of all moderator and model-wise tests |
| `17_plots_main.R` | Main figures (performance and associations) |
| `18_plots_mean_age_effect_modelwise.R` | Model-specific age effects on performance |
| `19_plots_within_modality.R` | Within-modality correlation plots |
| `20_fig4b_heatmap.py` | Fig 4B: heatmap of model-pair associations (Python) |
| `21_plots_descriptives.py` | Fig 1: cohort descriptives (Python) |

**Moderators** (scripts 01 and 03; see Methods and SM 1.7 for details)

| Moderator | Outcome | Levels / definition |
|---|---|---|
| Age | performance, associations | sample mean age (continuous; for associations: mean age at MRI) |
| Age range overlap | performance | % of the test-sample age range covered by the model's training age range (continuous) |
| Tissue overlap | performance (epigenetic) | test tissue matches training tissue (yes / no; multi-tissue clocks count as matching) |
| Tissue of training sample | performance (epigenetic) | cortex, blood, multi-tissue, buccal, cord blood |
| Tissue of test sample | performance (epigenetic) | blood, saliva, buccal, cord blood, placenta |
| Processing level (minimal vs extracted) | performance (brain) | minimally processed vs extracted imaging-derived features |
| Clock generation | performance (epigenetic) | first-generation gestational, first, second, fourth generation (MAEw); additionally third generation (DunedinPACE) for Pearson's r |
| Clock generation (2 levels) | associations | first vs next generation |
| Clock generation (4 levels) | associations | first, second, third (DunedinPACE), fourth generation |
| Saliva / buccal involvement | associations | saliva (or buccal) tissue in the test sample or the clock's training data (yes / no) |
| Cortex training tissue | associations | clock trained (partly) on cortical tissue (yes / no) |
| Time interval (revision, script 08) | associations | mean difference between age at MRI and age at DNAm assessment (continuous) |
| Proportion of female participants (revision, script 13) | performance (MAEw), associations | % female per cohort and timepoint (continuous) |

## 5. Reproduction: results, tables and figures

<!-- TO UPDATE when SM table and figure numbers are final -->

| Manuscript result | Script(s) | Main output (in `results/`) |
|---|---|---|
| Epigenetic clock performance, MAEw (Results; SM 2.1–2.4; Table S9) | 01 | `perf_epi_wMAE_noDunedinDNAmTL_*_final/` |
| Brain age model performance, MAEw (Results; SM 2.5–2.8; Table S11) | 01 | `perf_brain_wMAE_noDunedin_*_final/` |
| Performance, Pearson's r (SM 2.1–2.8; Tables S10, S12) | 03 | `perf_*_pearson_*_final/` |
| Performance, unweighted MAE (SM 2.1, 2.5) | 04 | `perf_*_MAE_*_final/` |
| Restricted subsample, performance (SM 2.1, 2.5) | 02 | `perf_*_wMAE_PURE_*_final/` |
| Within-modality correlations (SM 2.9) | 05 | `within_modality/` |
| Epigenetic PAR – brain PAR associations and moderators (Results; SM 2.10, 2.11; Table S13) | 01 | `assoc_model2_PAR_*_final/` |
| Associations under alternative covariate adjustment (SM 2.12) | 06, 07 | `assoc_model1/`, `assoc_model3/` |
| Associations in restricted subsample (SM 2.13) | 01 | `assoc_model2_PAR_*_pure*final/` |
| Time interval moderator | 08 | `assoc_model2_PAR_global_agediff_final/` |
| Age effects on unweighted MAE | 09, 10 | `revision/R2.1/` |
| Clock as fixed effect | 11 | `revision/R2.3/` |
| Omnibus tests with intercept | 12 | `revision/R2.5/` |
| Proportion female moderator | 13 | `revision/R3.2/` |
| Recalibrated offset, slope and MAE (SM recalibration tables) | 14, 15 | `revision/R4.3/` (`R4.3_table_*_recalibration.csv`) |
| FDR-corrected p-values (all moderators, model-wise results) | 16 | `revision/fdr/` |
| All figures | 17–20 | see table below |

**Figures**

All figure files are in `results/plots/` unless stated otherwise.

| Figure | Content | Script | File |
|---|---|---|---|
| 1 | Cohort descriptives (mean age ± SD, sample size, % female) | 21 | `plot_cohort_descriptives.png` |
| 2A | Epigenetic clocks: model-wise and pooled MAEw | 17 | `plot_perf_epi_modelwise_plus_overall_meta_raw_combined_cohortcolors.png` |
| 2B | Epigenetic clocks: model-wise and pooled Pearson's r | 17 | `plot_perf_epi_modelwise_plus_overall_pearson_combined_cohortcolors.png` |
| 2C | Epigenetic clocks: MAEw by age bin | 17 | `plot_violin_perf_epi_age_brackets_morebrackets_k.png` |
| 2D | Epigenetic clocks: model-specific age effects on MAEw | 18 | `plot_epi_modelwise_age_slopes_plus_pooled_by_generation.png` |
| 3A | Brain age models: model-wise and pooled MAEw | 17 | `plot_perf_brain_modelwise_plus_overall_meta_raw_combined_cohortcolors.png` |
| 3B | Brain age models: model-wise and pooled Pearson's r | 17 | `plot_perf_brain_modelwise_plus_overall_pearson_combined_cohortcolors.png` |
| 3C | Brain age models: MAEw by age bin | 17 | `plot_violin_perf_brain_age_wMAE_brackets_morebrackets_k.png` |
| 3D | Brain age models: model-specific age effects on MAEw | 18 | `plot_brain_wMAE_modelwise_age_slopes_plus_pooled_by_generation.png` |
| 4A | Epigenetic PAR – brain PAR associations per model pair | 17 | `plot_assoc_pairwise_with_margins_fixed_5cols_wide.png` |
| 4B | Heatmap of model-pair associations | 20 | `plot_assoc_pairwise_heatmap.png` |
| 4C | Associations by brain age bin | 17 | `plot_violin_assoc_brain_mean_age_brackets_modelcombi_meta_new.png` |
| SM S1 | Epigenetic clocks: Pearson's r by age bin | 17 | `plot_violin_perf_epi_age_pearson_brackets.png` |
| SM S2 | Within-modality correlations of epigenetic PARs | 19 | `results/within_modality/plots/corrplot_modelcombination_r_clustered_circle_epi-PAR.png` |
| SM S3 | Within-modality correlations of brain age PARs | 19 | `results/within_modality/plots/corrplot_modelcombination_r_clustered_circle_brain-PAR.png` |

<!-- TO UPDATE: figure numbers once final -->

Further figure files in `results/plots/` (e.g. alternative layouts and versions without gestational clocks or based on unweighted MAE) are produced by the same scripts but are not shown in the paper.

*Claude (Anthropic) was used to help streamline the scripts for reproducibility and to map outputs to the reported results.*

## 6. License

This code is released under the [MIT License](LICENSE).

## 7. Paper

This code accompanies:

<!-- TO UPDATE with the published reference and DOI once available -->
Staginnus, M., et al. Epigenetic and brain age across development: Performance and associations in the MIND consortium. *bioRxiv* (2026). https://doi.org/10.64898/2026.07.24.739162

## 8. Contact

For questions or problems, please open an issue on this GitHub repository or contact Esther Walton (E.Walton@bath.ac.uk).
