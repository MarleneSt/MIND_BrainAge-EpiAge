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

**Figures only.** To recreate the figures without running the analyses, download `figure_inputs.zip` from the repository's Releases page, unzip it in this folder, and run `Rscript run_all.R 17` (runs scripts 17–21; about 3 minutes). The zip contains the analysis outputs the figure scripts read.

**Run time.** About 3 hours on a standard desktop computer (Apple M4 Pro, 24 GB RAM). The longest scripts are `05` (≈ 85 min), `01` (≈ 50 min) and `07` (≈ 40 min); all other scripts take less than 5 minutes each. Most of the run time comes from the leave-one-cohort-out sensitivity analyses, which refit each model once per cohort (12–15 times); the main results do not depend on them.

## 4. Instructions for use

**Folder structure**

```
data/        cohort-level input files (one folder per cohort, as returned by the cohort-level scripts)
             + moderator tables (moderators_*_R.csv, moderators_*_models.csv)
             + cohort_age_sd.csv (SD of age per cohort x timepoint; used by script 14)
             + descriptives_table_20260625.csv (cohort descriptives; used by script 21)
results/     all outputs (created by the scripts)
archive/     superseded or one-off scripts (not part of the pipeline)
```

**Scripts**

| Script | Content |
|---|---|
| `01_main_meta_wMAE.R` | Reads and cleans all cohort-level files; performance (MAEw) meta-analyses for epigenetic clocks and brain age models, including moderators (age, age range overlap, tissue, processing level (minimal vs extracted), clock generation; see Methods and SM Table S1.7.1); main epigenetic PAR – brain PAR association meta-analyses and moderators; restricted subsample of associations |
| `02_meta_pure_performance.R` | MAEw in the restricted subsample (full age range and tissue overlap) |
| `03_meta_correlation.R` | Performance based on Pearson's r (correlation of predicted with chronological age), including moderators |
| `04_meta_unweightedMAE.R` | Performance based on unweighted MAE |
| `05_meta_within_modality.R` | Within-modality correlations (between epigenetic clocks; between brain age models) |
| `06_extracting_assoc_mod1_mod3.R` | Reads cohort-level association results without covariates (model 1) and with cell type adjustment (model 3) |
| `07_meta_assoc_model1_and_3.R` | Meta-analyses of models 1 and 3 (sensitivity analyses for covariate adjustment) |
| `08_R1_time_interval.R` | Time interval between DNAm and MRI assessment as moderator of the associations (revision) |
| `09_R2.1_age_effect_unweighted_MAE.R` | Pooled age effect on unweighted MAE, adjusted for test-sample age range (sensitivity analysis) (revision) |
| `10_R2.1_per_clock_age_effects.R` | As script 09 but per model; shown alongside the original MAEw estimates (revision) |
| `11_R2.3_clock_fixed_effect.R` | Age range and tissue overlap moderators with epigenetic clock as fixed effect (revision) |
| `12_R2.5_meta_with_intercept.R` | Omnibus tests with intercept for differences between models (MAEw), model pairs and clock generations (associations) (revision) |
| `13_R3.2_sex_moderator.R` | Proportion of female participants as moderator (revision) |
| `14_R4.3_recalibrated_metrics.R` | Offset, slope and MAE after per-cohort linear recalibration (revision) |
| `15_R4.3_rank_correlation.R` | Rank correlation between MAEw and recalibrated MAE (revision) |
| `16_fdr.R` | FDR correction of all moderator and model-wise tests |
| `17_plots_main.R` | Figures of performance and associations (Figs 2A–C, 3A–C, 4A, 4C; SM Fig S1) |
| `18_plots_mean_age_effect_modelwise.R` | Figures of model-specific age effects on performance (Figs 2D, 3D) |
| `19_plots_within_modality.R` | Figures of within-modality correlations (SM Figs S2, S3) |
| `20_fig4b_heatmap.py` | Figure of model pair associations as heatmap (Fig 4B) |
| `21_plots_descriptives.py` | Figure of cohort descriptives (Fig 1) |

## 5. Reproduction: results, tables and figures

<!-- TO UPDATE when SM table and figure numbers are final -->

### Results and tables

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

### Figures

All figure files are saved in `results/plots/`, except SM Figures S2 and S3 (`results/within_modality/plots/`). To recreate the figures without running the analyses, see "Figures only" in section 3.

| Figure | Content | Script | File |
|---|---|---|---|
| 1 | Cohort descriptives | 21 | `plot_cohort_descriptives.png` |
| 2A | Epigenetic clocks: model-wise and pooled MAEw | 17 | `plot_perf_epi_modelwise_plus_overall_meta_raw_combined_cohortcolors.png` |
| 2B | Epigenetic clocks: model-wise and pooled Pearson's r | 17 | `plot_perf_epi_modelwise_plus_overall_pearson_combined_cohortcolors.png` |
| 2C | Epigenetic clocks: MAEw by age bin | 17 | `plot_violin_perf_epi_age_brackets_morebrackets_k.png` |
| 2D | Epigenetic clocks: model-specific age effects on MAEw | 18 | `plot_epi_modelwise_age_slopes_plus_pooled_by_generation.png` |
| 3A | Brain age models: model-wise and pooled MAEw | 17 | `plot_perf_brain_modelwise_plus_overall_meta_raw_combined_cohortcolors.png` |
| 3B | Brain age models: model-wise and pooled Pearson's r | 17 | `plot_perf_brain_modelwise_plus_overall_pearson_combined_cohortcolors.png` |
| 3C | Brain age models: MAEw by age bin | 17 | `plot_violin_perf_brain_age_wMAE_brackets_morebrackets_k.png` |
| 3D | Brain age models: model-specific age effects on MAEw | 18 | `plot_brain_wMAE_modelwise_age_slopes_plus_pooled_by_generation.png` |
| 4A | Epigenetic PAR – brain PAR associations per model pair | 17 | `plot_assoc_pairwise_with_margins_fixed_5cols_wide.png` |
| 4B | Heatmap of model pair associations | 20 | `plot_assoc_pairwise_heatmap.png` |
| 4C | Associations by brain age bin | 17 | `plot_violin_assoc_brain_mean_age_brackets_modelcombi_meta_new.png` |
| SM S1 | Epigenetic clocks: Pearson's r by age bin | 17 | `plot_violin_perf_epi_age_pearson_brackets.png` |
| SM S2 | Within-modality correlations of epigenetic PARs | 19 | `within_modality/plots/corrplot_modelcombination_r_clustered_circle_epi-PAR.png` |
| SM S3 | Within-modality correlations of brain age PARs | 19 | `within_modality/plots/corrplot_modelcombination_r_clustered_circle_brain-PAR.png` |

<!-- TO UPDATE: figure numbers once final -->

Further figure files in `results/plots/` (e.g. alternative layouts and versions without gestational clocks or based on unweighted MAE) are produced by the same scripts but are not shown in the paper.

*Claude (Anthropic) was used to help streamline the scripts for reproducibility and to map outputs to the reported results.*

## 6. License

This code is released under the [MIT License](LICENSE).
