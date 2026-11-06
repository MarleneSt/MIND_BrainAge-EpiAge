# MIND_BrainAge-EpiAge

This repository contains the code, documentation, workflows, and supporting materials used for the MIND Consortium BrainAge–EpiAge project.

Pre-registration: https://doi.org/10.17605/OSF.IO/6CSM3

## Reproducing the paper

All results, tables and figures of the paper can be reproduced from the shareable cohort-level outputs in [`7.meta_analysis/`](7.meta_analysis/) with one command (`Rscript run_all.R`, about 3 hours). To recreate only the figures, download `figure_inputs.zip` from the [latest release](https://github.com/MarleneSt/MIND_BrainAge-EpiAge/releases/latest) (about 3 minutes). See the [`7.meta_analysis` README](7.meta_analysis/README.md) for instructions and a mapping of each figure and table to its script.

## Getting Started

Get started by reading through **BrainAge-EpiAge_Overview.docx**. 

Overall, the repository contains:

- `1.epi.preprocessing/` - Epigenetic preprocessing & QC workflows
- `2.epi.postprocessing/` – Epigenetic postprocessing workflows
- `3.epi.ages/` – Generate a variety of epigenetic ages 
- `4.brain.mri.preprocessing/` – MRI preprocessing & QC workflows
- `5.brain.ages/` – Generate a variety of brain ages 
- `6.analyses/` – Downstream analyses (cohort level)
- `7.meta_analysis/` – Consortium-level meta-analyses; reproduces all results, tables and figures of the paper (see `7.meta_analysis/README.md`)
- `brain_epi/` – Shared scripts, models, and workflow resources

## Interactive Dashboard

An interactive dashboard for exploring model performance and brain–epigenetic age associations across development is available at:
https://epi-brain-age-dashboard.streamlit.app/

The dashboard allows researchers to interactively filter results by cohort, age group, and model, and to download summary statistics as CSV files. It is optimised for desktop use. Note: the app may take ~30 seconds to load if it has been inactive (standard behaviour on Streamlit Community Cloud).

Source code: https://github.com/VilteBaltra/epi-brain-dashboard

## Required External Files

The following files are not stored in this repository because they exceed 100 MB.

### other_files/mind_dbn.sif

Download:
https://drive.google.com/file/d/1qEuOc40mCWq5uGO-_UuzquFpp6As0vjE/view

Place in:
`other_files/`

### other_files/ENIGMA-brainage-local.zip

Download:
https://drive.google.com/file/d/1NshqwnP9Sqj08oORJU1OEP6FY5V2D3PF/view

Place in:
`other_files/`

### other_files/pyment-public.zip

Download:
https://drive.google.com/file/d/1N5GODe8NngUoKg4mb10YibguPo1wmejo/view

Place in:
`other_files/`

### brain_epi/mri/scripts_models/pybrainage/ExtraTreesModel.onnx

Download:
https://drive.google.com/file/d/1v_xeJBY-gpN74dVb9UpeWAjL8HEgsJPJ/view

Place in:
`brain_epi/mri/scripts_models/pybrainage/`

## Citation

The preprint for this project is available on bioRxiv: Staginnus, M., et al. (2026). *Epigenetic and brain age across development: Performance and associations in the MIND consortium*. bioRxiv. https://doi.org/10.64898/2026.07.24.739162

If you use this repository, please acknowledge the MIND Consortium and cite the preprint and: 
* The project pre-registration: Staginnus, M. et al. (2025). *Associations between epigenetic and brain age across development: Findings from the MIND consortium*. Open Science Framework. https://doi.org/10.17605/OSF.IO/6CSM3
* The MIND consortium paper: Schuurmans, I. K. et al. (2026). Consortium profile: the methylation, imaging and NeuroDevelopment (MIND) consortium. Molecular Psychiatry, 31(2), 1177-1189. https://doi.org/10.1038/s41380-025-03203-w

Additionally, this repository draws heavily on other open-access resources, including from the ENIGMA consortium and different brain age and epigenetic age models. Please ensure that you acknowledge the specific sources (papers, R packages, etc.) when used. We have provided links to papers and repositories where relevant. If you notice an omission, please don't hesitate to contact us. 

## Contact

For questions or problems, please open an issue on this GitHub repository or contact Esther Walton (E.Walton@bath.ac.uk).
