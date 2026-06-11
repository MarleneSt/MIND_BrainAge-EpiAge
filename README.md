# MIND_BrainAge-EpiAge

This repository contains the code, documentation, workflows, and supporting materials used for the MIND Consortium BrainAge–EpiAge project.

Pre-registration: https://doi.org/10.17605/OSF.IO/6CSM3

## Getting Started

Get started by reading through **BrainAge-EpiAge_Overview.docx**. 

Overall, the repsository contains:

- `1.epi.preprocessing/` - Epigenetic preprocessing & QC workflows
- `2.epi.postprocessing/` – Epigenetic postprocessing workflows
- `3.epi.ages/` – Generate a variety of epigenetic ages 
- `4.brain.mri.preprocessing/` – MRI preprocessing & QC workflows
- `5.brain.ages/` – Generate a variety of brain ages 
- `6.analyses/` – Downstream analyses (cohort level)
- `brain_epi/` – Shared scripts, models, and workflow resources

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

The publication for this project is currently in preparation. In the meantime, if you use this repository, please acknowledge the MIND Consortium and cite: 
* The project pre-registration: Staginnus, M. et al. (2025). *Associations between epigenetic and brain age across development: Findings from the MIND consortium*. Open Science Framework. https://doi.org/10.17605/OSF.IO/6CSM3
* The MIND consortium paper: Schuurmans, I. K. et al. (2026). Consortium profile: the methylation, imaging and NeuroDevelopment (MIND) consortium. Molecular Psychiatry, 31(2), 1177-1189. https://doi.org/10.1038/s41380-025-03203-w

Additionally, this repository draws heavily on other open-access resources, including from the ENIGMA consortium and different brain age and epigenetic age models. Please ensure that you acknowledge the specific sources (papers, R packages, etc.) when used. We have provided links to papers and repositories where relevant. If you notice an omission, please don't hesitate to contact us. 
