## This script obtains epigenetic age predictions for AdaptAge, DamAge and DunedinPace
# conda activate /campaign/VB-FM5HPC-001/Vilte/condaenvs/conda-marlene
# export PATH=/campaign/VB-FM5HPC-001/Vilte/condaenvs/conda-marlene/bin:$PATH
# install all packages using 'pip install package_name' (before loading python), e.g. pip install sys

import sys
import pandas as pd
import pyaging as pya
import pyreadr
import os

# !! set working directory to 'brain_epi/epi' folder 
wd="/campaign/VB-FM5HPC-001/Vilte/Projects/Marlene-epi-clocks/brain_epi/epi"
os.chdir(wd)

# path to imputed DNAm file(s)
meth_dir = os.path.join(wd, "postprocessed/imputed")

# path to output directory
output_dir = os.path.join(wd, "models")

# path to 'functions' folder
functions_path = os.path.join(wd, "scripts", "functions")

# add 'functions' folder to path
sys.path.insert(0, functions_path)
from run_adaptage_damage_dunedin_prediction import run_adaptage_damage_dunedin_prediction

# !! define your dataset names and paths (to postprocessed DNAm data)
# !! for dataset name, add your array and time point (e.g., "450k_cord")
datasets = [
    ("450k_cord", os.path.join(meth_dir, "DNAm.imputed_wins_cord_450k.rds")),
    ("450k_F7", os.path.join(meth_dir, "DNAm.imputed_wins_F7_450k.rds")),
    ("450k_15up", os.path.join(meth_dir, "DNAm.imputed_wins_up15_450k.rds")),
    ("epic_15up", os.path.join(meth_dir, "DNAm.imputed_wins_up15_EPIC.rds"))
]

# obtain predictions for each dataset
for name, path in datasets:
    output_file = os.path.join(output_dir, f"epi_ages_{name}_pred_pyaging.csv")
    run_adaptage_damage_dunedin_prediction(path, output_file)
    print(f"Output for {name} saved in: {output_file}")

## !! see section below if you experience issues with EPIC array !!
## for me, epic .rds file was too large to read in using the function above
# hence, I converted the epic file from a matrix to a data.frame (in R), which made it easier for python to read it
# please do the same if you experience issues with the epic array
R
epic_df <- readRDS("/campaign/VB-FM5HPC-001/Vilte/Projects/Marlene-epi-clocks/DNAm.imputed_wins_up15_EPIC.rds")
epic_df <- as.data.frame(epic_df)
class(epic_df)
# re-saving the EPIC DNAm as a dataframe (DNAm.imputed_wins_up15_EPIC_df.rds)
saveRDS(epic_df, "DNAm.imputed_wins_up15_EPIC_df.rds", version=2)

# rerunning for EPIC 
output_file=os.path.join(output_dir, "epi_ages_epic_15up_pred_pyaging.csv")
run_adaptage_damage_dunedin_prediction('DNAm.imputed_wins_up15_EPIC_df.rds', output_file)
# ran successfully
