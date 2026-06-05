## This script obtains epigenetic age predictions for AltumAge and Horvath
# conda activate /campaign/VB-FM5HPC-001/Vilte/Tools/condaenvs/conda-marlene
# export PATH=/campaign/VB-FM5HPC-001/Vilte/Tools/condaenvs/conda-marlene/bin:$PATH
# install all packages using 'pip install package_name' (before loading python), e.g. pip install sys
# if you have any problems with installation of pyaging check here: https://pyaging.readthedocs.io/en/latest/installation.html 

import sys
import pandas as pd
import pyaging as pya
import os

# !! set working directory to 'brain_epi/epi' folder 
wd="/campaign/VB-FM5HPC-001/Vilte/Projects/Marlene-epi-clocks/brain_epi/epi"
os.chdir(wd)

# path to bmiq directory
bmiq_dir = os.path.join(wd, "postprocessed/imputed/bmiq")

# path to output directory
output_dir = os.path.join(wd, "models")

# path to 'functions' folder
functions_path = os.path.join(wd, "scripts", "functions")

# add 'functions' folder to path
sys.path.insert(0, functions_path)
from run_bmiq_prediction import run_bmiq_prediction

# !! define your BMIQ normalised dataset name(s) [on the right side]
# !! for dataset name [on the left side], add your array and time point (e.g., "450k_cord")
datasets = [
    ("450k_cord", os.path.join(bmiq_dir, "DNAm.imputed_wins_cord_450k.bmiq.csv")),
    ("450k_F7", os.path.join(bmiq_dir, "DNAm.imputed_wins_F7_450k.bmiq.csv")),
    ("450k_15up", os.path.join(bmiq_dir,  "DNAm.imputed_wins_up15_450k.bmiq.csv")),
    ("epic_15up", os.path.join(bmiq_dir, "DNAm.imputed_wins_up15_EPIC.bmiq.csv"))
]

# obtain predictions for each dataset
for name, path in datasets:
    output_file = os.path.join(output_dir, f"epi_ages_{name}_pred_pyaging_bmiq.csv")
    run_bmiq_prediction(path, output_file)
    print(f"Output for {name} saved in: {output_file}")