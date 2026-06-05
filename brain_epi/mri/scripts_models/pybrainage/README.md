## PyBrainAge

This script provides an end-to-end workflow for the experimental PyBrainAge. The original PyBrainAge has two limitations:
 - The provided models require Python versions that are now end-of-life.
 - The documentation for generating the required input files does not match the included scripts.

This version updates the process and standardizes the input generation. It retains the original permissive open source license.

### Requirements

FreeSurfer installed with a valid license, and FREESURFER_HOME set in your environment.

Python dependencies (install with):

```bash
pip install -r requirements.txt
```

### Input data
Prepare one T1-weighted NIfTI scan per subject in a folder containing no other NIfTI files. Including the participant ID in the filename is recommended.
Example: the folder t1/ contains

```
cmemprage.nii.gz  
tmemprage.nii.gz  
```

### Processing steps

Run the three stages from the root folder (here t1):

```bash
python 1freesurfer.py ./t1
python 2fs2brainage.py ./t1
# Edit brainage.csv to replace the placeholder age (42) with the actual age
python 3predict_csv.py ./t1/brainage.csv
```

 - Stage 1: Runs FreeSurfer on each T1 image, storing results in the root folder.
 - Stage 2: Extracts the required brain region volumes and writes brainage.csv to the root folder. The second column contains a placeholder age (42) that you should replace with the actual value (see Brain-PAD details in next section).
 - Stage 3: Predicts brain age, saving results to PyBrainAge_Output.csv and a plot to PyBrainAge_Output.png.

## Brain-PAD
Typical Brain-age analyses involve calculating Brain-PAD (Brain Predicted Age Difference), also referred to as BrainGAP, BrainAge Delta, or similar variations. Brain-PAD is determined by subtracting chronological age from Brain-age ("Age" minus "Brain-age" columns, which is already calculated for you in the PyBrainAge_Output.csv file). This metric can be utilised to examine associations with health outcomes. For further detailed discussion, refer to the work of Cole and Franke (2017). For example, larger Brain-PADs (older-appearing brains) have been associated to an increased risk in a future diagnosis of dementia in memory clinic patients (Biondo et al, 2022).

Statistical analyses involving Brain-PAD (or Brain-age) values require minimising the regression-to-the-mean effect by including linear and non-linear terms of age as covariates (see de Lange and Cole, 2020). 

## Acknowledgements and References
This software was co-created: by Andre Marquand<sup>1</sup>, Saige Rutherford<sup>1</sup>, Ayodeji Ijishakin<sup>2</sup>, Francesca Biondo<sup>2</sup> & James Cole<sup>2</sup>.
<sup>1</sup> The Donders Institute, Netherlands; <sup>2</sup> University College London (UCL), United Kingdom. The port to modern Python was completed by Chris Rorden at the Universit of South Carolina.

### Caveats:

This model has not been peer reviewed. The FreeSurfer processing has changed a bit since the original model. The update to modern Python has not been extensively tested or validated against the orignal predictions.

### Citations: 
We have not used this software in published work yet. However, if you need a reference for this software, feel free to cite this page or the Rutherford et al paper in the list of references below.

### References
Rutherford, Fraza, Dinga, Kia, Wolfers, Zabihi... Beckmann and Marquand (2022). Charting brain growth and aging at high spatial precision. eLife. [DOI](https://doi.org/10.7554/eLife.72904)

Cole and Franke (2017). Predicting age using neuroimaging: innovative brain ageing biomarkers. Trends in Neurosciences. [DOI](https://doi.org/10.1016/j.tins.2017.10.001)

Biondo, Jewell, Pritchard, Aarsland, Steves, Mueller and Cole (2022). Brain-age is associated with progression to dementia in memory clinic patients. Neuroimage: Clinical. [DOI](https://doi.org/10.1016/j.nicl.2022.103175)

de Lange & Cole (2020). Commentary: Correction procedures in brain-age prediction. Neuroimage: Clinical [DOI](https://doi.org/10.1016/j.nicl.2020.102229)