#!/bin/bash

# This script has been adapted from the extract.sh script and the instruction on how to extract subcortical ROIs from the ENIGMA consortium https://enigma.ini.usc.edu/protocols/imaging-protocols/

echo 'SUBJID,Left.Lateral.Ventricle,Left.Inf.Lat.Vent,Left.Cerebellum.White.Matter,Left.Cerebellum.Cortex,Left.Thalamus.Proper,Left.Caudate,Left.Putamen,Left.Pallidum,X3rd.Ventricle,X4th.Ventricle,Brain.Stem,Left.Hippocampus,Left.Amygdala,Left.Accumbens.area,Right.Lateral.Ventricle,Right.Inf.Lat.Vent,Right.Cerebellum.White.Matter,Right.Cerebellum.Cortex,Right.Thalamus.Proper,Right.Caudate,Right.Putamen,Right.Pallidum,Right.Hippocampus,Right.Amygdala,Right.Accumbens.area,CC_Posterior,CC_Mid_Posterior,CC_Central,CC_Mid_Anterior,CC_Anterior,lhCortexVol,rhCortexVol,lhCorticalWhiteMatterVol,rhCorticalWhiteMatterVol,SubCortGrayVol,TotalGrayVol,SupraTentorialVol,EstimatedTotalIntraCranialVol' > glasser_cerebellum_subcortical.csv

for subj_id in `ls -d sub*`; do 

printf "%s,"  "${subj_id}" >> glasser_cerebellum_subcortical.csv

for x in Left-Lateral-Ventricle Left-Inf-Lat-Vent Left-Cerebellum-White-Matter Left-Cerebellum-Cortex Left-Thalamus-Proper Left-Caudate Left-Putamen Left-Pallidum 3rd-Ventricle 4th-Ventricle Brain-Stem Left-Hippocampus Left-Amygdala Left-Accumbens-area Right-Lateral-Ventricle Right-Inf-Lat-Vent Right-Cerebellum-White-Matter Right-Cerebellum-Cortex Right-Thalamus-Proper Right-Caudate Right-Putamen Right-Pallidum Right-Hippocampus Right-Amygdala Right-Accumbens-area CC_Posterior CC_Mid_Posterior CC_Central CC_Mid_Anterior CC_Anterior; do

printf "%s," `grep  ${x} ${subj_id}/stats/aseg.stats | awk '{print $4}'` >> glasser_cerebellum_subcortical.csv

done

printf "%s," `cat ${subj_id}/stats/aseg.stats | grep lhCortexVol | awk -F, '{print $4}'` >> glasser_cerebellum_subcortical.csv
printf "%s," `cat ${subj_id}/stats/aseg.stats | grep rhCortexVol | awk -F, '{print $4}'` >> glasser_cerebellum_subcortical.csv
printf "%s," `cat ${subj_id}/stats/aseg.stats | grep lhCerebralWhiteMatterVol | awk -F, '{print $4}'` >> glasser_cerebellum_subcortical.csv
printf "%s," `cat ${subj_id}/stats/aseg.stats | grep rhCerebralWhiteMatterVol | awk -F, '{print $4}'` >> glasser_cerebellum_subcortical.csv
printf "%s," `cat ${subj_id}/stats/aseg.stats | grep SubCortGrayVol | awk -F, '{print $4}'` >> glasser_cerebellum_subcortical.csv
printf "%s," `cat ${subj_id}/stats/aseg.stats | grep TotalGrayVol | awk -F, '{print $4}'` >> glasser_cerebellum_subcortical.csv
# as there are multiple instances starting with SupraTentorialVol, the grep command here had to be altered 
printf "%s," `cat ${subj_id}/stats/aseg.stats | grep -w SupraTentorialVol, | awk -F, '{print $4}'` >> glasser_cerebellum_subcortical.csv
printf "%s" `cat ${subj_id}/stats/aseg.stats | grep IntraCranialVol | awk -F, '{print $4}'` >> glasser_cerebellum_subcortical.csv

echo "" >> glasser_cerebellum_subcortical.csv

done
