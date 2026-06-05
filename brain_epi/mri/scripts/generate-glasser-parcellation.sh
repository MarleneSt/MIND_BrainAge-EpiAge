#!/bin/sh

# this script is based on instructions/code from this github page https://github.com/ftadel/IntrAnat/issues/7 written by user ftadel

# change the output_data path to point to your SUBJECTS_DIR/folder with FreeSurfer pre-processed participants
Output_data=/brain_epi/mri/output

# Change to the directory where the participant folders are located
cd "$Output_data" || exit

# Loop over the directories starting with 'sub', this requires changing in line with the prefix of your subject folders 
for i in sub*/;
do
    # Remove the trailing slash to just get the directory name
    i=${i%/}

    echo "${i}"

    # Project annot files from fsaverage to subject (left and right hemispheres separately)
    mri_surf2surf --srcsubject fsaverage --trgsubject "${i}" --hemi lh --sval-annot "${Output_data}/lh.HCP-MMP1.annot" --tval "${Output_data}/${i}/label/lh.HCP-MMP1.annot"
    mri_surf2surf --srcsubject fsaverage --trgsubject "${i}" --hemi rh --sval-annot "${Output_data}/rh.HCP-MMP1.annot" --tval "${Output_data}/${i}/label/rh.HCP-MMP1.annot"

    # Convert from .annot to aseg.mgz volume (left hemi = original labels ; right hemi = 200 + original labels)
    mri_surf2volseg --o "${Output_data}/${i}/mri/HCP-MMP1+aseg.mgz" --label-cortex --i "${Output_data}/${i}/mri/aseg.mgz" --threads 1 --lh-annot "${Output_data}/${i}/label/lh.HCP-MMP1.annot" 0 --lh-cortex-mask "${Output_data}/${i}/label/lh.cortex.label" --lh-white "${Output_data}/${i}/surf/lh.white" --lh-pial "${Output_data}/${i}/surf/lh.pial" --rh-annot "${Output_data}/${i}/label/rh.HCP-MMP1.annot" 200 --rh-cortex-mask "${Output_data}/${i}/label/rh.cortex.label" --rh-white "${Output_data}/${i}/surf/rh.white" --rh-pial "${Output_data}/${i}/surf/rh.pial"

    # Constrain to cortical ribbon (and remove the rest of the ASEG stuff)
    mri_concat --combine --i "${Output_data}/${i}/mri/lh.ribbon.mgz" --i "${Output_data}/${i}/mri/rh.ribbon.mgz" --o "${Output_data}/${i}/mri/ribbon_lr.mgz"
    mri_mask "${Output_data}/${i}/mri/HCP-MMP1+aseg.mgz" "${Output_data}/${i}/mri/ribbon_lr.mgz" "${Output_data}/${i}/mri/HCP-MMP1.mgz"

    # Save stats
    mris_anatomical_stats -a "${Output_data}/${i}/label/lh.HCP-MMP1.annot" -f "${Output_data}/${i}/stats/lh.glasser.stats" "${i}" lh
    mris_anatomical_stats -a "${Output_data}/${i}/label/rh.HCP-MMP1.annot" -f "${Output_data}/${i}/stats/rh.glasser.stats" "${i}" rh

done

