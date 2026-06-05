#!/bin/bash

#could potentially load freesurfer here? 

# set path to your freesurfer output folder
SUB_dir=/brain_epi/mri/output

# set path to where the extracted features should be saved
OUTPUT_dir=/brain_epi/mri/extracted_features

cd $SUB_dir

# This particular version of the script will loop through the participants in the folder with directories that start with "sub"
    # and create a table from ?.stats files for all of those participants
    #focused on the Destrieux and aseg atlases 
    # change the subject base name below in line with your participant names

export SUBJECTS_DIR=$PWD

list="`ls -d sub*`" #change this in line with the base of your subject names

asegstats2table --subjects $list --meas volume --skip --tablefile $OUTPUT_dir/aseg.stats.txt
aparcstats2table --subjects $list --hemi lh --parc aparc.a2009s --meas thickness --skip --tablefile $OUTPUT_dir/lh.aparcstats.a2009s.txt
aparcstats2table --subjects $list --hemi rh --parc aparc.a2009s --meas thickness --skip --tablefile $OUTPUT_dir/rh.aparcstats.a2009s.txt
