#!/bin/sh
# define the folder where the outputs from the cross-sectional stream are saved

location=/brain_epi/mri/output

# the script will now change directory to the input folder

cd ${location}

# RUN THE BASE STEP:
# the script will loop through each file listed in the called index file (index_base),
# please make sure that you created an index file
# ensure that your subjects with multiple timepoints have the same base name (e.g., subj1)
# but different suffixes (e.g., _TP1, _TP2)
# Your index file should only list the base name
# whereas the suffix can be specified in the script below (see _TP1 and _TP2)

# If you have participants with differnt numbers of timepoints, run these in separate scripts
# as all timepoints called in the script need to be available 
# e.g., one for all subjects with just TP1, one for those with TP1 and TP2 etc. 

for i in `cat *index_base*`
do
echo ${i}
recon-all -base ${i} -tp ${i}_TP1 -tp ${i}_TP2 -all
tail -1 /${location}/${i}/scripts/recon-all-status.log >> ${location}/progress_base.txt
done
