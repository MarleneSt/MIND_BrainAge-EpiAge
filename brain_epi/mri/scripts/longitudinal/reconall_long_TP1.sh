#!/bin/sh
# define the folder where the outputs from the cross-sectional stream & base processing stepare saved

location=/brain_epi/mri/output

# the script will now change directory to the input folder

cd ${location}

# RUN THE LONG STEP:
# the script will loop through each file listed in the called index file (index_base),
# please make sure that you created an index file 
# ensure that your subjects with multiple timepoints have the same base name (e.g., subj1)
# but different suffixes (e.g., _TP1, _TP2)
# Your index file should only list the base name (consistent with the output of the BASE Step)
# whereas the suffix for the specific timepoint (see _TP1)

# This script should be run separately by time point (i.e., one script for TP1, one for TP2 etc.)

for i in `cat *index_base*`
do
echo ${i}
recon-all -long ${i}_TP1 ${i} -all
tail -1 ${location}/${i}/scripts/recon-all-status.log >> ${location}/progress_long.txt
done
