import os
import shutil
import gzip
import json

# Input folder containing T1-weighted scans
input_folder = '/brain_epi/mri/raw'

# Output BIDS directory
bids_folder = '/brain_epi/mri/raw_BIDS'

# Ensure output directory exists
os.makedirs(bids_folder, exist_ok=True)

# Create dataset_description.json
dataset_description = {
    "Name": "My Neuroimaging Dataset",
    "BIDSVersion": "1.8.0",
    "License": "CC0",
    "Authors": ["Your Name"],
    "Acknowledgements": "",
    "HowToAcknowledge": "",
    "Funding": ["Your Funding Source"],
    "ReferencesAndLinks": []
}
with open(os.path.join(bids_folder, 'dataset_description.json'), 'w') as f:
    json.dump(dataset_description, f, indent=4)

# Process files
for file in os.listdir(input_folder):
    if file.startswith('sub-') and file.endswith('_T1w.nii'): #please define here the prefix of your subject names as well as the suffix
        # Extract participant ID (e.g., sub-02 from sub-02_T1w.nii)
        participant_id = file.split('_')[0]
        
        # Create participant-specific anat directory
        participant_folder = os.path.join(bids_folder, participant_id, "anat")
        os.makedirs(participant_folder, exist_ok=True)
        
        # Define BIDS-compatible filename (e.g., sub-02_T1w.nii.gz)
        bids_filename = f"{participant_id}_T1w.nii.gz"
        
        # Copy or compress the file into the BIDS folder
        input_file = os.path.join(input_folder, file)
        output_file = os.path.join(participant_folder, bids_filename)
        
        if file.endswith('.nii'):  # Compress to .nii.gz
            with open(input_file, 'rb') as f_in:
                with gzip.open(output_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            shutil.copy(input_file, output_file)

print(f"BIDS dataset created successfully at {bids_folder}")
