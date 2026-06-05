#!/usr/bin/env python

import os
import argparse
import subprocess
import numpy as np
from tqdm import tqdm
from threading import Lock
from threading import Thread
import nibabel as nib
import matplotlib.pyplot as plt
import pandas as pd

# ---------------------- Argument parsing ---------------------- #
parser = argparse.ArgumentParser(description="Script for Pyment preprocessing with FreeSurfer and FSL.")
parser.add_argument("--raw_dir", required=True, help="Folder containing raw .nii or .nii.gz images")
parser.add_argument("--output_dir", required=True, help="Destination folder for FreeSurfer+FSL output (can be prepopulated with FS output)")
parser.add_argument("--labels_csv", required=False, help="Optional CSV file with subject IDs (SUBJID)")
parser.add_argument("--template", required=True, help="Path to MNI152_T1_1mm_brain.nii.gz template")
parser.add_argument("--processing", required=True, help="Path to Pyment freesurfer_and_fsl.sh script")
parser.add_argument('--qc_dir', required=False, help="Path to save pre-processing log")
parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads (default = 4)")
args = parser.parse_args()

log_file = os.path.join(args.qc_dir, "pyment_preprocessing_log.txt")
log_lock = Lock()
def log_message(msg):
    from datetime import datetime
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    full_msg = f'[{timestamp}] {msg}'
    # print suppressed
    with log_lock:
        with open(log_file, "a") as f:
            f.write(msg + "\n")


raw_folder = args.raw_dir
freesurfer_fsl_folder = args.output_dir
MNI152_TEMPLATE = args.template
preprocessing_script = args.processing
NUM_THREADS = args.threads

print("🧠 Starting preprocessing with Pyment CLI...")
print(f"➡️  Raw images folder: {raw_folder}")
print(f"➡️  Output folder: {freesurfer_fsl_folder}")
print(f"➡️  Template: {MNI152_TEMPLATE}")
print(f"➡️  Processing script: {preprocessing_script}")
print(f"➡️  Threads: {NUM_THREADS}")

# ---------------------- Check paths ---------------------- #
if not os.path.isfile(preprocessing_script):
    raise ValueError(f"❌ Unable to find FreeSurfer+FSL preprocessing script at: {preprocessing_script}")

if not os.path.isdir(freesurfer_fsl_folder):
    os.makedirs(freesurfer_fsl_folder)
    print(f"📁 Created output directory: {freesurfer_fsl_folder}")

# ---------------------- Threaded worker ---------------------- #
class Worker(Thread):
    def __init__(self, files, source, destination):
        super().__init__()
        self.files = files
        self.source = source
        self.destination = destination

    def run(self):
        thread_log_prefix = f"[Thread {self.name}]"

        for filename in tqdm(self.files, desc=f"Thread {self.name}", position=int(self.name), leave=False):
            if not (filename.endswith('.nii') or filename.endswith('.nii.gz')):
                log_message(f"{thread_log_prefix} Skipping non-NIfTI file: {filename}")
                continue

            subject = filename.replace('.nii.gz', '').replace('.nii', '')
            dest = os.path.join(self.destination, subject)

            brainmask_mgz = os.path.join(dest, 'mri', 'brainmask.mgz')
            if os.path.exists(brainmask_mgz):
                log_message(f"{thread_log_prefix} ✅ {subject}: 🧠 brainmask.mgz found — skipping FreeSurfer autorecon and starting at mri_convert")
                ret = subprocess.call([
                    "bash", preprocessing_script,
                    "-f", os.path.join(self.source, filename),
                    "-d", dest,
                    "-t", MNI152_TEMPLATE
                ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            else:
                log_message(f"{thread_log_prefix} 🚧 {subject}: brainmask.mgz not found, running full FreeSurfer and FSL pipeline")
                ret = subprocess.call([
                    "bash", preprocessing_script,
                    "-f", os.path.join(self.source, filename),
                    "-d", dest,
                    "-t", MNI152_TEMPLATE
                ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            if ret != 0:
                log_message(f"{thread_log_prefix} ❌ {subject}: Error during processing")
            else:
                log_message(f"{thread_log_prefix} ✅ {subject}: Preprocessing completed.")


# ---------------------- Run threads ---------------------- #
#  Identify subjects to process #
if args.labels_csv:
    df_labels = pd.read_csv(args.labels_csv)
    subject_ids = df_labels['SUBJID'].tolist()
    images = []
    for subj in subject_ids:
        for ext in ['.nii.gz', '.nii']:
            candidate_path = os.path.join(raw_folder, f"{subj}{ext}")
            if os.path.isfile(candidate_path):
                images.append(f"{subj}{ext}")
                break
else:
    # fallback: include all NIfTI images in the raw folder
    images = [f for f in os.listdir(raw_folder) if f.endswith('.nii') or f.endswith('.nii.gz')]


print(f"🔍 Found {len(images)} NIfTI images to process.")


batches = np.array_split(images, NUM_THREADS)

threads = [Worker(batch, raw_folder, freesurfer_fsl_folder) for batch in batches]
for i, worker in enumerate(threads):
    worker.name = str(i)
    worker.start()

for worker in threads:
    worker.join()


# ✅ Summarize success/failure
success_count = 0
fail_count = 0
with open(log_file) as f:
    for line in f:
        if "Preprocessing completed" in line:
            success_count += 1
        elif "Error during processing" in line:
            fail_count += 1

summary = f"✅ Finished preprocessing {success_count + fail_count} images. {fail_count} failed."
print(summary)
with open(log_file, "a") as f:
    f.write(summary + "\n")

print("✅ Preprocessing complete.")