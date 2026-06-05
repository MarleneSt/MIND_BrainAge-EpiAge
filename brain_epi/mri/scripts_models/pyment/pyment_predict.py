import os
import argparse
import nibabel as nib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from pyment.models import RegressionSFCN

# -----------------------------
# Argument Parsing
# -----------------------------
parser = argparse.ArgumentParser(description="Run brain age prediction using Pyment")
parser.add_argument('--output_dir', required=True, help="Path to subject folders with cropped.nii.gz (usually FS output folders)")
parser.add_argument('--output_csv', required=True, help="Path to save predictions CSV")
parser.add_argument('--qc_dir', required=False, help="Path to save QC images & predict log")
parser.add_argument('--labels_csv', required=False, help="Optional CSV file with subject IDs (SUBJID)")
args = parser.parse_args()

# -----------------------------
# Helper function for logging
# -----------------------------
def log(message):
    print(message)
    if args.qc_dir:
        os.makedirs(args.qc_dir, exist_ok=True)
        with open(os.path.join(args.qc_dir, "pyment_predict_log.txt"), "a") as log_file:
            log_file.write(message + "\n")

# -----------------------------
# Model Setup
# -----------------------------
WEIGHTS = 'brain-age-2022'
MIN_AGE = 3
MAX_AGE = 95

model = RegressionSFCN(weights=WEIGHTS, prediction_range=(MIN_AGE, MAX_AGE))

# -----------------------------
# Identify Subject IDs
# -----------------------------
if args.labels_csv:
    df_labels = pd.read_csv(args.labels_csv)
    subject_ids = df_labels['SUBJID'].tolist()
else:
    subject_ids = sorted(
    d for d in os.listdir(args.output_dir)
    if os.path.isdir(os.path.join(args.output_dir, d))
    and os.path.isfile(os.path.join(args.output_dir, d, 'mri', 'cropped.nii.gz'))
)

log(f"\U0001F9E0 Found {len(subject_ids)} subject(s) to process.")

# -----------------------------
# Predict Ages
# -----------------------------
predictions = []

for subj in tqdm(subject_ids, desc="Predicting"):
    cropped_path = os.path.join(args.output_dir, subj, 'mri', 'cropped.nii.gz')
    if not os.path.exists(cropped_path):
        log(f"⚠️ Skipping {subj}: cropped.nii.gz not found.")
        continue

    try:
        img = nib.load(cropped_path).get_fdata()
        img = np.expand_dims(img, axis=0)  # add channel dimension
        pred = model.predict(img, verbose=0)[0]
        pred = model.postprocess(pred)
        predictions.append({'SUBJID': subj, 'Pyment_predage': pred})
    except Exception as e:
        log(f"❌ Failed prediction for {subj}: {e}")

# -----------------------------
# Save CSV
# -----------------------------
pred_df = pd.DataFrame(predictions)
pred_df.to_csv(args.output_csv, index=False)
log(f"✅ Saved predictions to {args.output_csv}")

# -----------------------------
# QC Image Generation
# -----------------------------
if args.qc_dir:
    os.makedirs(args.qc_dir, exist_ok=True)
    log(f"🖼️  Saving QC images to: {args.qc_dir}")

    for entry in tqdm(predictions, desc="Creating QC images"):
        subj = entry['SUBJID']
        pred = entry['Pyment_predage']
        cropped_path = os.path.join(args.output_dir, subj, 'mri', 'cropped.nii.gz')

        try:
            img = nib.load(cropped_path).get_fdata()
            mid_slice = img.shape[2] // 2
            plt.imshow(img[:, :, mid_slice], cmap='gray')
            plt.axis('off')
            plt.title(f"{subj} | Predicted: {pred:.2f}")
            plt.savefig(os.path.join(args.qc_dir, f"{subj}.png"), bbox_inches='tight')
            plt.close()
        except Exception as e:
            log(f"⚠️ Could not create QC for {subj}: {e}")

log(f"✅ Finished brain age prediction for {len(predictions)} images")
