"""
dbn_antspynet.py

Run DeepBrainNet brain-age prediction on all T1 images in a folder
and save results to a CSV. The function should be run from within the
pre-specified conda environment. 

Example:
    python dbn_antspynet.py \
        --raw_dir /brain_epi/mri/raw \
        --output_csv /brain_epi/mri/models/DBN/DBN_predictions.csv \
        --subject_list_txt /brain_epi/mri/subject_list.txt \
        --num_workers 4

Notes:
- subject_list_txt = optional should contain IDs matching filenames WITHOUT .nii/.nii.gz
"""

import os

# ---------------------------------------------------------------------
# Environment settings BEFORE TensorFlow is loaded
# ---------------------------------------------------------------------
os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VTK_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

import sys
import random
import argparse
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import tensorflow as tf

# Force TensorFlow CPU-only
tf.config.set_visible_devices([], "GPU")

try:
    tf.config.experimental.enable_op_determinism(True)
except Exception:
    pass

tf.config.threading.set_intra_op_parallelism_threads(1)
tf.config.threading.set_inter_op_parallelism_threads(1)

# Global seeds (main process)
random.seed(0)
np.random.seed(0)
tf.random.set_seed(0)

import ants
import antspynet
import pandas as pd
from tqdm import tqdm


# -----------------------------
# Argument Parsing (Pyment-style)
# -----------------------------
parser = argparse.ArgumentParser(
    description="Run DeepBrainNet (ANTsPyNet) brain-age prediction on a folder of T1 NIfTIs"
)
parser.add_argument("--raw_dir", required=True, help="Folder containing .nii / .nii.gz images")
parser.add_argument("--output_csv", required=True, help="Path to output CSV file")
parser.add_argument(
    "--subject_list_txt",
    required=False,
    help=("Optional flag: Text file with one subject ID per line, matching filename WITHOUT .nii / .nii.gz "
          '(e.g. "sub-18_T1w")')
)
parser.add_argument(
    "--num_workers",
    type=int,
    default=int(os.environ.get("DBN_NUM_WORKERS", "4")),
    help="Number of worker processes (default: DBN_NUM_WORKERS env var or 4)"
)
parser.add_argument(
    "--overwrite",
    action="store_true",
    help="Optional flag: Overwrite output CSV instead of appending"
)
parser.add_argument(
    "--skip_existing",
    action="store_true",
    help="Optional flag: If output CSV exists, skip subjects already present in it"
)
args = parser.parse_args()


def remove_nifti_suffix(filename: str) -> str:
    if filename.endswith(".nii.gz"):
        return filename[:-7]
    elif filename.endswith(".nii"):
        return filename[:-4]
    return filename


def worker_init():
    """Runs once per worker process."""
    random.seed(0)
    np.random.seed(0)
    tf.random.set_seed(0)


def process_one_image(task):
    image_folder, filename = task
    subject_id = remove_nifti_suffix(filename)
    image_path = os.path.join(image_folder, filename)

    try:
        image = ants.image_read(image_path)
        deep = antspynet.utilities.brain_age(
            image,
            do_preprocessing=True,
            number_of_simulations=0
        )
        predicted_age = float(deep["predicted_age"])
        return {"SUBJID": subject_id, "DBN_predage": predicted_age}

    except Exception as e:
        print(f"Error processing {subject_id}: {e}")
        return None


def main():
    image_folder = args.raw_dir
    output_csv_path = args.output_csv
    subject_list_path = args.subject_list_txt

    if not os.path.isdir(image_folder):
        raise NotADirectoryError(f"Input folder does not exist or is not a directory: {image_folder}")

    out_dir = os.path.dirname(output_csv_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # Load optional subject list
    subject_ids_to_include = None
    if subject_list_path:
        if not os.path.exists(subject_list_path):
            raise FileNotFoundError(f"Subject list file not found: {subject_list_path}")
        with open(subject_list_path, "r") as f:
            subject_ids_to_include = {line.strip() for line in f if line.strip()}
        print(f"Loaded {len(subject_ids_to_include)} subject IDs from {subject_list_path}")

    # List NIfTIs
    image_files = [
        f for f in os.listdir(image_folder)
        if (f.endswith(".nii") or f.endswith(".nii.gz"))
        and (subject_ids_to_include is None or remove_nifti_suffix(f) in subject_ids_to_include)
    ]
    image_files = sorted(image_files)

    if not image_files:
        print("No matching .nii/.nii.gz files found in:", image_folder)
        sys.exit(0)

    # Overwrite handling
    if args.overwrite and os.path.exists(output_csv_path):
        os.remove(output_csv_path)

    # Skip-existing handling
    if args.skip_existing and os.path.exists(output_csv_path):
        try:
            existing = pd.read_csv(output_csv_path)
            if "SUBJID" in existing.columns:
                done = set(existing["SUBJID"].astype(str).tolist())
                before = len(image_files)
                image_files = [f for f in image_files if remove_nifti_suffix(f) not in done]
                print(f"Skipping {before - len(image_files)} subject(s) already in CSV.")
        except Exception as e:
            print(f"Warning: could not read existing CSV for --skip_existing ({e}). Continuing without skipping.")

    if not image_files:
        print("Nothing to do (all subjects already present).")
        sys.exit(0)

    print(f"Found {len(image_files)} image(s) to process.")
    print(f"Using {args.num_workers} worker process(es).")

    tasks = [(image_folder, f) for f in image_files]
    rows = []

    with ProcessPoolExecutor(max_workers=args.num_workers, initializer=worker_init) as executor:
        for result in tqdm(
            executor.map(process_one_image, tasks),
            total=len(tasks),
            desc="Processing Images",
            unit="file",
        ):
            if result is not None:
                rows.append(result)

    if not rows:
        print("No successful predictions were produced.")
        sys.exit(0)

    df = pd.DataFrame(rows)

    # If file exists and not overwriting: append without header
    if os.path.exists(output_csv_path) and not args.overwrite:
        df.to_csv(output_csv_path, mode="a", header=False, index=False)
    else:
        df.to_csv(output_csv_path, index=False)

    print(f"Predicted ages have been saved to {output_csv_path}")


if __name__ == "__main__":
    try:
        mp.set_start_method("spawn")
    except RuntimeError:
        pass
    main()
