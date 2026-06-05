#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import time

def report_elapsed_time(start_time, comment=""):
    elapsed = int(time.time() - start_time)
    if comment:
        print(f"{comment} completed in {elapsed} seconds.")
    return elapsed

def find_niftis(folder):
    out = []
    for name in sorted(os.listdir(folder)):
        if name.startswith('.'):
            continue
        if name.endswith('.nii.gz'):
            out.append(os.path.join(folder, name))
        elif name.endswith('.nii'):
            # avoid adding .nii if a .nii.gz with same stem also exists
            gz = os.path.join(folder, name + '.gz')
            if not os.path.exists(gz):
                out.append(os.path.join(folder, name))
    return out

def strip_ext(path):
    base = os.path.basename(path)
    if base.endswith('.nii.gz'):
        return base[:-7]
    if base.endswith('.nii'):
        return base[:-4]
    return os.path.splitext(base)[0]

def run(cmd, env):
    print(">>", " ".join(cmd))
    subprocess.run(cmd, check=True, env=env)

def main():
    ap = argparse.ArgumentParser(
        description="Batch FreeSurfer recon-all on all NIfTI files in a folder."
    )
    ap.add_argument("t1_dir", help="Folder containing T1 NIfTI files (*.nii or *.nii.gz)")
    ap.add_argument("--overwrite", action="store_true",
                    help="If a subject folder already exists, remove it and re-run")
    ap.add_argument("--openmp", type=int, default=max(os.cpu_count() - 1, 1),
                    help="Threads for recon-all (-openmp). Default: CPUs-1")
    args = ap.parse_args()
    start_time = time.time()
    t1_dir = os.path.abspath(args.t1_dir)
    if not os.path.isdir(t1_dir):
        print(f"ERROR: Not a directory: {t1_dir}")
        sys.exit(1)

    # Export SUBJECTS_DIR for this run
    env = os.environ.copy()
    env["SUBJECTS_DIR"] = t1_dir
    print(f"SUBJECTS_DIR set to: {t1_dir}")

    # Sanity check: recon-all available?
    try:
        subprocess.run(["recon-all", "-version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, env=env)
    except Exception as e:
        print("ERROR: 'recon-all' not found or FreeSurfer not sourced. Source SetUpFreeSurfer.sh first.")
        print(e)
        sys.exit(1)

    niftis = find_niftis(t1_dir)
    if not niftis:
        print("No .nii or .nii.gz files found.")
        sys.exit(0)

    for nifti in niftis:
        baseID = strip_ext(nifti)
        subj_dir = os.path.join(t1_dir, baseID)
        if os.path.isdir(subj_dir):
            if args.overwrite:
                print(f"Subject {baseID} already exists — removing for overwrite...")
                subprocess.run(["rm", "-rf", subj_dir], check=True)
            else:
                print(f"Subject {baseID} already exists — skipping (use --overwrite to re-run).")
                continue

        cmd = ["recon-all", "-i", nifti, "-s", baseID, "-all", "-openmp", str(args.openmp)]
        try:
            run(cmd, env)
        except subprocess.CalledProcessError as e:
            print(f"ERROR: recon-all failed for {baseID}")
            print(e)
            # continue to next file
            continue

    report_elapsed_time(start_time, "done")

if __name__ == "__main__":
    main()
