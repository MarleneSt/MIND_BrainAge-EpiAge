#!/usr/bin/env python3
import os
import sys
import csv
import tempfile
import subprocess

def run_cmd(cmd, env):
    res = subprocess.run(cmd, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"Command failed ({res.returncode}): {' '.join(cmd)}\n{res.stderr}")
    return res

def read_two_rows_csv(path):
    with open(path, newline="") as f:
        rows = list(csv.reader(f))
    if len(rows) < 2:
        raise ValueError(f"{path}: expected at least two rows")
    return rows[0], rows[1]

def subjects_in(root):
    subs = []
    for name in sorted(os.listdir(root)):
        p = os.path.join(root, name)
        if os.path.isdir(p) and os.path.isdir(os.path.join(p, "stats")):
            subs.append(name)
    return subs

def normalize_fs_name_to_template(name: str) -> str:
    # Match FreeSurfer label quirks to template naming
    name = name.replace("_G_and_S_", "_G&S_")
    name = name.replace("_and_", "&")
    return name

def build_feature_map_from_tables(subj, lh_csv, rh_csv, seg_csv):
    lh_header, lh_vals = read_two_rows_csv(lh_csv)
    rh_header, rh_vals = read_two_rows_csv(rh_csv)
    seg_header, seg_vals = read_two_rows_csv(seg_csv)

    fmap = {}

    # LH aparc thickness
    for k, v in zip(lh_header[1:], lh_vals[1:]):
        fmap[normalize_fs_name_to_template(k)] = v

    # RH aparc thickness
    for k, v in zip(rh_header[1:], rh_vals[1:]):
        fmap[normalize_fs_name_to_template(k)] = v

    # Aseg volumes (+ thalamus proper aliases)
    for k, v in zip(seg_header[1:], seg_vals[1:]):
        fmap[k] = v
        if k == "Left-Thalamus":
            fmap["Left-Thalamus-Proper"] = v
        elif k == "Right-Thalamus":
            fmap["Right-Thalamus-Proper"] = v

    return fmap

def main():
    if len(sys.argv) != 2:
        print("Usage: fs2brainage.py /path/to/freesurfer/subjects")
        sys.exit(1)

    subjects_dir = os.path.abspath(sys.argv[1])
    if not os.path.isdir(subjects_dir):
        print(f"Subjects dir not found: {subjects_dir}")
        sys.exit(1)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    template_path = os.path.join(script_dir, "ROIS_input_template.txt")
    if not os.path.isfile(template_path):
        print(f"Template not found: {template_path}")
        sys.exit(1)

    # Load template (one item per line)
    with open(template_path) as f:
        template = [line.strip() for line in f if line.strip()]

    if len(template) < 2 or template[0] != "ID" or template[1] != "Age":
        print("Template must start with two lines: 'ID' then 'Age'.")
        sys.exit(1)

    env = os.environ.copy()
    env["SUBJECTS_DIR"] = subjects_dir

    subjects = subjects_in(subjects_dir) 
    if not subjects:
        print("No subjects with a 'stats' folder found.")
        sys.exit(1)
        
    # Exclude one subject, change made by Marlene to remove a subjecz from the process that was causing it to fail 
    exclude = {"15821A", "9165A"}
    subjects = [s for s in subjects if s not in exclude]

    out_csv = os.path.join(subjects_dir, "brainage.csv")
    with open(out_csv, "w", newline="") as f_out:
        w = csv.writer(f_out)
        w.writerow(template)

        for subj in subjects:
            print(f"Processing {subj}...")
            stats_dir = os.path.join(subjects_dir, subj, "stats")
            need = [
                os.path.join(stats_dir, "lh.aparc.a2009s.stats"),
                os.path.join(stats_dir, "rh.aparc.a2009s.stats"),
                os.path.join(stats_dir, "aseg.stats"),
            ]
            if not all(os.path.exists(p) for p in need):
                print(f"  Skipping {subj}: missing one of lh/rh/aseg stats files")
                continue

            with tempfile.TemporaryDirectory() as tmp:
                lh_csv = os.path.join(tmp, "lh.txt")
                rh_csv = os.path.join(tmp, "rh.txt")
                seg_csv = os.path.join(tmp, "seg.txt")

                run_cmd([
                    "aparcstats2table", "--subjects", subj,
                    "--hemi", "lh", "--meas", "thickness",
                    "--parc", "aparc.a2009s", "--delimiter", "comma",
                    "--tablefile", lh_csv
                ], env)

                run_cmd([
                    "aparcstats2table", "--subjects", subj,
                    "--hemi", "rh", "--meas", "thickness",
                    "--parc", "aparc.a2009s", "--delimiter", "comma",
                    "--tablefile", rh_csv
                ], env)

                run_cmd([
                    "asegstats2table", "--subjects", subj,
                    "--delimiter", "comma",
                    "--tablefile", seg_csv
                ], env)

                fmap = build_feature_map_from_tables(subj, lh_csv, rh_csv, seg_csv)

            row = []
            missing = []
            for i, key in enumerate(template):
                if i == 0:      # ID
                    row.append(subj)
                elif i == 1:    # Age
                    row.append("42")
                else:
                    if key in fmap:
                        row.append(fmap[key])
                    else:
                        missing.append(key)
                        row.append("")

            if missing:
                raise RuntimeError(
                    f"{subj}: missing {len(missing)} ROI(s) required by template. "
                    f"Examples: {', '.join(missing[:5])}"
                )

            w.writerow(row)

    print(f"Wrote {out_csv}")

if __name__ == "__main__":
    main()
