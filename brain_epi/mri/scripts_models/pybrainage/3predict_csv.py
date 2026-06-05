#!/usr/bin/env python3
import os, sys, json, warnings
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import onnxruntime as ort

MODEL_NAME = "ExtraTreesModel.onnx"
FEATS_JSON = "features.json"
OUT_CSV = "PyBrainAge_Output.csv"
OUT_PNG = "PyBrainAge_Output.png"

def die(msg, code=1):
    print(f"ERROR: {msg}"); sys.exit(code)

def normalize_feature_names(cols):
    out = []
    for c in cols:
        c2 = c.replace("_and_", "&")
        if "Left-Thalamus-Proper" not in c2: c2 = c2.replace("Left-Thalamus","Left-Thalamus-Proper")
        if "Right-Thalamus-Proper" not in c2: c2 = c2.replace("Right-Thalamus","Right-Thalamus-Proper")
        out.append(c2)
    return out

def main():
    if len(sys.argv) != 2:
        print("Usage: python predict_csv.py /path/to/freesurfer.csv"); sys.exit(1)
    csv_path = os.path.abspath(sys.argv[1])
    if not os.path.isfile(csv_path): die(f"Input CSV not found: {csv_path}")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    model_path = os.path.join(script_dir, MODEL_NAME)
    feats_path = os.path.join(script_dir, FEATS_JSON)
    if not os.path.isfile(model_path): die(f"Missing model: {model_path}")

    # Load expected feature order
    feats = None
    if os.path.isfile(feats_path):
        try:
            feats = json.load(open(feats_path))["feature_names_in_"]
        except Exception as e:
            warnings.warn(f"Could not read {feats_path}: {e}")

    # ONNX session (CPU)
    try:
        sess = ort.InferenceSession(model_path, providers=["CPUExecutionProvider"])
        in_name = sess.get_inputs()[0].name
        out_name = sess.get_outputs()[0].name
    except Exception as e:
        die(f"Failed to load ONNX model: {e}")

    # Load CSV
    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        die(f"Failed to read CSV: {csv_path}: {e}")
    if df.shape[1] < 3: die("CSV must have at least 3 columns: ID, Age, features")

    ids = df.iloc[:, 0].astype(str)
    ages = pd.to_numeric(df.iloc[:, 1], errors="coerce")
    feat_df = df.iloc[:, 2:].copy()
    feat_df.columns = normalize_feature_names(feat_df.columns)

    # Enforce training order
    if feats:
        missing = [c for c in feats if c not in feat_df.columns]
        if missing:
            die(f"Missing features required by model: {len(missing)} (e.g., {', '.join(missing[:8])})")
        extra = [c for c in feat_df.columns if c not in feats]
        if extra:
            warnings.warn(f"Ignoring {len(extra)} extra columns (e.g., {', '.join(extra[:8])})")
        feat_df = feat_df[feats]

    # Validate numeric
    feat_df = feat_df.apply(pd.to_numeric, errors="coerce")
    bad = feat_df.isna().any(axis=1) | np.isinf(feat_df.values).any(axis=1) | ages.isna()
    if bad.any():
        bad_ids = ids[bad].tolist()
        die(f"Non-numeric or NaN/Inf values found for {len(bad_ids)} row(s), e.g., {bad_ids[:5]}")

    # IMPORTANT: Do NOT scale here (scaler is in the ONNX graph)
    X = feat_df.to_numpy(dtype=np.float32, copy=False)

    try:
        y_pred = sess.run([out_name], {in_name: X})[0].ravel()
    except Exception as e:
        die(f"ONNX inference failed: {e}")

    out = pd.DataFrame({"ID": ids.values,
                        "Age": ages.values.astype(float),
                        "BrainAge": y_pred.astype(float)})
    out["BrainPAD"] = out["BrainAge"] - out["Age"]

    out_dir = os.path.dirname(csv_path)
    out_csv = os.path.join(out_dir, OUT_CSV)
    out_png = os.path.join(out_dir, OUT_PNG)

    out.to_csv(out_csv, index=False)
    print(f"Wrote {out_csv} with {len(out)} rows")

    try:
        plt.figure(figsize=(6, 6))
        plt.scatter(out["Age"], out["BrainAge"], s=18)
        lim_min = float(np.floor(min(out["Age"].min(), out["BrainAge"].min())/5)*5)
        lim_max = float(np.ceil (max(out["Age"].max(), out["BrainAge"].max())/5)*5)
        plt.plot([lim_min, lim_max], [lim_min, lim_max])
        plt.xlim(lim_min, lim_max); plt.ylim(lim_min, lim_max)
        plt.xlabel("Chronological Age"); plt.ylabel("Predicted Brain Age")
        plt.title("Brain Age vs Chronological Age")
        plt.tight_layout(); plt.savefig(out_png, dpi=200); plt.close()
        print(f"Wrote {out_png}")
    except Exception as e:
        warnings.warn(f"Failed to save plot: {e}")

if __name__ == "__main__":
    main()
