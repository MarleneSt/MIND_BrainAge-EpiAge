"""
Fig 4B heatmap — Brain-PAR x Epi-PAR pooled associations
Requires: pip install plotly kaleido pandas numpy scipy
Run from the Reproducibilty folder (python3 20_fig4b_heatmap.py),
after 01_main_meta_wMAE.R; reads its model-wise estimates and writes to results/plots
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

# ── 1. Load data ───────────────────────────────────────────────────────────────
# CHANGED: input = deduplicated model-wise estimates (GenR duplicates removed; k = 2448)
data = pd.read_csv("results/assoc_model2_PAR_modelwise_final/modelwise_estimates.csv")   # [REPRODUCIBILITY] relative path to output of 01_main_meta_wMAE.R (was: "modelwise_estimates.csv")
data = data.rename(columns={"ci.lb": "ci_lb", "ci.ub": "ci_ub"})   # beta already named "beta"
data = data[(data["brain_model"] != "Pooled") & (data["epi_model"] != "Pooled")]

# ── 2. Pivot to matrix (rows = brain models, columns = epi clocks) ─────────────
beta  = data.pivot_table(index="brain_model", columns="epi_model", values="beta",  aggfunc="first")
ci_lb = data.pivot_table(index="brain_model", columns="epi_model", values="ci_lb", aggfunc="first")
ci_ub = data.pivot_table(index="brain_model", columns="epi_model", values="ci_ub", aggfunc="first")

# ── 3. Hierarchical clustering to order rows and columns ───────────────────────
z = beta.fillna(0).values
row_order = leaves_list(linkage(pdist(z,   metric="correlation"), method="average"))
col_order = leaves_list(linkage(pdist(z.T, metric="correlation"), method="average"))

brain_models = [beta.index[i]   for i in row_order]
epi_clocks   = [beta.columns[i] for i in col_order]

beta  = beta.reindex(index=brain_models,  columns=epi_clocks)
ci_lb = ci_lb.reindex(index=brain_models, columns=epi_clocks)
ci_ub = ci_ub.reindex(index=brain_models, columns=epi_clocks)

# ── 4. Mark significant cells (95% CI excludes zero) ──────────────────────────
significant = (ci_lb.values > 0) | (ci_ub.values < 0)

# Show beta value only in significant cells, empty string otherwise
cell_labels = np.where(significant, np.vectorize(lambda v: f"{v:.2f}")(beta.values), "")
cell_labels[np.isnan(beta.values)] = ""

# ── 5. Build heatmap ───────────────────────────────────────────────────────────
abs_max = float(np.nanmax(np.abs(beta.values)))

fig = go.Figure(go.Heatmap(
    z          = beta.values,
    x          = epi_clocks,
    y          = brain_models,
    customdata = np.stack([ci_lb.values, ci_ub.values], axis=-1),
    colorscale = "RdBu_r",
    zmid=0, zmin=-abs_max, zmax=abs_max,
    colorbar   = dict(title="β", thickness=14, len=0.8),
    text       = cell_labels,
    texttemplate = "%{text}",
    textfont   = {"size": 9},
    hovertemplate = (
        "<b>Brain:</b> %{y}<br>"
        "<b>Epi:</b> %{x}<br>"
        "<b>β</b> = %{z:.3f} [%{customdata[0]:.3f}, %{customdata[1]:.3f}]"
        "<extra></extra>"
    ),
))

fig.update_layout(
    title  = "", #keeping empty
    width  = max(1400, 120 * len(epi_clocks) + 260),
    height = max(500,   90 * len(brain_models) + 300),
    margin = dict(l=130, r=80, t=60, b=200),
    xaxis  = dict(tickangle=-45, side="bottom", title="Epigenetic clock"),
    yaxis  = dict(autorange="reversed", title="Brain age model"),
)

# ── 6. Save as 600 DPI PNG ─────────────────────────────────────────────────────
fig.write_image("results/plots/plot_assoc_pairwise_heatmap.png", scale=600/96, format="png")   # [REPRODUCIBILITY] saved to results/plots with a descriptive name (was: "fig4b_heatmap_dedup.png")
print("Saved: results/plots/plot_assoc_pairwise_heatmap.png")
