import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

df = pd.read_csv('descriptives_table_20260625.csv')

# rename cohorts
df["cohort"] = df["cohort"].replace("K2H_childhood", "Kids2Health Childhood")
df["cohort"] = df["cohort"].replace("K2H_infancy", "Kids2Health Infancy")
df["cohort"] = df["cohort"].replace("NICAP_T1", "NICAP")
df["cohort"] = df["cohort"].replace("Oregon_ADHD-1000", "Oregon ADHD-1000")
df["cohort"] = df["cohort"].replace("UCI_ECHO", "UCI ECHO")


# 2. separate "mean (sd)" values
def parse_mean_sd(val):
    if pd.isna(val) or val is None:
        return np.nan, np.nan
    val_str = str(val).strip()
    if "(" in val_str:
        parts = val_str.split("(")
        mean = float(parts[0].strip())
        sd = float(parts[1].replace(")", "").strip())
        return mean, sd
    else:
        try:
            return float(val_str), 0.0
        except ValueError:
            return np.nan, np.nan

# 3. extracts percentage from "234 (31%)"
def parse_sex_pct(val):
    if pd.isna(val) or val is None:
        return np.nan
    val_str = str(val).strip()
    if '(' in val_str:
        parts = val_str.split('(')
        pct_str = parts[1].replace(')', '').replace('%', '').strip()
        return float(pct_str)
    else:
        try:
            return float(val_str.replace('%', ''))
        except ValueError:
            return np.nan
        

df[["Age_brain_mean", "Age_brain_sd"]] = df["AGE_brain"].apply(
    lambda x: pd.Series(parse_mean_sd(x))
)
df[["Age_epi_mean", "Age_epi_sd"]] = df["AGE_epi"].apply(
    lambda x: pd.Series(parse_mean_sd(x))
)
df['sex_pct_female'] = df['sex'].apply(parse_sex_pct)

# Female percentage range for each cohort
female_range_map = (
    df.groupby("cohort")["sex_pct_female"]
    .agg(["min", "max"])
    .apply(
        lambda r: (
            f"{r['min']:.0f}%"
            if r["min"] == r["max"]
            else f"{r['min']:.0f}–{r['max']:.0f}%"
        ),
        axis=1,
    )
    .to_dict()
)

# 3. pick metrics per row
def resolve_model_stats(row):
    m = row["model"]
    if m == "brain":
        return row["Age_brain_mean"], row["Age_brain_sd"]
    elif m == "epi":
        return row["Age_epi_mean"], row["Age_epi_sd"]
    else:  
        mean = np.nanmean([row["Age_brain_mean"], row["Age_epi_mean"]])
        sd = np.nanmean([row["Age_brain_sd"], row["Age_epi_sd"]])
        return mean, sd


df[["mean_age", "sd_age"]] = df.apply(
    lambda r: pd.Series(resolve_model_stats(r)), axis=1
)
df = df.dropna(subset=["mean_age"])

model_map = {
    "epi": '#A9C7E8',
    "brain": '#B7D9BE',
    "overlap_brain_epi": '#E7B8AE'
}

df["modeln"] = df["model"].map(model_map)

# 4. Generate Forest Plot Coordinates with row offsets
unique_cohorts = sorted(df["cohort"].unique(), reverse=True)
model_offsets = {"brain": 0.5, "epi": 0.00, "overlap_brain_epi": -0.5}

plot_data = []
ytick_positions = []
ytick_labels = []

for block_idx, cohort in enumerate(unique_cohorts):
    cohort_center_y = block_idx * 1.5
    ytick_positions.append(cohort_center_y)
    ytick_labels.append(cohort)

    cohort_rows = df[df["cohort"] == cohort]
    for idx, row in cohort_rows.iterrows():
        offset = model_offsets.get(row["model"], 0.00)
        exact_y = cohort_center_y + offset

        plot_data.append(
            {
                "cohort": cohort,
                "mean": row["mean_age"],
                "sd": row["sd_age"],
                "y": exact_y,
                "size": row["total_n"],
                "color": row["modeln"], #
                "model": row["model"],
            }
        )

df_plot = pd.DataFrame(plot_data)

# 5. plot Settings 
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Helvetica", "Arial"]
sns.set_theme(style="white")
fig, ax = plt.subplots(figsize=(11, 9))

# Draw baseline separation rows
for base_y in ytick_positions:
    #ax.axhline(y=base_y, color="#e5e5e5", linestyle="-", linewidth=1.2, zorder=1)
    ax.axhline(y=base_y-0.7, color="#e5e5e5", linestyle="-", linewidth=1.2, zorder=1)
    #ax.axhline(
    #    y=base_y + 0.25, color="#f5f5f5", linestyle="--", linewidth=0.8, zorder=1
    #)
    #ax.axhline(
    #    y=base_y - 0.25, color="#f5f5f5", linestyle="--", linewidth=0.8, zorder=1
    #)



# 6. Plot Horizontal Error Bars representing Standard Deviation
ax.errorbar(
    x=df_plot["mean"],
    y=df_plot["y"],
    xerr=df_plot["sd"],
    fmt="none",  
    ecolor="#555555",
    elinewidth=1.2,
    capsize=3,  
    zorder=2,
)

# 7. Scatter Primary Data Bubbles
scatter = ax.scatter(
    x=df_plot["mean"],
    y=df_plot["y"],
    s=df_plot["size"] * 0.35,
    c=df_plot["color"],
    cmap="RdYlBu_r",
    alpha=0.9,
    edgecolors="#1a1a1a",
    linewidths=0.75,
    zorder=3,
    vmin=0,
    vmax=100,
)

# Label sub-rows next to the Y-axis

label_map = {"brain": "brain", "epi": "epi", "overlap_brain_epi": "overlap"}

color_map = {
    "epi": "#1F4E79",
    "brain": "#2E6B3A",
    "overlap_brain_epi": "#A23B2A"
}


for _, row in df_plot.iterrows():
    ax.text(
        x=-1.5,
        y=row["y"],
        s=label_map[row["model"]],
        fontsize=8.5,
        color=color_map[row["model"]], ##555555
        ha="right",
        va="center",
        #fontstyle="italic",
    )

# Axis adjustments
ax.set_yticks(ytick_positions)
ax.set_yticklabels(
    ytick_labels, fontsize=16, fontweight="bold", color="#1a1a1a"
)
ax.set_xlabel(
    "Age (Years, Mean ± SD)", fontsize=14, fontweight="bold", labelpad=12
)
ax.set_xlim(-5, 23)
ax.set_ylim(-0.75, len(unique_cohorts) * 1.5 - 0.25)
sns.despine(ax=ax, left=True, trim=True)

for cohort, y in zip(ytick_labels, ytick_positions):
    ax.text(
        x=-4,   # position between cohort and model labels
        y=y,
        s=female_range_map.get(cohort, ""),
        fontsize=10,
        color="#444444",
        ha="center",
        va="center",
    )
    

    
# 8. Legends
#cbar_ax = fig.add_axes([0.18, 0.02, 0.30, 0.02])
#cbar = fig.colorbar(scatter, cax=cbar_ax, orientation="horizontal")
#cbar.set_label("Sex Ratio (% Female)", fontsize=10, fontweight="bold")
#cbar.ax.tick_params(labelsize=9)

sizes = [100, 500, 1000, 3000]
legend_markers = [
    plt.scatter(
        [], [], s=sz * 0.35, marker="o", c="#bfbfbf", edgecolors="#2b2b2b", alpha=0.8
    )
    for sz in sizes
]
ax.legend(
    legend_markers,
    [str(sz) for sz in sizes],
    title="Sample Size ($n$)",
    title_fontsize=12,
    fontsize=12,
    loc="lower right",
    bbox_to_anchor=(1.2, 0.4),
    scatterpoints=1,
    labelspacing=1.4,
    frameon=True,
    facecolor="#ffffff",
    edgecolor="#dcdcdc",
    borderpad=1.6, 
)

header_y = len(unique_cohorts) * 1.5 - 0.2

ax.text(
    x=-5.5,
    y=header_y,
    s="Cohort",
    fontsize=9,
    fontweight="bold",
    ha="right",
)

ax.text(
    x=-4,
    y=header_y,
    s="% Female",
    fontsize=9,
    fontweight="bold",
    ha="center",
)

ax.text(
    x=-2,
    y=header_y,
    s="Model",
    fontsize=9,
    fontweight="bold",
    ha="center",
)

plt.suptitle(
    "Cohort descriptives",
    fontsize=20,
    fontweight="bold",
    y=0.96,
    x=0.53,
)
plt.subplots_adjust(bottom=0.15, left=0.22)

# Save high resolution
plt.savefig("Descriptives_Fig1_20260625.png", dpi=300, bbox_inches="tight")
plt.close()
