#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# exec(open("/beegfs/data/CostaLab/Spatial_ATAC_Analysis_KuppeLab/results/SpatialClones/Kidney-Spatial-multi-modal-project/04_Mitochondrial_Clone_Identification/00602_MSS_Trajectory.py").read())

"""
We load the MSS and trajectory values and correlate them using the Spearman correlation.


Input
-----
The MSS_per_spot.csv files written by 00500, one per sample, with one row per spot:

    (index)          spot ID
    Sample           sample name as used in 00500 (one library per tissue section)
    MSS              heteroplasmy-weighted MLC score sum of the spot
    mt_depth         mitochondrial sequencing depth of the spot (mgatk depthTable)
    MSS_per_depth    not used; recomputed here from MSS and mt_depth
    MSS_resid_depth  not used; recomputed here per sample
    Annotation       aligned cell type

The donor of each sample is set in sample_to_donor. Samples without an entry are their
own donor.

Stratified means that values per expressed as percentile, separate for each cell type.
So, a stratified MSS of 0.5 means in the middle of the observed values for this cell
type. This makes the values more comparable. 

"""

import os
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
from scipy.stats import pearsonr, spearmanr


color_map = {
    "ATL_DTL"           : "#1f77b4",
    "CD"                : "#ff7f0e",
    "DCT"               : "#2ca02c",
    "Endo"              : "black",
    "Fibroblast"        : "#d62728",
    "Immune"            : "#9467bd",
    "Injured PTs"       : "#8c564b",
    "Intercalated_Cell" : "#e377c2",
    "PT-S1"             : "#7f7f7f",
    "PT-S2"             : "#bcbd22",
    "PT-S3"             : "#17becf",
    "Pod"               : "#393b79",
    "Principal_Cell"    : "#ff7f0e",
    "TAL"               : "#8c6d31",
    "Unannotated"       : "#843c39",
    "VSMC_Pericytes"    : "#7b4173",
}


print("Settings for the cell type comparison.", flush = True)
# Folder that contains one sub-folder per sample with MSS_per_spot.csv (output of 00500).
vaf_use = 0.10
alts_threshold = 2
low_coverage_threshold = 10
input_directory = f"MinVAF_{vaf_use:.2f}_MinAltReads_{alts_threshold}_Pruning/"
output_directory = os.path.join(input_directory, "MSS_Analysis",)
mss_file_name = "MSS_per_spot.csv"
os.makedirs(output_directory, exist_ok = True)

# Samples included in this analysis (names as in 00500).
samples_to_analyse = ["57480", "57480_flowgel", "58152", "61968", "63055", "63079", "63409", "63511_reseq", "71116", "72351_reseq", "73266", "76992_2", "78454_flowgel"]
# samples_to_analyse = ["57480", "58152", "61968", "63055", "63079",]# "63511_reseq", "71116", "72351_reseq", "73266", "76992_2",]

# Donor per sample. Samples that are not listed are treated as their own donor.
# Example: if Sample1 and Sample1_flowgel come from the same patient, map both to "Donor1".
sample_to_donor = {
    "57480":         "57480",
    "57480_flowgel": "57480",
    "58152":         "58152",
    "61968":         "61968",
    "63055":         "63055",
    "63079":         "63079",
    "63409":         "63409",
    "63511":         "63511",
    "63511_reseq":   "63511",
    "71116":         "71116_flowgel",
    "72351_reseq":   "72351_flowgel",
    "73266":         "73266_flowgel",
    "76992_2":       "76992_flowgel",
    "78454_flowgel": "78454",
}

# Condition per sample (optional, only for plots and the donor table). Leave empty to skip.
sample_to_condition = {
    "57480":   "Fibrotic",
    "57480_flowgel": "Fibrotic",
    "58152":   "Healthy",
    "61968":   "Healthy",
    "63055":   "Healthy",
    "63079":   "Healthy",
    "63409":   "Healthy",
    "63511":   "Healthy",
    "71116":   "Fibrotic",
    "72351":   "Healthy",
    "73266":   "Healthy",
    "76992_2": "Healthy",
    #"78454":   "Fibrotic",
    "63511_reseq": "Healthy",
    "72351_reseq": "Healthy",
    "58393": "Healthy",
    "63212": "Healthy",
    "72330_2": "Healthy",
    "76997_2": "Fibrotic",
    "62791_2": "Fibrotic",
    "78454_flowgel": "Fibrotic",}

column_mss       = "MSS"
column_depth     = "mt_depth"
column_sample    = "Sample"
column_donor     = "donor"
column_cell_type = "Annotation"
column_condition = "condition"


print("Loading the MSS values of all samples.", flush = True)
mss_tables_per_sample = []
for sample_name in samples_to_analyse:
    mss_file_path = os.path.join(input_directory, sample_name, mss_file_name)
    if not os.path.isfile(mss_file_path):
        raise FileNotFoundError(f"Missing MSS file for {sample_name}: {mss_file_path}")
    # end if statement
    
    mss_table = pd.read_csv(mss_file_path, index_col = 0)
    if column_sample not in mss_table.columns:
        raise KeyError(f"{mss_file_path} has no column {column_sample!r}. Rerun 00500 so it saves the sample.")
    # end if statement
    if (mss_table[column_sample].astype(str) != sample_name).any():
        raise ValueError(f"{mss_file_path} contains spots of another sample than {sample_name}.")
    # end if statement
    
    mss_table[column_donor] = sample_to_donor.get(sample_name, sample_name)
    if sample_to_condition:
        mss_table[column_condition] = sample_to_condition.get(sample_name, np.nan)
    # end if statement
    print(f"{sample_name}: {len(mss_table)} spots, donor {mss_table[column_donor].iloc[0]}", flush = True)
    mss_tables_per_sample.append(mss_table)
# end for loop

combined = pd.concat(mss_tables_per_sample, axis = 0, verify_integrity = True)
combined["log1p_mt_depth"] = np.log1p(combined["mt_depth"])


print("Loading the trajectory.", flush=True)
meta_data_mega = pd.read_csv("/beegfs/data/CostaLab/Spatial_ATAC_Analysis_KuppeLab/results/analysis/ATAC/tables/20mm/ArchR/with_peak_ref/kidney7_scMEGA_metadata_annotation09.csv", index_col=0)
# Replace the # with an _.
meta_data_mega.index = meta_data_mega.index.str.replace("#", "_")
df_merged = pd.merge(meta_data_mega, meta_data_mega, left_index=True, right_index=True, how="outer")

df_merged = combined[["Sample", "mt_depth", "log1p_mt_depth", "MSS", "Annotation"]].join(meta_data_mega[["Clusters2", "Trajectory"]], how="outer")
df_merged = df_merged.dropna()


print("Showing the correlation between depth and MSS before filtering.", flush=True)
# spearmanr(df_merged["MSS"], df_merged["mt_depth"])
x = df_merged.mt_depth.values
y = df_merged.MSS.values
s_lin, i_lin = np.polyfit(x, y, 1)

celltypes = df_merged["Clusters2"].value_counts().index
fig, ax = plt.subplots(1, 1, figsize = (5,4))
for celltype_use in celltypes:
    spots_use = df_merged["Clusters2"].to_numpy() == celltype_use
    color = color_map.get(celltype_use, "#999999")
    ax.scatter(x[spots_use], y[spots_use], s=4, alpha=0.4, color=color, label=celltype_use)

xs = np.linspace(x.min(), x.max(), 100)
ax.plot(xs, s_lin*xs + i_lin, "r-", lw = 1)
ax.set_xlabel("mt_depth")
ax.set_ylabel("MSS")
ax.set_title(f"Raw depth (slope={s_lin:.3f})")
legend_use = ax.legend(title="Cell Type", bbox_to_anchor=(1.01, 1), loc="upper left", markerscale=2,
                       frameon=False, title_fontproperties={"weight": "bold", "size": 10}, fontsize=8)
for handle in legend_use.legend_handles:
    handle.set_alpha(1)

fig.tight_layout()
fig.savefig(f"{output_directory}/Scatter_Depth_vs_MSS_raw.png", dpi = 150)
plt.close(fig)


print("Removing spots with less than 20 reads.", flush=True)
df_merged = df_merged[(df_merged.mt_depth >= 20)]


print("Showing the correlation between depth and MSS after filtering.", flush=True)
# spearmanr(df_merged["MSS"], df_merged["mt_depth"])
x = df_merged.mt_depth.values
y = df_merged.MSS.values
s_lin, i_lin = np.polyfit(x, y, 1)

fig, ax = plt.subplots(1, 1, figsize = (5,4))
for celltype_use in celltypes:
    spots_use = df_merged["Clusters2"].to_numpy() == celltype_use
    color = color_map.get(celltype_use, "#999999")
    ax.scatter(x[spots_use], y[spots_use], s=4, alpha=0.4, color=color, label=celltype_use)

xs = np.linspace(x.min(), x.max(), 100)
ax.plot(xs, s_lin*xs + i_lin, "r-", lw = 1)
ax.set_xlabel("mt_depth")
ax.set_ylabel("MSS")
ax.set_title(f"Filtered Depth (slope={s_lin:.3f})")
legend_use = ax.legend(title="Cell Type", bbox_to_anchor=(1.01, 1), loc="upper left", markerscale=2,
                       frameon=False, title_fontproperties={"weight": "bold", "size": 10}, fontsize=8)
for handle in legend_use.legend_handles:
    handle.set_alpha(1)

fig.tight_layout()
fig.savefig(f"{output_directory}/Scatter_Depth_vs_MSS.png", dpi = 150)
plt.close(fig)


print("Performing the comparison.", flush=True)
results = []
rho, p_value = spearmanr(df_merged["MSS"], df_merged["Trajectory"])
result_to_add = pd.DataFrame({"Comparison": ["Global"], "Rho": [rho], "P": [p_value]}, index=["Global"])
results.append(result_to_add)


# celltypes = df_merged.Clusters2.unique()
for celltype_use in celltypes:
    print(celltype_use, flush=True)
    data_use = df_merged[(df_merged.Clusters2 == celltype_use)]
    rho, p_value = spearmanr(data_use["MSS"], data_use["Trajectory"])
    result_to_add = pd.DataFrame({"Comparison": [celltype_use], "Rho": [rho], "P": [p_value]}, index=[celltype_use])
    results.append(result_to_add)


# Percentile ranks computed separately within each cell type (average ranks for ties)
df_merged["MSS_percentile_within_celltype"] = (
    df_merged.groupby("Clusters2", observed=True)["MSS"].rank(pct=True)
)
df_merged["Trajectory_percentile_within_celltype"] = (
    df_merged.groupby("Clusters2", observed=True)["Trajectory"].rank(pct=True)
)
# Pearson on within-type percentile ranks = stratified Spearman correlation
stratified_rho, stratified_p_value = pearsonr(
    df_merged["MSS_percentile_within_celltype"],
    df_merged["Trajectory_percentile_within_celltype"],
)
result_to_add = pd.DataFrame(
    {"Comparison": ["Stratified"], "Rho": [stratified_rho], "P": [stratified_p_value]},
    index=["Stratified"],
)
results.append(result_to_add)
results = pd.concat(results, axis=0, verify_integrity=True)

is_celltype_row = results["Comparison"].isin(celltypes)
results["P_adjusted"] = np.nan
results.loc[is_celltype_row, "P_adjusted"] = multipletests(results.loc[is_celltype_row, "P"], method="fdr_by")[1]
results.to_csv(f"{output_directory}/Correlation_MSS_Trajectory.csv")


print("Showing the correlation between Trajectory and MSS.", flush=True)
x = df_merged.Trajectory.values
y = df_merged.MSS.values
s_lin, i_lin = np.polyfit(x, y, 1)

fig, ax = plt.subplots(1, 1, figsize = (5,4))
for celltype_use in celltypes:
    spots_use = df_merged["Clusters2"].to_numpy() == celltype_use
    color = color_map.get(celltype_use, "#999999")
    ax.scatter(x[spots_use], y[spots_use], s=4, alpha=0.4, color=color, label=celltype_use)

xs = np.linspace(x.min(), x.max(), 100)
ax.plot(xs, s_lin*xs + i_lin, "r-", lw = 1)
ax.set_xlabel("Trajectory")
ax.set_ylabel("MSS")
ax.set_title(f"linear (slope={s_lin:.3f})")
legend_use = ax.legend(title="Cell Type", bbox_to_anchor=(1.01, 1), loc="upper left", markerscale=2,
                       frameon=False, title_fontproperties={"weight": "bold", "size": 10}, fontsize=8)
for handle in legend_use.legend_handles:
    handle.set_alpha(1)

fig.tight_layout()
fig.savefig(f"{output_directory}/Scatter_Trajectory_vs_MSS.png", dpi = 150)
plt.close(fig)


print("Showing the correlation between Trajectory and MSS as a facet plot.", flush=True)
x = df_merged.Trajectory.values
y = df_merged.MSS.values

celltypes_number = len(celltypes)
ncols = int(np.ceil(np.sqrt(celltypes_number)))
nrows = int(np.ceil(celltypes_number / ncols))

fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols, 3 * nrows), sharex=True, sharey=True, squeeze=False)
axes = axes.ravel()
xs = np.linspace(x.min(), x.max(), 100)

handles = []
for ax, celltype_use in zip(axes, celltypes):
    ax.set_box_aspect(1)
    spots_use = df_merged["Clusters2"].to_numpy() == celltype_use
    color = color_map.get(celltype_use, "#999999")
    
    h = ax.scatter(x[spots_use], y[spots_use], s=4, alpha=0.4, color=color, label=celltype_use, rasterized=True)
    handles.append(h)
    
    # Per-cell-type linear fit (needs at least two distinct x values)
    s_lin, i_lin = np.polyfit(x[spots_use], y[spots_use], 1)
    ax.plot(xs, s_lin * xs + i_lin, "r-", lw=1)
    ax.set_title(f"{celltype_use}\nslope={s_lin:.3f}", fontsize=9)

# Remove empty panels in the last row.
for ax in axes[celltypes_number:]:
    ax.set_visible(False)

fig.supxlabel("Trajectory")
fig.supylabel("MSS")

leg = fig.legend(handles=handles, title="Cell Type", loc="center left",
                 bbox_to_anchor=(1.0, 0.5), markerscale=2, frameon=False,
                 title_fontproperties={"weight": "bold", "size": 10}, fontsize=8)
for h in leg.legend_handles:
    h.set_alpha(1)

fig.tight_layout()
fig.savefig(f"{output_directory}/Scatter_Trajectory_vs_MSS_facet.png", dpi=150, bbox_inches="tight")
plt.close(fig)


print("Showing the correlation between depth and MSS after filtering as a facet plot.", flush=True)
x = df_merged.mt_depth.values
y = df_merged.MSS.values

celltypes_number = len(celltypes)
ncols = int(np.ceil(np.sqrt(celltypes_number)))
nrows = int(np.ceil(celltypes_number / ncols))

fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols, 3 * nrows), sharex=True, sharey=True, squeeze=False)
axes = axes.ravel()
xs = np.linspace(x.min(), x.max(), 100)

handles = []
for ax, celltype_use in zip(axes, celltypes):
    ax.set_box_aspect(1)
    spots_use = df_merged["Clusters2"].to_numpy() == celltype_use
    color = color_map.get(celltype_use, "#999999")
    
    h = ax.scatter(x[spots_use], y[spots_use], s=4, alpha=0.4, color=color, label=celltype_use, rasterized=True)
    handles.append(h)
    
    # Per-cell-type linear fit (needs at least two distinct x values)
    s_lin, i_lin = np.polyfit(x[spots_use], y[spots_use], 1)
    ax.plot(xs, s_lin * xs + i_lin, "r-", lw=1)
    ax.set_title(f"{celltype_use}\nslope={s_lin:.3f}", fontsize=9)

# Remove empty panels in the last row.
for ax in axes[celltypes_number:]:
    ax.set_visible(False)

fig.supxlabel("mt_depth")
fig.supylabel("MSS")

leg = fig.legend(handles=handles, title="Cell Type", loc="center left",
                 bbox_to_anchor=(1.0, 0.5), markerscale=2, frameon=False,
                 title_fontproperties={"weight": "bold", "size": 10}, fontsize=8)
for h in leg.legend_handles:
    h.set_alpha(1)

fig.tight_layout()
fig.savefig(f"{output_directory}/Scatter_Depth_vs_MSS_facet.png", dpi=150,
            bbox_inches="tight")
plt.close(fig)


test1 = df_merged[df_merged.MSS > 0.75]
test2 = df_merged[df_merged.MSS < 0.75]
