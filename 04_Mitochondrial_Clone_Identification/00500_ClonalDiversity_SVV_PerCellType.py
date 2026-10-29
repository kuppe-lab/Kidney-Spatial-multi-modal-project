#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 27.10.2026

@author: grasshoff

Calculate the Effective Number of Species per cell type and sample to determine the 
diversity driven by the clones.
Then, compare the clonal diversity per cell type.
We consider all cell types as cell types of interest.
For each cell type of interest the differences in ENS between it and the median ENS of
all other cell types are calculated. This is done per sample, which allows us to do 
this in paired manner.
These differences are then compared using a paried Wilcoxon test.

"""

print("Packages.", flush=True)
import numpy as np
import pandas as pd
import scanpy as sc
import os
import matplotlib.pyplot as plt
from scipy.stats import false_discovery_control, wilcoxon

sc.set_figure_params(dpi_save=300, figsize=[10, 10])


def significance_stars(p):
    if np.isnan(p):
        return ""
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    
    return ""


def shannon_entropy(counts):
    """Shannon entropy in nats (natural log)."""
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()
    if total <= 0:
        return 0.0
    
    p = counts / total
    p = p[p > 0]
    if len(p) == 0:
        return 0.0
    
    return -float(np.sum(p * np.log(p)))


def benjamini_hochberg(p_values):
    """
    Benjamini-Hochberg false discovery rate (FDR) correction.
    NaNs are skipped and not used in the correction as an actual test.
    This function returns the adjusted p-values aligned with the input.
    NaNs are kept in the output. So, original NaNs are identified.
    """
    p_values = np.asarray(p_values, dtype = float)
    adjusted_p_values = np.full(len(p_values), np.nan)
    is_finite_p_value = ~np.isnan(p_values)
    if is_finite_p_value.sum() == 0:
        return adjusted_p_values
    
    adjusted_p_values[is_finite_p_value] = false_discovery_control(p_values[is_finite_p_value], method = "bh")
    return adjusted_p_values


def paired_wilcoxon_safe(deltas):
    """
    Two-sided Wilcoxon signed-rank on a 1-D array of differences.
    Returns (W, p, n_pairs, median_delta).
    """
    deltas_use = np.asarray(deltas, dtype=float)
    deltas_use = deltas_use[~np.isnan(deltas_use)]
    n = len(deltas_use)
    if n < 2:
        return (np.nan, np.nan, n, float(np.median(deltas_use)) if n else np.nan)
    
    nonzero = deltas_use[deltas_use != 0.0]
    if len(nonzero) < 2:
        return (np.nan, np.nan, n, float(np.median(deltas_use)))
    
    try:
        result = wilcoxon(deltas_use, alternative="two-sided", zero_method="wilcox", correction=False)
        w_effect_size = float(result.statistic)
        pvalue = float(result.pvalue)
    except ValueError:
        w_effect_size = np.nan
        pvalue = np.nan
    
    return (w_effect_size, pvalue, n, float(np.median(deltas_use)))


clone_definitions = {
    "57480": {
        "Clone_01": ["chrM_8614_T>C", "chrM_11435_G>A"],
        "Clone_02": ["chrM_8440_A>G", "chrM_12776_T>C"],
        "Clone_03": ["chrM_10640_T>C"],
        "Clone_04": ["chrM_14438_G>A"],
        "Clone_05": ["chrM_14123_T>C"],
        "Clone_06": ["chrM_3315_G>A"],
        "Clone_07": ["chrM_9100_A>C"],
        "Clone_08": ["chrM_13928_G>A"],
    },
    "63055": {
        "Clone_09": ["chrM_709_G>A", "chrM_16192_C>T", "chrM_16306_C>T"],
        "Clone_10": ["chrM_51_T>C", "chrM_2395_A>G"],
        "Clone_11": ["chrM_1051_A>C"],
        "Clone_12": ["chrM_3715_G>A"],
        "Clone_13": ["chrM_2118_T>C"],
        },
    "58152": {
        "Clone_14": ["chrM_16266_C>T"],
        "Clone_15": ["chrM_13691_A>G"],
        "Clone_16": ["chrM_16180_A>G"],
        },
    "61968": {
        "Clone_17": ["chrM_14180_T>C"],
        "Clone_18": ["chrM_15192_T>C"],
        "Clone_19": ["chrM_15132_T>C"],
        "Clone_20": ["chrM_12084_C>T"],
        },
    "63079": {
        "Clone_21": ["chrM_556_A>G"],
        "Clone_22": ["chrM_6505_T>C"],
        "Clone_23": ["chrM_13552_G>A"],
        },
    "63409": {
        "Clone_24": ["chrM_556_A>G"],
        "Clone_25": ["chrM_13198_G>A"],
        },
    "63511_reseq": {
        "Clone_26": ["chrM_16311_T>C"],
        "Clone_27": ["chrM_10775_G>A"],
        },
    "71116": {
        "Clone_28": ["chrM_531_T>C"],
        "Clone_29": ["chrM_2635_G>C"],
        },
    "72351_reseq": {
        "Clone_30": ["chrM_14470_T>C"],
        "Clone_31": ["chrM_16007_A>G"],
        },
    "73266": {
        "Clone_32": ["chrM_15432_C>T", "chrM_16188_C>T"],
        "Clone_33": ["chrM_3424_G>A"],
        "Clone_34": ["chrM_15612_G>A"],
        },
    "76992_2": {
        "Clone_35": ["chrM_422_T>C", "chrM_16349_A>G"],
        "Clone_36": ["chrM_15746_A>G"],
        },
#     "63212": {
#         "Clone_39": ["chrM_981_C>T"],
#         "Clone_40": ["chrM_12902_T>C"],
#         "Clone_41": ["chrM_15309_T>C"],
#         "Clone_42": ["chrM_11506_T>C"],
#         },
#     "72330_2": {
#         "Clone_43": ["chrM_3608_G>A"],
#         "Clone_44": ["chrM_64_C>A"],
#         },
    "57480_flowgel": {
        "Clone_45": ["chrM_2963_A>C"],
        "Clone_46": ["chrM_7869_T>C"],
    },
    "78454_flowgel": {
        "Clone_47": ["chrM_13957_G>A"],
        "Clone_48": ["chrM_15734_G>A"],
    },
}


print("Variables.", flush=True)
jitter_seed = 0 # We want to replicate the plots exactly.
alts_threshold = 2
vaf_use = 0.10
min_clonal_spots = 1
design_matrix_path = "DesignMatrix_VisiumHD.csv"
design_matrix = pd.read_csv(design_matrix_path, index_col=0)
order_by_median = True
samples = ["57480", "58152", "61968", "63055", "63079", "63409", "63511_reseq", "71116", "72351_reseq", "73266", "76992_2", "57480_flowgel", "78454_flowgel"]
fig_height = 5
dpi = 300
condition_colors = {"Non-Fibrotic": "#1f77b4", "Fibrotic": "#ff7f0e",}


input = f"MinVAF_{vaf_use:.2f}_MinAltReads_{alts_threshold}_Pruning/"
output = f"MinVAF_{vaf_use:.2f}_MinAltReads_{alts_threshold}_Pruning/ClonalDiversity_SVV_PerCellType/"
os.makedirs(output, exist_ok=True)


print("Validating the clone definitions.", flush=True)
# "WT" (wild type) labels spots without any tracked variant and "+" joins the
# names of several clones found in the same spot. So, neither may be used in a
# clone name. We check this here and produce a warning.
for sample_use, clones_of_sample in clone_definitions.items():
    seen_svv_sets = {}
    for clone_name, clone_svvs in clones_of_sample.items():
        if clone_name == "WT" or "+" in clone_name:
            raise ValueError(f"Clone {clone_name!r} in {sample_use!r}: a clone name may not be W' or contain +.")
        
        if len(clone_svvs) == 0:
            raise ValueError(f"Clone {clone_name!r} in {sample_use!r} has no spatially variable variants.")
        
        svv_set = frozenset(clone_svvs)
        if svv_set in seen_svv_sets:
            raise ValueError(f"Clones {seen_svv_sets[svv_set]!r} and {clone_name!r} in {sample_use!r} have the same SVVs.")
        
        seen_svv_sets[svv_set] = clone_name

ignored_samples = sorted(set(clone_definitions) - set(samples))
if ignored_samples:
    print(f"Clones of samples not in samples are ignored: {ignored_samples}", flush=True)


print("Calculating the Effective Number of Species per sample and cell type.", flush=True)
# Carrier spots are spots with at least one tracked SVV.
records = []
for sample_use in samples:
    print(f"Sample: {sample_use}", flush=True)
    if sample_use not in clone_definitions:
        print(f"No clone definitions for {sample_use}; skipping.", flush=True)
        continue
    
    adata = sc.read_h5ad(f"{input}{sample_use}/adata_pruned.h5ad")
    
    # The spot annotation was added in 00200.
    number_of_annotated_spots = int(adata.obs["Annotation"].notna().sum())
    print(f"Annotated spots: {number_of_annotated_spots} of {adata.n_obs}", flush=True)
    if number_of_annotated_spots == 0:
        print(f"WARNING: no spot of {sample_use} has an annotation in the AnnData object.", flush=True)
        continue
    
    adata = adata[adata.obs["Annotation"].notna()].copy()
    
    # A spot carries a clone if at least one of the clone's SVVs is detected.
    clones_of_sample   = clone_definitions[sample_use]
    clone_names        = sorted(clones_of_sample)
    spot_carries_clone = np.zeros((adata.n_obs, len(clone_names)), dtype=bool)
    for clone_position, clone_name in enumerate(clone_names):
        measured_svvs = [svv for svv in clones_of_sample[clone_name] if svv in adata.var_names]
        missing_svvs  = [svv for svv in clones_of_sample[clone_name] if svv not in adata.var_names]
        if missing_svvs:
            print(f"WARNING: SVVs of {clone_name} missing from {sample_use}: {missing_svvs}", flush=True)
        
        if len(measured_svvs) == 0:
            continue
        
        svv_values = adata[:, measured_svvs].X
        svv_values = svv_values.toarray() if hasattr(svv_values, "toarray") else svv_values
        spot_carries_clone[:, clone_position] = (svv_values > 0).any(axis=1)
    
    # Spots carrying several clones get a combined label, e.g. "Clone_A+Clone_B".
    # Spots without any clone are labelled "WT".
    clone_label_per_spot = []
    for carried_clones in spot_carries_clone:
        carried_clone_names = [clone_name for clone_name, is_carried in zip(clone_names, carried_clones) if is_carried]
        clone_label_per_spot.append("+".join(carried_clone_names) if carried_clone_names else "WT")
    
    clone_label_per_spot = pd.Series(clone_label_per_spot, index=adata.obs_names)
    is_carrier_spot      = clone_label_per_spot != "WT"
    annotation_of_spot   = adata.obs["Annotation"].astype(str)
    
    for cell_type in sorted(annotation_of_spot.unique()):
        is_celltype_spot = annotation_of_spot == cell_type
        carrier_labels   = clone_label_per_spot[is_celltype_spot & is_carrier_spot]
        if len(carrier_labels) > 0:
            ens_value = float(np.exp(shannon_entropy(carrier_labels.value_counts().values)))
        else:
            ens_value = np.nan
        
        records.append({
            "sample":                     sample_use,
            "group":                      design_matrix.loc[sample_use, "condition"],
            "cell_type":                  cell_type,
            "n_spots_celltype":           int(is_celltype_spot.sum()),
            "n_clonal_spots":             len(carrier_labels),
            "n_spots_composite_celltype": int(carrier_labels.str.contains("+", regex=False).sum()),
            "n_clones":                   carrier_labels.nunique(),
            "ENS":                        ens_value,
        })


print("Thresholding.", flush=True)
records_df = pd.DataFrame(records)
if records_df.empty:
    raise ValueError("No records. Check the clone definitions and the annotation transfer.")

# A (sample x cell type) entry is only used if it has at least min_svv_spots carrier spots.
records_df["passed_threshold"] = records_df["n_clonal_spots"] >= min_clonal_spots
kept_records = records_df[records_df["passed_threshold"]]
if kept_records.empty:
    raise ValueError(f"No (sample x cell type) entry has at least {min_clonal_spots} carrier spots.")


print("Paired Wilcoxon test: cell type vs median of the other cell types.", flush=True)
# Rows are samples (patients), columns are cell types, values are ENS.
ens_wide     = kept_records.pivot(index="sample", columns="cell_type", values="ENS")
delta_tables = []
stats_rows   = []
for cell_type in sorted(ens_wide.columns):
    paired_values = pd.DataFrame({
        "ENS":               ens_wide[cell_type],
        "ENS_others_median": ens_wide.drop(columns=[cell_type]).median(axis=1, skipna=True),
    }).dropna()
    paired_values["ENS_delta"] = paired_values["ENS"] - paired_values["ENS_others_median"]
    paired_values["cell_type"] = cell_type
    delta_tables.append(paired_values.reset_index())
    
    w_statistic, p_value, number_of_pairs, median_delta = paired_wilcoxon_safe(paired_values["ENS_delta"].values)
    stats_rows.append({
        "cell_type":       cell_type,
        "n_pairs":         number_of_pairs,
        "median_interest": paired_values["ENS"].median(),
        "median_rest":     paired_values["ENS_others_median"].median(),
        "median_delta":    median_delta,
        "W":               w_statistic,
        "P":               p_value,
    })

stats_df = pd.DataFrame(stats_rows)
stats_df["P_FDR"] = benjamini_hochberg(stats_df["P"].values)
stats_df["sig"]   = [significance_stars(p_adjusted) for p_adjusted in stats_df["P_FDR"]]
stats_df.to_csv(f"{output}diversity_celltype_stats.csv", index=False)

# The per-patient median of the other cell types and the delta are added to the
# values table. They stay NaN for entries that did not pass the threshold.
delta_df   = pd.concat(delta_tables)[["sample", "cell_type", "ENS_others_median", "ENS_delta"]]
records_df = records_df.merge(delta_df, on=["sample", "cell_type"], how="left")
records_df.to_csv(f"{output}diversity_per_celltype_values.csv", index=False)


print("Plotting.", flush=True)
random_generator  = np.random.default_rng(jitter_seed)
plot_records      = records_df[records_df["passed_threshold"]]
stars_of_celltype = stats_df.set_index("cell_type")["sig"]
plot_settings = [
    # (plotted column, y-axis label, title, output file)
    ("ENS",
     "ENS",
     "Effective Number of Species per Cell Type",
     f"{output}celltype_diversity_ENS"),
    ("ENS_delta",
     "Delta ENS",
     "Differences in Effective Number of Species per Cell Type",
     f"{output}celltype_diversity_delta_ENS"),
]

for value_column, y_axis_label, plot_title, plot_path in plot_settings:
    cell_types_sorted = sorted(plot_records["cell_type"].unique())
    if order_by_median:
        median_per_celltype = plot_records.groupby("cell_type")[value_column].median()
        cell_types_sorted   = median_per_celltype.sort_values(ascending=False).index.tolist()
    
    values_per_celltype = [plot_records.loc[plot_records["cell_type"] == cell_type, value_column].dropna().values
                           for cell_type in cell_types_sorted]
    all_values = np.concatenate(values_per_celltype)
    if len(all_values) == 0:
        print(f"No values to plot for {value_column}.", flush=True)
        continue
        
    figure, axis = plt.subplots(figsize=(max(1.2 * len(cell_types_sorted) + 2, 8), fig_height))
    
    if value_column == "ENS_delta":
        # Delta = 0 is the null hypothesis of the Wilcoxon test.
        # We assume there is no difference, so the deltas should all be 0.
        axis.axhline(0.0, color="grey", linestyle="--", linewidth=1, zorder=1)
    
    # Boxes only for cell types with at least two patients.
    box_positions = [position for position, values in enumerate(values_per_celltype) if len(values) >= 2]
    if box_positions:
        boxes = axis.boxplot([values_per_celltype[position] for position in box_positions],
                             positions    = box_positions,
                             widths       = 0.5,
                             patch_artist = True,
                             showfliers   = False,
                             medianprops  = dict(color="black"),
                             whiskerprops = dict(color="black"),
                             capprops     = dict(color="black"))
        for box in boxes["boxes"]:
            box.set_facecolor("lightgrey")
            box.set_alpha(0.4)
            box.set_edgecolor("black")
    
    value_range = max(all_values.max() - all_values.min(), 1e-9)
    for position, cell_type in enumerate(cell_types_sorted):
        celltype_records = plot_records[plot_records["cell_type"] == cell_type]
        for group_name, group_color in condition_colors.items():
            group_values = celltype_records.loc[celltype_records["group"] == group_name, value_column].dropna().values
            x_positions  = random_generator.normal(loc=position, scale=0.06, size=len(group_values))
            axis.plot(x_positions, group_values, "o", color=group_color, markeredgecolor="black", markersize=7, zorder=3)
        
        # The stars come from the same Wilcoxon test in both plots.
        star = stars_of_celltype.get(cell_type, "")
        if star and len(values_per_celltype[position]) > 0:
            axis.text(position, values_per_celltype[position].max() + 0.04 * value_range, star,
                      ha="center", va="bottom", fontsize=12, fontweight="bold")
    
    axis.set_xticks(range(len(cell_types_sorted)))
    # axis.set_xticklabels([f"{cell_type}\n(n={len(values)})" for cell_type, values in zip(cell_types_sorted, values_per_celltype)],
    #                      fontsize=8, rotation=30, ha="right")
    axis.set_xticklabels(cell_types_sorted,
                         fontsize=14, fontweight="bold", rotation=45, ha="right")
    axis.set_xlim(-0.5, len(cell_types_sorted) - 0.5)
    axis.set_ylabel(y_axis_label, fontsize=14, fontweight="bold")
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    legend_handles = [plt.Line2D([0], [0], marker="o", linestyle="", color=group_color, markeredgecolor="black", markersize=12, label=group_name)
                      for group_name, group_color in condition_colors.items()]
    axis.legend(handles=legend_handles, loc="best", frameon=False, fontsize=12)
    axis.set_title(plot_title, fontsize=16, fontweight="bold")
    plt.tight_layout()
    plt.savefig(f"{plot_path}.png", dpi=dpi, bbox_inches="tight")
    plt.savefig(f"{plot_path}.pdf", dpi=dpi, bbox_inches="tight")
    plt.close()
