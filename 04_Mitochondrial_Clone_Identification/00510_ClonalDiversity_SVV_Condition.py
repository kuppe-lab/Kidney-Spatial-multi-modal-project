#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Martin Grasshoff
Environment: spatial_sigurd_py_env

Calculate the Shannon entropy and the Effective Number of Species (ENS) per sample
from the sample-wise clone definitions of SVVs.
A spot is part of a clone if at least one of the SVVs of a clone is present.
In case a spot is part of multiple clones, the clones are combined (Clone1+Clone2).
This should not happen, but is added as a safety measure.
Spots that are not part of any clone are set to by wildtype (WT) and not used.

Then, the two conditions are compared using the two-sided Mann-Whitney U test.

"""

print("Packages.", flush=True)
import numpy as np
import pandas as pd
import scanpy as sc
import os
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import warnings

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


clone_definitions = {
    "57480": {
        "8614_T_C_11435_G_A": ["chrM_8614_T>C", "chrM_11435_G>A"],
        "8440_A_G_12776_T_C": ["chrM_8440_A>G", "chrM_12776_T>C"],
        "10640_T_C": ["chrM_10640_T>C"],
        "14438_G_A": ["chrM_14438_G>A"],
        "14123_T_C": ["chrM_14123_T>C"],
        "3315_G_A": ["chrM_3315_G>A"],
        "9100_A_C": ["chrM_9100_A>C"],
        "13928_G_A": ["chrM_13928_G>A"],
    },
    "57480_flowgel": {
        "7869_T_C": ["chrM_7869_T>C"],
        "13369_T_C_15758_A_G": ["chrM_13369_T>C", "chrM_15758_A>G"],
    },
    "58152": {
        "16266_C_T": ["chrM_16266_C>T"],
        "13691_A_G": ["chrM_13691_A>G"],
        "16180_A_G": ["chrM_16180_A>G"],
    },
    "61968": {
        "14180_T_C": ["chrM_14180_T>C"],
        "15192_T_C": ["chrM_15192_T>C"],
        "15132_T_C": ["chrM_15132_T>C"],
        "12084_C_T": ["chrM_12084_C>T"],
    },
    "63055": {
        "709_G_A_16192_C_T_16306_C_T": ["chrM_709_G>A", "chrM_16192_C>T", "chrM_16306_C>T"],
        "51_T_C_2395_A_G": ["chrM_51_T>C", "chrM_2395_A>G"],
        "1051_A_C": ["chrM_1051_A>C"],
        "3715_G_A": ["chrM_3715_G>A"],
        "2118_T_C": ["chrM_2118_T>C"],
    },
    "63079": {
        "556_A_G": ["chrM_556_A>G"],
        "6505_T_C": ["chrM_6505_T>C"],
        "13552_G_A": ["chrM_13552_G>A"],
    },
    "63409": {
        "556_A_G": ["chrM_556_A>G"],
        "13198_G_A": ["chrM_13198_G>A"],
    },
    "63511_reseq": {
        "16311_T_C": ["chrM_16311_T>C"],
        "10775_G_A": ["chrM_10775_G>A"],
    },
    "71116": {
        "531_T_C": ["chrM_531_T>C"],
        "2635_G_C": ["chrM_2635_G>C"],
    },
    "72351_reseq": {
        "14470_T_C": ["chrM_14470_T>C"],
        "16007_A_G": ["chrM_16007_A>G"],
    },
    "73266": {
        "15612_G_A": ["chrM_15612_G>A"],
        "3424_G_A": ["chrM_3424_G>A"],
        "4849_G_A_16188_C_T": ["chrM_4849_G>A", "chrM_16188_C>T"],
    },
    "76992_2": {
        "3337_G_A": ["chrM_3337_G>A"],
        "15746_A_G": ["chrM_15746_A>G"],
    },
        "78454": {
        "63_T_G_65_T_G_16362_T_C": ["chrM_63_T>G", "chrM_65_T>G", "chrM_16362_T>C"],
        "564_G_A": ["chrM_564_G>A"],
    },
    "78454_flowgel": {
        "15734_G_A": ["chrM_15734_G>A"],
        "556_A_G": ["chrM_556_A>G"],
    },
}


print("Variables.", flush=True)
jitter_seed = 0 # We want to replicate the plots exactly.
alts_threshold = 2
vaf_use = 0.10
design_matrix_path = "DesignMatrix_VisiumHD.csv"
design_matrix = pd.read_csv(design_matrix_path, index_col=0)
samples = ["57480", "57480_flowgel", "58152", "61968", "63055", "63079", "63409", "63511_reseq", "71116", "72351_reseq", "73266", "76992_2", "78454_flowgel"]
metrics          = ["shannon", "ENS"]
groups_compared  = ["Fibrotic", "Non-Fibrotic"]
fig_height       = 5
dpi              = 300
condition_colors = {"Non-Fibrotic": "#1f77b4", "Fibrotic": "#ff7f0e",}
main_title       = "SVV diversity: Fibrotic vs Non-Fibrotic"


input = f"MinVAF_{vaf_use:.2f}_MinAltReads_{alts_threshold}_Pruning/"
output = f"MinVAF_{vaf_use:.2f}_MinAltReads_{alts_threshold}_Pruning/Diversity_SVV_ManualCombinations/"
output_prefix = f"{output}diversity_comparison_svv_manualcombinations"
os.makedirs(output, exist_ok=True)


print("Validating the clone definitions.", flush=True)
# "WT" labels spots without any tracked variant and "+" joins the names of several
# clones found in the same spot. So, neither may be used in a clone name.
for sample_use, clones_of_sample in clone_definitions.items():
    seen_svv_sets = {}
    for clone_name, clone_svvs in clones_of_sample.items():
        if clone_name == "WT" or "+" in clone_name:
            raise ValueError(f"Clone {clone_name!r} in {sample_use!r}: a clone name may not be WT or contain +.")
        
        if len(clone_svvs) == 0:
            raise ValueError(f"Clone {clone_name!r} in {sample_use!r} has no spatially variable variants.")
        
        svv_set = frozenset(clone_svvs)
        if svv_set in seen_svv_sets:
            raise ValueError(f"Clones {seen_svv_sets[svv_set]!r} and {clone_name!r} in {sample_use!r} have the same SVVs.")
        
        seen_svv_sets[svv_set] = clone_name

ignored_samples = sorted(set(clone_definitions) - set(samples))
if ignored_samples:
    print(f"Clones of samples not in samples are ignored: {ignored_samples}", flush=True)


print("Calculating the diversity per sample.", flush=True)
# Carrier spots are spots with at least one tracked SVV.
records           = []
per_clone_records = []
for sample_use in samples:
    print(f"Sample: {sample_use}", flush=True)
    if sample_use not in clone_definitions:
        print(f"No clone definitions for {sample_use}; skipping.", flush=True)
        continue
    
    adata = sc.read_h5ad(f"{input}{sample_use}/adata_pruned.h5ad")
    group_of_sample = design_matrix.loc[sample_use, "condition"]
    
    # A spot carries a clone if at least one of the clone's SVVs is detected.
    clones_of_sample   = clone_definitions[sample_use]
    clone_names        = sorted(clones_of_sample)
    spot_carries_clone = np.zeros((adata.n_obs, len(clone_names)), dtype=bool)
    for clone_position, clone_name in enumerate(clone_names):
        measured_svvs = [svv for svv in clones_of_sample[clone_name] if svv in adata.var_names]
        missing_svvs  = [svv for svv in clones_of_sample[clone_name] if svv not in adata.var_names]
        if missing_svvs:
            warnings.warn(f"SVVs of {clone_name} missing from {sample_use}: {missing_svvs}")
        
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
    carrier_labels       = clone_label_per_spot[clone_label_per_spot != "WT"]
    
    # Solo spots carry only this clone, composite spots carry this clone and at least one other.
    number_of_clones_per_spot = spot_carries_clone.sum(axis=1)
    for clone_position, clone_name in enumerate(clone_names):
        n_spots_any  = int(spot_carries_clone[:, clone_position].sum())
        n_spots_solo = int((spot_carries_clone[:, clone_position] & (number_of_clones_per_spot == 1)).sum())
        per_clone_records.append({
            "sample":            sample_use,
            "group":             group_of_sample,
            "clone":             clone_name,
            "svvs":              "+".join(sorted(clones_of_sample[clone_name])),
            "n_spots_solo":      n_spots_solo,
            "n_spots_any":       n_spots_any,
            "n_spots_composite": n_spots_any - n_spots_solo,
        })
    
    if len(carrier_labels) > 0:
        shannon_value = shannon_entropy(carrier_labels.value_counts().values)
        ens_value     = float(np.exp(shannon_value))
    else:
        warnings.warn(f"No carrier spots in {sample_use}.")
        shannon_value = np.nan
        ens_value     = np.nan
    
    records.append({
        "sample":            sample_use,
        "group":             group_of_sample,
        "n_spots":           adata.n_obs,
        "n_spots_wt":        adata.n_obs - len(carrier_labels),
        "n_clonal_spots":    len(carrier_labels),
        "n_spots_composite": int(carrier_labels.str.contains("+", regex=False).sum()),
        "n_clones":          carrier_labels.nunique(),
        "shannon":           shannon_value,
        "ENS":               ens_value,
    })

records_df = pd.DataFrame(records)
if records_df.empty:
    raise ValueError("No records. Check the clone definitions.")

records_df.to_csv(f"{output_prefix}_values.csv", index=False)
pd.DataFrame(per_clone_records).to_csv(f"{output_prefix}_per_clone.csv", index=False)


print("Mann-Whitney U test.", flush=True)
# Samples without carrier spots (NaN) are not used in the test.
stats_rows = []
for metric in metrics:
    values_of_group = [records_df.loc[records_df["group"] == group_name, metric].dropna().values
                       for group_name in groups_compared]
    if min(len(group_values) for group_values in values_of_group) < 1:
        print(f"Not enough samples for {metric}; skipping the test.", flush=True)
        u_statistic, p_value, rank_biserial = np.nan, np.nan, np.nan
    else:
        u_statistic, p_value = mannwhitneyu(values_of_group[0], values_of_group[1], alternative="two-sided")
        # Rank-biserial correlation as effect size, ranging from -1 to 1.
        rank_biserial = 2 * u_statistic / (len(values_of_group[0]) * len(values_of_group[1])) - 1
    
    stats_rows.append({
        "metric":   metric,
        "group_a":  groups_compared[0],
        "group_b":  groups_compared[1],
        "n_a":      len(values_of_group[0]),
        "n_b":      len(values_of_group[1]),
        "median_a": np.median(values_of_group[0]) if len(values_of_group[0]) else np.nan,
        "median_b": np.median(values_of_group[1]) if len(values_of_group[1]) else np.nan,
        "U":        float(u_statistic),
        "P":        float(p_value),
        "r":        float(rank_biserial),
    })

stats_df = pd.DataFrame(stats_rows)
stats_df["sig"] = [significance_stars(p_value) for p_value in stats_df["P"]]
stats_df.to_csv(f"{output_prefix}_stats.csv", index=False)


print("Plotting.", flush=True)
random_generator = np.random.default_rng(jitter_seed)
stats_of_metric  = stats_df.set_index("metric")
figure, axes = plt.subplots(1, len(metrics), figsize=(4 * len(metrics), fig_height), squeeze=False)
for axis, metric in zip(axes[0], metrics):
    values_of_group = [records_df.loc[records_df["group"] == group_name, metric].dropna().values
                       for group_name in groups_compared]
    if min(len(group_values) for group_values in values_of_group) < 1:
        axis.text(0.5, 0.5, f"Not enough samples\nfor {metric}", ha="center", va="center", transform=axis.transAxes)
        axis.set_xticks([])
        axis.set_yticks([])
        axis.set_title(metric, fontsize=16, fontweight="bold")
        continue
    
    boxes = axis.boxplot(values_of_group,
                         positions    = [0, 1],
                         widths       = 0.5,
                         patch_artist = True,
                         showfliers   = False,
                         medianprops  = dict(color="black"),
                         whiskerprops = dict(color="black"),
                         capprops     = dict(color="black"))
    for box, group_name in zip(boxes["boxes"], groups_compared):
        box.set_facecolor(condition_colors[group_name])
        box.set_alpha(0.4)
        box.set_edgecolor("black")
    
    # Each sample is one point. This matters with few samples.
    for position, (group_name, group_values) in enumerate(zip(groups_compared, values_of_group)):
        x_positions = random_generator.normal(loc=position, scale=0.05, size=len(group_values))
        axis.plot(x_positions, group_values, "o", color=condition_colors[group_name], markeredgecolor="black", markersize=7, zorder=3)
    
    # Significance bracket between the two groups.
    all_values     = np.concatenate(values_of_group)
    value_range    = max(all_values.max() - all_values.min(), 1e-9)
    bracket_height = all_values.max() + 0.08 * value_range
    bracket_tail   = 0.03 * value_range
    axis.plot([0, 0, 1, 1], [bracket_height - bracket_tail, bracket_height, bracket_height, bracket_height - bracket_tail],
              linewidth=1, color="black")
    axis.text(0.5, bracket_height, stats_of_metric.loc[metric, "sig"], ha="center", va="bottom", fontsize=12, fontweight="bold")
    axis.set_ylim(top=bracket_height + 0.15 * value_range)
    
    axis.set_xticks([0, 1])
    axis.set_xticklabels([f"{group_name}\n(n={len(group_values)})" for group_name, group_values in zip(groups_compared, values_of_group)],
                         fontsize=14, fontweight="bold")
    axis.set_ylabel(metric, fontsize=14, fontweight="bold")
    axis.set_title(f"{metric} (p = {stats_of_metric.loc[metric, 'P']:.3f})", fontsize=16, fontweight="bold")
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)

figure.suptitle(main_title, fontsize=14, fontweight="bold")
plt.tight_layout()
figure.subplots_adjust(top=0.85) # So the title does not overlap with the panels.
plt.savefig(f"{output_prefix}.png", dpi=dpi, bbox_inches="tight")
plt.savefig(f"{output_prefix}.pdf", dpi=dpi, bbox_inches="tight")
plt.close()