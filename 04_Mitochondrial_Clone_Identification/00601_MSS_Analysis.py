#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Comparison of the mitochondrial local constraint (MLC) score sum (MSS) between
cell types, corrected for mitochondrial sequencing depth.

Reference cell types reach the minimum number of spots in every sample.
They form the baseline. In each sample, deviations_from_reference_average compares the
depth-adjusted effect of every cell type with the median effect of the reference cell
types. For a reference cell type, the baseline is the median of all other reference cell
types, excluding itself. For a non-reference cell type, it is the median of all reference cell
types.

"""

import os
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests


print("Settings for the cell type comparison.", flush=True)
vaf_use = 0.10
alts_threshold = 2
input_directory = f"/beegfs/data/CostaLab/Spatial_ATAC_Analysis_KuppeLab/results/SpatialClones/MinVAF_{vaf_use:.2f}_MinAltReads_{alts_threshold}_Pruning/"

# Samples included in this analysis.
samples_to_analyse = ["57480", "58152", "61968", "63055", "63079", "63409", "63511_reseq", "71116", "72351_reseq", "73266", "76992_2", "57480_flowgel", "78454_flowgel"]
design_matrix_path = "DesignMatrix_VisiumHD.csv"
design_matrix = pd.read_csv(design_matrix_path, index_col=0)
sample_to_condition = design_matrix["condition"].to_dict()
sample_to_patient = design_matrix["donor"].to_dict()
samples_without_donor = [sample_name for sample_name in samples_to_analyse if pd.isna(sample_to_patient.get(sample_name))]
if samples_without_donor:
    raise ValueError(f"No donor in {design_matrix_path} for: {samples_without_donor}")

column_mss       = "MSS"
column_depth     = "mt_depth"
column_sample    = "Sample"
column_patient   = "patient"
column_cell_type = "Annotation"
column_condition = "condition"

condition_colors = {"Non-Fibrotic": "#1f77b4", "Fibrotic": "#ff7f0e",}

# A cell type enters the model of a sample only with at least this many spots.
# Reference cell types reach this minimum in every sample. They form the baseline.
minimum_spots_per_cell_type_in_sample = 20
minimum_patients_for_test = 3
output_directory_cell_type = os.path.join(input_directory, "MSS_Analysis", "CellType_DepthModel_raw_ols_linear")
mss_file_name = "MSS_per_spot.csv"
os.makedirs(output_directory_cell_type, exist_ok=True)


print("Loading the MSS values of all samples.", flush=True)
mss_tables_per_sample = []
for sample_name in samples_to_analyse:
    mss_file_path = os.path.join(input_directory, sample_name, mss_file_name)
    if not os.path.isfile(mss_file_path):
        raise FileNotFoundError(f"Missing MSS file for {sample_name}: {mss_file_path}")
    
    mss_table = pd.read_csv(mss_file_path, index_col=0)
    if column_sample not in mss_table.columns:
        raise KeyError(f"{mss_file_path} has no column {column_sample!r}.")
    
    if (mss_table[column_sample].astype(str) != sample_name).any():
        raise ValueError(f"{mss_file_path} contains spots of another sample than {sample_name}.")
    
    mss_table[column_patient] = sample_to_patient.get(sample_name, sample_name)
    if sample_to_condition:
        mss_table[column_condition] = sample_to_condition.get(sample_name, np.nan)
    
    print(f"{sample_name}: {len(mss_table)} spots, patient {mss_table[column_patient].iloc[0]}", flush=True)
    mss_tables_per_sample.append(mss_table)

combined = pd.concat(mss_tables_per_sample, axis=0, verify_integrity=True)


print("Checking the input data frame.", flush=True)
required_columns = [column_mss, column_depth, column_sample, column_patient, column_cell_type]
missing_columns  = [column_name for column_name in required_columns if column_name not in combined.columns]
if missing_columns:
    raise KeyError(f"combined misses the columns {missing_columns}. Available columns: {list(combined.columns)}.")


print("Removing spots with coverage less than 20.", flush=True)
combined = combined[combined[column_depth] >= 20]
print(f"Loaded {len(combined)} spots from {combined[column_sample].nunique()} samples and {combined[column_patient].nunique()} patients.", flush=True)

condition_available = column_condition in combined.columns
columns_to_keep = required_columns + ([column_condition] if condition_available else [])

spot_table = combined[columns_to_keep].copy()
spot_table = spot_table.rename(columns = {column_mss: "mss",
                                          column_depth: "mt_depth",
                                          column_sample: "sample",
                                          column_patient: "patient",
                                          column_cell_type: "cell_type"})
if condition_available:
    spot_table = spot_table.rename(columns = {column_condition: "condition"})

spot_table = spot_table.dropna(subset = ["mss", "mt_depth", "sample", "patient", "cell_type"]).copy()

spot_table["sample"]       = spot_table["sample"].astype(str)
spot_table["patient"]      = spot_table["patient"].astype(str)
spot_table["cell_type"]    = spot_table["cell_type"].astype(str)
spot_table["log_mt_depth"] = np.log(spot_table["mt_depth"])

patients_per_sample = spot_table.groupby("sample")["patient"].nunique()
if (patients_per_sample > 1).any():
    raise ValueError(f"These samples are assigned to more than one patient: {list(patients_per_sample[patients_per_sample > 1].index)}.")

sample_to_patient = spot_table.groupby("sample")["patient"].first()
patient_to_condition = spot_table.groupby("patient")["condition"].first() if condition_available else None


print("Defining the reference cell types.", flush=True)
spots_per_cell_type_and_sample = pd.crosstab(spot_table["cell_type"], spot_table["sample"])
spots_per_cell_type_and_sample.to_csv(os.path.join(output_directory_cell_type, "spots_per_celltype_and_sample_raw_ols.csv"))

reference_cell_types = spots_per_cell_type_and_sample.index[
    (spots_per_cell_type_and_sample >= minimum_spots_per_cell_type_in_sample).all(axis=1)
].tolist()
non_reference_cell_types = [cell_type for cell_type in spots_per_cell_type_and_sample.index if cell_type not in reference_cell_types]

# If we only have 1 reference cell type, we are using a too strict threshold.
# Throw an error in that case, but this should not happen with reasonable values.
if len(reference_cell_types) < 2:
    raise ValueError("Fewer than two reference cell types. Lower minimum_spots_per_cell_type_in_sample or merge cell types.")


print("Calculating the cell type effects per sample.", flush=True)
deviation_rows = []
for sample_name, sample_spots in spot_table.groupby("sample"):
    spots_per_cell_type = sample_spots["cell_type"].value_counts()
    cell_types_in_sample = spots_per_cell_type[spots_per_cell_type >= minimum_spots_per_cell_type_in_sample].index
    sample_spots = sample_spots[sample_spots["cell_type"].isin(cell_types_in_sample)].copy()
    
    # MSS ~ cell type + log depth. The first cell type (in alphabetical order) is
    # the reference and gets 0. The reference cancels out when the deviations from
    # the reference median are taken.
    sample_spots["cell_type_label"] = pd.Categorical(sample_spots["cell_type"])
    cell_type_categories = sample_spots["cell_type_label"].cat.categories
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fitted_model = smf.ols("mss ~ C(cell_type_label) + log_mt_depth", data=sample_spots).fit()
    
    # Coefficients are the MSS value. Dividing by the mean MSS of the spots in the model makes
    # them relative, so samples with different MSS levels are comparable.
    mean_mss_in_model = sample_spots["mss"].mean()
    # Checking if the mean MSS is 0 or less. That should not happen.
    # In this case we have a serious quality issue and detect nearly no varaints.
    # So, we throw an error.
    if mean_mss_in_model <= 0:
        raise ValueError(f"Sample {sample_name}: all spots in the model have MSS 0; relative effects are undefined.")
    
    effects_in_sample = pd.Series(0.0, index=cell_type_categories)
    for cell_type in cell_type_categories[1:]:
        effects_in_sample[cell_type] = float(fitted_model.params[f"C(cell_type_label)[T.{cell_type}]"])
    
    effects_in_sample = effects_in_sample / mean_mss_in_model
    
    # Deviation from the median of the other reference cell types (non-reference: all reference cell types).
    reference_effects = effects_in_sample[effects_in_sample.index.isin(reference_cell_types)]
    for cell_type, effect in effects_in_sample.items():
        baseline_effects = reference_effects.drop(cell_type, errors="ignore")
        if len(baseline_effects) == 0:
            continue
        
        deviation_rows.append({"sample": sample_name,
                               "patient": sample_to_patient[sample_name],
                               "cell_type": cell_type,
                               "is_reference_cell_type": cell_type in reference_cell_types,
                               "effect": effect,
                               "deviation_from_reference_median": effect - baseline_effects.median()})

deviation_per_sample = pd.DataFrame(deviation_rows)
deviation_per_sample.to_csv(os.path.join(output_directory_cell_type, "celltype_deviation_per_sample_raw_ols.csv"), index=False)


print("Averaging the deviations per patient and then testing.", flush=True)
# If we have multiple samples from the same patient, we average the values.
# This way, we do not count a patient twice.
deviation_per_patient = (
    deviation_per_sample.groupby(["patient", "cell_type", "is_reference_cell_type"])["deviation_from_reference_median"]
                        .agg(deviation_from_reference_median="mean", n_samples="size")
                        .reset_index()
)
if condition_available:
    deviation_per_patient["condition"] = deviation_per_patient["patient"].map(patient_to_condition)

deviation_per_patient.to_csv(os.path.join(output_directory_cell_type, "celltype_deviation_per_patient_raw_ols.csv"), index=False)

summary_rows = []
for cell_type, patient_values in deviation_per_patient.groupby("cell_type"):
    deviations = patient_values["deviation_from_reference_median"].values
    if len(deviations) >= minimum_patients_for_test and not np.all(deviations == 0):
        p_value = wilcoxon(deviations).pvalue
    else:
        p_value = np.nan
    
    summary_rows.append({"cell_type": cell_type,
                         "is_reference_cell_type": bool(patient_values["is_reference_cell_type"].iloc[0]),
                         "n_patients": len(deviations),
                         "median_deviation": np.median(deviations),
                         "mean_deviation": np.mean(deviations),
                         "fraction_patients_higher": np.mean(deviations > 0),
                         "P": p_value})

summary_table = pd.DataFrame(summary_rows)

# Benjamini-Hochberg correction over all tested cell types.
summary_table["P_adj"] = np.nan
rows_to_correct = summary_table["P"].notna()
if rows_to_correct.any():
    summary_table.loc[rows_to_correct, "P_adj"] = multipletests(summary_table.loc[rows_to_correct, "P"], method="fdr_bh")[1]

# The deviations are relative differences (0.25 means 25% higher).
summary_table["percent_difference_median"] = summary_table["median_deviation"] * 100
summary_table = summary_table.sort_values("median_deviation", ascending=False)
summary_table.to_csv(os.path.join(output_directory_cell_type, "celltype_comparison_summary_raw_ols.csv"), index=False)


print("Plotting the per-patient deviations.", flush=True)
random_generator = np.random.default_rng(0)
cell_type_order = summary_table["cell_type"].tolist()

figure, axis = plt.subplots(figsize = (max(6, 0.7 * len(cell_type_order)), 5))
for position, cell_type in enumerate(cell_type_order):
    cell_type_patients = deviation_per_patient[deviation_per_patient["cell_type"] == cell_type]
    jitter = random_generator.uniform(-0.15, 0.15, len(cell_type_patients))
    if condition_available:
        point_colors = cell_type_patients["condition"].map(condition_colors).fillna("black")
    else:
        point_colors = "black"
    
    axis.scatter(position + jitter, cell_type_patients["deviation_from_reference_median"],
                 c=point_colors, edgecolor="black", s=40, zorder=3)
    axis.hlines(np.median(cell_type_patients["deviation_from_reference_median"]), position - 0.3, position + 0.3,
                color="black", linewidth=2, zorder=4)

# Significance stars are placed after all points are drawn, so they sit above the data.
# This also means that the legend cannot see them and you should confirm that the
# legend does not cover anything important.
lower_limit, upper_limit = axis.get_ylim()
axis.set_ylim(lower_limit, upper_limit + 0.08 * (upper_limit - lower_limit))
for position, cell_type in enumerate(cell_type_order):
    adjusted_p_value = summary_table.loc[summary_table["cell_type"] == cell_type, "P_adj"].iloc[0]
    if pd.notna(adjusted_p_value) and adjusted_p_value < 0.05:
        axis.text(position, upper_limit, "*", ha="center", va="bottom", fontsize=14)

axis.axhline(0, color="grey", linewidth=0.8, linestyle="--", zorder=0)
axis.set_xticks(range(len(cell_type_order)))
axis.set_xticklabels(cell_type_order, fontsize=14, fontweight="bold", rotation=45, ha="right")
axis.set_ylabel("MSS Difference", fontsize=14, fontweight="bold")
axis.set_title("MSS Differences per Cell Type", fontsize=16, fontweight="bold")
legend = axis.legend(handles=[plt.Line2D([0], [0], marker="o", linestyle="",
                                         markerfacecolor=color, markeredgecolor="black", label=condition_name)
                              for condition_name, color in condition_colors.items()],
                     title="Condition", loc="upper left", bbox_to_anchor=(1.01, 1), prop={"weight": "bold", "size": 14})
legend.get_title().set_fontweight("bold")
legend.get_title().set_fontsize(14)

figure.tight_layout()
figure.savefig(os.path.join(output_directory_cell_type, "celltype_deviation_raw_ols_joint_model.png"), dpi=300, bbox_inches="tight")
figure.savefig(os.path.join(output_directory_cell_type, "celltype_deviation_raw_ols_joint_model.pdf"), dpi=300, bbox_inches="tight")
plt.close(figure)
