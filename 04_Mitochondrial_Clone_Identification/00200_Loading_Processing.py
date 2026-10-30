#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 26.10.2026

@author: Martin Grasshoff
Environment: spatial_sigurd_py_env

We load the data per sample and determine the SVVs.

Execute the script like this:
python 00200_Loading_Processing.py <interger_index>

The script needs as input:
- DesignMatrix_VisiumHD.csv
  - The files indicated in the design matrix are of course also needed.
- kidney7_marker_peak_integrate_scATAC_res0.9_annoted2_meta.csv
- kidney7_15mm_with_ref_2_meta_data.csv
- chrM_refAllele.txt # The mitochondrial reference alleles

Outputs are:
- variants per coverage plot
- strand correlation plot
- variants of interest CSV for further downstream analysis
- spatially variable variants derived from the VOIs
- pruned AnnData object
- pruned AnnData object with clonal information
- plots of the SVVs in space

"""

print("Loading functions.", flush=True)
import spatial_sigurd_py
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import math
import os
import numpy as np
import pandas as pd
import gc
import seaborn as sns
import anndata as ad
import sys
import scanpy as sc
import scipy.sparse as sp
ad.settings.allow_write_nullable_strings = True

def plot_variants_per_coverage(coverage_df, sample_use, output_path):
    """
    Parameters
    ----------
    coverage_df : pd.DataFrame
        DataFrame with at least a "avg_coverage" column (one row per variant).
    sample_use : str
        Sample name, used for title and output filename.
    output_path : str or Path
        Output directory where the plot will be saved.
    """
    
    print("We get the number of variants with at least X coverage.", flush=True)
    # Build coverage threshold summary
    thresholds = list(range(5, 91, 5))
    variants = [
        (coverage_df["avg_coverage"] >= cov).sum()
        for cov in thresholds
    ]
    
    variants_per_coverage = pd.DataFrame({
        "Coverage": thresholds,
        "Variants": variants
    })
    
    print("We plot the number of variants retained per coverage.", flush=True)
    
    # Plot
    sns.set_theme(style = "whitegrid")
    plt.figure(figsize=(6, 4))
    ax = sns.lineplot(data = variants_per_coverage, x="Coverage", y="Variants", marker="o")
    ax.set_xlabel("Coverage")
    ax.set_ylabel("Variants")
    ax.set_title(f"Variants {sample_use}")
    
    # Save plot
    output_file = f"{output_path}/Variants_Per_Coverage_{sample_use}.png"
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()


print("Parameters.", flush=True)
alts_threshold = 2
minimum_fw_rev_reads = 2
low_coverage_threshold = 10
radius_factor_use = 1.15
vaf_use = 0.10
genomic_reference_path = "positions/chrM_refAllele.txt"
design_matrix_path = "DesignMatrix_VisiumHD.csv"
design_matrix = pd.read_csv(design_matrix_path)
index_use = int(sys.argv[1])
# index_use = 8
meta_data_09_path = "kidney7_marker_peak_integrate_scATAC_res0.9_annoted2_meta.csv"
meta_data_25_path = "kidney7_integrate_scRNA_res2.5_annotation.csv"
meta_data_flowgel_path = "kidney7_15mm_with_ref_2_meta_data.csv"
sample_use = design_matrix.loc[index_use,"patient"]
positions_path = design_matrix.loc[index_use,"positions"]
coverage_use = design_matrix.loc[index_use,"coverage"]
metadata_sample_name = design_matrix.loc[index_use,"metadata_name"]
output_directory = f"MinVAF_{vaf_use:.2f}_MinAltReads_{alts_threshold}_Pruning/"
output_directory_sample = f"{output_directory}{sample_use}"
print("Using:", flush=True)
print("Sample:   " + sample_use, flush=True)
print("Coverage: " + str(coverage_use), flush=True)
print("Radius:   " + str(radius_factor_use), flush=True)


# We make the necessary directory.
if not os.path.isdir(f"{output_directory_sample}/Variant_On_Spatial/"):
    os.makedirs(f"{output_directory_sample}/Variant_On_Spatial/")


print("Setting scanpy output folder.", flush=True)
sc.settings.figdir = f"{output_directory_sample}"
sc.set_figure_params(dpi_save=150, figsize=[10,10])


print("Loading the data.", flush=True)
adata = spatial_sigurd_py.readwrite.LoadingMGATK_typewise_sparse(
    design_matrix_path = design_matrix_path,
    patient = sample_use,
    reference_path = genomic_reference_path,
    samples_path = None,
    barcodes_path = None,
    patient_column = "patient",
    type_use = "ATAC",
    chromosome_prefix = "chrM",
    min_cells = 2,
    cells_include = None,
    cells_exclude = None,
    verbose = True)


print("Loading the positions.", flush=True)
positions = pd.read_csv(positions_path, header = None)
positions[0] = sample_use + "_" + positions[0].astype(str) + "-1"
positions = positions[positions[1] == 1]
positions = positions[[0, 5, 4]]
positions.columns=["CellID", "X", "Y"]
positions = positions[positions["CellID"].isin(adata.obs_names)]


print("Loading the meta data.", flush=True)
meta_data_flowgel = pd.read_csv(meta_data_flowgel_path, index_col=0)
meta_data_flowgel["Sample"] = meta_data_flowgel["Sample"].astype(str).replace(r"Diseased_|Healthy_", "", regex=True)
meta_data_flowgel["Sample"] = meta_data_flowgel["Sample"].astype(str) + "_flowgel"
barcodes = meta_data_flowgel.index.str.replace(r".*_", "", regex=True)
meta_data_flowgel.index = meta_data_flowgel["Sample"].to_numpy() + "_" + barcodes.to_numpy()
meta_data_flowgel["Annotation"] = meta_data_flowgel["annotation"].str.replace(r"-|/| ", "_", regex=True)
meta_data_flowgel = meta_data_flowgel[["Sample", "Annotation"]]

meta_data_09 = pd.read_csv(meta_data_09_path, index_col=0)
meta_data_09.index = meta_data_09.index.str.replace(r"^P", "", regex=True)
meta_data_09["Sample"] = meta_data_09["Sample"].str.replace(r"^P", "", regex=True)
meta_data_09["Annotation"] = meta_data_09["Clusters2"].str.replace(r"-|/| ", "_", regex=True)
meta_data_09 = meta_data_09[["Sample", "Annotation"]]

meta_data = pd.concat([meta_data_flowgel, meta_data_09], axis=0, join="outer", verify_integrity=True)


print("Aligning the cell type names.", flush=True)
meta_data["Annotation"] = meta_data["Annotation"].str.replace(r"ATL_DTL", "ATL", regex=True)
meta_data["Annotation"] = meta_data["Annotation"].str.replace(r"PC_Endothelial|AVR_Endothelial|Endothelial|Endo", "Endothelial Cells", regex=True)
meta_data["Annotation"] = meta_data["Annotation"].str.replace(r"FIB|Fibroblast", "Fibroblasts", regex=True)
meta_data["Annotation"] = meta_data["Annotation"].str.replace(r"PV_Cap_fibroblasts", "PV-Cap Fibroblasts", regex=True)
meta_data["Annotation"] = meta_data["Annotation"].str.replace(r"Immune cells|Immune_cells|Immune", "Immune Cells", regex=True)
meta_data["Annotation"] = meta_data["Annotation"].str.replace(r"Pod|POD", "Podocytes", regex=True)
meta_data["Annotation"] = meta_data["Annotation"].str.replace(r"PT_S\d(_S\d)?", "PT", regex=True)
meta_data["Annotation"] = meta_data["Annotation"].str.replace(r"Injured_PTs", "Injured-PT", regex=True)
meta_data["Annotation"] = meta_data["Annotation"].str.replace(r"VSMCs?_Pericytes?", "VSMC-Pericytes", regex=True)


print("Adding the cell type annotation.", flush=True)
metadata_spot_ids   = metadata_sample_name + adata.obs_names.str[len(sample_use):]
annotation_per_spot = pd.Series(meta_data["Annotation"].reindex(metadata_spot_ids).values, index=adata.obs_names)
adata.obs["Annotation"] = annotation_per_spot


print("Filtering the data.", flush=True)
cells_include = adata.obs_names[metadata_spot_ids.isin(meta_data.index) & adata.obs_names.isin(positions["CellID"])]
adata = spatial_sigurd_py.pp.Filtering(
    adata = adata,
    cells_include = cells_include,
    cells_exclude = None,
    fraction_threshold = vaf_use,
    keep_refs = True,
    alts_threshold = alts_threshold,
    min_cells_per_variant = None,
    min_variants_per_cell = None,
    verbose = True,)

gc.collect()

# We add the positions to the annData object. Since we have already removed all
# spots that are not in positions and the meta data, we are not adding spots
# that should have been removed. This would happen with reindex since we would
# introduce NAs.
positions = positions.set_index("CellID").reindex(adata.obs_names)
adata.obsm["spatial"] = positions


print("Performing the pruning.", flush=True)
adata_pruned = spatial_sigurd_py.pp.Pruning(
    adata = adata,
    spatial_key = "spatial",
    x_col = "X",
    y_col = "Y",
    radius_factor = radius_factor_use,
    drop_zero_variants = True,
    drop_zero_cells = False,
    copy = True,
    verbose = True,)


print("Saving the output.", flush=True)
adata_pruned.write_h5ad(f"{output_directory_sample}/adata_pruned.h5ad", compression="gzip", compression_opts=9)
# adata_pruned = sc.read_h5ad(f"{output_directory_sample}/adata_pruned.h5ad")


print("Plotting the coverage distribution.", flush=True)
vals = np.asarray(adata_pruned.var["avg_coverage"].values, dtype = float)
# Remove NaN or inf values
vals = vals[np.isfinite(vals)]
n = vals.size
if n == 0:
    raise ValueError("No valid numeric data found in 'adata_pruned.var['avg_coverage']'.")

# Freedman‑Diaconis bin width for the histogram
q1, q3 = np.percentile(vals, [25, 75])
iqr = q3 - q1
bin_width = 2 * iqr / (n ** (1.0 / 3.0))

# If IQR == 0 fall back to 1 bin and basically generate a bar plot.
if bin_width <= 0:
    bins = 1
else:
    bins = max(1, math.ceil((vals.max() - vals.min()) / bin_width))


# Plot
fig, ax = plt.subplots(figsize=(6, 4))
ax.hist(vals, bins=bins)
ax.set_xlabel("Average coverage")
ax.set_ylabel("Count")
ax.set_title("Histogram of average coverage")

# Save
out_path = f"{output_directory_sample}/{sample_use}_avg_coverage_hist.png"
fig.tight_layout()
fig.savefig(out_path, dpi=150)
plt.show()
plt.close(fig)


print("We get the number of variants with at least X coverage.", flush=True)
plot_variants_per_coverage(coverage_df = adata_pruned.var, 
                           sample_use = sample_use, 
                           output_path = f"{output_directory_sample}/")


print("Selecting variants.", flush=True)
vois = spatial_sigurd_py.tl.Variant_Selection_VMR(
    adata = adata_pruned,
    stabilize_variance = True,
    low_coverage_threshold = low_coverage_threshold,
    minimum_fw_rev_reads = minimum_fw_rev_reads,
    verbose = True,
)


print("Saving variants.", flush=True)
vois.to_csv(f"{output_directory_sample}/VOIs_All.csv")
# vois = pd.read_csv(f"{output_directory_sample}/VOIs_All.csv", index_col=0)


print("We save the matrices as CSVs.", flush = True)
alts_matrix = adata_pruned.layers["alt_reads"]
alts_df = pd.DataFrame(alts_matrix.toarray(), index = adata_pruned.obs_names, columns = adata_pruned.var_names)
alts_df.to_csv(f"{output}{sample_use}/alts_matrix.csv")

refs_matrix = adata_pruned.layers["ref_reads"]
refs_df = pd.DataFrame(refs_matrix.toarray(), index = adata_pruned.obs_names, columns = adata_pruned.var_names)
refs_df.to_csv(f"{output}{sample_use}/refs_matrix.csv")

coverage_matrix = adata_pruned.layers["coverage"]
coverage_df = pd.DataFrame(coverage_matrix.toarray(), index = adata_pruned.obs_names, columns = adata_pruned.var_names)
coverage_df.to_csv(f"{output}{sample_use}/coverage_matrix.csv")

fraction_matrix = adata_pruned.layers["fraction"]
fraction_df = pd.DataFrame(fraction_matrix.toarray(), index = adata_pruned.obs_names, columns = adata_pruned.var_names)
fraction_df.to_csv(f"{output}{sample_use}/fraction_matrix.csv")


print("Generating a VOIs scatter plot.", flush=True)
vois["vmr_log10_neg"] = -np.log10(vois["vmr"])
vois_filtered = vois.query("mean_coverage > @coverage_use")
mask_red = (vois_filtered["vmr_log10_neg"] > 2) & (vois_filtered["strand_correlation"] > 0.65)
vois_filtered["point_color"] = ["red" if r else "gray" for r in mask_red]
fig, ax = plt.subplots(figsize=(8, 6))

ax.scatter(vois_filtered["vmr_log10_neg"], vois_filtered["strand_correlation"], c=vois_filtered["point_color"], alpha=0.8, edgecolor="k", label="Variants")
ax.axhline(y = 0.65, color = "black", linestyle = "--", linewidth = 1.5, label="y = 0.65")
ax.axvline(x = 2, color = "black", linestyle = "--", linewidth = 1.5, label="x = 2")

ax.set_xlabel("vmr_log10_neg")
ax.set_ylabel("strand_correlation")
ax.set_title("Strand correlation vs. VMR (log10, negative)")
ax.legend()
ax.grid(True, which = "both", ls = ":", alpha=0.5)

fig.savefig(f"{output_directory_sample}/strand_corr_plot.png", dpi=150, bbox_inches="tight")


print("We get the spatial clones.", flush=True)
vois_filtered = vois.query("strand_correlation > 0.65 and mean_coverage > @coverage_use and cells_detected >= 5")
adata_pruned_subset = adata_pruned[:,vois_filtered.variant].copy()
adata_clones, variant_df, svg_list, variant_df_enhanced = spatial_sigurd_py.tl.Spatial_Clone_Detection(
    adata = adata_pruned_subset,
    spatial_key = "spatial",
    filter_by_min_spots = 5,
    normalize = True,
    neighbors_for_ratio = 8,
    neighbors_for_ratio_absolute = True,
    detect_svg_normalize_lap = False,
    detect_svg_filter_peaks = True,
    detect_svg_S = 1,
    detect_svg_cal_pval = True,
    p_value_threshold_svvs = 0.05,
    perform_low_pass_enhancement = False,
    low_pass_enhancement_ratio_low_freq = 15,
    low_pass_enhancement_c=0.003,
    low_pass_enhancement_normalize_lap = False,
    perform_ftu = False,
    ftu_resolution=(0.7, 1.8, 0.1),
    ftu_weight_by_freq = False,
    ftu_normalize_lap = False,
    neighbors_for_ftu = 8,
    ftu_random_state = 0,
    ftu_on_error = "warn",
    verbose = True,
)

maximum_number_of_spots = math.ceil(math.sqrt(adata_pruned_subset.shape[0]))
variant_df = variant_df[variant_df.n_cells < maximum_number_of_spots]
mask = variant_df.cutoff_gft_score.astype(bool)
mask &= variant_df.fdr < 0.05
svg_list = variant_df.index[mask].tolist()


print("Plotting the SVGs.", flush=True)
plot_svgs = svg_list[:8]
sc.pl.spatial(adata_clones, color=plot_svgs, size=1.6, cmap="magma", use_raw=False, spot_size=25, save="_SVGs_Top8.png")


print("Saving the output.", flush=True)
adata_clones.write_h5ad(f"{output_directory_sample}/adata_spagft.h5ad", compression="gzip", compression_opts=9)
# adata_clones = sc.read_h5ad(f"{output_directory_sample}/adata_spagft.h5ad")


# This is useful if the data is reloaded again.
# print("Selecting the SVVs from the data again.", flush=True)
# variant_df = adata_clones.var.copy()
# variant_df = variant_df[variant_df.n_cells < maximum_number_of_spots]
# mask = variant_df.cutoff_gft_score.astype(bool)
# mask &= variant_df.fdr < 0.05
# svg_list = variant_df.index[mask].tolist()


print("Saving the SVGs specifically.", flush = True)
pd.DataFrame(svg_list, columns=["SVV"]).to_csv(f"{output_directory_sample}/SpatiallyVariableVariants.csv", index=False)


print("We plot the variants on the spatial plot using the original data.", flush=True)
coords = adata_clones.obsm["spatial"]
x = coords[:, 0]
y = coords[:, 1]
cmap = ListedColormap(["#B4B4B4", "#FF6879"]) # 0 as grey and 1 as light red
for variant_use in svg_list:
    print(variant_use, flush=True)
    if "X_original" in adata_clones.layers:
        count_values = adata_clones[:, variant_use].layers["X_original"]
    else:
        count_values = adata_clones[:, variant_use].X
    
    if sp.issparse(count_values):
        count_values = count_values.toarray().ravel()
    else:
        count_values = np.array(count_values).ravel()
    
    # Binarizing the values
    binary = (count_values > 0).astype(int)
    order = np.argsort(count_values)
    
    # Plotting and saving
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(x[order], y[order], c=binary[order], cmap=cmap, s=5)
    ax.invert_yaxis()
    ax.axis("off")
    ax.set_title(variant_use, pad=4)
    ax.margins(x=0.01, y=0.01)
    
    filename = f"{output_directory_sample}/Variant_On_Spatial/{variant_use}_spatial.png"
    plt.savefig(filename, dpi=150, bbox_inches="tight", pad_inches=0)
    plt.close()
