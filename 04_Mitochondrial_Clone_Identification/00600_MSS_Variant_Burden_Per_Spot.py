#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Martin Grasshoff
Environment: spatial_sigurd_py_env

We get the variant burden per spot. We also use the mtDNA locus
constraint (MLC) score sum (MSS).

Then, we plot the value per spot on the spatial plot.

"""

print("Packages.", flush=True)
import anndata as ad
import numpy as np
import pandas as pd
import os
import sys


print("Variables.", flush=True)
alts_threshold = 2
vaf_use = 0.10
design_matrix_path = "DesignMatrix_VisiumHD.csv"
design_matrix = pd.read_csv(design_matrix_path)
index_use = int(sys.argv[1])
# index_use = 10
sample_use = design_matrix.loc[index_use,"sample"]
coverage_use = design_matrix.loc[index_use,"coverage"]
mt_depth_path = design_matrix.loc[index_use,"depth_table"]
path_mlc_values = "MLC_supplementary_dataset_7.tsv"
input_directory = f"MinVAF_{vaf_use:.2f}_MinAltReads_{alts_threshold}_Pruning/"
output_directory = f"MinVAF_{vaf_use:.2f}_MinAltReads_{alts_threshold}_Pruning/"
output_directory_sample = f"{output_directory}{sample_use}"
print("Using:", flush=True)
print("Input:    " + input_directory, flush=True)
print("Output:   " + output_directory, flush=True)
print("Sample:   " + sample_use, flush=True)
print("Coverage: " + str(coverage_use), flush=True)
os.makedirs(output_directory_sample, exist_ok=True)


print("Loading the variants of interest.", flush=True)
vois = pd.read_csv(f"{input_directory}{sample_use}/VOIs_All.csv", index_col=0)
vois = vois.query("strand_correlation > 0.65 and vmr > 0.01 and mean_coverage > @coverage_use")
vois = vois["variant"].values


print("Loading the AnnData object from 00200.", flush=True)
adata_pruned = ad.read_h5ad(f"{input_directory}{sample_use}/adata_pruned.h5ad")
# A variant that is listed twice counts once.
# This should not happen, but since this can happen when supplying a hand-written
# VOI list, we better check.
vois = pd.unique(vois)
# Spots x variants. The MSS uses the allele fractions of all variants of interest.
adata_vois = adata_pruned[:, vois]


print("Loading the mitochondrial sequencing depth.", flush=True)
mt_depth = pd.read_csv(mt_depth_path, sep="\t", header=None, names=["CellID", "mt_depth"])
mt_depth.index = sample_use + "_" + mt_depth["CellID"].astype(str)
mt_depth = mt_depth["mt_depth"]


print("Adding the cell type annotation from 00200.", flush=True)
# 00200 only keeps in-tissue spots, so all spots of the AnnData object are used.
annotation_per_spot = adata_vois.obs["Annotation"]
print(f"Annotated spots: {annotation_per_spot.notna().sum()} of {adata_vois.n_obs}", flush=True)
if annotation_per_spot.notna().sum() == 0:
    raise ValueError(f"No spot of {sample_use} has an annotation in adata_pruned.h5ad.")

annotation_per_spot = annotation_per_spot.dropna()
annotated_spots = annotation_per_spot.index


print("Calculating the MSS per spot.", flush=True)
mlc_values = pd.read_table(path_mlc_values, sep = "\t")
mlc_variant_names = ("chrM_" + mlc_values["Position"].astype(str) + "_" + mlc_values["Reference"].astype(str)
                     + ">" + mlc_values["Alternate"].astype(str))
mlc_per_variant = pd.Series(mlc_values["MLC_score"].values, index=mlc_variant_names)
# Variants without an MLC score contribute 0.
mlc_per_variant_aligned = mlc_per_variant.reindex(adata_vois.var_names).fillna(0)
fraction_annotated_spots = adata_vois[annotated_spots].layers["fraction"]
mss_per_spot = pd.Series(np.asarray(fraction_annotated_spots @ mlc_per_variant_aligned.to_numpy()).ravel(), index=annotated_spots)


print("Saving the MSS per spot.", flush=True)
mss_table = pd.DataFrame({"Sample":     sample_use,
                          "MSS":        mss_per_spot,
                          "mt_depth":   mt_depth.reindex(annotated_spots),
                          "Annotation": annotation_per_spot},
                         index = annotated_spots)
mss_table.to_csv(f"{output_directory_sample}/MSS_per_spot.csv")
print(f"Saved {len(mss_table)} spots to {output_directory_sample}/MSS_per_spot.csv", flush=True)
