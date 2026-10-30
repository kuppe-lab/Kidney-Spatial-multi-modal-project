#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 13.5.2026

@author: grasshoff

Environment: spatial_sigurd_py_env

We add the clones to the meta data for each sample.
Each clone is written in the following way:
sample_id: clone_name: svv1, svv2

Output is a CSV file with the clonal information and the annotation.
"""

print("Packages.", flush=True)
import numpy as np
import pandas as pd
import scanpy as sc
import warnings

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
    "78454_flowgel": {
        "15734_G_A": ["chrM_15734_G>A"],
        "556_A_G": ["chrM_556_A>G"],
    },
}


# Name of the new obs column.
clone_column_name = "clone_manual"

# Output filename in {input}{sample}/.
output_filename = "meta_data_with_clones.csv"


print("Variables.", flush=True)
vaf_use = 0.10
alts_threshold = 2
# Sample we want to actually use.
samples = ["57480", "57480_flowgel", "58152", "61968", "63055", "63079", "63409", "63511_reseq", "71116", "72351_reseq", "73266", "76992_2", "78454", "78454_flowgel"]
input_directory = f"MinVAF_{vaf_use:.2f}_MinAltReads_{alts_threshold}_Pruning/"


# Validate sample-wise clone definitions.
forbidden_names      = {"WT"}
forbidden_characters = {"+"}

if len(clone_definitions) == 0:
    raise ValueError("clone_definitions is empty.")

samples_without_definitions = sorted(set(samples) - set(clone_definitions.keys()))
definitions_without_samples = sorted(set(clone_definitions.keys()) - set(samples))
if samples_without_definitions:
    warnings.warn(f"No clone definitions for {samples_without_definitions}, these samples are skipped.")

if definitions_without_samples:
    warnings.warn(f"{definitions_without_samples} are in clone_definitions but not in samples, these definitions are ignored.")

sample_clone_signatures = {}
sample_tracked_svvs     = {}
for sample_id, clones in clone_definitions.items():
    if not isinstance(clones, dict) or len(clones) == 0:
        raise ValueError(f"Sample '{sample_id}': clone definitions must be a non-empty dictionary of clone_name: [svvs].")
    
    clone_signatures        = {}
    clone_name_by_signature = {}
    for clone_name, clone_svvs in clones.items():
        if clone_name in forbidden_names:
            raise ValueError(f"Sample '{sample_id}': clone name '{clone_name}' is forbidden.")
        
        for forbidden_character in forbidden_characters:
            if forbidden_character in clone_name:
                raise ValueError(f"Sample '{sample_id}': clone name '{clone_name}' contains the forbidden character '{forbidden_character}'.")
        
        if not isinstance(clone_svvs, (list, tuple, set, frozenset)) or len(clone_svvs) == 0:
            raise ValueError(f"Sample '{sample_id}': clone '{clone_name}' needs a non-empty list of SVVs.")
        
        clone_signature = frozenset(clone_svvs)
        if clone_signature in clone_name_by_signature:
            raise ValueError(f"Sample '{sample_id}': clones '{clone_name_by_signature[clone_signature]}' and '{clone_name}' have identical SVV sets.")
        
        clone_name_by_signature[clone_signature] = clone_name
        clone_signatures[clone_name]             = clone_signature
    
    clone_names_for_check = list(clone_signatures.keys())
    for first_index in range(len(clone_names_for_check)):
        for second_index in range(first_index + 1, len(clone_names_for_check)):
            first_clone_name  = clone_names_for_check[first_index]
            second_clone_name = clone_names_for_check[second_index]
            shared_svvs       = clone_signatures[first_clone_name] & clone_signatures[second_clone_name]
            if shared_svvs:
                warnings.warn(f"NOTICE: sample '{sample_id}': clones '{first_clone_name}' and '{second_clone_name}' share {sorted(shared_svvs)}, spots carrying these get combined labels.")
    
    sample_clone_signatures[sample_id] = clone_signatures
    sample_tracked_svvs[sample_id]     = sorted({svv for clone_svvs in clones.values() for svv in clone_svvs})


# We check if the SVVs are all present in the data. An error is thrown if any one is 
# missing. We only load the object in backed mode, which only loads the obs and var
# information. So, this is quicker than loading the entire object.
print("Checking SVV presence.", flush=True)
missing_svvs = []
for sample_use in samples:
    if sample_use not in clone_definitions:
        continue
    
    input_path   = f"{input_directory}{sample_use}/adata_pruned.h5ad"
    adata_backed = sc.read_h5ad(input_path, backed="r")
    variant_names_in_sample = set(adata_backed.var_names)
    adata_backed.file.close()
    
    for svv in sample_tracked_svvs[sample_use]:
        if svv not in variant_names_in_sample:
            missing_svvs.append(f"'{svv}' in sample '{sample_use}'")

if missing_svvs:
    raise ValueError(f"SVVs missing from the AnnData object: {', '.join(missing_svvs)}.")


# Sample-wise loop: assign the labels and save the meta data.
number_of_samples_done    = 0
number_of_samples_skipped = 0
for sample_use in samples:
    print(f"Sample: {sample_use}", flush=True)
    
    if sample_use not in clone_definitions:
        # No clone definition for this sample.
        # No message, since we already have tested for this above.
        number_of_samples_skipped += 1
        continue
    
    input_path  = f"{input_directory}{sample_use}/adata_pruned.h5ad"
    output_path = f"{input_directory}{sample_use}/{sample_use}_{output_filename}"
    
    print(f"Loading: {input_path}", flush=True)
    adata = sc.read_h5ad(input_path)
    
    tracked_svvs     = sample_tracked_svvs[sample_use]
    clone_signatures = sample_clone_signatures[sample_use]
    clone_names      = list(clone_definitions[sample_use].keys())
    
    svv_matrix   = adata[:, tracked_svvs].X
    svv_matrix   = svv_matrix.toarray() if hasattr(svv_matrix, "toarray") else svv_matrix
    svv_detected = svv_matrix > 0
    
    number_of_spots = svv_detected.shape[0]
    clone_labels    = np.empty(number_of_spots, dtype=object)
    for spot_index in range(number_of_spots):
        spot_svvs = frozenset(
            tracked_svvs[svv_index]
            for svv_index in range(len(tracked_svvs))
            if svv_detected[spot_index, svv_index]
        )
        if len(spot_svvs) == 0:
            clone_labels[spot_index] = "WT"
            continue
        
        matching_clone_names     = [clone_name for clone_name in clone_names if clone_signatures[clone_name] & spot_svvs]
        clone_labels[spot_index] = "+".join(sorted(matching_clone_names))
    
    # Convert to category with a set order. 
    # First, the defined clones (as in the dictionary).
    # Second, the combined labels that were observed.
    # Third, the wildtype clones.
    # We do not drop clones if they are empty, just in case we are working with a 
    # subset AnnData object.
    observed_labels  = pd.unique(pd.Series(clone_labels))
    composite_labels = sorted([label for label in observed_labels if label != "WT" and "+" in label])
    categories       = list(clone_names) + composite_labels
    if "WT" in observed_labels:
        categories.append("WT")
    
    adata.obs[clone_column_name] = pd.Categorical(
        clone_labels,
        categories = categories,
        ordered    = False,
    )
    
    # Quick label-wise summary while running to verify the definitions.
    summary = (
        pd.Series(clone_labels)
        .value_counts()
        .rename_axis(clone_column_name)
        .reset_index(name = "number_of_spots")
    )
    print(f"Assignment summary ({number_of_spots} spots total):", flush=True)
    for _, summary_row in summary.iterrows():
        print(f"{summary_row[clone_column_name]:40s} {int(summary_row['number_of_spots']):>8d}", flush=True)
    
    meta_data = adata.obs
    meta_data.to_csv(output_path)
    number_of_samples_done += 1

print(f"Done. {number_of_samples_done} samples annotated, {number_of_samples_skipped} samples skipped.", flush=True)
