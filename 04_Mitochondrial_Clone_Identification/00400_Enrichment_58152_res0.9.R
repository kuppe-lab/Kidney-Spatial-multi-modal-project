"""
@author: Martin Grasshoff
Environment: Spatial_ATAC_MT_Genotyping_Kuppe_Lab

Performing the enrichment analysis for 58152
"""
print("Libraries.")
suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(tibble)
  library(ComplexHeatmap)
})


print("Significance legend.")
legend_sig <- ComplexHeatmap::Legend(
  title      = "Significance",
  pch        = c("*", "**", "***"),
  type       = "points",
  labels     = c("< 0.05", "< 0.01", "< 0.001"),
  direction  = "vertical",
  title_gp   = grid::gpar(fontsize = 12, fontface = "plain"),
  labels_gp  = grid::gpar(fontsize = 12),
  border     = FALSE,
  background = "transparent",
)


print("Variables.")
sample_use <- "58152"
input <- "MinVAF_0.10_MinAltReads_2_Pruning/"
output <- "MinVAF_0.10_MinAltReads_2_Pruning/"
input_use <- paste0(input, sample_use, "/")
output_use <- paste0(output, sample_use, "/Enrichment_kidney7_res0.9_Variant/")
vois <- read.csv(paste0(input_use, "SpatiallyVariableVariants.csv"))[,1]
vois <- gsub(">", "_", vois)
dir.create(output_use, showWarnings = FALSE, recursive = TRUE)


print("We load the meta data.")
meta_data <- read.table("kidney7_marker_peak_integrate_scATAC_res0.9_annoted2_meta.csv", sep = ",", header = TRUE, row.names = 1)
meta_data[["Sample"]] <- gsub("^P", "", meta_data[["Sample"]])
rownames(meta_data) <- gsub("^P", "", rownames(meta_data))
meta_data[["CellType"]] <- meta_data[["Clusters2"]]


print("We load the fraction.")
fraction <- read.table(paste0(input, sample_use, "/fraction_matrix.csv"), sep = ",", row.names = 1, header = TRUE)
colnames(fraction) <- gsub("\\.", "_", colnames(fraction))
fraction[fraction > 0] <- 1
meta_data <- meta_data[rownames(meta_data) %in% rownames(fraction), ]
fraction <- fraction[rownames(meta_data), ]
meta_data[, vois] <- fraction[, vois]
for(i in seq_along(vois)){
  voi_use <- vois[i]
  meta_data[, voi_use] <- ifelse(meta_data[, voi_use] == 1, "Alt", "Ref")
}

print("We perform the enrichment test.")
celltypes <- sort(unique(meta_data[["CellType"]]))
results <- c()
cells_per_variant_and_celltype <- matrix(NA, nrow = length(vois), ncol = length(celltypes), dimnames = list(vois, celltypes))
for(i in seq_along(vois)){
  voi_use <- vois[i]
  print(voi_use)
  for(j in 1:length(celltypes)){
    celltype_use <- celltypes[j]
    print(celltype_use)
    count_matrix <- matrix(NA, nrow = 2, ncol = 2, dimnames = list(c(voi_use, "Other"), c(celltype_use, "Other")))
    count_matrix[voi_use, celltype_use] <- nrow(meta_data[meta_data[, voi_use] == "Alt" & meta_data[, "CellType"] == celltype_use, ])
    count_matrix["Other", celltype_use] <- nrow(meta_data[meta_data[, voi_use] != "Alt" & meta_data[, "CellType"] == celltype_use, ])
    count_matrix[voi_use, "Other"] <- nrow(meta_data[meta_data[, voi_use] == "Alt" & meta_data[, "CellType"] != celltype_use, ])
    count_matrix["Other", "Other"] <- nrow(meta_data[meta_data[, voi_use] != "Alt" & meta_data[, "CellType"] != celltype_use, ])
    res <- fisher.test(count_matrix)
    res <- data.frame(VOI = voi_use, CellType = celltype_use, OddsRatio = res$estimate, P = res$p.value)
    results <- rbind(results, res)
    cells_per_variant_and_celltype[voi_use, celltype_use] <- count_matrix[1,1]
  }
}
results[["P_adj"]] <- p.adjust(results[["P"]], method = "fdr")
results[["P_adj_log10"]] <- -log10(results[["P_adj"]])
write.table(results, paste0(output_use, "/Enrichment_Variant_CellType_Enhanced_kidney7_res0.9_Merged.csv"), sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(cells_per_variant_and_celltype, paste0(output_use, "/cells_per_variant_and_celltype_enhanced_kidney7_res0.9_Merged.csv"), sep = ",", col.names = TRUE, row.names = TRUE, quote = FALSE)
