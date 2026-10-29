# We get the clonal cells from all samples and perform the enrichment.
# We only check clonal cells versus cell type.

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
  title_gp   = grid::gpar(fontsize = 12, fontface = "bold"),
  labels_gp  = grid::gpar(fontsize = 12),
  border     = FALSE,
  background = "transparent",
)


print("Variables.")
input <- "MinVAF_0.10_MinAltReads_2_Pruning/"
output <- "MinVAF_0.10_MinAltReads_2_Pruning/"

# The samples are inconsistently named. I need a list to keep track of them.
samples <- data.frame(Name = c("57480", "63055", "63079", "63511_reseq", "72351_reseq", "57480_flowgel", "78454_flowgel", "71116", "76992_2", "73266", "63409", "58152", "61968"),
                      Name_In_MetaData = c("57480", "63055", "63079", "63511_flowgel", "72351_flowgel", "57480_flowgel", "78454_flowgel", "71116_flowgel", "76992_flowgel", "73266_flowgel", "63409", "58152", "61968"))
vois_list <- list()
for(i in 1:nrow(samples)){
  sample_use <- samples[i, "Name"]
  name_metadata_use <- samples[i, "Name_In_MetaData"]
  input_use <- paste0(input, sample_use, "/")
  vois_use <- read.table(paste0(input_use, "SpatiallyVariableVariants.csv"), header = TRUE, sep = ",")[["SVV"]]
  vois_use <- gsub(">", "_", vois_use)
  vois_list[[name_metadata_use]] <- vois_use
}
input_use <- paste0(input, sample_use, "/")
output_use <- paste0(output, "Enrichment_Variant_Merged_Samples/")
dir.create(output_use, showWarnings = FALSE, recursive = TRUE)


print("We load the meta data.")
meta_data_flowgel_new <- read.table("meta_data/kidney7_15mm_with_ref_2/meta_data.csv", sep = ",", header = TRUE, row.names = 1)
meta_data_flowgel_new[["Sample"]] <- gsub("Diseased_|Healthy_", "", meta_data_flowgel_new[["Sample"]])
meta_data_flowgel_new[["Sample"]] <- paste0(meta_data_flowgel_new[["Sample"]], "_flowgel")
rownames(meta_data_flowgel_new) <- paste(meta_data_flowgel_new[["Sample"]], gsub(".*_", "", rownames(meta_data_flowgel_new)), sep = "_")
meta_data_flowgel_new[["annotation"]] <- gsub("-|/| ", "_", meta_data_flowgel_new[["annotation"]])
meta_data_flowgel_new[["CellType"]] <- meta_data_flowgel_new[["annotation"]]
meta_data_flowgel_new <- meta_data_flowgel_new[, c("Sample", "CellType")]

meta_data_09 <- read.table("meta_data/kidney7_marker_peak_integrate_scATAC_res0.9_annoted2_meta.csv", sep = ",", header = TRUE, row.names = 1)
meta_data_09[["Sample"]] <- gsub("^P", "", meta_data_09[["Sample"]])
rownames(meta_data_09) <- gsub("^P", "", rownames(meta_data_09))
meta_data_09[["CellType"]] <- meta_data_09[["Clusters2"]]
meta_data_09 <- meta_data_09[, c("Sample", "CellType")]


print("We align the cell type names.")
meta_data <- rbind(meta_data_flowgel_new, meta_data_09)
meta_data[["Clonal"]] <- NA
meta_data$CellType <- gsub("ATL_DTL", "ATL", meta_data$CellType)
meta_data$CellType <- gsub("PC_Endothelial|AVR_Endothelial|Endothelial", "Endo", meta_data$CellType)
meta_data$CellType <- gsub("Fibroblast", "FIB", meta_data$CellType)
meta_data$CellType <- gsub("Immune cells", "Immune", meta_data$CellType)
meta_data$CellType <- gsub("Injured PTs", "Injured_PTs", meta_data$CellType)
meta_data$CellType <- gsub("Pod", "POD", meta_data$CellType)
meta_data$CellType <- gsub("PT_S1|PT_S1_S2|PT_S2|PT_S2_S3|PT_S3", "PT", meta_data$CellType)
meta_data$CellType <- gsub("VSMC_Pericytes", "VSMC_Pericyte", meta_data$CellType)


print("We load the fraction and add the clonal information to the meta data.")
for(i in 1:nrow(samples)){
  print(paste0(i, " out of ", nrow(samples)))
  sample_use <- samples[i, "Name"]
  name_metadata_use <- samples[i, "Name_In_MetaData"]
  vois_use <- vois_list[[name_metadata_use]]
  fraction <- read.table(paste0(input, sample_use, "/fraction_matrix.csv"), sep = ",", row.names = 1, header = TRUE)
  colnames(fraction) <- gsub("\\.", "_", colnames(fraction))
  fraction <- fraction[, vois_use, drop = FALSE]
  rownames(fraction) <- gsub(sample_use, name_metadata_use, rownames(fraction))
  fraction <- fraction[rownames(fraction) %in% rownames(meta_data), ]
  fraction[fraction > 0] <- 1
  fraction <- rowSums(fraction)
  fraction <- fraction > 0
  meta_data[names(fraction), "Clonal"] <- fraction
}
meta_data <- meta_data[!is.na(meta_data$Clonal), ]
meta_data$Clonal <- ifelse(meta_data$Clonal, "Clonal", "WT")
rm(fraction)


print("We save the clonal information.")
write.table(meta_data, paste0(output_use, "Clonal_Information.csv"), row.names = TRUE, col.names = TRUE, sep = ",", quote = FALSE)


print("We perform the enrichment test.")
result_enrichment_function <- EnrichmentFisherTest(InputDataFrame = meta_data, Trait1 = "Clonal", Trait2 = "CellType", Output_Path = output_use, Prefix = "Merged_Samples", remove_always_NA = FALSE, do_plot = TRUE, save_results = TRUE,
                                                   pdf_or_png = "both", fig_width = 9, fig_height = 7, fig_units = "in")


celltypes <- sort(unique(meta_data[["CellType"]]))
results <- c()
cells_per_variant_and_celltype <- matrix(NA, nrow = 1, ncol = length(celltypes), dimnames = list("Clonal", celltypes))
for(j in 1:length(celltypes)){
  celltype_use <- celltypes[j]
  print(celltype_use)
  count_matrix <- matrix(NA, nrow = 2, ncol = 2, dimnames = list(c("Clonal", "Other"), c(celltype_use, "Other")))
  count_matrix["Clonal", celltype_use] <- nrow(meta_data[meta_data[, "Clonal"] == "Clonal" & meta_data[, "CellType"] == celltype_use, ])
  count_matrix["Other",  celltype_use] <- nrow(meta_data[meta_data[, "Clonal"] != "Clonal" & meta_data[, "CellType"] == celltype_use, ])
  count_matrix["Clonal", "Other"]      <- nrow(meta_data[meta_data[, "Clonal"] == "Clonal" & meta_data[, "CellType"] != celltype_use, ])
  count_matrix["Other",  "Other"]      <- nrow(meta_data[meta_data[, "Clonal"] != "Clonal" & meta_data[, "CellType"] != celltype_use, ])
  res <- fisher.test(count_matrix)
  res <- data.frame(VOI = "Clonal", CellType = celltype_use, OddsRatio = res$estimate, P = res$p.value)
  results <- rbind(results, res)
  cells_per_variant_and_celltype["Clonal", celltype_use] <- count_matrix[1,1]
}
results[["P_adj"]] <- p.adjust(results[["P"]], method = "fdr")
results[["P_adj_log10"]] <- -log10(results[["P_adj"]])
write.table(results, paste0(output_use, "/Enrichment_Variant_CellType.csv"), sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(cells_per_variant_and_celltype, paste0(output_use, "/cells_per_variant_and_celltype.csv"), sep = ",", col.names = TRUE, row.names = TRUE, quote = FALSE)


print("We generate a heatmap for the new annotation.")
# We set the P value to 0 and the log10 value to 1, if the OddsRatio is smaller than 1.
results[["P"]] <- ifelse(results[["OddsRatio"]] < 1, 0, results[["P"]])
results[["P_adj"]] <- ifelse(results[["OddsRatio"]] < 1, 0, results[["P_adj"]])
results[["P_adj_log10"]] <- ifelse(results[["OddsRatio"]] < 1, 0, results[["P_adj_log10"]])

p <- results
rownames(p) <- NULL
p <- p |>
  dplyr::select(VOI, CellType, P_adj_log10) |>
  tidyr::pivot_wider(names_from = CellType, values_from = P_adj_log10) |>
  tibble::column_to_rownames("VOI") |>
  as.matrix()
scores <- results
rownames(scores) <- NULL
scores <- scores |>
  dplyr::select(VOI, CellType, OddsRatio) |>
  tidyr::pivot_wider(names_from = CellType, values_from = OddsRatio) |>
  tibble::column_to_rownames("VOI") |>
  as.matrix()

score_max <- ceiling(max(p, na.rm = TRUE))
if(score_max <= 0) score_max <- 1
col_fun <- circlize::colorRamp2(c(0, score_max), c("white", "red"))
p[p < -log10(0.05)] <- NA
rownames(p) <- gsub("chrM_", "", rownames(p))
rownames(p) <- sub("_(?=.$)", ">", rownames(p), perl = TRUE)
rownames(p) <- gsub("_", "", rownames(p))


legend_P <- ComplexHeatmap::Legend(
  title     = "-log10(P)",
  col_fun   = circlize::colorRamp2(c(0, ceiling(max(p, na.rm = TRUE))), c("white", "red")),
  at        = pretty(c(0, ceiling(max(p, na.rm = TRUE)))),
  labels    = pretty(c(0, ceiling(max(p, na.rm = TRUE)))),
  title_gp  = grid::gpar(fontsize = 12, fontface = "bold"),
  labels_gp = grid::gpar(fontsize = 12),
  border    = FALSE,
  background = "transparent",
)

p1 <- ComplexHeatmap::Heatmap(p, name = "-log10(P)",
  column_title = "Merged Samples",
  rect_gp = grid::gpar(col = "black", lwd = 1),
  row_title = "",
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  col = col_fun,
  row_names_gp = grid::gpar(fontsize = 15),
  column_names_gp = grid::gpar(fontsize = 15),
  column_names_rot = 45,
  column_title_gp = grid::gpar(fontsize = 20, fontface = "bold"),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_heatmap_legend = FALSE,
  bottom_annotation = NULL,
  heatmap_legend_param = list(legend_height = unit(4, "cm")),
  row_names_max_width = unit(12, "cm"),
  cell_fun = function(j, i, x, y, w, h, fill){
    if(is.na(p[i,j])){
      txt <- ""
    } else if(p[i,j] > -log10(0.001)){
      txt <- "***"
    } else if(p[i,j] > -log10(0.01)){
      txt <- "**"
    } else if(p[i,j] > -log10(0.05)){
      txt <- "*"
    } else{
      txt <- ""
    }
    grid::grid.text(txt, x, y, gp = grid::gpar(fontface = "bold"))
  }
)
png(paste0(output_use, "/Heatmap_Annotation.png"), width = 10, height = 10, units = "in", res = 600, type = "cairo")
ComplexHeatmap::draw(p1, padding = unit(c(2,3,2,2), "mm"), heatmap_legend_side = "right", annotation_legend_list = list(legend_sig, legend_P))
dev.off()
pdf(paste0(output_use, "/Heatmap_Annotation.pdf"), width = 10, height = 10)
ComplexHeatmap::draw(p1, padding = unit(c(2,3,2,2), "mm"), heatmap_legend_side = "right", annotation_legend_list = list(legend_sig, legend_P))
dev.off()


na_check <- colSums(!is.na(p)) > 0
p_subset <- p[, na_check, drop = FALSE]
p1 <- ComplexHeatmap::Heatmap(p_subset, name = "-log10(P)",
  column_title = "Merged Samples",
  rect_gp = grid::gpar(col = "black", lwd = 1),
  row_title = "",
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  col = col_fun,
  row_names_gp = grid::gpar(fontsize = 15),
  column_names_gp = grid::gpar(fontsize = 15),
  column_names_rot = 45,
  column_title_gp = grid::gpar(fontsize = 20, fontface = "bold"),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_heatmap_legend = FALSE,
  bottom_annotation = NULL,
  heatmap_legend_param = list(legend_height = unit(4, "cm")),
  row_names_max_width = unit(12, "cm"),
  cell_fun = function(j, i, x, y, w, h, fill){
    if(is.na(p_subset[i,j])){
      txt <- ""
    } else if(p_subset[i,j] > -log10(0.001)){
      txt <- "***"
    } else if(p_subset[i,j] > -log10(0.01)){
      txt <- "**"
    } else if(p_subset[i,j] > -log10(0.05)){
      txt <- "*"
    } else{
      txt <- ""
    }
    grid::grid.text(txt, x, y, gp = grid::gpar(fontface = "bold"))
  }
)
png(paste0(output_use, "/Heatmap_Annotation_Smaller_Without_Legend.png"), width = 4, height = 5, units = "in", res = 600, type = "cairo")
ComplexHeatmap::draw(p1, padding = unit(c(2,4,2,2), "mm"), heatmap_legend_side = "right")
dev.off()
pdf(paste0(output_use, "/Heatmap_Annotation_Smaller_Without_Legend.pdf"), width = 4, height = 5)
ComplexHeatmap::draw(p1, padding = unit(c(2,4,2,2), "mm"), heatmap_legend_side = "right")
dev.off()

# With legend.
p1 <- ComplexHeatmap::Heatmap(p_subset, name = "-log10(P)",
  column_title = "Combined Samples",
  rect_gp = grid::gpar(col = "black", lwd = 1),
  row_title = "",
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  col = col_fun,
  row_names_gp = grid::gpar(fontsize = 15, fontface = "bold"),
  column_names_gp = grid::gpar(fontsize = 15, fontface = "bold"),
  column_names_rot = 45,
  column_title_gp = grid::gpar(fontsize = 20, fontface = "bold"),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_heatmap_legend = FALSE,
  bottom_annotation = NULL,
  heatmap_legend_param = list(legend_height = unit(4, "cm")),
  row_names_max_width = unit(12, "cm"),
  cell_fun = function(j, i, x, y, w, h, fill){
    if(is.na(p_subset[i,j])){
      txt <- ""
    } else if(p_subset[i,j] > -log10(0.001)){
      txt <- "***"
    } else if(p_subset[i,j] > -log10(0.01)){
      txt <- "**"
    } else if(p_subset[i,j] > -log10(0.05)){
      txt <- "*"
    } else{
      txt <- ""
    }
    grid::grid.text(txt, x, y, gp = grid::gpar(fontface = "bold"))
  }
)
png(paste0(output_use, "/Heatmap_Annotation_Smaller.png"), width = 4, height = 3.5, units = "in", res = 600, type = "cairo")
ComplexHeatmap::draw(p1, padding = unit(c(2,4,2,3), "mm"), heatmap_legend_side = "right", annotation_legend_list = list(legend_sig, legend_P))
dev.off()
pdf(paste0(output_use, "/Heatmap_Annotation_Smaller.pdf"), width = 4, height = 3.5)
ComplexHeatmap::draw(p1, padding = unit(c(2,4,2,3), "mm"), heatmap_legend_side = "right", annotation_legend_list = list(legend_sig, legend_P))
dev.off()
