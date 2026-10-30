"""
@author: Martin Grasshoff
Environment: Spatial_ATAC_MT_Genotyping_Kuppe_Lab

Performing the enrichment analysis for 73266
"""
print("Libraries.")
suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(tibble)
  library(grid)
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
input <- "MinVAF_0.10_MinAltReads_2_Pruning/"
output <- "MinVAF_0.10_MinAltReads_2_Pruning/"
sample_use <- "73266"
input_use <- paste0(input, sample_use, "/")
output_use <- paste0(output, sample_use, "/Enrichment_Variant_PT_CellTypes_Separate/")
vois <- c("chrM_15612_G_A", "chrM_3424_G_A", "chrM_4849_G_A", "chrM_16188_C_T")
dir.create(output_use, showWarnings = FALSE, recursive = TRUE)


print("We load the meta data.")
meta_data <- read.table("/beegfs/data/CostaLab/Spatial_ATAC_Analysis_KuppeLab/data/processed/scATAC/r_objects/kidney7_15mm_with_ref_2/meta_data.csv", sep = ",", header = TRUE, row.names = 1)
meta_data[["Sample"]] <- gsub("Diseased_|Healthy_", "", meta_data[["Sample"]])
rownames(meta_data) <- gsub("Diseased_|Healthy_", "", rownames(meta_data))
meta_data[["annotation"]] <- gsub("-|/| ", "_", meta_data[["annotation"]])
meta_data[["CellType"]] <- meta_data[["annotation"]]

print("We load the fraction.")
fraction <- read.table(paste0(input, sample_use, "/fraction_matrix.csv"), sep = ",", row.names = 1, header = TRUE)
colnames(fraction) <- gsub("\\.", "_", colnames(fraction))
rownames(fraction) <- gsub("_flowgel", "", rownames(fraction))
fraction[fraction > 0] <- 1
meta_data <- meta_data[rownames(meta_data) %in% rownames(fraction), ]
fraction <- fraction[rownames(meta_data), ]
meta_data[, vois] <- fraction[, vois]
for(i in seq_along(vois)){
  voi_use <- vois[i]
  meta_data[, voi_use] <- ifelse(meta_data[, voi_use] == 1, "Alt", "Ref")
}


print("Merging the variants.")
meta_data[["Merged1"]] <- paste0(meta_data[, "chrM_4849_G_A"], "_", meta_data[, "chrM_16188_C_T"])
meta_data[["Merged1"]] <- ifelse(meta_data[["Merged1"]] %in% c("Alt_Alt", "Alt_Ref", "Ref_Alt"), "Alt", meta_data[["Merged1"]])
meta_data[["Merged1"]] <- ifelse(meta_data[["Merged1"]] == "Ref_Ref", "Ref", meta_data[["Merged1"]])

meta_data[, c("chrM_4849_G_A", "chrM_16188_C_T")] <- NULL


print("We perform the enrichment test.")
vois <- c("chrM_15612_G_A", "chrM_3424_G_A", "Merged1")
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
  title_gp  = grid::gpar(fontsize = 12, fontface = "plain"),
  labels_gp = grid::gpar(fontsize = 12),
  border    = FALSE,
  background = "transparent",
)
rownames(p)[3] <- "chrM_4849_G_A+\chrM_16188_C_T"
p1 <- ComplexHeatmap::Heatmap(p, name = "-log10(P)",
  column_title = sample_use,
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
  heatmap_legend_param = list(legend_height = grid::unit(4, "cm")),
  row_names_max_width = grid::unit(12, "cm"),
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
ComplexHeatmap::draw(p1, padding = grid::unit(c(2,3,2,2), "mm"), heatmap_legend_side = "right", annotation_legend_list = list(legend_sig, legend_P))
dev.off()
pdf(paste0(output_use, "/Heatmap_Annotation.pdf"), width = 10, height = 10)
ComplexHeatmap::draw(p1, padding = grid::unit(c(2,3,2,2), "mm"), heatmap_legend_side = "right", annotation_legend_list = list(legend_sig, legend_P))
dev.off()


na_check <- colSums(!is.na(p)) > 0
p_subset <- p[, na_check, drop = FALSE]
p1 <- ComplexHeatmap::Heatmap(p_subset, name = "-log10(P)",
  column_title = sample_use,
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
  heatmap_legend_param = list(legend_height = grid::unit(4, "cm")),
  row_names_max_width = grid::unit(12, "cm"),
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
ComplexHeatmap::draw(p1, padding = grid::unit(c(2,4,2,2), "mm"), heatmap_legend_side = "right")
dev.off()
pdf(paste0(output_use, "/Heatmap_Annotation_Smaller_Without_Legend.pdf"), width = 4, height = 5)
ComplexHeatmap::draw(p1, padding = grid::unit(c(2,4,2,2), "mm"), heatmap_legend_side = "right")
dev.off()

# With legend.
p1 <- ComplexHeatmap::Heatmap(p_subset, name = "-log10(P)",
  column_title = sample_use,
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
  heatmap_legend_param = list(legend_height = grid::unit(4, "cm")),
  row_names_max_width = grid::unit(12, "cm"),
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
png(paste0(output_use, "/Heatmap_Annotation_Smaller.png"), width = 4, height = 5, units = "in", res = 600, type = "cairo")
ComplexHeatmap::draw(p1, padding = grid::unit(c(2,4,2,3), "mm"), heatmap_legend_side = "right", annotation_legend_list = list(legend_sig, legend_P))
dev.off()
pdf(paste0(output_use, "/Heatmap_Annotation_Smaller.pdf"), width = 4, height = 5)
ComplexHeatmap::draw(p1, padding = grid::unit(c(2,4,2,3), "mm"), heatmap_legend_side = "right", annotation_legend_list = list(legend_sig, legend_P))
dev.off()
