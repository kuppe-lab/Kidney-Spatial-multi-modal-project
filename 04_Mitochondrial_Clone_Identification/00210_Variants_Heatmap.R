"""
@author: Martin Grasshoff
Environment: Spatial_ATAC_MT_Genotyping_Kuppe_Lab

We plot the mtVars as a heatmap.

Execute the script like this:
python 00200_Loading_Processing.py <interger_index>

The script needs as input:
- the SVVs from 00200
- the fraction matrix
- positions

Outputs is a heatmap of the SVVs.

"""

print("Libraries.")
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(grid)
})


print("Variables.")
samples <- c("57480", "57480_flowgel", "58152", "61968", "63055", "63079", "63409", "63511_reseq", "71116", "72351_reseq", "73266", "76992_2", "78454", "78454_flowgel")
positions_paths <- read.csv("DesignMatrix_VisiumHD.csv")
positions_paths <- setNames(positions_paths$positions, positions_paths$sample)

for(sample_use in samples){
  print(sample_use)
  input <- paste0("MinVAF_0.10_MinAltReads_2_Pruning/", sample_use, "/")
  output <- paste0("MinVAF_0.10_MinAltReads_2_Pruning/", sample_use, "/")
  position_use <- positions_paths[[sample_use]]
  vois <- read.table(paste0(input, "SpatiallyVariableVariants.csv"), header=TRUE, sep=",")[, 1]
  vois <- gsub(">", "_", vois)

  print("We load the enhanced FTU results and other meta data.")
  positions <- read.table(position_use, header=FALSE, row.names=1, sep=",")
  rownames(positions) <- paste0("P", sample_use, "_", rownames(positions), "-1")
  positions <- positions[positions[, 1] == 1, ]
  positions <- positions[, c(5, 4)]
  colnames(positions) <- c("X", "Y")

  print("We load the fraction.")
  fraction <- read.table(paste0(input, "/fraction_matrix.csv"), sep=",", row.names=1, header=TRUE)
  rownames(fraction) <- paste0("P", rownames(fraction))
  colnames(fraction) <- gsub("\\.", "_", colnames(fraction))
  spots <- intersect(rownames(positions), rownames(fraction))
  fraction <- fraction[spots, ]


  print("Generating the heatmap for all VOIs.")
  fraction_use <- t(fraction[, vois, drop = FALSE])
  fraction_use <- fraction_use[, colSums(fraction_use) > 0, drop = FALSE]
  binary_matrix <- fraction_use > 0

  # We order the variants. First, the most numerous variant, then the rest in 
  # decreasing order.
  spots_per_variant <- rowSums(binary_matrix)
  variant_order <- order(-spots_per_variant, rownames(binary_matrix))
  binary_matrix <- binary_matrix[variant_order, , drop = FALSE]
  fraction_use <- fraction_use[variant_order, , drop = FALSE]

  # Second, we order the spots by their binary values.
  sort_keys_per_variant <- list()
  for (variant_index in seq_len(nrow(binary_matrix))) {
    sort_keys_per_variant[[variant_index]] <- -as.integer(binary_matrix[variant_index, ])
  }
  spot_order <- do.call(order, sort_keys_per_variant)

  fraction_reordered <- fraction_use[, spot_order, drop = FALSE]
  fraction_reordered[fraction_reordered > 0] <- 1 # We set all values to one to improve the plot.

  hm <- ComplexHeatmap::Heatmap(fraction_reordered,
    name = "VAF",
    row_title_gp = grid::gpar(fontsize=16, fontface="bold"),
    row_names_gp = grid::gpar(fontsize=16, fontface="bold"),
    row_names_max_width = unit(32, "cm"),
    column_title = "MT Variants",
    column_title_gp = grid::gpar(fontsize=16, fontface="bold"),
    column_names_gp = grid::gpar(fontsize=16, fontface="bold"),
    col = circlize::colorRamp2(seq(0, max(fraction_reordered, na.rm=TRUE), length.out=9), c("#FCFCFC", "#FFEDB0", "#FFDF5F", "#FEC510", "#FA8E24", "#F14C2B", "#DA2828", "#BE2222", "#A31D1D")),
    show_column_dend = FALSE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    heatmap_legend_param = list(border="#000000", grid_height=unit(10, "mm"), direction="horizontal"),
    # bottom_annotation = ha,
    border=TRUE,
    use_raster = TRUE,
    show_heatmap_legend = FALSE,)
  png(paste0(output, "Heatmap_FTU_reordered_Enhanced.png"), width=7, height=6, units="in", res=300, type="cairo")
  ComplexHeatmap::draw(hm, padding = unit(c(2, 2, 2, 5), "mm"))
  dev.off()
  pdf(paste0(output, "Heatmap_FTU_reordered_Enhanced.pdf"), width=7, height=6)
  ComplexHeatmap::draw(hm, padding = unit(c(2, 2, 2, 5), "mm"))
  dev.off()
}
