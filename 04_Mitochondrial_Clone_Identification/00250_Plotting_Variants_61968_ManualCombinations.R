"""
Plotting the clonal definitions for 61968
Environment: Spatial_ATAC_MT_Genotyping_Kuppe_Lab
"""
print("Libraries.")
suppressPackageStartupMessages({
  library(ggplot2)
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
colours_voi <- c("14180_T_C" = "#1f77b4",
                 "15192_T_C" = "#9467bd",
                 "15132_T_C" = "#8c564b",
                 "12084_C_T" = "#7f7f7f",
                 "WT" = "#d3d3d3")
break_names <- c("14180T>C", "15192T>C", "15132T>C", "12084C>T", "WT")


print("We get the used sample.")
sample_use <- "61968"
condition_use <- "Non-Fibrotic"
print(sample_use)
input_use <- paste0(input, sample_use, "/")
output_use <- paste0(output, sample_use, "/")
dir.create(output_use, showWarnings = FALSE, recursive = TRUE)


print("We load the meta data.")
meta_data <- read.table(paste0(input_use, sample_use, "_meta_data_with_clones.csv"), sep = ",", header = TRUE, row.names = 1)


print("Adding positions.")
positions_paths <- read.csv("DesignMatrix_VisiumHD.csv")
positions_paths <- setNames(positions_paths$positions, positions_paths$sample)
position_use <- positions_paths[[sample_use]]
positions <- read.table(position_use, header=FALSE, row.names=1, sep=",")

rownames(positions) <- paste0(sample_use, "_", rownames(positions), "-1")
positions <- positions[positions[, 1] == 1,, drop = FALSE]
positions <- positions[, c(5, 4)]
colnames(positions) <- c("X", "Y")
positions <- positions[rownames(positions) %in% rownames(meta_data),, drop = FALSE]
meta_data <- meta_data[rownames(meta_data) %in% rownames(positions),, drop = FALSE]
positions <- positions[rownames(meta_data),, drop = FALSE]
meta_data <- cbind(meta_data, positions)


print("Assigning colors.")
meta_data[["plot_colour"]] <- colours_voi[meta_data[["clone_manual"]]]
meta_data <- rbind(subset(meta_data, clone_manual == "WT"), subset(meta_data, clone_manual != "WT"))

p <- ggplot2::ggplot(meta_data, ggplot2::aes(x = X, y = -Y, colour = plot_colour)) +
  ggplot2::geom_point() +
  ggplot2::ggtitle(paste0(condition_use, " ", sample_use)) +
  ggplot2::scale_colour_identity(
    name   = "Variants", 
    guide  = "legend",
    labels = break_names,
    breaks = unname(colours_voi)) +
  ggplot2::theme_void() +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 6))) +
  ggplot2::theme(plot.title          = ggplot2::element_text(size = 22, face = "bold", hjust = 0.5),
                 plot.title.position = "plot",
                 legend.title        = ggplot2::element_text(size = 20, face = "bold"),
                 legend.text         = ggplot2::element_text(size = 16, face = "bold"))
ggplot2::ggsave(paste0(output_use, "VOIs_scatter_ManualCombinations.png"), p, width = 9.25, height = 8, dpi = 300)
ggplot2::ggsave(paste0(output_use, "VOIs_scatter_ManualCombinations.pdf"), p, width = 9.25, height = 8, dpi = 300)
p <- p + ggplot2::theme_void() + ggplot2::theme(plot.title = ggplot2::element_blank(), legend.position = "none")
ggplot2::ggsave(paste0(output_use, "VOIs_scatter_ManualCombinations_void.png"), p, width = 8, height = 8, dpi = 300)
ggplot2::ggsave(paste0(output_use, "VOIs_scatter_ManualCombinations_void.pdf"), p, width = 8, height = 8, dpi = 300)
