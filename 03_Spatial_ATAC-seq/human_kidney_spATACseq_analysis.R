##### Load required libraries -------------------------------------------------
library(Seurat)
library(Signac)
library(harmony)
library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(GenomeInfoDb)
library(Matrix)
library(dplyr)
library(Seurat)
library(patchwork)
library(Signac)
library(parallel)
library(ArchR)
library(ggsci)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SeuratDisk)
library(tidyverse)
library(purrr)
library(cowplot)
library(tidyverse)
library(purrr)
library(MOJITOO)
library(scMEGA)
library(SeuratDisk)
library(tidyverse)
library(purrr)
library(JASPAR2020)
library(TFBSTools)
library(magrittr)
library(igraph)
library(ggraph)

set.seed(1)

##### default settings ---------------------------------------------------------
addArchRGenome("hg38")
addArchRThreads(threads = 30)
options(device = function() png(type = "cairo"))

#### parameter setting
res_type <- "15mm"
peak_type <- "with_peak_ref"

##### create arrow file --------------------------------------------------------
samples_names <- c("73266", "63511", "72351", "76992", "71116", "57480", "78454")
names <- c("Healthy_73266", "Healthy_63511", "Healthy_72351", "Healthy_76992",
           "Diseased_71116", "Diseased_57480", "Diseased_78454"
          )

# -------------------
# 1. Use single-cell ATAC-seq regions from reference Muto, Y., Wilson, P.C.,
#    Ledru, N. et al. Single cell transcriptional and chromatin accessibility
#    profiling redefine cellular heterogeneity in the adult human kidney. Nat
#    Commun 12, 2190 (2021). https://doi.org/10.1038/s41467-021-22368-w
# -------------------
get_reference_peaks <- function() {
        fragments <- "1945850f-c96a-4e48-92b7-f84808c5461c-fragment.tsv.bgz"

        frags_obj <- CreateFragmentObject(path = fragments)
        gr.peaks <- CallPeaks(object = frags_obj, macs2.path = "/data/uu762107/.conda/envs/my-pipeline/bin/macs2")

        peak_counts <- FeatureMatrix(
          fragments = frags_obj,
          features = gr.peaks,
          cells = NULL # will auto-detect barcodes
        )
        write.table(rownames(mat), file = "peaks.tsv",
                    quote = FALSE, row.names = FALSE, col.names = FALSE)
}

make_archR_object <- function() {
	input_files <- c()
        for (i in 1:length(names)) {
		input_files <- append(input_files, c(file.path(frags_paths[samples_names[i]])))
        }
        names(input_files) <- names

        position_files <- c()
        for (i in 1:length(names)) {
		position_files <- append(position_files, c(file.path(pos_paths[samples_names[i]], "tissue_positions_list.csv")))
        }
        names(position_files) <- names

        genomeAnnotation <- createGenomeAnnotation(genome = "BSgenome.Hsapiens.UCSC.hg38")

        arrow_files <- createArrowFiles(
                                        inputFiles = input_files,
                                        sampleNames = names(input_files),
                                        addTileMat = FALSE,
                                        addGeneScoreMat = FALSE,
                                        minTSS = 0,
                                        minFrags = 0,
					maxFrags = Inf,
                                        force = TRUE,
                                        genomeAnnotation = genomeAnnotation
        )
        #### create an arch project ---------------------------------------------------
        project1 <- ArchRProject(
                                 ArrowFiles = arrow_files,
                                 genomeAnnotation = genomeAnnotation,
                                 outputDirectory = "ArchRProject",
                                 copyArrows = FALSE
        )

	saveArchRProject(ArchRProj = project1, outputDirectory = paste0(result_path, "human_kidney_", res_type, "_without_filtration"), load = TRUE)	
	return(project1)
}

filter_off_tissue_pixels <- function(project1) {
       position_files <- c()
        for (i in 1:length(names)) {
                position_files <- append(position_files, c(file.path(pos_paths[samples_names[i]], "tissue_positions_list.csv")))
        }
        names(position_files) <- names

        all_ontissue <- c()
        for (i in seq_along(position_files)) {
            positions <- read.csv(position_files[i], header = FALSE)
            positions$V1 <- paste(names[i], "#", positions$V1, "-1", sep = "")
            on_tissue <- positions$V1 [which(positions$V2 == 1)]
            all_ontissue <- c(all_ontissue, on_tissue)
        }

        project1 <- project1[project1$cellNames %in% all_ontissue]
        return(project1)
}

plot_qc <- function() {

        metadata <- as.data.frame(project1@cellColData)
        colours <- ArchR:::paletteDiscrete( values = unique(metadata$Sample))
        metadata$Sample <- factor(metadata$Sample, levels = names)

        P1 <- ggplot(metadata, aes(x = Sample, y = TSSEnrichment, fill = Sample)) +
		geom_violin(trim = TRUE, color = "white", show.legend = F) +
      		geom_boxplot(width = 0.1, position = position_dodge(0.9), show.legend = F,
			     fill = "white", outlier.size = 0, outlier.stroke = 0) +
		scale_fill_manual(values = colours) +
      		theme_cowplot() +
	    	theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
      		ylab("TSSEnrichment")

        P2 <- ggplot(metadata, aes(x = Sample, y = log10(nFrags), fill = Sample)) +
                geom_violin(trim = TRUE, color = "white", show.legend = F) +
                geom_boxplot(width = 0.1, position = position_dodge(0.9), show.legend = F,
			     fill = "white", outlier.size = 0, outlier.stroke = 0) +
                scale_fill_manual(values = colours) +
                theme_cowplot() +
                scale_y_continuous(limits = c(3, 5)) +
                theme(axis.text.x = element_text(colour = "black", family = "Times", size = 10)) +
                ylab(expression(Log[10] * "(" * nFragment * ")")) + xlab("")

        ggsave(paste0(figure_path, "/figures/pdf/", res_type,"/ArchR/", peak_type, "/samples_qc.pdf"),
	       do.call(plot_grid, c(list(P1, P2), ncol = 1, align = "v")), width = 8, height = 5)

        gg <- ArchR:::ggPoint(
		    	      x = pmin(log10(metadata$nFrags), 5) + rnorm(length(metadata$nFrags), sd = 0.00001),
		    	      y = metadata$TSSEnrichment + rnorm(length(metadata$nFrags), sd = 0.00001),
		    	      colorDensity = TRUE,
		    	      xlim = c(2.5, 5),
		    	      ylim = c(0, 20),
		    	      baseSize = 6,
		    	      continuousSet = "sambaNight",
		    	      xlabel = "Log 10 (Unique Fragments)",
		    	      ylabel = "TSS Enrichment",
		    	      rastr = TRUE) +
          geom_hline(yintercept = 4, lty = "dashed", size = 0.25) +
          geom_vline(xintercept = log10(1000), lty = "dashed", size = 0.25)

        ggsave(paste0(figure_path, "/figures/pdf/", res_type, "/ArchR/", peak_type, "/totalcell_qc.pdf"), gg, width = 3, height = 3)
}

plot_cor_depth <- function() {

        embed <- project1@reducedDims$IterativeLSI$matSVD
        counts <- subset(project1@cellColData, select = "nFrags")
        embed <- embed[rownames(x = counts), ]
        n = 10
        n <- Signac:::SetIfNull(x = n, y = ncol(x = embed))
        embed <- embed[, seq_len(length.out = n)]
        counts$nFrags <- as.numeric(counts$nFrags)
        depth.cor <- as.data.frame(cor(x = embed, y = counts$nFrags))
        depth.cor$Component <- rownames(depth.cor)
        colnames(depth.cor)[1] <- "correaltion"
        depth.cor$Component <- factor(depth.cor$Component, levels = depth.cor$Component)

        P3 <- ggplot(depth.cor, aes(Component, correaltion)) +
		geom_point() +
                ylab("Correlation") +
                ylim(c(-1, 1)) +
                theme_light() +
                ggtitle("Correlation between depth and LSI")

        ggsave(paste0(figure_path, "/figures/pdf/", res_type, "/ArchR/", peak_type, "/correaltion_with_depth.pdf"), P3, width = 4, height = 4)
}

plot_embeddings <- function(group = "Clusters") {
	P4 <- plotEmbedding(ArchRProj = project1, colorBy = "cellColData", name = group, embedding = "UMAP")
        ggsave(paste0(figure_path, "/figures/pdf/", res_type, "/ArchR/", peak_type, "/archR_umap_", group, "_nofilterDoublet.pdf"), P4, width = 5, height = 5)
}

plot_tss_enrichment <- function() {
        df <- as.data.frame(getCellColData(project1, select = c(
								"nFrags",
								"TSSEnrichment",
								"Sample"
								)))
        samples <- unique(df$Sample)

        df$log10_nFrags <- log10(df$nFrags)

        plots <- lapply(samples, function(s) {

                df_sub <- df[df$Sample == s, ]

                p <- ArchR:::ggPoint(
		   		     x = df_sub$log10_nFrags,
		   		     y = df_sub$TSSEnrichment,
		   		     labelMeans = TRUE,
		   		     colorDensity = TRUE,
		   		     continuousSet = "sambaNight",
		   		     xlabel = "Log10 Unique Fragments",
		   		     ylabel = "TSS Enrichment",
		   		     xlim = c(log10(500), quantile(df$log10_nFrags, probs = 0.99)),
		   		     ylim = c(0, quantile(df$TSSEnrichment, probs = 0.99))
		     		     ) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

                ggsave(paste0(figure_path, "/figures/pdf/", res_type, "/ArchR/", peak_type, "/archR_tss_enrichment", s, ".pdf"), p, width = 5, height = 5)
            })


        p <- plotFragmentSizes(ArchRProj = project1)
        ggsave(paste0(figure_path, "/figures/pdf/", res_type, "/ArchR/", peak_type, "/archR_plotFragmentSizes_profile.pdf"), p, width = 5, height = 5)


        p <- plotTSSEnrichment(ArchRProj = project1)
        ggsave(paste0(figure_path, "/figures/pdf/", res_type, "/ArchR/", peak_type, "/archR_tss_enrichment_profile.pdf"), p, width = 5, height = 5)

        cellCounts <- table(project1$Sample)
        cellCounts

        cellCountsDF <- as.data.frame(cellCounts)
        colnames(cellCountsDF) <- c("Sample", "nCells")

        p <- ggplot(cellCountsDF, aes(x = Sample, y = nCells)) +
      		geom_col() +
	    	theme_minimal() +
      		labs(title = "Number of cells per sample",
		         y = "Cells",
	     		 x = "Sample"
	       		 ) +
		theme(axis.text.x = element_text(angle = 45, hjust = 1))
        ggsave(paste0(figure_path, "/figures/pdf/", res_type, "/ArchR/", peak_type, "/archR_number_of_cells.pdf"), p, width = 5, height = 5)
}


### -------------------
### 2. Load data
### -------------------
project1 <- make_archR_object()

project1 <- loadArchRProject(paste0(result_path, "kidney7_15mm_with_ref"))
project1 <- filter_off_tissue_pixels(project1)


plot_qc()

plot_tss_enrichment()

project1 <- subsetArchRProject(
                              ArchRProj = project1,
                              cells = project1$cellNames[project1$TSSEnrichment >= 2 & project1$nFrags >= 1000 & project1$nFrags <= 10000000],
                              force = TRUE
)
project1

ido_peaks <- import(paste0(result_path, "peaks.bed"), format = "bed")

project1 <- addPeakSet(project1, peakSet = ido_peaks, force = TRUE)
project1 <- addPeakMatrix(project1)

project1 <- addGeneScoreMatrix(project1, force = TRUE)
project1


### -------------------
### 3. Dimension reduction
### -------------------
project1 <- addIterativeLSI(
                            ArchRProj = project1,
                            useMatrix = "PeakMatrix",
                            name = "IterativeLSI",
                            iterations = 4,
                            clusterParams = list(
                                                 resolution = c(0.1, 0.2, 0.4),
                                                 sampleCells = 50000,
                                                 n.start = 10
                                            ),
                            varFeatures = 125000,
                            dimsToUse = 1:30,
                            force = TRUE,
                            seed = 1
)
project1

plot_cor_depth()

### -------------------
### 4. Batch correction
### -------------------
project1 <- addHarmony(
                       ArchRProj = project1,
                       reducedDims = "IterativeLSI",
                       name = "Harmony",
                       groupBy = "Sample",
		       dimsToUse = 2:30,
                       force = TRUE,
                       seed = 1
                       )

## -------------------
## 5. Clustering
## -------------------
project1 <- addClusters(
                        input = project1,
                        reducedDims = "Harmony",
                        method = "Seurat",
                        name = "Clusters",
                        resolution = 2,
			sampleCells = 70000,
                        force = TRUE,
                        seed = 1
                        )

### -------------------
### 6. Embedding
### -------------------
project1 <- addUMAP(
                    ArchRProj = project1,
                    reducedDims = "Harmony",
                    name = "UMAP",
                    nNeighbors = 30,
                    minDist = 0.2,
                    metric = "cosine",
                    force = TRUE,
                    seed = 1
                    )

plot_embeddings("Sample")
plot_embeddings("Clusters")

### -------------------
### 7. Imputation with MAGIC
### -------------------
project1 <- addImputeWeights(
          ArchRProj = project1,
          reducedDims = "IterativeLSI",
          threads = getArchRThreads(),
          seed = 1,
          verbose = TRUE,
          logFile = createLogFile("addImputeWeights")
)

validGroups <- names(table(project1$Clusters))[table(project1$Clusters) > 6]

### -------------------
### 8. Cell type annotation
### -------------------
### Identify the marker genes for each cluster
markersGS <- getMarkerFeatures(
  ArchRProj = project1,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  testMethod = "wilcoxon",
  useGroups = validGroups,
  maxCells = 3000,
  k = 200
)

heatmapGenes <- plotMarkerHeatmap(
            seMarker = markersGS,
            cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
            transpose = TRUE
)

pdf(paste0(figure_path, "/figures/pdf/", res_type, "/ArchR/", peak_type, "/Gene-Marker-Heatmap.pdf"), width = 9, height = 6)
print(heatmapGenes)
dev.off()

fdr_thr <- 0.05
fc_thr <- 0.001

markers <- getMarkers(markersGS, cutOff = paste0("FDR <= ", fdr_thr," & Log2FC >= ", fc_thr))

markerGenes <- data.frame()
for (name in names(markers)){
    if (dim(markers[[name]])[1] > 0) {
       df <- markers[[name]]
       df$cluster <- name
       if(nrow(markerGenes) == 0){
           markerGenes <- df
       }
       else{
           markerGenes <- rbind(markerGenes, df)
       }
    }
}

rownames(markerGenes) <- 1:dim(markerGenes)[1]

top_clusters <- as.data.frame(markerGenes) %>% group_by(cluster) %>% top_n(n = 20, wt = Log2FC)
write.csv(top_clusters, paste0(figure_path, "/tables/", res_type, "/ArchR/", peak_type, "/top20_clusters_celltypes.csv"))

markersPeaks <- getMarkerFeatures(
       ArchRProj = project1,
       useMatrix = "PeakMatrix",
       groupBy = "Clusters",
       useGroups = validGroups,
       bias = c("TSSEnrichment", "log10(nFrags)"),
       testMethod = "binomial",
       binarize = TRUE
)

heatmapPeaks <- plotMarkerHeatmap(
       seMarker = markersPeaks,
       cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
       transpose = TRUE
)

pdf(paste0(figure_path, "/figures/pdf/", res_type, "/ArchR/", peak_type, "/Peak-Marker-Heatmap.pdf"), width = 9, height = 6)
print(heatmapPeaks)
dev.off()

fdr_thr <- 0.05
fc_thr <- 0.001

markers <- getMarkers(markersPeaks, cutOff = paste0("FDR <= ", fdr_thr," & Log2FC >= ", fc_thr))

markerPeaks <- data.frame()
for (name in names(markers)){
    if (dim(markers[[name]])[1] > 0) {
       df <- markers[[name]]
       df$cluster <- name
       if(nrow(markerPeaks) == 0){
           markerPeaks <- df
       }
       else{
           markerPeaks <- rbind(markerPeaks, df)
       }
    }
}

rownames(markerPeaks) <- 1:dim(markerPeaks)[1]

top_clusters <- as.data.frame(markerPeaks) %>% group_by(cluster) %>% top_n(n = 20, wt = Log2FC)
write.csv(top_clusters, paste0(figure_path, "/tables/", res_type, "/ArchR/", peak_type, "/kidney7_top20_peaks_celltypes", res, ".csv"))

new.col <- c("C1" = "PT-S2-S3",
             "C2" = "AVR Endothelial",
             "C3" = "Fibroblast",
             "C4" = "PC Endothelial",
             "C5" = "Endothelial",
             "C6" = "ATL",
             "C7" = "ATL",
             "C8" = "VSMC/Pericyte",
             "C9" = "TAL",
             "C10"= "CD",
             "C11"= "DCT",
             "C12"= "Injured PTs",
             "C13"= "PT-S1-S2",
             "C14"= "PT-S2-S3",
             "C15"= "PT-S1",
             "C16"= "Immune",
             "C17"= "PV/Cap fibroblasts",
             "C18"= "Pod",
             "C19"= "VSMC/Pericyte",
             "C20"= "Immune",
             "C21"= "Fibroblast"
             )

project1$annotation <- mapLabels(project1$Clusters, newLabels = new.col)

disc_pal <- c(
    "Unannotated"        = "#E0E0E0",
    "PT-S1-S2"           = "#006D2C",
    "PT-S2-S3"           = "#31A354",
    "PT-S1"              = "#74C476",
    "Injured PTs"        = "#FFD700",
    "ATL"                = "#9ECAE1",
    "TAL"                = "#4292C6",
    "DCT"                = "#084594",
    "CD"                 = "#6A51A3",
    "Pod"                = "#E31A1C",
    "Endothelial"        = "#FF7F00",
    "PC Endothelial"     = "#FDAE6B",
    "AVR Endothelial"    = "#D94801",
    "VSMC/Pericyte"      = "#FDBF6F",
    "Fibroblast"         = "#8C564B",
    "PV/Cap fibroblasts" = "#C49C94",
    "Immune"             = "#F781BF",
    "other"              = "#f0f4f9"
    )

P5 <- plotEmbedding(
    ArchRProj = project1,
    embedding = "UMAP",
    colorBy = "cellColData",
    name = "annotation",
    pal = disc_pal,
    sampleCells = NULL,
    plotAs = "points",
    labelMeans = FALSE
    )

ggsave(paste0(figure_path, "/figures/pdf/", res_type, "/ArchR/", peak_type, "/archR_umap_celltypes.pdf"), P5,
       width = 5, height = 5)

markersGS <- getMarkerFeatures(
  ArchRProj = project1,
  useMatrix = "GeneScoreMatrix",
  groupBy = "annotation",
  testMethod = "wilcoxon",
  maxCells = 3000,
  k = 200
)

log2fc <- assay(markersGS, "Log2FC")

top2 <- lapply(seq_len(ncol(log2fc)), function(i) {
    x <- log2fc[, i]
    idx <- order(x, decreasing = TRUE, na.last = NA)[seq_len(min(2, sum(!is.na(x))))]

    data.frame(
        Feature = rownames(markersGS)[idx],
        Log2FC = x[idx],
        Column = colnames(markersGS)[i],
        row.names = NULL
    )
})

top2 <- do.call(rbind, top2)

subsetSE <- markersGS[which(rowData(markersGS)$name %in% top2$Feature),]

p <- plotMarkerHeatmap(seMarker = subsetSE,
                      cutOff = "FDR <= 1 & Log2FC >= 0",
                      labelMarkers = top2$Feature
                      )

pdf(paste0(figure_path, "/figures/pdf/", res_type, "/ArchR/", peak_type, "/Gene-Marker-Heatmap_selected_genes_celltypes.pdf"), width = 6, height = 9)
print(p)
dev.off()

###--------------------
### 9. Convert ArchR to Seurat
###--------------------
# Peak matrix (binary accessibility matrix)
peak_matrix <- getMatrixFromProject(project1, useMatrix = "PeakMatrix")
peak_matrix

# Cell metadata
cell_metadata <- getCellColData(project1)

# Get peak annotations
peak_gr <- rowRanges(peak_matrix)

# Convert to sparse matrix
peak_counts <- assay(peak_matrix)

# Row names as peak strings
peak_names <- paste0(seqnames(peak_gr), ":", start(peak_gr), "-", end(peak_gr))
rownames(peak_counts) <- peak_names

# Cell metadata to data.frame
cell_metadata <- as.data.frame(cell_metadata)

chrom_assay <- CreateChromatinAssay(
  counts = peak_counts,
  sep = c(":", "-"),
  ranges = peak_gr
)

obj <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = cell_metadata
)

###--------------------
### 10. Compute TF activity
###--------------------
JASPAR <- JASPAR()
JASPARConnect <- RSQLite::dbConnect(RSQLite::SQLite(), db(JASPAR))

opts <- list(
   species = c("Homo sapiens", "Mus musculus"),
   collection = "CORE",
   tax_group = "vertebrates",
   all_versions = FALSE
)

pfm <- TFBSTools::getMatrixSet(JASPARConnect,  opts)

# add motif information
obj <- AddMotifs(
	       	 object = obj,
	       	 genome = BSgenome.Hsapiens.UCSC.hg38,
	       	 pfm = pfm,
	       	 assay = "peaks"
)

obj <- RunChromVAR(
                   object = obj,
                   genome = BSgenome.Hsapiens.UCSC.hg38,
                   assay = "peaks"
)

###--------------------
### 11. Trajectory analysis
###--------------------
obj <- setDimRed(obj, getReducedDims(project1, "Harmony"), "atacHarmony", key="lsi")
embedd <- Embeddings(obj[["atacHarmony"]])

obj <- RunDiffusionMap(obj,
 		       reduction = "atacHarmony"
  		       )
dm <- DiffusionMap(as.matrix(embedd))

obj <- AddTrajectory(object = obj,
                     trajectory = c("PT-S1", "PT-S2", "PT-S3", "Injured PTs"),
                     group.by = "annotation",
                     reduction = "dm",
                     dims = 1:3,
                     use.all = FALSE)

# we only plot the cells that are in this trajectory
obj <- obj[, !is.na(obj$Trajectory)]

p <- TrajectoryPlot(object = obj,
                    reduction = "dm",
                    continuousSet = "blueYellow",
                    size = 1,
                    addArrow = FALSE)

pdf(paste0(result_path, "/figures/pdf/", res_type, "/ArchR/", peak_type, "/scMEGA_trajectory_analysis.pdf"), height = 5, width = 5)
print(p)
dev.off()

p <- DimPlot(object = obj,
            group.by = "annotation",
            reduction = "dm",
	    label = TRUE) + NoLegend()

pdf(paste0(result_path, "/figures/pdf/", res_type, "/ArchR/", peak_type, "/PT_scMEGA_trajectory_analysis.pdf"), height = 5, width = 5)
print(p)
dev.off()

atac.assay = "peaks"
tf.assay = "chromvar"
trajectory.name = "Trajectory"
groupEvery = 1

trajMM <- GetTrajectory(
    obj,
    assay = tf.assay,
    trajectory.name = trajectory.name,
    groupEvery = groupEvery,
    slot = "data",
    smoothWindow = 7,
    log2Norm = FALSE
  )

rownames(trajMM) <- obj@assays[[atac.assay]]@motifs@motif.names

ht <- TrajectoryHeatmap(
    trajMM,
    pal = paletteContinuous(set = "solarExtra"),
    varCutOff = 0,
    maxFeatures = nrow(trajMM),
    limits = c(-2, 2),
    labelRows = TRUE,
    labelTop = 50,
    name = "TF activity"
  )

pdf(paste0(result_path, "/figures/pdf/", res_type, "/ArchR/", peak_type, "/PT_scMEGA_heatmap.pdf"), height = 7, width = 5)
print(ht)
dev.off()

for(tf.use in rownames(trajMM)){
        df.tf.activity <- assay(trajMM)
        df.tf.activity <- t(scale(t(df.tf.activity)))
        df.tf.activity <- as.data.frame(df.tf.activity)
        colnames(df.tf.activity) <- seq(0, 100, groupEvery)[-1]
        df.tf.activity$tf <- toupper(rownames(df.tf.activity))
        df.tf.activity <- tidyr::pivot_longer(df.tf.activity, -tf,
                                            names_to = "pseudotime",
                                            values_to = "value")
        df.tf.activity$pseudotime <- as.numeric(df.tf.activity$pseudotime)

        df.tf.activity$data <- "TF activity"
        df.plot <- subset(df.tf.activity, tf == tf.use)

        p <- ggplot(df.plot, aes(x = pseudotime, y = value, color = data)) +
        geom_smooth(method = "loess",se = FALSE) +
        ggtitle(tf.use) +
        cowplot::theme_cowplot() +
        ylab("") +
        theme(legend.title = element_blank())

        ggsave(p, filename = paste0(result_path, "/figures/pdf/", res_type, "/ArchR/", peak_type, "/TFs/trajectory_", tf.use, ".pdf"),
              bg = "transparent", width = 4, height = 3)
}


