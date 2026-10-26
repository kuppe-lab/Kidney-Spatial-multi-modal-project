The Flex-based single cell RNA-seq raw fastq files were processed with the Costalab scrna_seurat_pipeline (https://github.com/CostaLab/scrna_seurat_pipeline.git). We ran the pipeline in a local singularity container. The parameters are shown in run_example.sh.

After finishing the above pipeline, we first annotate the scRNA-seq results and then integrate the data with Spatial ATAC-seq data for annotation.
