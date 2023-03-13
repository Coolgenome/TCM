#--------------------------------------------------------------
# filename : extract.R
# Date : 2022-05-06
# contributor : Yanshuo Chu
# function: extract
#--------------------------------------------------------------

print('<==== extract.R ====>')
rm(list=ls())

suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(ggplot2)
})

figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE144649/5_extractT")
if(!dir.exists(figurePath)){
    dir.create(figurePath, recursive = T)
}
setwd(figurePath)


seuratObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE144649/3_pca/nPC_50/UMAP_dist_0.1_nneighbor_50/GSE144649_UMAP_dist_0.1_nneighbor_50_CLUSTER_res_0.3/cluster.rds")

Idents(seuratObj) <- seuratObj$seurat_clusters

CD4_Clusters <- c( 0, 2, 4 )
CD8_Clusters <- c( 1, 3, 5, 17, 18 )

CD4_Obj <- subset(seuratObj, idents = CD4_Clusters)

saveRDS(CD4_Obj, file.path(getwd(), paste0('CD4.rds')))

CD8_Obj <- subset(seuratObj, idents = CD8_Clusters)

saveRDS(CD8_Obj, file.path(getwd(), paste0('CD8.rds')))
