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

figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE173351/6_extractT")
if(!dir.exists(figurePath)){
    dir.create(figurePath, recursive = T)
}
setwd(figurePath)

seuratObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE173351/4_harmony/nPC_30/UMAP_dist_0.1_nneighbor_50/GSE173351_UMAP_dist_0.1_nneighbor_50_CLUSTER_res_0.3/cluster.rds")

Idents(seuratObj) <- seuratObj$seurat_clusters

## mait_markers <- c("LTB", "KLRB1", "IL7R", "GZMK", "TRAV1-2", "SLC4A10")
## pdf(file.path(getwd(), "mait_bubbleplot.pdf"))
## DotPlot(seuratObj, features = mait_markers)
## dev.off()

CD4_Clusters <- c( 0, 5, 13, 11)
CD8_Clusters <- c( 1, 3, 4, 8)

## 2 MAIT
## 6 NKT
## 10 PROLIFERATIVE
## 12 firboblast
## 7 DC
## 9 B
## 3 CD4+CD8+

CD4_Obj <- subset(seuratObj, idents = CD4_Clusters)
saveRDS(CD4_Obj, file.path(getwd(), paste0('CD4.rds')))

CD8_Obj <- subset(seuratObj, idents = CD8_Clusters)
saveRDS(CD8_Obj, file.path(getwd(), paste0('CD8.rds')))
