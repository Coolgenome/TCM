#--------------------------------------------------------------
# filename : extract.R
# Date : 2022-09-01
# contributor : Yanshuo Chu
# function: extract
#--------------------------------------------------------------

print('<==== extract.R ====>')

library(optparse)
library(tidyverse)
library(Seurat)
library(GEOquery)

seuratObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/b1_harmony/2_run_harmony/nPC_30/UMAP_dist_0.1_nneighbor_50/GSE169246_UMAP_dist_0.1_nneighbor_50_CLUSTER_res_0.3/cluster.rds")

TCellClusters <- c(0, 1, 2, 3, 7, 9, 13, 16, 20)

Idents(seuratObj) <- seuratObj$seurat_clusters

subObj <- subset(seuratObj, idents = TCellClusters)

figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/subT1_extract_from_b2/1_extract/outs")
if(!dir.exists(figurePath)){
    dir.create(figurePath, recursive = T)
}
setwd(figurePath)
saveRDS(subObj, file.path(getwd(), paste0('subObj', "_", Sys.Date(), '.rds')))
