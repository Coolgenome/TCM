#'--------------------------------------------------------------
#' filename : extractData.R
#' Date : 2021-07-18
#' contributor : Yanshuo Chu
#' function: extractData
#'--------------------------------------------------------------

print('<==== extractData.R ====>')

library(optparse)
library(tidyverse)
library(Seurat)

NKGDT_PATH <- "/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/validate/NKGDT_V5/nPC_10/UMAP_dist_0.1_nneighbor_35/p1NKGDT_V4_10_UMAP_dist_0.1_nneighbor_35_CLUSTER_res_0.3/cluster.rds"

seuratObj <- readRDS(NKGDT_PATH)


Idents(seuratObj) <- seuratObj$seurat_clusters


subSeuratObj <- subset(seuratObj, idents = c(0, 2, 5), invert = T)

saveRDS(subSeuratObj, file.path("/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/validate/NKTMAIT_V6/data.rds"))
