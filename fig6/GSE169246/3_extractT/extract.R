#--------------------------------------------------------------
# filename : extract.R
# Date : 2022-03-04
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

TCellsFolder <- "/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/TCells"
if(!dir.exists(TCellsFolder)){
    dir.create(TCellsFolder)
}
setwd(TCellsFolder)

seuratObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/merged/nPC_30/UMAP_dist_0.1_nneighbor_50/GSE169246_UMAP_dist_0.1_nneighbor_50_CLUSTER_res_0.3/cluster.rds")

seuratObj$sample <- stringr::str_extract(pattern = "(?<=[ATGC]{16}.).*", Cells(seuratObj))
Idents(seuratObj) <- seuratObj$sample
pdf(file.path(getwd(), "orig_ident.pdf"))
DimPlot(seuratObj, label = F) +
  theme(legend.position = "none")
dev.off()


Idents(seuratObj) <- seuratObj$seurat_clusters
seuratObj_P <- subset(seuratObj, idents = 12)
saveRDS(seuratObj_P, file.path(paste0("P", "_", Sys.Date(), ".rds")))

TCellClusters <- c(0,1,2,8,14,16,21,23)
seuratObj_T <- subset(seuratObj, idents = TCellClusters)
saveRDS(seuratObj_T, file.path(paste0('T', "_", Sys.Date(), '.rds')))

CD4EXP <- seuratObj@assays$RNA@data["CD4",]
CD8EXP <- matrixStats::colMaxs(as.matrix(seuratObj@assays$RNA@data[c("CD8A", "CD8B"), ]))
seuratObj_T$cell.type <- "CD4"
seuratObj_T$cell.type[CD4EXP < CD8EXP] <- "CD8"

Idents(seuratObj_T) <- seuratObj_T$cell.type
pdf(file.path(getwd(), "T_bubbleplot.pdf"))
DotPlot(seuratObj_T, features = c("CD4", "CD8A", "CD8B", "FOXP3"))
dev.off()

CD4SeuratObj <- subset(seuratObj_T, idents = "CD4")
saveRDS(CD4SeuratObj, file.path(paste0("CD4T", "_", Sys.Date(), ".rds")))

CD8SeuratObj <- subset(seuratObj_T, idents = "CD8")
saveRDS(CD8SeuratObj, file.path(paste0("CD8T", "_", Sys.Date(), ".rds")))

## B plasma:
## 5
## 6
## 9
## 13
## 20

## Monocyte:
## 4
## 22
## 7
## 10
## 15

## unknown:
##   17
## 19
## 20
## 22

## mast:
## 18

## proliferative T:
## 12

## proliferative B:
## 13

## Treg:
##   8

## NK/NKT:
##   11
##   3


