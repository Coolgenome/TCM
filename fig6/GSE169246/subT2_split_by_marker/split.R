#--------------------------------------------------------------
# filename : split.R
# Date : 2022-09-01
# contributor : Yanshuo Chu
# function: split
#--------------------------------------------------------------

print('<==== split.R ====>')

rm(list=ls())
library(tidyverse)
library(Seurat)

figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/subT2_split_by_marker/outs")
if(!dir.exists(figurePath)){
    dir.create(figurePath, recursive = T)
}
setwd(figurePath)

seuratObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/subT1_extract_from_b2/1_extract/outs/subObj_2022-09-01.rds")

seuratObj <- Seurat::FindVariableFeatures(seuratObj)

hgs <- VariableFeatures(seuratObj)
hgs <- union(hgs, c("CD3D", "CD3E", "CD4", "CD8A", "CD8B"))

seuratObj <- ScaleData(seuratObj, features = hgs)

targetMatrix <- seuratObj@assays$RNA@scale.data[c("CD4", "CD8A", "CD8B"),]
targetMatrix <- as.matrix(targetMatrix)
targetMatrix[2,] <- apply(targetMatrix[2:3,], 2, max)
targetMatrix <- targetMatrix[1:2,]
targetMatrix <- t(targetMatrix)


targetTibble <- as_tibble(targetMatrix) %>%
    mutate(barcode = rownames(targetMatrix))

threhold <- 0.1527
targetTibble$CellType <- "Else"
targetTibble$CellType[targetTibble$CD4 - targetTibble$CD8A > threhold] <- "CD4"
targetTibble$CellType[targetTibble$CD4 - targetTibble$CD8A < -threhold] <- "CD8"
table(targetTibble$CellType)

CD4Obj <- subset(seuratObj, cells = targetTibble$barcode[targetTibble$CellType == "CD4"])
CD8Obj <- subset(seuratObj, cells = targetTibble$barcode[targetTibble$CellType == "CD8"])

## saveRDS(CD4Obj, file.path(getwd(), paste0('CD4Obj', "_", Sys.Date(), '.rds')))
saveRDS(CD8Obj, file.path(getwd(), paste0('CD8Obj', "_", Sys.Date(), '.rds')))
