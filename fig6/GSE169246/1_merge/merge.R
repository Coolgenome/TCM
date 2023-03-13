#--------------------------------------------------------------
# filename : merge.R
# Date : 2022-02-16
# contributor : Yanshuo Chu
# function: merge
#--------------------------------------------------------------

print('<==== merge.R ====>')
rm(list=ls())

library(data.table)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(harmony)

## seuratObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/merged/merged.obj")

data <- Read10X(data.dir = "/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/data/GSE169246/raw/RNA", gene.column = 1)
seuratObj = CreateSeuratObject(counts = data)

seuratObj <- seuratObj %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = 100, verbose = FALSE)

figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result", "GSE169246", "merged")
if(!dir.exists(figurePath)){
    dir.create(figurePath, recursive = T)
}
setwd(figurePath)
saveRDS(seuratObj, "merged.obj")

