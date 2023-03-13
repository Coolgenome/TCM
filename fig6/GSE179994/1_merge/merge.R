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

rcs <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/data/GSE179994/raw/GSE179994_all.Tcell.rawCounts.rds")
meta.data <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/data/GSE179994/raw/GSE179994_Tcell.metadata.tsv") %>%
    replace_na(list(celltype = 'NA', cluster = 'NA'))
rownames(meta.data) <- meta.data$cellid
seuratObj <- CreateSeuratObject(counts = rcs, meta.data = as.data.frame(meta.data))

seuratObj <- seuratObj %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(npcs = 20, verbose = FALSE)

figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result", "GSE179994", "merged")
if(!dir.exists(figurePath)){
    dir.create(figurePath, recursive = T)
}
setwd(figurePath)

saveRDS(seuratObj, "merged.obj")


