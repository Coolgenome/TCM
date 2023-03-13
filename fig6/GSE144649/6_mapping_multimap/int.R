#--------------------------------------------------------------
# filename : int.R
# Date : 2022-04-29
# contributor : Yanshuo Chu
# function: int
#--------------------------------------------------------------

print('<==== int.R ====>')

suppressMessages({
    library(optparse)
    library(tidyverse)
    library(Seurat)
    library(SeuratObject)
    library(cowplot)
    library(MultiMap)
})

CD4PredictedT <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE179994/MappingResult_MultiMap/CD4/predictedT_2022-04-29.rds")
