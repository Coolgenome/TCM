#--------------------------------------------------------------
# filename : inject.R
# Date : 2022-08-31
# contributor : Yanshuo Chu
# function: inject
#--------------------------------------------------------------

print('<==== inject.R ====>')

library(optparse)
library(tidyverse)
library(Seurat)
library(GEOquery)

seuratObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/merged/merged.obj")

TotalAntiPDL1ChemoPatients <- c("P019",
                                "P010",
                                "P012",
                                "P007",
                                "P017",
                                "P001",
                                "P002",
                                "P014",
                                "P004",
                                "P005",
                                "P016")

TotalChemoPatients <- c("P022",
                        "P011",
                        "P020",
                        "P008",
                        "P013",
                        "P025",
                        "P018",
                        "P023",
                        "P024",
                        "P003",
                        "P028")

AllResponsePatients <- c("P019",
                         "P010",
                         "P012",
                         "P007",
                         "P022",
                         "P011",
                         "P020",
                         "P008",
                         "P013")

seuratObj@meta.data$Sample <- stringr::str_extract(Cells(seuratObj), "(?<=^.{10,20}\\.).+")
seuratObj@meta.data$Patient <- stringr::str_extract(seuratObj@meta.data$Sample, "P\\d+")
seuratObj@meta.data$Tissue <- stringr::str_extract(seuratObj@meta.data$Sample, "\\w$")
seuratObj@meta.data$TumorTreatment <- stringr::str_extract(seuratObj@meta.data$Sample, "^[a-zA-Z]+")
seuratObj@meta.data$isResponse <- "NR"
seuratObj@meta.data$isResponse[seuratObj@meta.data$Patient %in% AllResponsePatients] <- "R"
seuratObj@meta.data$isResponse[seuratObj@meta.data$Patient == "P028"] <- "-"
seuratObj@meta.data$TreatmentType <- "PDL1+Chemo"
seuratObj@meta.data$TreatmentType[seuratObj@meta.data$Patient %in% TotalChemoPatients] <- "Chemo"
seuratObj@meta.data$group <- paste0(seuratObj@meta.data$TreatmentType, "-", seuratObj@meta.data$TumorTreatment, "-", seuratObj@meta.data$isResponse)

seuratObj$batch <- seuratObj$Patient

figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/b1_harmony/1_injectBatchinfo/outs")
if(!dir.exists(figurePath)){
    dir.create(figurePath, recursive = T)
}
setwd(figurePath)
saveRDS(seuratObj, file.path(paste0('harmony_input.rds')))

