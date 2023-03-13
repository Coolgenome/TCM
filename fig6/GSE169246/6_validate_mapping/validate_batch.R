#--------------------------------------------------------------
# filename : validate.R
# Date : 2022-08-31
# contributor : Yanshuo Chu
# function: validate
#--------------------------------------------------------------

print('<==== validate.R ====>')

rm(list=ls())

library(Seurat)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggpubr)

figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/6_validate_mapping/outs")
if(!dir.exists(figurePath)){
    dir.create(figurePath, recursive = T)
}
setwd(figurePath)


seuratObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/merged/nPC_30/UMAP_dist_0.1_nneighbor_50/GSE169246_UMAP_dist_0.1_nneighbor_50_CLUSTER_res_0.3/cluster.rds")

## Idents(seuratObj) <- seuratObj$orig.ident
## pdf(file.path(getwd(), "orig_ident.pdf"))
## DimPlot(seuratObj, label =T)
## dev.off()

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

Idents(seuratObj) <- seuratObj$Sample
pdf(file.path(getwd(), "sample.pdf"), width = 15)
DimPlot(seuratObj, label =T)
dev.off()

Idents(seuratObj) <- seuratObj$Patient
pdf(file.path(getwd(), "patient.pdf"))
DimPlot(seuratObj, label =T)
dev.off()

Idents(seuratObj) <- seuratObj$Tissue
pdf(file.path(getwd(), "tissue.pdf"))
DimPlot(seuratObj, label =T)
dev.off()

Idents(seuratObj) <- seuratObj$TreatmentType
pdf(file.path(getwd(), "treatment_type.pdf"))
DimPlot(seuratObj, label =T)
dev.off()
