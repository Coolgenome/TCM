#--------------------------------------------------------------
# filename : stat.R
# Date : 2022-03-04
# contributor : Yanshuo Chu
# function: stat
#--------------------------------------------------------------

print('<==== stat.R ====>')
rm(list=ls())

suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(ggplot2)
    library(stringr)
    library(pheatmap)
})

CD4Obj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/TCells/forMapping/CD4.rds")
CD8Obj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/TCells/forMapping/CD8.rds")

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

CD8_meta <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/MappingResult_filter/CD8/meta_2022-05-25.tsv")
CD8_meta$Sample <- stringr::str_extract(CD8_meta$ID, "(?<=^.{10,20}\\.).+")
CD8_meta$Patient <- stringr::str_extract(CD8_meta$Sample, "P\\d+")
CD8_meta$Tissue <- stringr::str_extract(CD8_meta$Sample, "\\w$")
CD8_meta$TumorTreatment <- stringr::str_extract(CD8_meta$Sample, "^[a-zA-Z]+")
CD8_meta$isResponse <- "NR"
CD8_meta$isResponse[CD8_meta$Patient %in% AllResponsePatients] <- "R"
CD8_meta$isResponse[CD8_meta$Patient == "P028"] <- "-"
CD8_meta$TreatmentType <- "PDL1+Chemo"
CD8_meta$TreatmentType[CD8_meta$Patient %in% TotalChemoPatients] <- "Chemo"
CD8_meta$group <- paste0(CD8_meta$TreatmentType, "-", CD8_meta$TumorTreatment, "-", CD8_meta$isResponse)

CD8Obj@meta.data$Sample <- stringr::str_extract(Cells(CD8Obj), "(?<=^.{10,20}\\.).+")
CD8Obj@meta.data$Patient <- stringr::str_extract(CD8Obj@meta.data$Sample, "P\\d+")
CD8Obj@meta.data$Tissue <- stringr::str_extract(CD8Obj@meta.data$Sample, "\\w$")
CD8Obj@meta.data$TumorTreatment <- stringr::str_extract(CD8Obj@meta.data$Sample, "^[a-zA-Z]+")
CD8Obj@meta.data$isResponse <- "NR"
CD8Obj@meta.data$isResponse[CD8Obj@meta.data$Patient %in% AllResponsePatients] <- "R"
CD8Obj@meta.data$isResponse[CD8Obj@meta.data$Patient == "P028"] <- "-"
CD8Obj@meta.data$TreatmentType <- "PDL1+Chemo"
CD8Obj@meta.data$TreatmentType[CD8Obj@meta.data$Patient %in% TotalChemoPatients] <- "Chemo"
CD8Obj@meta.data$group <- paste0(CD8Obj@meta.data$TreatmentType, "-", CD8Obj@meta.data$TumorTreatment, "-", CD8Obj@meta.data$isResponse)


CD4_meta <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/MappingResult_filter/CD4/meta_2022-05-25.tsv")
CD4_meta$Sample <- stringr::str_extract(CD4_meta$ID, "(?<=^.{10,20}\\.).+")
CD4_meta$Patient <- stringr::str_extract(CD4_meta$Sample, "P\\d+")
CD4_meta$Tissue <- stringr::str_extract(CD4_meta$Sample, "\\w$")
CD4_meta$TumorTreatment <- stringr::str_extract(CD4_meta$Sample, "^[a-zA-Z]+")
CD4_meta$isResponse <- "NR"
CD4_meta$isResponse[CD4_meta$Patient %in% AllResponsePatients] <- "R"
CD4_meta$isResponse[CD4_meta$Patient == "P028"] <- "-"
CD4_meta$TreatmentType <- "PDL1+Chemo"
CD4_meta$TreatmentType[CD4_meta$Patient %in% TotalChemoPatients] <- "Chemo"
CD4_meta$group <- paste0(CD4_meta$TreatmentType, "-", CD4_meta$TumorTreatment, "-", CD4_meta$isResponse)

CD4Obj@meta.data$Sample <- stringr::str_extract(Cells(CD4Obj), "(?<=^.{10,20}\\.).+")
CD4Obj@meta.data$Patient <- stringr::str_extract(CD4Obj@meta.data$Sample, "P\\d+")
CD4Obj@meta.data$Tissue <- stringr::str_extract(CD4Obj@meta.data$Sample, "\\w$")
CD4Obj@meta.data$TumorTreatment <- stringr::str_extract(CD4Obj@meta.data$Sample, "^[a-zA-Z]+")
CD4Obj@meta.data$isResponse <- "NR"
CD4Obj@meta.data$isResponse[CD4Obj@meta.data$Patient %in% AllResponsePatients] <- "R"
CD4Obj@meta.data$isResponse[CD4Obj@meta.data$Patient == "P028"] <- "-"
CD4Obj@meta.data$TreatmentType <- "PDL1+Chemo"
CD4Obj@meta.data$TreatmentType[CD4Obj@meta.data$Patient %in% TotalChemoPatients] <- "Chemo"
CD4Obj@meta.data$group <- paste0(CD4Obj@meta.data$TreatmentType, "-", CD4Obj@meta.data$TumorTreatment, "-", CD4Obj@meta.data$isResponse)

allMD <- list(CD8 = CD8_meta, CD4 = CD4_meta)
allObj <- list(CD8 = CD8Obj, CD4 = CD4Obj)

allT <- bind_rows(CD4_meta, CD8_meta)
allTSampleCN <- allT %>%
    group_by(Patient, Sample) %>%
    count

targetTreatments = c("Pre", "Post")
targetCellSet = c("AllT", "Split")
targetTreatmentType = c("PDL1+Chemo", "Chemo")


figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/5_stat_filter/outs4groups")
if(!dir.exists(figurePath)){
    dir.create(figurePath, recursive = T)
}
setwd(figurePath)

# CD8 fraction ################################################################

hs_genes <- c("HSPA1A", "HSPA1B")

for(oi in 1:length(allObj)){
    tempName <- names(allObj[oi])
    tempObj <- allObj[[oi]]
    Idents(tempObj) <- tempObj$group
    for(hg in hs_genes){
        pdf(file.path(getwd(), paste0(tempName, "_", hg, "_violinplot.pdf")), width = 10, height = 4)
        print(VlnPlot(tempObj, features = hg, pt.size = 0) + geom_boxplot(color = "black", fill = NA, outlier.shape = NA, width = 0.2))
        dev.off()
    }
}

