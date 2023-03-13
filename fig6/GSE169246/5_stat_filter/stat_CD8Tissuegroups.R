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
CD8_meta$TreatmentType <- "PDL1+Chemo"
CD8_meta$TreatmentType[CD8_meta$Patient %in% TotalChemoPatients] <- "Chemo"

CD4_meta <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/MappingResult_filter/CD4/meta_2022-05-25.tsv")
CD4_meta$Sample <- stringr::str_extract(CD4_meta$ID, "(?<=^.{10,20}\\.).+")
CD4_meta$Patient <- stringr::str_extract(CD4_meta$Sample, "P\\d+")
CD4_meta$Tissue <- stringr::str_extract(CD4_meta$Sample, "\\w$")
CD4_meta$TumorTreatment <- stringr::str_extract(CD4_meta$Sample, "^[a-zA-Z]+")
CD4_meta$isResponse <- "NR"
CD4_meta$isResponse[CD4_meta$Patient %in% AllResponsePatients] <- "R"
CD4_meta$TreatmentType <- "PDL1+Chemo"
CD4_meta$TreatmentType[CD4_meta$Patient %in% TotalChemoPatients] <- "Chemo"

allMD <- list(CD8 = CD8_meta)

allT <- bind_rows(CD8_meta)
allTSampleCN <- allT %>%
    group_by(Patient, Sample) %>%
    count

targetTreatments = c("Pre", "Post")
targetCellSet = c("AllT", "Split")
targetTreatmentType = c("PDL1+Chemo", "Chemo")

figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/5_stat_filter/outs_tissuetype")
if(!dir.exists(figurePath)){
    dir.create(figurePath, recursive = T)
}
setwd(figurePath)


for(ami in 1:length(allMD)){
    tempName <- names(allMD[ami])
    tempMD <- allMD[[ami]] %>%
        mutate(group = Tissue)

    ClusterT <- tempMD %>%
        group_by(predicted.celltype) %>%
        dplyr::count()
    N <- dim(tempMD)[1]

    matRoe <- matrix(rep(0, length(unique(tempMD$predicted.celltype)) * 2),
                     nrow = length(unique(tempMD$predicted.celltype)),
                     ncol = 2)
    matCell <- matrix(rep(0, length(unique(tempMD$predicted.celltype)) * 2),
                      nrow = length(unique(tempMD$predicted.celltype)),
                      ncol = 2)
    rownames(matRoe) <- sort(unique(tempMD$predicted.celltype))
    colnames(matRoe) <- unique(tempMD$group)
    rownames(matCell) <- sort(unique(tempMD$predicted.celltype))
    colnames(matCell) <- unique(tempMD$group)

    ggplotdata <- c()
    for(tempGroup in unique(tempMD$group)){
        targetTcellMD <- tempMD %>% filter(group == tempGroup)
        M <- dim(targetTcellMD)[1]
        for(gsci in unique(tempMD$predicted.celltype)){
            k <- ClusterT$n[ClusterT$predicted.celltype == gsci]
            if(gsci %in% targetTcellMD$predicted.celltype){
                n <- dim(targetTcellMD %>% filter(predicted.celltype == gsci))[1]
            } else{
                n <- 0
            }
            tempData <- tibble(group = tempGroup,
                               Cluster = gsci,
                               Roe = (n/M) / (k/N))
            ggplotdata <- bind_rows(ggplotdata, tempData)
            tempRowName <- as.character(gsci)
            tempColName <- tempGroup
            matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
            matCell[tempRowName, tempColName] <- n
        }
    }
    matRoe[matRoe>10] = 10
    my.breaks <- seq(0, 1.99, by=0.01)
    my.colors <- c(
        colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)/2),
        colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks)/2))
    my.breaks.extra <- seq(2, 10, by = (10 - 2)/99)
    my.colors.extra <- colorRampPalette(colors = c("#FD9AA0", "#550000"))(length(my.breaks.extra))
    my.breaks <- c(my.breaks, my.breaks.extra)
    my.colors <- c(my.colors, my.colors.extra)
    widthi = ncol(matRoe)/12 + 1.8
    heighti = nrow(matRoe)/7 + 2
    print(rownames(matRoe))
    pdf(file.path(figurePath, paste0(tempName, "_2grouptissue_response_roe_heatmap.pdf")), width = widthi, height = heighti)
    print(pheatmap(matRoe,
                   color = my.colors,
                   breaks = my.breaks,
                   show_colnames = T,
                   show_rownames = T,
                   cluster_cols = F,
                   cluster_rows = F,
                   border_color = F))
    dev.off()
    pdf(file.path(figurePath, paste0(tempName, "_2grouptissue_response_cell_heatmap.pdf")), width = widthi * 1.5, height = heighti)
    print(pheatmap(matCell,
                   display_numbers = T,
                   show_colnames = T,
                   show_rownames = T,
                   cluster_cols = T,
                   cluster_rows = F,
                   number_format = "%.0f",
                   border_color = F))
    dev.off()
}
