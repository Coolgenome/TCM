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


CD8QueryPath <- "/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE179994/TCells/ForMapping/CD8.rds"
CD4QueryPath <- "/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE179994/TCells/ForMapping/CD4.rds"

CD8PredictedTPath <- "/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE179994/MappingResult_MultiMap/CD8/predictedT_2022-04-29.rds"
CD4PredictedTPath <- "/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE179994/MappingResult_MultiMap/CD4/predictedT_2022-04-29.rds"

CD8QueryObj <- readRDS(CD8QueryPath)
CD4QueryObj <- readRDS(CD4QueryPath)

CD8PredictedT <- readRDS(CD8PredictedTPath)
CD4PredictedT <- readRDS(CD4PredictedTPath)

clinicT <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/data/GSE179994/ClinicData.txt") %>%
    filter(`Treatment Hx` %in% c("On treatment"))
noResponseSamples <- clinicT %>%
    filter(Response == "No") %>%
    pull(`Sample Name`)
ResponseSamples <- clinicT %>%
    filter(Response == "Yes") %>%
    pull(`Sample Name`)

CD8_meta <- as_tibble(CD8QueryObj@meta.data, rownames = NA) %>%
    mutate(predicted.celltype = CD8PredictedT$MostFreq,
           barcode = Cells(CD8QueryObj))
CD8_meta$isResponse <- "No"
CD8_meta$isResponse[CD8_meta$sample %in% ResponseSamples] <- "Yes"

CD4_meta <- as_tibble(CD4QueryObj@meta.data, rownames = NA) %>%
    mutate(predicted.celltype = CD4PredictedT$MostFreq,
           barcode = Cells(CD4QueryObj))
CD4_meta$isResponse <- "No"
CD4_meta$isResponse[CD4_meta$sample %in% ResponseSamples] <- "Yes"

allMD <- list(CD8 = CD8_meta, CD4 = CD4_meta)

for(ami in 1:length(allMD)){
    tempName <- names(allMD[ami])
    tempMD <- allMD[[ami]] %>%
        rename(Sample = sample) %>%
        rename(Patient = patient)
    TotalSampleCellNum <- tempMD %>%
        group_by(Patient, Sample) %>%
        count()

    totalT <- c()
    for(tempCluster in unique(tempMD$predicted.celltype)){
        TNR <- tempMD %>%
            filter(predicted.celltype == tempCluster) %>%
            group_by(Patient, Sample) %>%
            count()
        TNR$Frac <- 0.0
        TNR$Frac <- TNR$n / TotalSampleCellNum$n[match(TNR$Sample, TotalSampleCellNum$Sample)]
        TNR$isResponse <- "No"
        TNR$isResponse[TNR$Sample %in% ResponseSamples] <- "Yes"
        TNR$cluster <- tempCluster
        totalT <- bind_rows(totalT, TNR)
    }
    g <- totalT %>%
        ggstatsplot::grouped_ggbetweenstats(
                         data = .,
                         x = isResponse,
                         y = Frac,
                         grouping.var = cluster,
                         xlab = "",
                         ylab = "Sample fraction",
                         ## pairwise.display = "aiwl", # display only significant pairwise comparisons
                         p.adjust.method = "fdr", # adjust p-values for multiple tests using this method
                         ggtheme = theme_classic(),
                         package = "ggsci",
                         palette = "default_jco",
                         plotgrid.args = list(ncol = 1))
    figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE179994/MappingResult_MultiMap", tempName)
    if(!dir.exists(figurePath)){
        dir.create(figurePath, recursive = T)
    }
    ggsave(file.path(figurePath, paste0(tempName, "response_bar.pdf")), g, width = 200, height = 1200, units = "mm")

    ClusterT <- tempMD %>%
        group_by(predicted.celltype) %>%
        dplyr::count()
    N <- dim(tempMD)[1]
    YesMD <- tempMD %>%
        filter(isResponse == "Yes") %>%
        rename(Group = isResponse)
    NoMD <- tempMD %>%
        filter(isResponse == "No") %>%
        rename(Group = isResponse)
    totalTargetTcellMD <- bind_rows(YesMD, NoMD)
    matRoe <- matrix(rep(0, length(unique(tempMD$predicted.celltype)) * 2),
                     nrow = length(unique(tempMD$predicted.celltype)),
                     ncol = 2)
    matCell <- matrix(rep(0, length(unique(tempMD$predicted.celltype)) * 2),
                     nrow = length(unique(tempMD$predicted.celltype)),
                     ncol = 2)
    rownames(matRoe) <- sort(unique(tempMD$predicted.celltype))
    colnames(matRoe) <- c("Yes", "No")
    rownames(matCell) <- sort(unique(tempMD$predicted.celltype))
    colnames(matCell) <- c("Yes", "No")
    ggplotdata <- c()
    for(tempGroup in c("Yes", "No")){
        targetTcellMD <- totalTargetTcellMD %>% filter(Group == tempGroup)
        M <- dim(targetTcellMD)[1]
        for(gsci in unique(tempMD$predicted.celltype)){
            k <- ClusterT$n[ClusterT$predicted.celltype == gsci]
            if(gsci %in% targetTcellMD$predicted.celltype){
                n <- dim(targetTcellMD %>% filter(predicted.celltype == gsci))[1]
            } else{
                n <- 0
            }
            tempData <- tibble(isResponse = tempGroup,
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
    ## colnames(matRoe) <- paste0(colnames(matRoe), " (", colSums(matCell)[colnames(matRoe)], ")")
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
    pdf(file.path(figurePath, paste0(tempName, "_response_roe_heatmap.pdf")), width = widthi, height = heighti)
    print(pheatmap(matRoe,
                   color = my.colors,
                   breaks = my.breaks,
                   show_colnames = T,
                   show_rownames = T,
                   cluster_cols = F,
                   cluster_rows = F,
                   border_color = F))
    dev.off()
    pdf(file.path(figurePath, paste0(tempName, "_response_cell_heatmap.pdf")), width = widthi * 1.5, height = heighti)
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
