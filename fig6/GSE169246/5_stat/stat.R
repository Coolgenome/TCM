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

CD8_meta <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/MappingResult/CD8/meta_2022-03-04.tsv")
CD8_meta$Sample <- stringr::str_extract(CD8_meta$ID, "(?<=^.{10,20}\\.).+")
CD8_meta$Patient <- stringr::str_extract(CD8_meta$Sample, "P\\d+")
CD8_meta$Tissue <- stringr::str_extract(CD8_meta$Sample, "\\w$")
CD8_meta$TumorTreatment <- stringr::str_extract(CD8_meta$Sample, "^[a-zA-Z]+")
CD8_meta$isResponse <- "No"
CD8_meta$isResponse[CD8_meta$Patient %in% AllResponsePatients] <- "Yes"
CD8_meta$TreatmentType <- "PDL1+Chemo"
CD8_meta$TreatmentType[CD8_meta$Patient %in% TotalChemoPatients] <- "Chemo"

CD4_meta <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/MappingResult/CD4/meta_2022-03-04.tsv")
CD4_meta$Sample <- stringr::str_extract(CD4_meta$ID, "(?<=^.{10,20}\\.).+")
CD4_meta$Patient <- stringr::str_extract(CD4_meta$Sample, "P\\d+")
CD4_meta$Tissue <- stringr::str_extract(CD4_meta$Sample, "\\w$")
CD4_meta$TumorTreatment <- stringr::str_extract(CD4_meta$Sample, "^[a-zA-Z]+")
CD4_meta$isResponse <- "No"
CD4_meta$isResponse[CD4_meta$Patient %in% AllResponsePatients] <- "Yes"
CD4_meta$TreatmentType <- "PDL1+Chemo"
CD4_meta$TreatmentType[CD4_meta$Patient %in% TotalChemoPatients] <- "Chemo"

allMD <- list(CD8 = CD8_meta, CD4 = CD4_meta)

allT <- bind_rows(CD4_meta, CD8_meta)
allTSampleCN <- allT %>%
    group_by(Patient, Sample) %>%
    count

targetTreatments = c("Pre", "Post")
targetCellSet = c("AllT", "Split")
targetTreatmentType = c("PDL1+Chemo", "Chemo")

for(tTTT in targetTreatmentType){
    for(tCS in targetCellSet){
        for(tT in targetTreatments){
            for(ami in 1:length(allMD)){
                tempName <- names(allMD[ami])
                tempMD <- allMD[[ami]] %>%
                    filter(Tissue == "t") %>%
                    filter(TumorTreatment == tT) %>%
                    filter(TreatmentType == tTTT)
                if(tCS == "Split"){
                    TotalSampleCellNum <- tempMD %>%
                        group_by(Patient, Sample) %>%
                        count()
                }else{
                    TotalSampleCellNum <- allTSampleCN %>%
                        filter(Sample %in% tempMD$Sample)
                }
                totalT <- c()
                for(tempCluster in unique(tempMD$predicted.celltype)){
                    TNR <- tempMD %>%
                        filter(predicted.celltype == tempCluster) %>%
                        group_by(Patient, Sample) %>%
                        count()
                    TNR$Frac <- 0.0
                    TNR$Frac <- TNR$n / TotalSampleCellNum$n[match(TNR$Sample, TotalSampleCellNum$Sample)]
                    TNR$isResponse <- "No"
                    TNR$isResponse[TNR$Patient %in% AllResponsePatients] <- "Yes"
                    TNR$cluster <- tempCluster
                    totalT <- bind_rows(totalT, TNR)
                    missedSamples <- setdiff(
                        unique(TotalSampleCellNum$Sample),
                        unique(TNR$Sample))
                    if(length(missedSamples) >0){
                        Patient = TotalSampleCellNum$Patient[match(missedSamples, TotalSampleCellNum$Sample)]
                        isResponse = rep("No", length(missedSamples))
                        isResponse[Patient %in% AllResponsePatients] = "Yes"
                        tempTNR <- tibble(
                            Patient = Patient,
                            Sample = missedSamples,
                            Frac = 0,
                            isResponse = isResponse,
                            cluster = tempCluster
                        )
                        totalT <- bind_rows(totalT, tempTNR)
                    }
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
                figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/MappingResult",
                                        tempName,
                                        tTTT)
                if(!dir.exists(figurePath)){
                    dir.create(figurePath, recursive = T)
                }
                ggsave(file.path(figurePath, paste0(tempName, "_CellSet-", tTTT, "-", tCS, "_Treatment-", tT,"_response_bar.pdf")), g, width = 200, height = 1200, units = "mm")
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
                pdf(file.path(figurePath, paste0(tempName, "_CellSet-", tTTT, "-", tCS, "_Treatment-", tT, "_response_roe_heatmap.pdf")), width = widthi, height = heighti)
                print(pheatmap(matRoe,
                               color = my.colors,
                               breaks = my.breaks,
                               show_colnames = T,
                               show_rownames = T,
                               cluster_cols = F,
                               cluster_rows = F,
                               border_color = F))
                dev.off()
                pdf(file.path(figurePath, paste0(tempName, "_CellSet-", tTTT, "-", tCS, "_Treatment-", tT, "_response_cell_heatmap.pdf")), width = widthi * 1.5, height = heighti)
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
        }
    }
}
