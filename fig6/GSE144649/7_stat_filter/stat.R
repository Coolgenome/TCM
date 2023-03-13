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

clinic_info <- tibble(
    patient = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "NC1", "NC2",
                "NC3", "NC4", "NC5", "NC6"),
    treatment = c("lpilimumab", "lpilimumab", "lpilimumab+Nivolumab",
                  "lpilimumab+Nivolumab", "lpilimumab+Nivolumab",
                  "lpilimumab+Nivolumab", "lpilimumab+Nivolumab",
                  "lpilimumab+Nivolumab", "lpilimumab+Nivolumab", "none",
                  "lpilimumab", "lpilimumab+Nivolumab", "lpilimumab+Nivolumab",
                  "lpilimumab+Nivolumab"),
    diagnosis = c("colitis", "colitis", "colitis", "entero-colitis",
                  "entero-colitis", "entero-colitis", "colitis", "colitis",
                  "enteritis", "normal", "enteritis", "normal", "enteritis",
                  "normal" ),
    response = c("PD", "PR", "PR", "PD", "SD", "PR", "PD", "PD", "PR", "*", "PD",
                 "PR", "PD", "PD"))


CD4PredictedObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE144649/6_mapping_filter/CD4/querySeuratObj_2022-05-25.rds")
CD4PredictedObj$patient <- stringr::str_extract(CD4PredictedObj$orig.ident, "^.*(?=-.+$)")
CD4_meta <- as_tibble(CD4PredictedObj@meta.data, rownames = NA) %>%
    filter(patient %in% clinic_info$patient)
CD4_meta$treatment <- clinic_info$treatment[match(CD4_meta$patient, clinic_info$patient)]
CD4_meta$diagnosis <- clinic_info$diagnosis[match(CD4_meta$patient, clinic_info$patient)]
CD4_meta$response <- clinic_info$response[match(CD4_meta$patient, clinic_info$patient)]
CD4_meta$dataset <- "CD4"
CD4_meta$sample <- CD4_meta$orig.ident
CD4_meta <- CD4_meta %>% filter(response != "*")
CD4_meta$isResponse <- "Yes"
CD4_meta$isResponse[CD4_meta$response %in% c("PD", "SD")] <- "No"
    
CD8PredictedObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE144649/6_mapping_filter/CD8/querySeuratObj_2022-05-25.rds")
CD8PredictedObj$patient <- stringr::str_extract(CD8PredictedObj$orig.ident, "^.*(?=-.+$)")
CD8_meta <- as_tibble(CD8PredictedObj@meta.data, rownames = NA) %>%
    filter(patient %in% clinic_info$patient)
CD8_meta$treatment <- clinic_info$treatment[match(CD8_meta$patient, clinic_info$patient)]
CD8_meta$diagnosis <- clinic_info$diagnosis[match(CD8_meta$patient, clinic_info$patient)]
CD8_meta$response <- clinic_info$response[match(CD8_meta$patient, clinic_info$patient)]
CD8_meta$dataset <- "CD8"
CD8_meta$sample <- CD8_meta$orig.ident
CD8_meta <- CD8_meta %>% filter(response != "*")
CD8_meta$isResponse <- "Yes"
CD8_meta$isResponse[CD8_meta$response %in% c("PD", "SD")] <- "No"

allMD <- list(CD8 = CD8_meta, CD4 = CD4_meta)
allT <- bind_rows(CD4_meta, CD8_meta)

allTSampleCN <- allT %>%
    group_by(patient, sample) %>%
    count

targetTreatments = c("lpilimumab+Nivolumab", "lpilimumab")
targetCellSet = c("AllT", "Split")
targetTreatmentType = c("lpilimumab+Nivolumab", "lpilimumab")

for(tCS in targetCellSet){
    for(tT in targetTreatments){
        for(ami in 1:length(allMD)){
            tempName <- names(allMD[ami])
            tempMD <- allMD[[ami]] %>% filter(treatment == tT)
            if(tCS == "Split"){
                TotalSampleCellNum <- tempMD %>%
                    group_by(patient, sample) %>%
                    count()
            }else{
                TotalSampleCellNum <- allTSampleCN %>%
                    filter(sample %in% tempMD$sample)
            }
            totalT <- c()
            for(tempCluster in unique(tempMD$predicted.id)){
                TNR <- tempMD %>%
                    filter(predicted.id == tempCluster) %>%
                    group_by(patient, sample, isResponse) %>%
                    count()
                TNR$Frac <- 0.0
                TNR$Frac <- TNR$n / TotalSampleCellNum$n[match(TNR$sample, TotalSampleCellNum$sample)]

                TNR$cluster <- tempCluster
                totalT <- bind_rows(totalT, TNR)
                missedSamples <- setdiff(
                    unique(TotalSampleCellNum$sample),
                    unique(TNR$sample))

                if(length(missedSamples) >0){
                    Patient = TotalSampleCellNum$patient[match(missedSamples, TotalSampleCellNum$sample)]
                    isResponse = rep("No", length(missedSamples))
                    isResponse[Patient %in% clinic_info$patient[clinic_info$response %in% c("CR", "PR")]] = "Yes"
                    tempTNR <- tibble(
                        patient = Patient,
                        sample = missedSamples,
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
            figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE144649/7_MappingResult_filter",
                                    tempName)
            if(!dir.exists(figurePath)){
                dir.create(figurePath, recursive = T)
            }
            ggsave(file.path(figurePath, paste0(tempName, "-", tCS, "_Treatment-", tT,"_response_bar.pdf")), g, width = 200, height = 1200, units = "mm")
            ClusterT <- tempMD %>%
                group_by(predicted.id) %>%
                dplyr::count()
            N <- dim(tempMD)[1]
            YesMD <- tempMD %>%
                filter(isResponse == "Yes") %>%
                rename(Group = isResponse)
            NoMD <- tempMD %>%
                filter(isResponse == "No") %>%
                rename(Group = isResponse)
            totalTargetTcellMD <- bind_rows(YesMD, NoMD)
            matRoe <- matrix(rep(0, length(unique(tempMD$predicted.id)) * 2),
                             nrow = length(unique(tempMD$predicted.id)),
                             ncol = 2)
            matCell <- matrix(rep(0, length(unique(tempMD$predicted.id)) * 2),
                              nrow = length(unique(tempMD$predicted.id)),
                              ncol = 2)
            rownames(matRoe) <- sort(unique(tempMD$predicted.id))
            colnames(matRoe) <- c("Yes", "No")
            rownames(matCell) <- sort(unique(tempMD$predicted.id))
            colnames(matCell) <- c("Yes", "No")
            ggplotdata <- c()
            for(tempGroup in c("Yes", "No")){
                targetTcellMD <- totalTargetTcellMD %>% filter(Group == tempGroup)
                M <- dim(targetTcellMD)[1]
                for(gsci in unique(tempMD$predicted.id)){
                    k <- ClusterT$n[ClusterT$predicted.id == gsci]
                    if(gsci %in% targetTcellMD$predicted.id){
                        n <- dim(targetTcellMD %>% filter(predicted.id == gsci))[1]
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
            pdf(file.path(figurePath, paste0(tempName, "-", tCS, "_Treatment-", tT, "_response_roe_heatmap.pdf")), width = widthi, height = heighti)
            print(pheatmap(matRoe,
                           color = my.colors,
                           breaks = my.breaks,
                           show_colnames = T,
                           show_rownames = T,
                           cluster_cols = F,
                           cluster_rows = F,
                           border_color = F))
            dev.off()
            pdf(file.path(figurePath, paste0(tempName, "-", tCS, "_Treatment-", tT, "_response_cell_heatmap.pdf")), width = widthi * 1.5, height = heighti)
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
