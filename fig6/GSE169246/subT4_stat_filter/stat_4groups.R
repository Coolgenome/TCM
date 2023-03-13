#--------------------------------------------------------------
# filename : stat_4groups.R
# Date : 2022-09-01
# contributor : Yanshuo Chu
# function: stat_4groups
#--------------------------------------------------------------

print('<==== stat_4groups.R ====>')

rm(list=ls())

suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(ggplot2)
    library(stringr)
    library(pheatmap)
})

CD4Obj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/subT2_split_by_marker/outs/CD4.rds")
CD8Obj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/subT2_split_by_marker/outs/CD8.rds")

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

CD8_meta <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/subT3_mapping_filter_split_by_marker/CD8/meta_2022-09-01.tsv")
CD8_meta$Sample <- stringr::str_extract(CD8_meta$ID, "(?<=^.{10,20}\\.).+")
CD8_meta$Patient <- stringr::str_extract(CD8_meta$Sample, "P\\d+")
CD8_meta$Tissue <- stringr::str_extract(CD8_meta$Sample, "\\w$")
CD8_meta$TumorTreatment <- stringr::str_extract(CD8_meta$Sample, "^[a-zA-Z]+")
CD8_meta$isResponse <- "NR"
CD8_meta$isResponse[CD8_meta$Patient %in% AllResponsePatients] <- "R"
CD8_meta$isResponse[CD8_meta$Patient %in% "P028"] <- "NA"
CD8_meta$TreatmentType <- "PDL1+Chemo"
CD8_meta$TreatmentType[CD8_meta$Patient %in% TotalChemoPatients] <- "Chemo"

CD4_meta <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/subT3_mapping_filter_split_by_marker/CD4/meta_2022-09-01.tsv")
CD4_meta$Sample <- stringr::str_extract(CD4_meta$ID, "(?<=^.{10,20}\\.).+")
CD4_meta$Patient <- stringr::str_extract(CD4_meta$Sample, "P\\d+")
CD4_meta$Tissue <- stringr::str_extract(CD4_meta$Sample, "\\w$")
CD4_meta$TumorTreatment <- stringr::str_extract(CD4_meta$Sample, "^[a-zA-Z]+")
CD4_meta$isResponse <- "NR"
CD4_meta$isResponse[CD4_meta$Patient %in% AllResponsePatients] <- "R"
CD4_meta$isResponse[CD4_meta$Patient %in% "P028"] <- "NA"
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

figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/subT4_stat_filter/outs")
if(!dir.exists(figurePath)){
    dir.create(figurePath, recursive = T)
}
setwd(figurePath)

for(tt in targetTreatmentType){
    for(ami in 1:length(allMD)){
        tempName <- names(allMD[ami])
        tempMD <- allMD[[ami]] %>%
            filter(Tissue == "t") %>%
            filter(TumorTreatment != "Prog") %>% 
            filter(isResponse != "NA") %>% 
            mutate(group = paste0(TumorTreatment, "_", isResponse)) %>%
            filter(TreatmentType == tt)

        ClusterT <- tempMD %>%
            group_by(predicted.celltype) %>%
            dplyr::count()
        N <- dim(tempMD)[1]

        matRoe <- matrix(rep(0, length(unique(tempMD$predicted.celltype)) * 2),
                         nrow = length(unique(tempMD$predicted.celltype)),
                         ncol = 4)
        matCell <- matrix(rep(0, length(unique(tempMD$predicted.celltype)) * 2),
                          nrow = length(unique(tempMD$predicted.celltype)),
                          ncol = 4)
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
        matRoe <- matRoe[, c("Pre_R", "Pre_NR", "Post_R", "Post_NR")]

        matRoe = matRoe[rowSums(matCell)>50,]
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
        pdf(file.path(figurePath, paste0("t_TreatmentType_", tt, "_", tempName, "_4group_response_roe_heatmap.pdf")), width = widthi, height = heighti)
        print(pheatmap(matRoe,
                       color = my.colors,
                       breaks = my.breaks,
                       show_colnames = T,
                       show_rownames = T,
                       cluster_cols = F,
                       cluster_rows = F,
                       border_color = F))
        dev.off()
        pdf(file.path(figurePath, paste0("t_TreatmentType_", tt, "_", tempName, "_4group_response_cell_heatmap.pdf")), width = widthi * 1.5, height = heighti)
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

# blood #######################################################################
for(tt in targetTreatmentType){
    for(ami in 1:length(allMD)){
        tempName <- names(allMD[ami])
        tempMD <- allMD[[ami]] %>%
            filter(Tissue == "b") %>%
            filter(TumorTreatment != "Prog") %>% 
            filter(isResponse != "NA") %>%
            mutate(group = paste0(TumorTreatment, "_", isResponse)) %>%
            filter(TreatmentType == tt)

        ClusterT <- tempMD %>%
            group_by(predicted.celltype) %>%
            dplyr::count()
        N <- dim(tempMD)[1]

        matRoe <- matrix(rep(0, length(unique(tempMD$predicted.celltype)) * 2),
                         nrow = length(unique(tempMD$predicted.celltype)),
                         ncol = 4)
        matCell <- matrix(rep(0, length(unique(tempMD$predicted.celltype)) * 2),
                          nrow = length(unique(tempMD$predicted.celltype)),
                          ncol = 4)
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
        matRoe <- matRoe[, c("Pre_R", "Pre_NR", "Post_R", "Post_NR")]

        matRoe = matRoe[rowSums(matCell)>50,]
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
        pdf(file.path(figurePath, paste0("b_TreatmentType_", tt, "_", tempName, "_4group_response_roe_heatmap.pdf")), width = widthi, height = heighti)
        print(pheatmap(matRoe,
                       color = my.colors,
                       breaks = my.breaks,
                       show_colnames = T,
                       show_rownames = T,
                       cluster_cols = F,
                       cluster_rows = F,
                       border_color = F))
        dev.off()
        pdf(file.path(figurePath, paste0("b_TreatmentType_", tt, "_", tempName, "_4group_response_cell_heatmap.pdf")), width = widthi * 1.5, height = heighti)
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





# CD8 fraction ################################################################

CD8C4Frac = CD8_meta %>%
    filter(Tissue == "t") %>% 
    group_by(TreatmentType, TumorTreatment, isResponse, predicted.celltype) %>%
    count %>%
    mutate(frac = n / dim(CD8_meta %>% filter(Tissue == "t"))[1]) %>%
    filter(predicted.celltype == 4)

g <- ggplot(CD8C4Frac %>% filter(TumorTreatment != "Prog")) +
    geom_bar(aes(x = TreatmentType, y = frac, fill = TreatmentType), stat = "identity") +
    facet_grid(TumorTreatment ~ isResponse) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("CD8C4Frac.pdf", g)
 
