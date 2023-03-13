
###############################################################################
#'                                  figure 6                                 '#
###############################################################################
rm(list=ls())
library(tidyverse)
library(Seurat)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(wesanderson)
library(viridis)
library(ggsci)
library(tidytext)
library(ggpubr)
library(cowplot)
library(facetscales)
library(latex2exp)
library(ggstatsplot)
set.seed(123)
CD8_path <- '/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/data/ClinicInfoBk/CD8_cluster_2021-07-22.rds'
CD4_path <- '/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/data/ClinicInfoBk/CD4_cluster_2021-07-22.rds'
NK_path <- '/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/data/ClinicInfoBk/NKT_cluster_2021-07-22.rds'
Proliferative_path <- '/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/data/ClinicInfoBk/Proliferative_cluster_2021-07-22.rds'
CD8_Obj <- readRDS(CD8_path)
CD4_Obj <- readRDS(CD4_path)
NK_Obj <- readRDS(NK_path)
Proliferative_Obj <- readRDS(Proliferative_path)
allPath <- list(CD8 = CD8_path, CD4 = CD4_path, NK = NK_path, Proliferative = Proliferative_path)
allObj <- list(CD8 = CD8_Obj, CD4 = CD4_Obj, NK = NK_Obj, Proliferative = Proliferative_Obj)
common_columns <- c()
for(oi in 1:length(allObj)){
    tempName <- names(allObj[oi])
    tempObj <- allObj[[oi]]
    tempObj$barcode <- Cells(tempObj)
    if(oi == 1){
        common_columns <- colnames(tempObj@meta.data)
    }else{
        common_columns <- intersect(common_columns, colnames(tempObj@meta.data))
    }
}
allMD <- c()
for(oi in 1:length(allObj)){
    tempName <- names(allObj[oi])
    tempObj <- allObj[[oi]]
    tempMD <- as_tibble(tempObj@meta.data, rownames = NA) %>%
        dplyr::select(all_of(common_columns)) %>% mutate(barcode = Cells(tempObj))
    allMD <- bind_rows(allMD, tempMD)
}
colorsForDataType <- c("#6DCCDD", "#EDCAE0", "#F494BE", "#F9B26C", "#A6ADCC", "#C4DA5D")
colorForClass1 <- c("#C4DA5D", "#6DCCDD", "#F494BE", "#EDCAE0")
dataTypeLevel <- c("CD4", "CD8", "NK", "NKT", "MAIT", "Proliferative")
figurePath <- "/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/figures4"
if(!dir.exists(figurePath)){
    dir.create(figurePath)
}

allMD$g_seurat_clusters[allMD$DataSet != "Proliferative"] <- paste0(allMD$DataSet[allMD$DataSet != "Proliferative"], "_c", allMD$seurat_clusters[allMD$DataSet != "Proliferative"])
allMD$g_seurat_clusters[allMD$DataSet == "Proliferative"] <- allMD$DataSet[allMD$DataSet == "Proliferative"]


ClusterT <- allMD %>% group_by(g_seurat_clusters) %>% dplyr::count()
N <- dim(allMD)[1]
RPs <- c("su001", "su002", "su003", "su004", "su009", "su012")
NRPs <- c("su005", "su006", "su007", "su008", "su010")
preRPs_pattern <- paste0(RPs, ".pre")
postRPs_pattern <- paste0(RPs, ".post")
preNRPs_pattern <- paste0(NRPs, ".pre")
postNRPs_pattern <- paste0(NRPs, ".post")
allPs <- unique(allMD %>% filter(batch == "BCC_Yost.rds") %>% pull(orig.ident))
preRPOIs <- grep(paste0(preRPs_pattern, collapse = "|"), allPs, value = T)
postRPOIs <- grep(paste0(postRPs_pattern, collapse = "|"), allPs, value = T)
preNRPOIs <- grep(paste0(preNRPs_pattern, collapse = "|"), allPs, value = T)
postNRPOIs <- grep(paste0(postNRPs_pattern, collapse = "|"), allPs, value = T)
preRTCellMD <- allMD %>% filter(orig.ident %in% preRPOIs) %>% mutate(tcellgroup = "R-pre")
postRTCellMD <- allMD %>% filter(orig.ident %in% postRPOIs) %>% mutate(tcellgroup = "R-post")
preNRTCellMD <- allMD %>% filter(orig.ident %in% preNRPOIs) %>% mutate(tcellgroup = "NR-pre")
postNRTCellMD <- allMD %>% filter(orig.ident %in% postNRPOIs) %>% mutate(tcellgroup = "NR-post")
totalTargetTcellMD <- bind_rows(preRTCellMD, postRTCellMD)
totalTargetTcellMD <- bind_rows(totalTargetTcellMD, preNRTCellMD)
totalTargetTcellMD <- bind_rows(totalTargetTcellMD, postNRTCellMD)
matRoe <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 4),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 4)
rownames(matRoe) <- sort(unique(allMD$g_seurat_clusters))
colnames(matRoe) <- c("R-pre", "R-post", "NR-pre", "NR-post")
matCell <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 4),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 4)
rownames(matCell) <- sort(unique(allMD$g_seurat_clusters))
colnames(matCell) <- c("R-pre", "R-post", "NR-pre", "NR-post")
ggplotdata <- c()
for(tempGroup in c("R-pre", "R-post", "NR-pre", "NR-post")){
    targetTcellMD <- totalTargetTcellMD %>% filter(tcellgroup == tempGroup)
    M <- dim(targetTcellMD)[1]
    for(gsci in unique(allMD$g_seurat_clusters)){
        k <- ClusterT$n[ClusterT$g_seurat_clusters == gsci]
        if(gsci %in% targetTcellMD$g_seurat_clusters){
            n <- dim(targetTcellMD %>% filter(g_seurat_clusters == gsci))[1]
        } else{
            n <- 0
        }
        tempData <- tibble(Group = tempGroup,
                           Cluster = gsci,
                           Roe = (n/M) / (k/N))
        ggplotdata <- bind_rows(ggplotdata, tempData)
        tempRowName <- gsci
        tempColName <- tempGroup
        matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
        matCell[tempRowName, tempColName] <- n
    }
}
matRoe[matRoe>10] = 10
colnames(matRoe) <- paste0(colnames(matRoe), " (", colSums(matCell)[colnames(matRoe)], ")")
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
pdf(file.path(figurePath, paste0("comparison1_roe_heatmap.pdf")), width = widthi, height = heighti)
print(pheatmap(matRoe,
               color = my.colors,
               breaks = my.breaks,
               ## annotation_row = cellType_col,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figurePath, paste0("comparison1_cell_heatmap.pdf")), width = widthi * 1.5, height = heighti)
print(pheatmap(matCell,
               display_numbers = T,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = T,
               cluster_rows = F,
               number_format = "%.0f",
               border_color = F))
dev.off()
colnames(matRoe) <- paste0("1_", colnames(matRoe))
saveRDS(matRoe, file.path(figurePath, paste0('matRoe1.rds')))


##################################################################################################################################################################
#' comparison 3: Metastatic tumor tissue         SKCM (all melanoma project) Compare by tissue source (organSite: Skin vs. lymph node vs. bowel vs. Adrenal gland) '#
##################################################################################################################################################################
library(stringr)
ClusterT <- allMD %>% group_by(g_seurat_clusters) %>% dplyr::count()
N <- dim(allMD)[1]
totalTargetTcellMD <- allMD %>% filter(CancerType == "SKCM")
matRoe <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) *
                        length(unique(totalTargetTcellMD$OrganSite))),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = length(unique(totalTargetTcellMD$OrganSite)))
matCell <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) *
                        length(unique(totalTargetTcellMD$OrganSite))),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = length(unique(totalTargetTcellMD$OrganSite)))
rownames(matRoe) <- sort(unique(allMD$g_seurat_clusters))
colnames(matRoe) <- unique(totalTargetTcellMD$OrganSite)
rownames(matCell) <- sort(unique(allMD$g_seurat_clusters))
colnames(matCell) <- unique(totalTargetTcellMD$OrganSite)
ggplotdata <- c()
for(tempOrganSite in unique(totalTargetTcellMD$OrganSite)){
    targetTcellMD <- totalTargetTcellMD %>% filter(OrganSite == tempOrganSite)
    M <- dim(targetTcellMD)[1]
    for(gsci in unique(allMD$g_seurat_clusters)){
        k <- ClusterT$n[ClusterT$g_seurat_clusters == gsci]
        if(gsci %in% targetTcellMD$g_seurat_clusters){
            n <- dim(targetTcellMD %>% filter(g_seurat_clusters == gsci))[1]
        } else{
            n <- 0
        }
        tempData <- tibble(OrganSite = tempOrganSite,
                           Cluster = gsci,
                           Roe = (n/M) / (k/N))
        ggplotdata <- bind_rows(ggplotdata, tempData)
        tempRowName <- gsci
        tempColName <- tempOrganSite
        matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
        matCell[tempRowName, tempColName] <- n
    }
}
matRoe[matRoe>10] = 10
colnames(matRoe) <- paste0(colnames(matRoe), " (", colSums(matCell)[colnames(matRoe)], ")")
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
pdf(file.path(figurePath, paste0("comparison3_roe_heatmap.pdf")), width = widthi, height = heighti)
print(pheatmap(matRoe,
               color = my.colors,
               breaks = my.breaks,
               ## annotation_col = cellType_col,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figurePath, paste0("comparison3_cell_heatmap.pdf")), width = widthi * 1.5, height = heighti)
print(pheatmap(matCell,
               display_numbers = T,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               number_format = "%.0f",
               border_color = F))
dev.off()
colnames(matRoe) <- paste0("12_", colnames(matRoe))
saveRDS(matRoe, file.path(figurePath, paste0('matRoe12.rds')))



library(stringr)
ClusterT <- allMD %>% group_by(g_seurat_clusters) %>% dplyr::count()
N <- dim(allMD)[1]
MCLBloodMD <- allMD %>% filter(OrganSite == "Blood") %>% filter(CancerType == "MCL") %>% mutate(Group = "MCL_Blood")
MCLBMMD <- allMD %>% filter(OrganSite == "BM") %>% filter(CancerType == "MCL") %>% mutate(Group = "MCL_BM")
HealthyBloodMD <- allMD %>% filter(OrganSite == "Blood") %>% filter(TissueType == "Healthy donor") %>% mutate(Group = "Healthy_Blood")
HealthyBMMD <- allMD %>% filter(OrganSite == "BM") %>% filter(TissueType == "Healthy donor") %>% mutate(Group = "Healthy_BM")
totalTargetTcellMD <- bind_rows(MCLBloodMD, MCLBMMD, HealthyBloodMD, HealthyBMMD)
matRoe <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 4),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 4)
matCell <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 4),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 4)
rownames(matRoe) <- sort(unique(allMD$g_seurat_clusters))
colnames(matRoe) <- c("MCL_Blood", "MCL_BM", "Healthy_Blood", "Healthy_BM")
rownames(matCell) <- sort(unique(allMD$g_seurat_clusters))
colnames(matCell) <- c("MCL_Blood", "MCL_BM", "Healthy_Blood", "Healthy_BM")
ggplotdata <- c()
for(tempGroup in c("MCL_Blood", "MCL_BM", "Healthy_Blood", "Healthy_BM")){
    targetTcellMD <- totalTargetTcellMD %>% filter(Group == tempGroup)
    M <- dim(targetTcellMD)[1]
    for(gsci in unique(allMD$g_seurat_clusters)){
        k <- ClusterT$n[ClusterT$g_seurat_clusters == gsci]
        if(gsci %in% targetTcellMD$g_seurat_clusters){
            n <- dim(targetTcellMD %>% filter(g_seurat_clusters == gsci))[1]
        } else{
            n <- 0
        }
        tempData <- tibble(Group = tempGroup,
                           Cluster = gsci,
                           Roe = (n/M) / (k/N))
        ggplotdata <- bind_rows(ggplotdata, tempData)
        tempRowName <- gsci
        tempColName <- tempGroup
        matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
        matCell[tempRowName, tempColName] <- n
    }
}
matRoe[matRoe>10] = 10
colnames(matRoe) <- paste0(colnames(matRoe), " (", colSums(matCell)[colnames(matRoe)], ")")
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
pdf(file.path(figurePath, paste0("comparison5_1_roe_heatmap.pdf")), width = widthi, height = heighti)
print(pheatmap(matRoe,
               color = my.colors,
               breaks = my.breaks,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figurePath, paste0("comparison5_1_cell_heatmap.pdf")), width = widthi * 1.5, height = heighti)
print(pheatmap(matCell,
               display_numbers = T,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = T,
               cluster_rows = F,
               number_format = "%.0f",
               border_color = F))
dev.off()
colnames(matRoe) <- paste0("16_", colnames(matRoe))
saveRDS(matRoe, file.path(figurePath, paste0('matRoe16.rds')))



library(stringr)
## allMD$g_seurat_clusters[allMD$DataSet != "Proliferative"] <- paste0(allMD$DataSet[allMD$DataSet != "Proliferative"], "_c", allMD$seurat_clusters[allMD$DataSet != "Proliferative"])
## allMD$g_seurat_clusters[allMD$DataSet == "Proliferative"] <- allMD$DataSet[allMD$DataSet == "Proliferative"]
ClusterT <- allMD %>% group_by(g_seurat_clusters) %>% dplyr::count()
N <- dim(allMD)[1]
LBCLMD <- allMD %>% filter(OrganSite == "LN") %>% filter(CancerType == "LBCL") %>% mutate(Group = "LBCL")
FLMD <- allMD %>% filter(OrganSite == "LN") %>% filter(CancerType == "FL") %>% mutate(Group = "FL")
RLNMD <- allMD %>% filter(OrganSite == "LN") %>% filter( TissueType %in% c("Healthy donor")) %>% mutate(Group = "Reactive_LN")
totalTargetTcellMD <- bind_rows(LBCLMD, FLMD, RLNMD)
matRoe <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 3),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 3)
matCell <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 3),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 3)
rownames(matRoe) <- sort(unique(allMD$g_seurat_clusters))
colnames(matRoe) <- c("LBCL", "FL", "Reactive_LN")
rownames(matCell) <- sort(unique(allMD$g_seurat_clusters))
colnames(matCell) <- c("LBCL", "FL", "Reactive_LN")
ggplotdata <- c()
for(tempGroup in c("LBCL", "FL", "Reactive_LN")){
    targetTcellMD <- totalTargetTcellMD %>% filter(Group == tempGroup)
    M <- dim(targetTcellMD)[1]
    for(gsci in unique(allMD$g_seurat_clusters)){
        k <- ClusterT$n[ClusterT$g_seurat_clusters == gsci]
        if(gsci %in% targetTcellMD$g_seurat_clusters){
            n <- dim(targetTcellMD %>% filter(g_seurat_clusters == gsci))[1]
        } else{
            n <- 0
        }
        tempData <- tibble(CancerType = tempGroup,
                           Cluster = gsci,
                           Roe = (n/M) / (k/N))
        ggplotdata <- bind_rows(ggplotdata, tempData)
        tempRowName <- gsci
        tempColName <- tempGroup
        matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
        matCell[tempRowName, tempColName] <- n
    }
}
matRoe[matRoe>10] = 10
colnames(matRoe) <- paste0(colnames(matRoe), " (", colSums(matCell)[colnames(matRoe)], ")")
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
pdf(file.path(figurePath, paste0("comparison5_2_roe_heatmap.pdf")), width = widthi, height = heighti)
print(pheatmap(matRoe,
               color = my.colors,
               breaks = my.breaks,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figurePath, paste0("comparison5_2_cell_heatmap.pdf")), width = widthi * 1.5, height = heighti)
print(pheatmap(matCell,
               display_numbers = T,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = T,
               cluster_rows = F,
               number_format = "%.0f",
               border_color = F))
dev.off()
colnames(matRoe) <- paste0("17_", colnames(matRoe))
saveRDS(matRoe, file.path(figurePath, paste0('matRoe17.rds')))


library(stringr)
ClusterT <- allMD %>% group_by(g_seurat_clusters) %>% dplyr::count()
N <- dim(allMD)[1]
AMLBMMD <- allMD %>% filter(OrganSite == "BM") %>% filter(CancerType == "AML") %>% mutate(Group = "AML_BM")
MCLBMMD <- allMD %>% filter(TissueType == "Healthy donor") %>% filter(OrganSite == "BM") %>% mutate(Group = "Healthy_BM")
totalTargetTcellMD <- bind_rows(AMLBMMD, MCLBMMD)
matRoe <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 2),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 2)
matCell <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 2),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 2)
rownames(matRoe) <- sort(unique(allMD$g_seurat_clusters))
colnames(matRoe) <- c("AML_BM", "Healthy_BM")
rownames(matCell) <- sort(unique(allMD$g_seurat_clusters))
colnames(matCell) <- c("AML_BM", "Healthy_BM")
ggplotdata <- c()
for(tempGroup in c("AML_BM", "Healthy_BM")){
    targetTcellMD <- totalTargetTcellMD %>% filter(Group == tempGroup)
    M <- dim(targetTcellMD)[1]
    for(gsci in unique(allMD$g_seurat_clusters)){
        k <- ClusterT$n[ClusterT$g_seurat_clusters == gsci]
        if(gsci %in% targetTcellMD$g_seurat_clusters){
            n <- dim(targetTcellMD %>% filter(g_seurat_clusters == gsci))[1]
        } else{
            n <- 0
        }
        tempData <- tibble(CancerType = tempGroup,
                           Cluster = gsci,
                           Roe = (n/M) / (k/N))
        ggplotdata <- bind_rows(ggplotdata, tempData)
        tempRowName <- gsci
        tempColName <- tempGroup
        matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
        matCell[tempRowName, tempColName] <- n
    }
}
matRoe[matRoe>10] = 10
colnames(matRoe) <- paste0(colnames(matRoe), " (", colSums(matCell)[colnames(matRoe)], ")")
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
pdf(file.path(figurePath, paste0("comparison5_3_roe_heatmap.pdf")), width = widthi, height = heighti)
print(pheatmap(matRoe,
               color = my.colors,
               breaks = my.breaks,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figurePath, paste0("comparison5_3_cell_heatmap.pdf")), width = widthi * 1.5, height = heighti)
print(pheatmap(matCell,
               display_numbers = T,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = T,
               cluster_rows = F,
               number_format = "%.0f",
               border_color = F))
dev.off()
colnames(matRoe) <- paste0("18_", colnames(matRoe))
saveRDS(matRoe, file.path(figurePath, paste0('matRoe18.rds')))

library(stringr)
ClusterT <- allMD %>% group_by(g_seurat_clusters) %>% dplyr::count()
N <- dim(allMD)[1]
totalTargetTcellMD <- allMD %>% filter(CancerType %in% c("Head and neck", "HNSC"))
matRoe <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) *
                        3),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 3)
rownames(matRoe) <- sort(unique(allMD$g_seurat_clusters))
colnames(matRoe) <- unique(totalTargetTcellMD$TissueType)
matCell <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) *
                        3),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 3)
rownames(matCell) <- sort(unique(allMD$g_seurat_clusters))
colnames(matCell) <- unique(totalTargetTcellMD$TissueType)
ggplotdata <- c()
for(tempCancerType in unique(totalTargetTcellMD$TissueType)){
    targetTcellMD <- totalTargetTcellMD %>% filter(TissueType == tempCancerType)
    M <- dim(targetTcellMD)[1]
    for(gsci in unique(allMD$g_seurat_clusters)){
        k <- ClusterT$n[ClusterT$g_seurat_clusters == gsci]
        if(gsci %in% targetTcellMD$g_seurat_clusters){
            n <- dim(targetTcellMD %>% filter(g_seurat_clusters == gsci))[1]
        } else{
            n <- 0
        }
        tempData <- tibble(CancerType = tempCancerType,
                           Cluster = gsci,
                           Roe = (n/M) / (k/N))
        ggplotdata <- bind_rows(ggplotdata, tempData)
        tempRowName <- gsci
        tempColName <- tempCancerType
        matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
        matCell[tempRowName, tempColName] <- n
    }
}
matRoe[matRoe>10] = 10
colnames(matRoe) <- paste0(colnames(matRoe), " (", colSums(matCell)[colnames(matRoe)], ")")
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
pdf(file.path(figurePath, paste0("comparison6_roe_heatmap.pdf")), width = widthi, height = heighti)
print(pheatmap(matRoe,
               color = my.colors,
               breaks = my.breaks,
               ## annotation_col = cellType_col,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figurePath, paste0("comparison6_cell_heatmap.pdf")), width = widthi * 1.5, height = heighti)
print(pheatmap(matCell,
               display_numbers = T,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = T,
               cluster_rows = F,
               number_format = "%.0f",
               border_color = F))
dev.off()
colnames(matRoe) <- paste0("19_", colnames(matRoe))
saveRDS(matRoe, file.path(figurePath, paste0('matRoe19.rds')))



library(stringr)
ClusterT <- allMD %>% group_by(g_seurat_clusters) %>% dplyr::count()
N <- dim(allMD)[1]
totalTargetTcellMD <- allMD %>% filter(CancerType %in% c("LUAD"))
matRoe <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) *
                        length(unique(totalTargetTcellMD$TissueType))),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = length(unique(totalTargetTcellMD$TissueType)))
rownames(matRoe) <- sort(unique(allMD$g_seurat_clusters))
colnames(matRoe) <- unique(totalTargetTcellMD$TissueType)
matCell <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) *
                        length(unique(totalTargetTcellMD$TissueType))),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = length(unique(totalTargetTcellMD$TissueType)))
rownames(matCell) <- sort(unique(allMD$g_seurat_clusters))
colnames(matCell) <- unique(totalTargetTcellMD$TissueType)
ggplotdata <- c()
for(tempCancerType in unique(totalTargetTcellMD$TissueType)){
    targetTcellMD <- totalTargetTcellMD %>% filter(TissueType == tempCancerType)
    if(tempCancerType == "Metastatic tumor tissue"){
        BrainMetTCellMD <- allMD %>% filter(Patient %in% c("Lung_GSE123902_LX255B", "Lung_GSE123902_LX681", "Lung_GSE123902_LX701")) %>% mutate(tcellgroup = "BrainMet")
        targetTcellMD <- targetTcellMD %>% filter(! barcode %in% BrainMetTCellMD$barcode)
    }
    M <- dim(targetTcellMD)[1]
    for(gsci in unique(allMD$g_seurat_clusters)){
        k <- ClusterT$n[ClusterT$g_seurat_clusters == gsci]
        if(gsci %in% targetTcellMD$g_seurat_clusters){
            n <- dim(targetTcellMD %>% filter(g_seurat_clusters == gsci))[1]
        } else{
            n <- 0
        }
        tempData <- tibble(CancerType = tempCancerType,
                           Cluster = gsci,
                           Roe = (n/M) / (k/N))
        ggplotdata <- bind_rows(ggplotdata, tempData)
        tempRowName <- gsci
        tempColName <- tempCancerType
        matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
        matCell[tempRowName, tempColName] <- n
    }
}
matRoe[matRoe>10] = 10
colnames(matRoe) <- paste0(colnames(matRoe), " (", colSums(matCell)[colnames(matRoe)], ")")
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
pdf(file.path(figurePath, paste0("comparison7_roe_heatmap.pdf")), width = widthi, height = heighti)
print(pheatmap(matRoe,
               color = my.colors,
               breaks = my.breaks,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figurePath, paste0("comparison7_cell_heatmap.pdf")), width = widthi * 1.5, height = heighti)
print(pheatmap(matCell,
               display_numbers = T,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = T,
               cluster_rows = F,
               number_format = "%.0f",
               border_color = F))
dev.off()
colnames(matRoe) <- paste0("7_", colnames(matRoe))
saveRDS(matRoe, file.path(figurePath, paste0('matRoe7.rds')))


test <- allMD %>% filter(batch == "Leukocyte10x.T.rds") %>% group_by(Sample, TissueType) %>% dplyr::count()

library(stringr)
ClusterT <- allMD %>% group_by(g_seurat_clusters) %>% dplyr::count()
N <- dim(allMD)[1]
totalTargetTcellMD <- allMD %>% filter(TissueType %in% c("Primary tumor tissue")) %>% filter(CancerType %in% c("LUAD", "LUSC"))
matRoe <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) *
                        length(unique(totalTargetTcellMD$CancerType))),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = length(unique(totalTargetTcellMD$CancerType)))
rownames(matRoe) <- sort(unique(allMD$g_seurat_clusters))
colnames(matRoe) <- unique(totalTargetTcellMD$CancerType)
matCell <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) *
                        length(unique(totalTargetTcellMD$CancerType))),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = length(unique(totalTargetTcellMD$CancerType)))
rownames(matCell) <- sort(unique(allMD$g_seurat_clusters))
colnames(matCell) <- unique(totalTargetTcellMD$CancerType)
ggplotdata <- c()
for(tempCancerType in unique(totalTargetTcellMD$CancerType)){
    targetTcellMD <- totalTargetTcellMD %>% filter(CancerType == tempCancerType)
    M <- dim(targetTcellMD)[1]
    for(gsci in unique(allMD$g_seurat_clusters)){
        k <- ClusterT$n[ClusterT$g_seurat_clusters == gsci]
        if(gsci %in% targetTcellMD$g_seurat_clusters){
            n <- dim(targetTcellMD %>% filter(g_seurat_clusters == gsci))[1]
        } else{
            n <- 0
        }
        tempData <- tibble(CancerType = tempCancerType,
                           Cluster = gsci,
                           Roe = (n/M) / (k/N))
        ggplotdata <- bind_rows(ggplotdata, tempData)
        tempRowName <- gsci
        tempColName <- tempCancerType
        matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
        matCell[tempRowName, tempColName] <- n
    }
}
matRoe[matRoe>10] = 10
colnames(matRoe) <- paste0(colnames(matRoe), " (", colSums(matCell)[colnames(matRoe)], ")")
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
pdf(file.path(figurePath, paste0("comparison8_roe_heatmap.pdf")), width = widthi, height = heighti)
print(pheatmap(matRoe,
               color = my.colors,
               breaks = my.breaks,
               ## annotation_col = cellType_col,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figurePath, paste0("comparison8_cell_heatmap.pdf")), width = widthi * 1.5, height = heighti)
print(pheatmap(matCell,
               display_numbers = T,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = T,
               cluster_rows = F,
               number_format = "%.0f",
               border_color = F))
dev.off()





allMD %>% filter(batch == "BCC_Yost.rds") %>% group_by(Patient) %>% dplyr::count()

#' mark each cluster name #####################################################
## allMD$g_seurat_clusters[allMD$DataSet != "Proliferative"] <- paste0(allMD$DataSet[allMD$DataSet != "Proliferative"], "_c", allMD$seurat_clusters[allMD$DataSet != "Proliferative"])
## allMD$g_seurat_clusters[allMD$DataSet == "Proliferative"] <- allMD$DataSet[allMD$DataSet == "Proliferative"]
ClusterT <- allMD %>% group_by(g_seurat_clusters) %>% dplyr::count()
N <- dim(allMD)[1]
ptTCellMD <- allMD %>% filter(Patient %in% c("BCC_Yost_bcc.su010")) %>% mutate(tcellgroup = "Pembrolizumab")
ctTCellMD <- allMD %>% filter(Patient %in% c("BCC_Yost_bcc.su004")) %>% mutate(tcellgroup = "Cemiplimab")
totalTargetTcellMD <- bind_rows(ptTCellMD, ctTCellMD)
matRoe <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 2),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 2)
rownames(matRoe) <- sort(unique(allMD$g_seurat_clusters))
colnames(matRoe) <- c("Pembrolizumab", "Cemiplimab")
matCell <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 2),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 2)
rownames(matCell) <- sort(unique(allMD$g_seurat_clusters))
colnames(matCell) <- c("Pembrolizumab", "Cemiplimab")
for(tempGroup in c("Pembrolizumab", "Cemiplimab")){
    targetTcellMD <- totalTargetTcellMD %>% filter(tcellgroup == tempGroup)
    M <- dim(targetTcellMD)[1]
    for(gsci in unique(allMD$g_seurat_clusters)){
        k <- ClusterT$n[ClusterT$g_seurat_clusters == gsci]
        if(gsci %in% targetTcellMD$g_seurat_clusters){
            n <- dim(targetTcellMD %>% filter(g_seurat_clusters == gsci))[1]
        } else{
            n <- 0
        }
        tempRowName <- gsci
        tempColName <- tempGroup
        matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
        matCell[tempRowName, tempColName] <- n
    }
}
matRoe[matRoe>10] <- 10
colnames(matRoe) <- paste0(colnames(matRoe), " (", colSums(matCell)[colnames(matRoe)], ")")
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
pdf(file.path(figurePath, paste0("BCC_Yost_Treatment_PembrolizumabVSCemiplimab_roe.pdf")), width = widthi, height = heighti)
print(pheatmap(matRoe,
               color = my.colors,
               breaks = my.breaks,
               ## annotation_row = cellType_col,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figurePath, paste0("BCC_Yost_Treatment_PembrolizumabVSCemiplimab_cell.pdf")), width = widthi, height = heighti)
print(pheatmap(matCell,
               color = my.colors,
               breaks = my.breaks,
               ## annotation_row = cellType_col,
               display_numbers = T,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               number_format = "%.0f",
               border_color = F))
dev.off()
colnames(matRoe) <- paste0("2_", colnames(matRoe))
saveRDS(matRoe, file.path(figurePath, paste0('matRoe2.rds')))


###############################################################################
#'                                Breast_GSE114725.rds                       '#
###############################################################################


allMD %>% filter(batch == "Breast_GSE114725.rds") %>% group_by(Patient) %>% dplyr::count()

#' mark each cluster name #####################################################
## allMD$g_seurat_clusters[allMD$DataSet != "Proliferative"] <- paste0(allMD$DataSet[allMD$DataSet != "Proliferative"], "_c", allMD$seurat_clusters[allMD$DataSet != "Proliferative"])
## allMD$g_seurat_clusters[allMD$DataSet == "Proliferative"] <- allMD$DataSet[allMD$DataSet == "Proliferative"]
ClusterT <- allMD %>% group_by(g_seurat_clusters) %>% dplyr::count()
N <- dim(allMD)[1]
Her2TCellMD <- allMD %>% filter(Patient %in% c("Breast_GSE114725_BC7")) %>% mutate(tcellgroup = "Her2")
TNBCTCellMD <- allMD %>% filter(Patient %in% c("Breast_GSE114725_BC3", "Breast_GSE114725_BC5")) %>% mutate(tcellgroup = "TNBC")
ERPRTCellMD <- allMD %>% filter(Patient %in% c("Breast_GSE114725_BC1", "Breast_GSE114725_BC4")) %>% mutate(tcellgroup = "ERPR")
ERTCellMD <- allMD %>% filter(Patient %in% c("Breast_GSE114725_BC2", "Breast_GSE114725_BC6")) %>% mutate(tcellgroup = "ER")
totalTargetTcellMD <- bind_rows(Her2TCellMD, TNBCTCellMD, ERPRTCellMD, ERTCellMD)
matRoe <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 4),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 4)
rownames(matRoe) <- sort(unique(allMD$g_seurat_clusters))
colnames(matRoe) <- c("Her2", "TNBC", "ERPR", "ER")
matCell <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 4),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 4)
rownames(matCell) <- sort(unique(allMD$g_seurat_clusters))
colnames(matCell) <- c("Her2", "TNBC", "ERPR", "ER")
for(tempGroup in c("Her2", "TNBC", "ERPR", "ER")){
    targetTcellMD <- totalTargetTcellMD %>% filter(tcellgroup == tempGroup)
    M <- dim(targetTcellMD)[1]
    for(gsci in unique(allMD$g_seurat_clusters)){
        k <- ClusterT$n[ClusterT$g_seurat_clusters == gsci]
        if(gsci %in% targetTcellMD$g_seurat_clusters){
            n <- dim(targetTcellMD %>% filter(g_seurat_clusters == gsci))[1]
        } else{
            n <- 0
        }
        tempRowName <- gsci
        tempColName <- tempGroup
        matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
        matCell[tempRowName, tempColName] <- n
    }
}
matRoe[matRoe>10] - 10
colnames(matRoe) <- paste0(colnames(matRoe), " (", colSums(matCell)[colnames(matRoe)], ")")
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
pdf(file.path(figurePath, paste0("Breast_GSE114725_Her2ERPR_roe.pdf")), width = widthi, height = heighti)
print(pheatmap(matRoe,
               color = my.colors,
               breaks = my.breaks,
               ## annotation_row = cellType_col,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figurePath, paste0("Breast_GSE114725_Her2ERPR_cell.pdf")), width = widthi, height = heighti)
print(pheatmap(matCell,
               color = my.colors,
               breaks = my.breaks,
               ## annotation_row = cellType_col,
               display_numbers = T,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               number_format = "%.0f",
               border_color = F))
dev.off()
colnames(matRoe) <- paste0("15_", colnames(matRoe))
saveRDS(matRoe, file.path(figurePath, paste0('matRoe15.rds')))

###############################################################################
#'                    HeadNeck_Puram.rds and HN_GSE139324.rds  #
###############################################################################

## both treatment naive

allMD %>% filter(batch == "HeadNeck_Puram.rds") %>% group_by(orig.ident, Patient, Sample, TissueType) %>% dplyr::count()

seuratObj <- readRDS("/rsrch3/home/genomic_med/ychu2/data/data_SJZ/T/All/HeadNeck_Puram.rds")


#' mark each cluster name #####################################################
## allMD$g_seurat_clusters[allMD$DataSet != "Proliferative"] <- paste0(allMD$DataSet[allMD$DataSet != "Proliferative"], "_c", allMD$seurat_clusters[allMD$DataSet != "Proliferative"])
## allMD$g_seurat_clusters[allMD$DataSet == "Proliferative"] <- allMD$DataSet[allMD$DataSet == "Proliferative"]
ClusterT <- allMD %>% group_by(g_seurat_clusters) %>% dplyr::count()
N <- dim(allMD)[1]
PrimaryTCellMD <- allMD %>% filter(batch %in% c("HeadNeck_Puram.rds")) %>% filter(TissueType == "Primary tumor tissue") %>% mutate(tcellgroup = "Primary")
MetastaticTCellMD <- allMD %>% filter(batch %in% c("HeadNeck_Puram.rds")) %>% filter(TissueType == "Metastatic tumor tissue") %>% mutate(tcellgroup = "Metastatic")
totalTargetTcellMD <- bind_rows(PrimaryTCellMD, MetastaticTCellMD)
matRoe <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 2),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 2)
rownames(matRoe) <- sort(unique(allMD$g_seurat_clusters))
colnames(matRoe) <- c("Primary", "Metastatic")
matCell <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 2),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 2)
rownames(matCell) <- sort(unique(allMD$g_seurat_clusters))
colnames(matCell) <- c("Primary", "Metastatic")
for(tempGroup in c("Primary", "Metastatic")){
    targetTcellMD <- totalTargetTcellMD %>% filter(tcellgroup == tempGroup)
    M <- dim(targetTcellMD)[1]
    for(gsci in unique(allMD$g_seurat_clusters)){
        k <- ClusterT$n[ClusterT$g_seurat_clusters == gsci]
        if(gsci %in% targetTcellMD$g_seurat_clusters){
            n <- dim(targetTcellMD %>% filter(g_seurat_clusters == gsci))[1]
        } else{
            n <- 0
        }
        tempRowName <- gsci
        tempColName <- tempGroup
        matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
        matCell[tempRowName, tempColName] <- n
    }
}
matRoe[matRoe > 2] <- 2
my.breaks <- seq(0, 2, by=0.01)
my.colors <- c(
    colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)/2),
    colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks)/2))
widthi = ncol(matRoe)/14 + 3
heighti = nrow(matRoe)/7 + 3
pdf(file.path(figurePath, paste0("HeadNeck_Puram_primaryVSmeta_roe.pdf")), width = widthi, height = heighti)
print(pheatmap(matRoe,
               color = my.colors,
               breaks = my.breaks,
               ## annotation_row = cellType_col,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figurePath, paste0("HeadNeck_Puram_primaryVSmeta_cell.pdf")), width = widthi, height = heighti)
print(pheatmap(matCell,
               color = my.colors,
               breaks = my.breaks,
               ## annotation_row = cellType_col,
               display_numbers = T,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               number_format = "%.0f",
               border_color = F))
dev.off()


###############################################################################
#'                                LUAD vs LUSC                               '#
###############################################################################

## test <- allMD %>% filter(TissueType %in% c("Primary tumor tissue")) %>% filter(CancerType %in% c("LUAD", "LUSC"))
## tt <- table(test$batch)
## paste0(names(tt), " (", tt, ")")

library(stringr)
## allMD$g_seurat_clusters[allMD$DataSet != "Proliferative"] <- paste0(allMD$DataSet[allMD$DataSet != "Proliferative"], "_c", allMD$seurat_clusters[allMD$DataSet != "Proliferative"])
## allMD$g_seurat_clusters[allMD$DataSet == "Proliferative"] <- allMD$DataSet[allMD$DataSet == "Proliferative"]
ClusterT <- allMD %>% group_by(g_seurat_clusters) %>% dplyr::count()
N <- dim(allMD)[1]
totalTargetTcellMD <- allMD %>% filter(TissueType %in% c("Primary tumor tissue")) %>% filter(CancerType %in% c("LUAD", "LUSC"))
matRoe <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) *
                        length(unique(totalTargetTcellMD$CancerType))),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = length(unique(totalTargetTcellMD$CancerType)))
rownames(matRoe) <- sort(unique(allMD$g_seurat_clusters))
colnames(matRoe) <- unique(totalTargetTcellMD$CancerType)
matCell <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) *
                        length(unique(totalTargetTcellMD$CancerType))),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = length(unique(totalTargetTcellMD$CancerType)))
rownames(matCell) <- sort(unique(allMD$g_seurat_clusters))
colnames(matCell) <- unique(totalTargetTcellMD$CancerType)
ggplotdata <- c()
for(tempCancerType in unique(totalTargetTcellMD$CancerType)){
    targetTcellMD <- totalTargetTcellMD %>% filter(CancerType == tempCancerType)
    M <- dim(targetTcellMD)[1]
    for(gsci in unique(allMD$g_seurat_clusters)){
        k <- ClusterT$n[ClusterT$g_seurat_clusters == gsci]
        if(gsci %in% targetTcellMD$g_seurat_clusters){
            n <- dim(targetTcellMD %>% filter(g_seurat_clusters == gsci))[1]
        } else{
            n <- 0
        }
        tempData <- tibble(CancerType = tempCancerType,
                           Cluster = gsci,
                           Roe = (n/M) / (k/N))
        ggplotdata <- bind_rows(ggplotdata, tempData)
        tempRowName <- gsci
        tempColName <- tempCancerType
        matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
        matCell[tempRowName, tempColName] <- n
    }
}
matRoe[matRoe>10] - 10
colnames(matRoe) <- paste0(colnames(matRoe), " (", colSums(matCell)[colnames(matRoe)], ")")
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
pdf(file.path(figurePath, paste0("LUADvsLUSC_roe_heatmap.pdf")), width = widthi, height = heighti)
print(pheatmap(matRoe,
               color = my.colors,
               breaks = my.breaks,
               ## annotation_col = cellType_col,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figurePath, paste0("LUADvsLUSC_cell_heatmap.pdf")), width = widthi * 1.5, height = heighti)
print(pheatmap(matCell,
               display_numbers = T,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = T,
               cluster_rows = F,
               number_format = "%.0f",
               border_color = F))
dev.off()
colnames(matRoe) <- paste0("3_", colnames(matRoe))
saveRDS(matRoe, file.path(figurePath, paste0('matRoe3.rds')))

###############################################################################
## Lung_GSE123902.rds
###############################################################################
allMD %>% filter(batch == "Lung_GSE123902.rds") %>% group_by(Patient, Sample, TissueType) %>% dplyr::count()


## allMD$g_seurat_clusters[allMD$DataSet != "Proliferative"] <- paste0(allMD$DataSet[allMD$DataSet != "Proliferative"], "_c", allMD$seurat_clusters[allMD$DataSet != "Proliferative"])
## allMD$g_seurat_clusters[allMD$DataSet == "Proliferative"] <- allMD$DataSet[allMD$DataSet == "Proliferative"]
ClusterT <- allMD %>% group_by(g_seurat_clusters) %>% dplyr::count()
N <- dim(allMD)[1]
BrainMetTCellMD <- allMD %>% filter(Patient %in% c("Lung_GSE123902_LX255B", "Lung_GSE123902_LX681", "Lung_GSE123902_LX701")) %>% mutate(tcellgroup = "BrainMet")
AdrenalMetTCellMD <- allMD %>% filter(Patient %in% c("Lung_GSE123902_LX699")) %>% mutate(tcellgroup = "AdrenalMet")
PrimaryTCellMD <- allMD %>% filter(Patient %in% c("Lung_GSE123902_LX661")) %>% mutate(tcellgroup = "Primary")
totalTargetTcellMD <- bind_rows(BrainMetTCellMD, AdrenalMetTCellMD, PrimaryTCellMD)
matRoe <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 3),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 3)
rownames(matRoe) <- sort(unique(allMD$g_seurat_clusters))
colnames(matRoe) <- c("BrainMet", "AdrenalMet", "Primary")
matCell <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 3),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 3)
rownames(matCell) <- sort(unique(allMD$g_seurat_clusters))
colnames(matCell) <- c("BrainMet", "AdrenalMet", "Primary")
for(tempGroup in c("BrainMet", "AdrenalMet", "Primary")){
    targetTcellMD <- totalTargetTcellMD %>% filter(tcellgroup == tempGroup)
    M <- dim(targetTcellMD)[1]
    for(gsci in unique(allMD$g_seurat_clusters)){
        k <- ClusterT$n[ClusterT$g_seurat_clusters == gsci]
        if(gsci %in% targetTcellMD$g_seurat_clusters){
            n <- dim(targetTcellMD %>% filter(g_seurat_clusters == gsci))[1]
        } else{
            n <- 0
        }
        tempRowName <- gsci
        tempColName <- tempGroup
        matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
        matCell[tempRowName, tempColName] <- n
    }
}
matRoe[matRoe>10] - 10
colnames(matRoe) <- paste0(colnames(matRoe), " (", colSums(matCell)[colnames(matRoe)], ")")
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
pdf(file.path(figurePath, paste0("Lung_GSE123902_BrainMetaVSPriamry_roe.pdf")), width = widthi, height = heighti)
print(pheatmap(matRoe,
               color = my.colors,
               breaks = my.breaks,
               ## annotation_row = cellType_col,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figurePath, paste0("Lung_GSE123902_BrainMetaVSPriamry_cell.pdf")), width = widthi, height = heighti)
print(pheatmap(matCell,
               color = my.colors,
               breaks = my.breaks,
               display_numbers = T,
               number_format = "%.0f",
               ## annotation_row = cellType_col,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
colnames(matRoe) <- paste0("5_", colnames(matRoe))
saveRDS(matRoe, file.path(figurePath, paste0('matRoe5.rds')))
## Stage 1, primary treatment naive,  adjacent normal vs smoker vs never smoker

## never smoker sample: LUAD-T-cells.rds_P5-T2-Epcam-neg
## smoker sample: Lung_GSE123902.rds_GSM3516663 Lung_GSE123902.rds_GSM3516667 Lung_GSE123902.rds_GSM3516672 LUAD-T-cells.rds_P2-T7-Epcam-neg
## adjacent normal sample, never: Lung_GSE123902.rds_GSM3516666 Lung_GSE123902.rds_GSM3516676
## adjacent normal sample, smoker: Lung_GSE123902.rds_GSM3516673 Lung_GSE123902.rds_GSM3516675
## allMD$g_seurat_clusters[allMD$DataSet != "Proliferative"] <- paste0(allMD$DataSet[allMD$DataSet != "Proliferative"], "_c", allMD$seurat_clusters[allMD$DataSet != "Proliferative"])
## allMD$g_seurat_clusters[allMD$DataSet == "Proliferative"] <- allMD$DataSet[allMD$DataSet == "Proliferative"]
ClusterT <- allMD %>% group_by(g_seurat_clusters) %>% dplyr::count()
N <- dim(allMD)[1]
T_N_TCellMD <- allMD %>% filter(Sample %in% c("LUAD-T-cells.rds_P5-T2-Epcam-neg")) %>% mutate(tcellgroup = "Tumor_Never")
T_S_TCellMD <- allMD %>% filter(Sample %in% c("Lung_GSE123902.rds_GSM3516663", "Lung_GSE123902.rds_GSM3516667", "Lung_GSE123902.rds_GSM3516672", "LUAD-T-cells.rds_P2-T7-Epcam-neg")) %>% mutate(tcellgroup = "Tumor_Smoker")
N_N_TCellMD <- allMD %>% filter(Sample %in% c("Lung_GSE123902.rds_GSM3516666", "Lung_GSE123902.rds_GSM3516676")) %>% mutate(tcellgroup = "AdjacentNormal_Never")
N_S_TCellMD <- allMD %>% filter(Sample %in% c("Lung_GSE123902.rds_GSM3516673" , "Lung_GSE123902.rds_GSM3516675")) %>% mutate(tcellgroup = "AdjacentNormal_Smoker")
totalTargetTcellMD <- bind_rows(T_N_TCellMD, T_S_TCellMD, N_N_TCellMD, N_S_TCellMD)
matRoe <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 4),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 4)
rownames(matRoe) <- sort(unique(allMD$g_seurat_clusters))
colnames(matRoe) <- c("Tumor_Never", "Tumor_Smoker", "AdjacentNormal_Never", "AdjacentNormal_Smoker")
matCell <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 4),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 4)
rownames(matCell) <- sort(unique(allMD$g_seurat_clusters))
colnames(matCell) <- c("Tumor_Never", "Tumor_Smoker", "AdjacentNormal_Never", "AdjacentNormal_Smoker")
for(tempGroup in c("Tumor_Never", "Tumor_Smoker", "AdjacentNormal_Never", "AdjacentNormal_Smoker")){
    targetTcellMD <- totalTargetTcellMD %>% filter(tcellgroup == tempGroup)
    M <- dim(targetTcellMD)[1]
    for(gsci in unique(allMD$g_seurat_clusters)){
        k <- ClusterT$n[ClusterT$g_seurat_clusters == gsci]
        if(gsci %in% targetTcellMD$g_seurat_clusters){
            n <- dim(targetTcellMD %>% filter(g_seurat_clusters == gsci))[1]
        } else{
            n <- 0
        }
        tempRowName <- gsci
        tempColName <- tempGroup
        matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
        matCell[tempRowName, tempColName] <- n
    }
}
matRoe[matRoe>10] - 10
colnames(matRoe) <- paste0(colnames(matRoe), " (", colSums(matCell)[colnames(matRoe)], ")")
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
pdf(file.path(figurePath, paste0("LUAD_smokerVSnever_roe.pdf")), width = widthi, height = heighti)
print(pheatmap(matRoe,
               color = my.colors,
               breaks = my.breaks,
               ## annotation_row = cellType_col,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figurePath, paste0("LUAD_smokerVSnever_cell.pdf")), width = widthi, height = heighti)
print(pheatmap(matCell,
               color = my.colors,
               breaks = my.breaks,
               display_numbers = T,
               number_format = "%.0f",
               ## annotation_row = cellType_col,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
colnames(matRoe) <- paste0("6_", colnames(matRoe))
saveRDS(matRoe, file.path(figurePath, paste0('matRoe6.rds')))

###############################################################################
## Lung_GSE127465.rds
###############################################################################
allMD %>% filter(batch == "Lung_GSE127465.rds") %>% group_by(Patient) %>% dplyr::count()

## allMD$g_seurat_clusters[allMD$DataSet != "Proliferative"] <- paste0(allMD$DataSet[allMD$DataSet != "Proliferative"], "_c", allMD$seurat_clusters[allMD$DataSet != "Proliferative"])
## allMD$g_seurat_clusters[allMD$DataSet == "Proliferative"] <- allMD$DataSet[allMD$DataSet == "Proliferative"]
ClusterT <- allMD %>% group_by(g_seurat_clusters) %>% dplyr::count()
N <- dim(allMD)[1]
SquamousTCellMD <- allMD %>% filter(Patient %in% c("Lung_GSE127465_p1")) %>% mutate(tcellgroup = "Squamous")
AdenoTCellMD <- allMD %>% filter(Patient %in% c("Lung_GSE127465_p4", "Lung_GSE127465_p5", "Lung_GSE127465_p6", "Lung_GSE127465_p7")) %>% mutate(tcellgroup = "Adeno")
totalTargetTcellMD <- bind_rows(SquamousTCellMD, AdenoTCellMD)
matRoe <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 2),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 2)
rownames(matRoe) <- sort(unique(allMD$g_seurat_clusters))
colnames(matRoe) <- c("Squamous", "Adeno")
matCell <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 2),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 2)
rownames(matCell) <- sort(unique(allMD$g_seurat_clusters))
colnames(matCell) <- c("Squamous", "Adeno")
for(tempGroup in c("Squamous", "Adeno")){
    targetTcellMD <- totalTargetTcellMD %>% filter(tcellgroup == tempGroup)
    M <- dim(targetTcellMD)[1]
    for(gsci in unique(allMD$g_seurat_clusters)){
        k <- ClusterT$n[ClusterT$g_seurat_clusters == gsci]
        if(gsci %in% targetTcellMD$g_seurat_clusters){
            n <- dim(targetTcellMD %>% filter(g_seurat_clusters == gsci))[1]
        } else{
            n <- 0
        }
        tempRowName <- gsci
        tempColName <- tempGroup
        matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
        matCell[tempRowName, tempColName] <- n
    }
}
matRoe[matRoe>10] - 10
colnames(matRoe) <- paste0(colnames(matRoe), " (", colSums(matCell)[colnames(matRoe)], ")")
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
pdf(file.path(figurePath, paste0("Lung_GSE127465_SquamousVSAdeno_roe.pdf")), width = widthi, height = heighti)
print(pheatmap(matRoe,
               color = my.colors,
               breaks = my.breaks,
               ## annotation_row = cellType_col,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figurePath, paste0("Lung_GSE127465_SquamousVSAdeno_cell.pdf")), width = widthi, height = heighti)
print(pheatmap(matCell,
               color = my.colors,
               breaks = my.breaks,
               display_numbers = T,
               number_format = "%.0f",
               ## annotation_row = cellType_col,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
colnames(matRoe) <- paste0("4_", colnames(matRoe))
saveRDS(matRoe, file.path(figurePath, paste0('matRoe4.rds')))

###############################################################################
#'                                     Liver_GSE125449
###############################################################################

test <- allMD %>% filter(batch == "Liver_GSE125449.rds") %>% group_by(Patient) %>% dplyr::count()



## allMD$g_seurat_clusters[allMD$DataSet != "Proliferative"] <- paste0(allMD$DataSet[allMD$DataSet != "Proliferative"], "_c", allMD$seurat_clusters[allMD$DataSet != "Proliferative"])
## allMD$g_seurat_clusters[allMD$DataSet == "Proliferative"] <- allMD$DataSet[allMD$DataSet == "Proliferative"]
ClusterT <- allMD %>% group_by(g_seurat_clusters) %>% dplyr::count()
N <- dim(allMD)[1]
iCCA_TCellMD <- allMD %>% filter(Patient %in% paste0("Liver_GSE125449_", c("S09_P04_LCP25", "S19_P11_LCP39", "S305_P06_LCP56", "S300_P02_LCP60", "S365_P22_LCP66"))) %>% mutate(tcellgroup = "iCCA")
HCC_TCellMD <- allMD %>% filter(Patient %in% paste0("Liver_GSE125449_", c("S02_P01_LCP21", "S07_P02_LCP28", "S12_P07_LCP30", "S351_P10_LCP34", "S364_P21_LCP65"))) %>% mutate(tcellgroup = "HCC")
totalTargetTcellMD <- bind_rows(iCCA_TCellMD, HCC_TCellMD)
matRoe <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 2),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 2)
rownames(matRoe) <- sort(unique(allMD$g_seurat_clusters))
colnames(matRoe) <- c("iCCA", "HCC")
matCell <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 2),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 2)
rownames(matCell) <- sort(unique(allMD$g_seurat_clusters))
colnames(matCell) <- c("iCCA", "HCC")
for(tempGroup in c("iCCA", "HCC")){
    targetTcellMD <- totalTargetTcellMD %>% filter(tcellgroup == tempGroup)
    M <- dim(targetTcellMD)[1]
    for(gsci in unique(allMD$g_seurat_clusters)){
        k <- ClusterT$n[ClusterT$g_seurat_clusters == gsci]
        if(gsci %in% targetTcellMD$g_seurat_clusters){
            n <- dim(targetTcellMD %>% filter(g_seurat_clusters == gsci))[1]
        } else{
            n <- 0
        }
        tempRowName <- gsci
        tempColName <- tempGroup
        matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
        matCell[tempRowName, tempColName] <- n
    }
}
matRoe[matRoe>10] = 10
colnames(matRoe) <- paste0(colnames(matRoe), " (", colSums(matCell)[colnames(matRoe)], ")")
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
pdf(file.path(figurePath, paste0("Lung_GSE125449_iCCAvsHCC_roe.pdf")), width = widthi, height = heighti)
print(pheatmap(matRoe,
               color = my.colors,
               breaks = my.breaks,
               ## annotation_row = cellType_col,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figurePath, paste0("Lung_GSE125449_iCCAvsHCC_cell.pdf")), width = widthi, height = heighti)
print(pheatmap(matCell,
               color = my.colors,
               breaks = my.breaks,
               display_numbers = T,
               number_format = "%.0f",
               ## annotation_row = cellType_col,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
colnames(matRoe) <- paste0("13_", colnames(matRoe))
saveRDS(matRoe, file.path(figurePath, paste0('matRoe13.rds')))




###############################################################################
#'                                 Mel 115978                                '#
###############################################################################


test <- allMD %>% filter(batch == "Melanoma_GSE115978.rds") %>% group_by(Patient) %>% dplyr::count()

## Treatment and Lesion type comparison
## pici: Post-ICI
## allMD$g_seurat_clusters[allMD$DataSet != "Proliferative"] <- paste0(allMD$DataSet[allMD$DataSet != "Proliferative"], "_c", allMD$seurat_clusters[allMD$DataSet != "Proliferative"])
## allMD$g_seurat_clusters[allMD$DataSet == "Proliferative"] <- allMD$DataSet[allMD$DataSet == "Proliferative"]
ClusterT <- allMD %>% group_by(g_seurat_clusters) %>% dplyr::count()
N <- dim(allMD)[1]
M_pici_TCellMD <- allMD %>% filter(Patient %in% paste0("Mel", c("60", "72", "74", "75", "78", "88", "94", "126", "110", "121.1", "106", "75", "98", "102"))) %>% mutate(tcellgroup = "M_pici")
M_none_TCellMD <- allMD %>% filter(Patient %in% paste0("Mel", c("71", "79", "80", "81", "82", "116", "103", "112"))) %>% mutate(tcellgroup = "M_none")
P_none_TCellMD <- allMD %>% filter(Patient %in% paste0("Mel", c("84", "129", "105"))) %>% mutate(tcellgroup = "P_none")
totalTargetTcellMD <- bind_rows(M_pici_TCellMD, M_none_TCellMD, P_none_TCellMD)
matRoe <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 3),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 3)
rownames(matRoe) <- sort(unique(allMD$g_seurat_clusters))
colnames(matRoe) <- c("M_pici", "M_none", "P_none")
matCell <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 3),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 3)
rownames(matCell) <- sort(unique(allMD$g_seurat_clusters))
colnames(matCell) <- c("M_pici", "M_none", "P_none")
for(tempGroup in c("M_pici", "M_none", "P_none")){
    targetTcellMD <- totalTargetTcellMD %>% filter(tcellgroup == tempGroup)
    M <- dim(targetTcellMD)[1]
    for(gsci in unique(allMD$g_seurat_clusters)){
        k <- ClusterT$n[ClusterT$g_seurat_clusters == gsci]
        if(gsci %in% targetTcellMD$g_seurat_clusters){
            n <- dim(targetTcellMD %>% filter(g_seurat_clusters == gsci))[1]
        } else{
            n <- 0
        }
        tempRowName <- gsci
        tempColName <- tempGroup
        matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
        matCell[tempRowName, tempColName] <- n
    }
}
matRoe[matRoe>10] - 10
colnames(matRoe) <- paste0(colnames(matRoe), " (", colSums(matCell)[colnames(matRoe)], ")")
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
pdf(file.path(figurePath, paste0("Melanoma_GSE115978_treatmentAndCancerType_roe.pdf")), width = widthi, height = heighti)
print(pheatmap(matRoe,
               color = my.colors,
               breaks = my.breaks,
               ## annotation_row = cellType_col,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figurePath, paste0("Melanoma_GSE115978_treatmentAndCancerType_cell.pdf")), width = widthi, height = heighti)
print(pheatmap(matCell,
               color = my.colors,
               breaks = my.breaks,
               display_numbers = T,
               number_format = "%.0f",
               ## annotation_row = cellType_col,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
colnames(matRoe) <- paste0("11_", colnames(matRoe))
saveRDS(matRoe, file.path(figurePath, paste0('matRoe11.rds')))

###############################################################################
#'                                   120575                                  '#
###############################################################################


test <- allMD %>% filter(batch == "Melanoma_GSE120575.rds") %>% group_by(Patient) %>% dplyr::count()


annotationT <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/data/annotation_melanoma_120575.txt")
## pre vs post
library(stringr)
## allMD$g_seurat_clusters[allMD$DataSet != "Proliferative"] <- paste0(allMD$DataSet[allMD$DataSet != "Proliferative"], "_c", allMD$seurat_clusters[allMD$DataSet != "Proliferative"])
## allMD$g_seurat_clusters[allMD$DataSet == "Proliferative"] <- allMD$DataSet[allMD$DataSet == "Proliferative"]
ClusterT <- allMD %>% group_by(g_seurat_clusters) %>% dplyr::count()
N <- dim(allMD)[1]
responderT <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/data/1-s2.0-S0092867418313941-mmc1.txt")
colnames(responderT) <- c("Patient", "Sample", "Status")
tempSampleNames <- annotationT %>% filter(Therapy == "PD1") %>% pull(samples)
responderTCellMD <- allMD %>% filter(batch == "Melanoma_GSE120575.rds") %>%
    filter(Sample %in% paste0("Melanoma_GSE120575.rds_", responderT %>% filter(Status == "R") %>% pull(Sample))) %>%
    filter(Sample %in% paste0("Melanoma_GSE120575.rds_", tempSampleNames)) %>%
    mutate(tcellgroup = "R")
nresponderTCellMD <- allMD %>% filter(batch == "Melanoma_GSE120575.rds") %>%
    filter(Sample %in% paste0("Melanoma_GSE120575.rds_", responderT %>% filter(Status == "NR") %>% pull(Sample))) %>%
    filter(Sample %in% paste0("Melanoma_GSE120575.rds_", tempSampleNames)) %>%
    mutate(tcellgroup = "NR")
totalTargetTcellMD <- bind_rows(responderTCellMD, nresponderTCellMD)
matRoe <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 4),
                  nrow = length(unique(allMD$g_seurat_clusters)),
                  ncol = 4)
rownames(matRoe) <- sort(unique(allMD$g_seurat_clusters))
colnames(matRoe) <- c(paste0("R-", c("Pre", "Post")), paste0("NR-", c("Pre", "Post")))
matCell <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 4),
                  nrow = length(unique(allMD$g_seurat_clusters)),
                  ncol = 4)
rownames(matCell) <- sort(unique(allMD$g_seurat_clusters))
colnames(matCell) <- c(paste0("R-", c("Pre", "Post")), paste0("NR-", c("Pre", "Post")))
ggplotdata <- c()
for(tempGroup in c("R", "NR")){
    targetTcellMD <- totalTargetTcellMD %>% filter(tcellgroup == tempGroup)
    for(tempSampleType in c("Pre", "Post")){
        tempSampleTargetTcellMD <- targetTcellMD[str_detect(targetTcellMD$Sample, tempSampleType),]
        if(tempSampleType == "Pre"){
            tempSampleTargetTcellMD <- tempSampleTargetTcellMD %>% filter(Patient %in% c("Melanoma_GSE120575_PP1", "Melanoma_GSE120575_PP6"))
        }
        M <- dim(tempSampleTargetTcellMD)[1]
        for(gsci in unique(allMD$g_seurat_clusters)){
            k <- ClusterT$n[ClusterT$g_seurat_clusters == gsci]
            if(gsci %in% tempSampleTargetTcellMD$g_seurat_clusters){
                n <- dim(tempSampleTargetTcellMD %>% filter(g_seurat_clusters == gsci))[1]
            } else{
                n <- 0
            }
            tempData <- tibble(Group = tempGroup,
                                SampleType = tempSampleType,
                                Cluster = gsci,
                                Roe = (n/M) / (k/N))
            ggplotdata <- bind_rows(ggplotdata, tempData)
            tempRowName <- gsci
            tempColName <- paste0(tempGroup, "-", tempSampleType)
            matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
            matCell[tempRowName, tempColName] <- n
        }
    }
}
matRoe[matRoe>10] = 10
colnames(matRoe) <- paste0(colnames(matRoe), " (", colSums(matCell)[colnames(matRoe)], ")")
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
pdf(file.path(figurePath, paste0("comparison2_PD1_roe_heatmap.pdf")), width = widthi, height = heighti)
print(pheatmap(matRoe,
                color = my.colors,
                breaks = my.breaks,
                ## annotation_col = cellType_col,
                show_colnames = T,
                show_rownames = T,
                cluster_cols = F,
                cluster_rows = F,
                border_color = F))
dev.off()
pdf(file.path(figurePath, paste0("comparison2_PD1_cell_heatmap.pdf")), width = widthi * 1.5, height = heighti)
print(pheatmap(matCell,
                display_numbers = T,
                show_colnames = T,
                show_rownames = T,
                cluster_cols = F,
                cluster_rows = F,
                number_format = "%.0f",
                border_color = F))
dev.off()




## pre vs post
library(stringr)
## allMD$g_seurat_clusters[allMD$DataSet != "Proliferative"] <- paste0(allMD$DataSet[allMD$DataSet != "Proliferative"], "_c", allMD$seurat_clusters[allMD$DataSet != "Proliferative"])
## allMD$g_seurat_clusters[allMD$DataSet == "Proliferative"] <- allMD$DataSet[allMD$DataSet == "Proliferative"]
ClusterT <- allMD %>% group_by(g_seurat_clusters) %>% dplyr::count()
N <- dim(allMD)[1]
responderT <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/data/1-s2.0-S0092867418313941-mmc1.txt")
colnames(responderT) <- c("Patient", "Sample", "Status")
responderTCellMD <- allMD %>% filter(batch == "Melanoma_GSE120575.rds") %>%
    filter(Sample %in% paste0("Melanoma_GSE120575.rds_", responderT %>% filter(Status == "R") %>% pull(Sample))) %>%
    mutate(tcellgroup = "R")
nresponderTCellMD <- allMD %>% filter(batch == "Melanoma_GSE120575.rds") %>%
    filter(Sample %in% paste0("Melanoma_GSE120575.rds_", responderT %>% filter(Status == "NR") %>% pull(Sample))) %>%
    mutate(tcellgroup = "NR")
totalTargetTcellMD <- bind_rows(responderTCellMD, nresponderTCellMD)
matRoe <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 4),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 4)
rownames(matRoe) <- sort(unique(allMD$g_seurat_clusters))
colnames(matRoe) <- c(paste0("R-", c("Pre", "Post")), paste0("NR-", c("Pre", "Post")))
matCell <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 4),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 4)
rownames(matCell) <- sort(unique(allMD$g_seurat_clusters))
colnames(matCell) <- c(paste0("R-", c("Pre", "Post")), paste0("NR-", c("Pre", "Post")))
ggplotdata <- c()
for(tempGroup in c("R", "NR")){
    targetTcellMD <- totalTargetTcellMD %>% filter(tcellgroup == tempGroup)
    for(tempSampleType in c("Pre", "Post")){
        tempSampleTargetTcellMD <- targetTcellMD[str_detect(targetTcellMD$Sample, tempSampleType),]
        if(tempSampleType == "Pre"){
            tempSampleTargetTcellMD <- tempSampleTargetTcellMD %>% filter(Patient %in% c("Melanoma_GSE120575_PP1", "Melanoma_GSE120575_PP6"))
        }
        M <- dim(tempSampleTargetTcellMD)[1]
        for(gsci in unique(allMD$g_seurat_clusters)){
            k <- ClusterT$n[ClusterT$g_seurat_clusters == gsci]
            if(gsci %in% tempSampleTargetTcellMD$g_seurat_clusters){
                n <- dim(tempSampleTargetTcellMD %>% filter(g_seurat_clusters == gsci))[1]
            } else{
                n <- 0
            }
            tempData <- tibble(Group = tempGroup,
                               SampleType = tempSampleType,
                               Cluster = gsci,
                               Roe = (n/M) / (k/N))
            ggplotdata <- bind_rows(ggplotdata, tempData)
            tempRowName <- gsci
            tempColName <- paste0(tempGroup, "-", tempSampleType)
            matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
            matCell[tempRowName, tempColName] <- n
        }
    }
}
matRoe[matRoe>10] = 10
colnames(matRoe) <- paste0(colnames(matRoe), " (", colSums(matCell)[colnames(matRoe)], ")")
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
pdf(file.path(figurePath, paste0("comparison2_roe_heatmap.pdf")), width = widthi, height = heighti)
print(pheatmap(matRoe,
               color = my.colors,
               breaks = my.breaks,
               ## annotation_col = cellType_col,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
## pdf(file.path(figurePath, paste0("comparison2_cell_heatmap.pdf")), width = widthi * 1.5, height = heighti)
## print(pheatmap(matCell,
##                display_numbers = T,
##                show_colnames = T,
##                show_rownames = T,
##                cluster_cols = T,
##                cluster_rows = F,
##                number_format = "%.0f",
##                border_color = F))
## dev.off()
colnames(matRoe) <- paste0("9_", colnames(matRoe))
saveRDS(matRoe, file.path(figurePath, paste0('matRoe9.rds')))



annotationT <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/data/annotation_melanoma_120575.txt")
annotationT$isPost <- str_detect(annotationT$samples, "Post.*")
annotationT <- annotationT %>% filter(isPost) %>% filter(! samples %in% c("Post_P1", "Post_P1_2", "Post_P6"))
PD1CTLA4_samples <- annotationT %>% filter(Therapy == "CTLA4+PD1") %>% pull(samples)
PD1_samples <- annotationT %>% filter(Therapy == "PD1") %>% pull(samples)
CTLA4_samples <- c("Pre_P1", "Pre_P6")
## pre vs post
library(stringr)
## allMD$g_seurat_clusters[allMD$DataSet != "Proliferative"] <- paste0(allMD$DataSet[allMD$DataSet != "Proliferative"], "_c", allMD$seurat_clusters[allMD$DataSet != "Proliferative"])
## allMD$g_seurat_clusters[allMD$DataSet == "Proliferative"] <- allMD$DataSet[allMD$DataSet == "Proliferative"]
ClusterT <- allMD %>% group_by(g_seurat_clusters) %>% dplyr::count()
N <- dim(allMD)[1]
PD1CTLA4_TCellMD <- allMD %>% filter(Sample %in% paste0("Melanoma_GSE120575.rds_", PD1CTLA4_samples)) %>% mutate(tcellgroup = "PD1CTLA4")
PD1_TCellMD <- allMD %>% filter(Sample %in% paste0("Melanoma_GSE120575.rds_", PD1_samples)) %>% mutate(tcellgroup = "PD1")
CTLA4_TCellMD <- allMD %>% filter(Sample %in% paste0("Melanoma_GSE120575.rds_", CTLA4_samples)) %>% mutate(tcellgroup = "CTLA4")
totalTargetTcellMD <- bind_rows(PD1CTLA4_TCellMD, PD1_TCellMD, CTLA4_TCellMD)
matRoe <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 3),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 3)
rownames(matRoe) <- sort(unique(allMD$g_seurat_clusters))
colnames(matRoe) <- c("PD1CTLA4", "PD1", "CTLA4")
matCell <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 3),
                  nrow = length(unique(allMD$g_seurat_clusters)),
                  ncol = 3)
rownames(matCell) <- sort(unique(allMD$g_seurat_clusters))
colnames(matCell) <- c("PD1CTLA4", "PD1", "CTLA4")
ggplotdata <- c()
for(tempGroup in c("PD1CTLA4", "PD1", "CTLA4")){
    targetTcellMD <- totalTargetTcellMD %>% filter(tcellgroup == tempGroup)
    M <- dim(targetTcellMD)[1]
    for(gsci in unique(allMD$g_seurat_clusters)){
        k <- ClusterT$n[ClusterT$g_seurat_clusters == gsci]
        if(gsci %in% targetTcellMD$g_seurat_clusters){
            n <- dim(targetTcellMD %>% filter(g_seurat_clusters == gsci))[1]
        } else{
            n <- 0
        }
        tempData <- tibble(Group = tempGroup,
                           Cluster = gsci,
                           Roe = (n/M) / (k/N))
        ggplotdata <- bind_rows(ggplotdata, tempData)
        tempRowName <- gsci
        tempColName <- tempGroup
        matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
        matCell[tempRowName, tempColName] <- n
    }
}
matRoe[matRoe>10] = 10
colnames(matRoe) <- paste0(colnames(matRoe), " (", colSums(matCell)[colnames(matRoe)], ")")
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
pdf(file.path(figurePath, paste0("Melanoma_GSE120575_treatment_roe_heatmap.pdf")), width = widthi, height = heighti)
print(pheatmap(matRoe,
               color = my.colors,
               breaks = my.breaks,
               ## annotation_col = cellType_col,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figurePath, paste0("Melanoma_GSE120575_treatment_cell_heatmap.pdf")), width = widthi * 1.5, height = heighti)
print(pheatmap(matCell,
               display_numbers = T,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = T,
               cluster_rows = F,
               number_format = "%.0f",
               border_color = F))
dev.off()
colnames(matRoe) <- paste0("10_", colnames(matRoe))
saveRDS(matRoe, file.path(figurePath, paste0('matRoe10.rds')))



test <- allMD %>% filter(batch == "Melanoma_GSE120575.rds") %>% group_by(Sample) %>% dplyr::count()

annotationT <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/data/annotation_melanoma_120575.txt")
library(stringr)
## allMD$g_seurat_clusters[allMD$DataSet != "Proliferative"] <- paste0(allMD$DataSet[allMD$DataSet != "Proliferative"], "_c", allMD$seurat_clusters[allMD$DataSet != "Proliferative"])
## allMD$g_seurat_clusters[allMD$DataSet == "Proliferative"] <- allMD$DataSet[allMD$DataSet == "Proliferative"]
ClusterT <- allMD %>% group_by(g_seurat_clusters) %>% dplyr::count()
N <- dim(allMD)[1]
NRAS_TCellMD <- allMD %>% filter(Sample %in% paste0("Melanoma_GSE120575.rds_Pre_", c("P12"))) %>% mutate(tcellgroup = "NRAS")
BRAF_TCellMD <- allMD %>% filter(Sample %in% paste0("Melanoma_GSE120575.rds_Pre_", c("P20", "P28", "P35"))) %>% mutate(tcellgroup = "BRAF")
totalTargetTcellMD <- bind_rows(NRAS_TCellMD, BRAF_TCellMD)
matRoe <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 2),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 2)
rownames(matRoe) <- sort(unique(allMD$g_seurat_clusters))
colnames(matRoe) <- c("NRAS", "BRAF")
matCell <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 2),
                  nrow = length(unique(allMD$g_seurat_clusters)),
                  ncol = 2)
rownames(matCell) <- sort(unique(allMD$g_seurat_clusters))
colnames(matCell) <- c("NRAS", "BRAF")
ggplotdata <- c()
for(tempGroup in c("NRAS", "BRAF")){
    targetTcellMD <- totalTargetTcellMD %>% filter(tcellgroup == tempGroup)
    M <- dim(targetTcellMD)[1]
    for(gsci in unique(allMD$g_seurat_clusters)){
        k <- ClusterT$n[ClusterT$g_seurat_clusters == gsci]
        if(gsci %in% targetTcellMD$g_seurat_clusters){
            n <- dim(targetTcellMD %>% filter(g_seurat_clusters == gsci))[1]
        } else{
            n <- 0
        }
        tempData <- tibble(Group = tempGroup,
                           Cluster = gsci,
                           Roe = (n/M) / (k/N))
        ggplotdata <- bind_rows(ggplotdata, tempData)
        tempRowName <- gsci
        tempColName <- tempGroup
        matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
        matCell[tempRowName, tempColName] <- n
    }
}
matRoe[matRoe>10] = 10
colnames(matRoe) <- paste0(colnames(matRoe), " (", colSums(matCell)[colnames(matRoe)], ")")
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
pdf(file.path(figurePath, paste0("Melanoma_GSE120575_NRASvsBRAF_roe_heatmap.pdf")), width = widthi, height = heighti)
print(pheatmap(matRoe,
               color = my.colors,
               breaks = my.breaks,
               ## annotation_col = cellType_col,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figurePath, paste0("Melanoma_GSE120575_NRASvsBRAF_cell_heatmap.pdf")), width = widthi * 1.5, height = heighti)
print(pheatmap(matCell,
               display_numbers = T,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = T,
               cluster_rows = F,
               number_format = "%.0f",
               border_color = F))
dev.off()
colnames(matRoe) <- paste0("8_", colnames(matRoe))
saveRDS(matRoe, file.path(figurePath, paste0('matRoe8.rds')))




###############################################################################
#'                                     colon                                 '#
###############################################################################

ClusterT <- allMD %>% group_by(g_seurat_clusters) %>% dplyr::count()
N <- dim(allMD)[1]
ReC_Patients <- paste0("Leukocyte_", c("P0309", "P0202", "P0411", "P1025", "P1026"))
ReC_TCellMD <- allMD %>%
    filter(batch == "Leukocyte10x.T.rds" | batch == "LeukocyteSmartseq2.T.rds") %>%
    filter(Patient %in% ReC_Patients & TissueType == "Primary tumor tissue") %>%
    mutate(tcellgroup = "ReC")
RCC_TCellMD <- allMD %>%
    filter(batch == "Leukocyte10x.T.rds" | batch == "LeukocyteSmartseq2.T.rds") %>%
    filter((!Patient %in% ReC_Patients) & TissueType == "Primary tumor tissue") %>%
    mutate(tcellgroup = "RCC")
totalTargetTcellMD <- bind_rows(ReC_TCellMD, RCC_TCellMD)
matRoe <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 2),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 2)
rownames(matRoe) <- sort(unique(allMD$g_seurat_clusters))
colnames(matRoe) <- c("ReC", "RCC")
matCell <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 2),
                  nrow = length(unique(allMD$g_seurat_clusters)),
                  ncol = 2)
rownames(matCell) <- sort(unique(allMD$g_seurat_clusters))
colnames(matCell) <- c("ReC", "RCC")
ggplotdata <- c()
for(tempGroup in c("ReC", "RCC")){
    targetTcellMD <- totalTargetTcellMD %>% filter(tcellgroup == tempGroup)
    M <- dim(targetTcellMD)[1]
    for(gsci in unique(allMD$g_seurat_clusters)){
        k <- ClusterT$n[ClusterT$g_seurat_clusters == gsci]
        if(gsci %in% targetTcellMD$g_seurat_clusters){
            n <- dim(targetTcellMD %>% filter(g_seurat_clusters == gsci))[1]
        } else{
            n <- 0
        }
        tempData <- tibble(Group = tempGroup,
                           Cluster = gsci,
                           Roe = (n/M) / (k/N))
        ggplotdata <- bind_rows(ggplotdata, tempData)
        tempRowName <- gsci
        tempColName <- tempGroup
        matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
        matCell[tempRowName, tempColName] <- n
    }
}
matRoe[matRoe>10] = 10
colnames(matRoe) <- paste0(colnames(matRoe), " (", colSums(matCell)[colnames(matRoe)], ")")
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
pdf(file.path(figurePath, paste0("colon_ReCvsRCC_roe_heatmap.pdf")), width = widthi, height = heighti)
print(pheatmap(matRoe,
               color = my.colors,
               breaks = my.breaks,
               ## annotation_col = cellType_col,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figurePath, paste0("colon_ReCvsRCC_cell_heatmap.pdf")), width = widthi * 1.5, height = heighti)
print(pheatmap(matCell,
               display_numbers = T,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = T,
               cluster_rows = F,
               number_format = "%.0f",
               border_color = F))
dev.off()
colnames(matRoe) <- paste0("14_", colnames(matRoe))
saveRDS(matRoe, file.path(figurePath, paste0('matRoe14.rds')))

## allMD$g_seurat_clusters[allMD$DataSet != "Proliferative"] <- paste0(allMD$DataSet[allMD$DataSet != "Proliferative"], "_c", allMD$seurat_clusters[allMD$DataSet != "Proliferative"])
## allMD$g_seurat_clusters[allMD$DataSet == "Proliferative"] <- allMD$DataSet[allMD$DataSet == "Proliferative"]
ClusterT <- allMD %>% group_by(g_seurat_clusters) %>% dplyr::count()
N <- dim(allMD)[1]
ReC_Patients <- paste0("Leukocyte_", c("P0309", "P0202", "P0411", "P1025", "P1026"))
Moderate_Patients <- paste0("Leukocyte_", c("P0309", "P0411", "P0720", "P0728", "P0413", "P0202", "P0323", "P0410", "P0613", "P1026", "P0123"))
LowOrModerate_Patients <- paste0("Leukocyte_", c("P0305", "P1212", "P1228", "P0825", "P0408"))
High_Patients <- paste0("Leukocyte_", c("P0104", "P1025"))
Moderate_TCellMD <- allMD %>% filter(Patient %in% Moderate_Patients & (! Patient %in% ReC_Patients) & TissueType == "Primary tumor tissue") %>% mutate(tcellgroup = "Moderate")
LowOrModerate_TCellMD <- allMD %>% filter(Patient %in% LowOrModerate_Patients & (! Patient %in% ReC_Patients) & TissueType == "Primary tumor tissue") %>% mutate(tcellgroup = "LowOrModerate")
High_TCellMD <- allMD %>% filter(Patient %in% High_Patients & (! Patient %in% ReC_Patients) & TissueType == "Primary tumor tissue") %>% mutate(tcellgroup = "High")
totalTargetTcellMD <- bind_rows(Moderate_TCellMD, LowOrModerate_TCellMD, High_TCellMD)
matRoe <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 3),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = 3)
rownames(matRoe) <- sort(unique(allMD$g_seurat_clusters))
colnames(matRoe) <- c("Moderate", "LowOrModerate", "High")
matCell <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) * 3),
                  nrow = length(unique(allMD$g_seurat_clusters)),
                  ncol = 3)
rownames(matCell) <- sort(unique(allMD$g_seurat_clusters))
colnames(matCell) <- c("Moderate", "LowOrModerate", "High")
ggplotdata <- c()
for(tempGroup in c("Moderate", "LowOrModerate", "High")){
    targetTcellMD <- totalTargetTcellMD %>% filter(tcellgroup == tempGroup)
    M <- dim(targetTcellMD)[1]
    for(gsci in unique(allMD$g_seurat_clusters)){
        k <- ClusterT$n[ClusterT$g_seurat_clusters == gsci]
        if(gsci %in% targetTcellMD$g_seurat_clusters){
            n <- dim(targetTcellMD %>% filter(g_seurat_clusters == gsci))[1]
        } else{
            n <- 0
        }
        tempData <- tibble(Group = tempGroup,
                           Cluster = gsci,
                           Roe = (n/M) / (k/N))
        ggplotdata <- bind_rows(ggplotdata, tempData)
        tempRowName <- gsci
        tempColName <- tempGroup
        matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
        matCell[tempRowName, tempColName] <- n
    }
}
matRoe[matRoe > 2] <- 2
my.breaks <- seq(0, 2, by=0.01)
my.colors <- c(
    colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)/2),
    colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks)/2))
widthi = ncol(matRoe)/9 + 3
heighti = nrow(matRoe)/9 + 3
pdf(file.path(figurePath, paste0("colon_differentiate_roe_heatmap.pdf")), width = widthi, height = heighti)
print(pheatmap(matRoe,
               color = my.colors,
               breaks = my.breaks,
               ## annotation_col = cellType_col,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = T,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figurePath, paste0("colon_differentiate_cell_heatmap.pdf")), width = widthi * 1.5, height = heighti)
print(pheatmap(matCell,
               display_numbers = T,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = T,
               cluster_rows = F,
               number_format = "%.0f",
               border_color = F))
dev.off()








###############################################################################
#'                         merge all matrix together                         '#
###############################################################################


library(stringr)
library(randomcoloR)
obj.list = list()
fileNames = c(paste0("matRoe", 1:19, ".rds"))
colNum = 0
totalMatrix <- c()
for (fi in seq_along(fileNames)) {
    fileFullPath <- file.path("/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/figures4/matrix2", fileNames[fi])
    tempMatrix <- readRDS(fileFullPath)
    tempMatrix <- tempMatrix[rownames(tempMatrix) != "CD8_c14",]
    colNum = colNum + length(colnames(tempMatrix))
    totalMatrix <- cbind(totalMatrix, tempMatrix)
}
colnames(totalMatrix)
group_number <- str_extract(colnames(totalMatrix), "^[0-9]+")
cell_number <- str_extract(colnames(totalMatrix), "(?<=\\()[0-9]+(?=\\))")
remove_groups <- c(2, 4, 16)
remove_locations <- c(12, 13, 53)
remain_locations <- (!(group_number %in% remove_groups)) & (as.numeric(cell_number) >= 200) & (!(1:dim(totalMatrix)[2] %in% remove_locations))
totalMatrix <- totalMatrix[,remain_locations]
my.breaks <- seq(0, 1.99, by=0.01)
my.colors <- c(
    colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)/2),
    colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks)/2))
my.breaks.extra <- seq(2, 10, by = (10 - 2)/99)
my.colors.extra <- colorRampPalette(colors = c("#FD9AA0", "#550000"))(length(my.breaks.extra))
my.breaks <- c(my.breaks, my.breaks.extra)
my.colors <- c(my.colors, my.colors.extra)
widthi = ncol(totalMatrix)/7 + 1.8
heighti = nrow(totalMatrix)/7 + 2
pdf(file.path(figurePath, paste0("totalmatrix_roe_heatmap.pdf")), width = widthi, height = heighti)
print(pheatmap(totalMatrix,
               color = my.colors,
               breaks = my.breaks,
               ## annotation_row = cellType_col,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = T,
               cluster_rows = T,
               border_color = F))
dev.off()
pdf(file.path(figurePath, paste0("totalmatrix_roe_heatmap_CFRT.pdf")), width = widthi, height = heighti)
print(pheatmap(totalMatrix,
               color = my.colors,
               breaks = my.breaks,
               ## annotation_row = cellType_col,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = T,
               border_color = F))
dev.off()
group_number <- str_extract(colnames(totalMatrix), "^[0-9]+")
lung_group <- c(3, 5, 6, 7)
melan_group <- c(8:12)
lungMatrix <- totalMatrix[,group_number %in% lung_group]
melanMatrix <- totalMatrix[,group_number %in% melan_group]
elseMatrix <- totalMatrix[, (! group_number %in% lung_group) & (! group_number %in% melan_group)]
allMatrix <- list(total = totalMatrix, lung = lungMatrix, melan = melanMatrix, other = elseMatrix)
annotation_info <- c(
    "BCC",
    "",
    "Primary lung tumor",
    "Lung GSE127465",
    "Lung GSE123902",
    "LUAD smoker info",
    "LUAD",
    "",
    "Melanoma GSE127505",
    "Melanoma_GSE120575_treatment_roe_heatmap",
    "Melanoma_GSE115978_treatmentAndCancerType_roe",
    "SKCM",
    "Liver_GSE125449",
    "Colon RCCvsReC",
    "Breast_GSE114725_Her2ERPR_roe",
    "",
    "Lymph node ",
    "AML",
    "Head and neck"
)
annotation_info_CancerType <- c(
    "BCC",
    "",
    "NSCLC",
    "NSCLC",
    "NSCLC",
    "NSCLC",
    "NSCLC",
    "",
    "SKCM",
    "SKCM",
    "SKCM",
    "SKCM",
    "HCC",
    "CRCA",
    "BRCA",
    "",
    "Lymphoma",
    "AML",
    "HNC"
)

for(mi in 1:length(allMatrix)){
    tempMatrixName <- names(allMatrix[mi])
    tempMatrix <- allMatrix[[mi]]

    group_number <- as.numeric(str_extract(colnames(tempMatrix), "^[0-9]+"))
    group_number[group_number == 5] <- 7
    bad_groups <- names(table(group_number)[table(group_number) < 2])
    if(length(bad_groups) > 0){
        tempMatrix <- tempMatrix[, !group_number %in% bad_groups]
        group_number <- as.numeric(str_extract(colnames(tempMatrix), "^[0-9]+"))
        group_number[group_number == 5] <- 7
    }

    tempMatrix <- tempMatrix[, order(group_number, decreasing = F)]
    group_number <- group_number[order(group_number, decreasing = F)]
    tempMatrixColNames <- str_extract(colnames(tempMatrix), "(?<=_).*")
    colnames(tempMatrix) <- tempMatrixColNames
    group_info <- annotation_info[group_number]
    cancer_info <- annotation_info_CancerType[group_number]
    col_annotation <- data.frame(ClinicGroup = group_info, CancerType = cancer_info)
    rownames(col_annotation) <- tempMatrixColNames
    newOrder <- with(col_annotation, order(CancerType, ClinicGroup))
    group_number <- group_number[newOrder]
    col_annotation <- col_annotation[newOrder,]
    tempMatrix <- tempMatrix[,newOrder]
    if(mi == 1){
        newOrder <- c(
            2, 1,
            3, 4, 5, 6,
            10, 9, 8, 7,
            12,11,
            14, 13,
            17,16,15,
            20,19,18,
            22,23,21,
            26,27,24,25,
            28,29,
            30,31,
            32,33,
            35,34,
            36,42)
        tempMatrix <- tempMatrix[,newOrder]
        col_annotation <- col_annotation[newOrder,]
        group_number <- group_number[newOrder]
    }
    gaps <- cumsum(table(group_number)[as.character(unique(group_number))])
    clinic_group_gaps_col <- gaps[1:(length(gaps)-1)]
    colorsForDataType <- c("#6DCCDD", "#EDCAE0", "#F494BE", "#F9B26C", "#A6ADCC", "#C4DA5D")
    a <- length(unique(col_annotation$CancerType))
    b <- length(unique(col_annotation$ClinicGroup))
    ## mycolors_CancerType <- colorRampPalette(colorsForDataType)(length(unique(col_annotation$CancerType)))
    total_colors <- distinctColorPalette(a + b)
    mycolors_CancerType <- total_colors[1:a]
    names(mycolors_CancerType) <- unique(cancer_info)
    ## mycolors_ClinicGroup <- colorRampPalette(colorsForDataType)(length(unique(col_annotation$ClinicGroup)))
    mycolors_ClinicGroup <- total_colors[(a+1):(a+b)]
    names(mycolors_ClinicGroup) <- unique(group_info)
    mycolors <- list(ClinicGroup = mycolors_ClinicGroup, CancerType = mycolors_CancerType)
    my.breaks <- seq(0, 1.99, by=0.01)
    my.colors <- c(
        colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)/2),
        colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks)/2))
    my.breaks.extra <- seq(2, 10, by = (10 - 2)/99)
    my.colors.extra <- colorRampPalette(colors = c("#FD9AA0", "#550000"))(length(my.breaks.extra))
    my.breaks <- c(my.breaks, my.breaks.extra)
    my.colors <- c(my.colors, my.colors.extra)
    widthi = ncol(tempMatrix)/7 + 8
    heighti = nrow(tempMatrix)/7 + 2
    pdf(file.path(figurePath, paste0(tempMatrixName, "_roe_heatmap_CFRT.pdf")), width = widthi, height = heighti)
    print(pheatmap(tempMatrix,
                   color = my.colors,
                   breaks = my.breaks,
                   annotation_col = col_annotation,
                   annotation_names_col = T,
                   annotation_colors = mycolors,
                   gaps_col = clinic_group_gaps_col,
                   show_colnames = T,
                   show_rownames = T,
                   cluster_cols = F,
                   cluster_rows = T,
                   border_color = F))
    dev.off()
}
