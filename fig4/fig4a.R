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

figurePath <- "/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/scripts/pipelines/step4_figs/FigsForPaper/fig4/output"
if(!dir.exists(figurePath)){
    dir.create(figurePath)
}



allMD$g_seurat_clusters <- factor(allMD$g_seurat_clusters,
                                  levels = c(paste0("CD8_c", 0:13),
                                             paste0("CD4_c", 0:11),
                                             "NKT_c0",
                                             "MAIT-like_c1",
                                             "MAIT_c2",
                                             "MAIT-like_c3",
                                             "NKT_c4",
                                             paste0("Proliferative_c", 0:7))
                                  )
DataSetL <- list(CD8 = "CD8", CD4 = "CD4", Innate = c("MAIT", "NKT", "MAIT-like"), Proliferative = "Proliferative")
allMD$DataSetCategory <- allMD$DataSet
allMD$DataSetCategory[allMD$DataSetCategory %in% c("MAIT", "NKT", "MAIT-like")] <- "Innate"
for(dsLi in 1:length(DataSetL)){
    tempDataSetName <- names(DataSetL[dsLi])
    tempDataSetIDs <- DataSetL[[dsLi]]
    ggplotdata <- allMD %>% filter(DataSet %in% tempDataSetIDs)
    clusterNum <- ggplotdata %>% group_by(g_seurat_clusters) %>% dplyr::count()
    gL <- list()
    for(tT in c("Healthy donor", "Uninvolved normal tissue", "Primary tumor tissue", "Metastatic tumor tissue")){
        tgd <- ggplotdata %>% filter(TissueType == tT) %>%
            group_by(g_seurat_clusters, CancerType, TissueType) %>%
            dplyr::count()
        tgd$frac <- 0
        tgd$frac <- tgd$n / clusterNum$n[match(tgd$g_seurat_clusters, clusterNum$g_seurat_clusters)]
        g <- tgd %>% ggplot() +
            geom_bar(aes(x = g_seurat_clusters, y = frac,
                         fill = CancerType), stat = "identity") +
            facet_grid(. ~ TissueType, scales = "free", space = "free") +
            theme_classic() +
            theme(
                text = element_text(size = 10),
                axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
                legend.direction = "vertical",
                strip.background = element_blank(),
                legend.key.size = unit(2.5, "mm")) +
            guides(fill = guide_legend(ncol = 1, override.aes = list(size = 0.05)))
        gL[[length(gL)+1]] <- g
    }
    g <- plot_grid(plotlist = gL, ncol = 2, nrow= 2, axis = "bt", align = 'vh')
    ggsave(file.path(figurePath, paste0(tempDataSetName, "_cluster-cancertype_fraction.pdf")), g,  width = 210, height = 297/2, units = "mm")
}
DataSetL <- list(CD8 = "CD8", CD4 = "CD4", Innate = c("MAIT", "NKT", "MAIT-like"), Proliferative = "Proliferative")
allMD$DataSetCategory <- allMD$DataSet
allMD$DataSetCategory[allMD$DataSetCategory %in% c("MAIT", "NKT", "MAIT-like")] <- "Innate"
clusterNum <- allMD %>% group_by(g_seurat_clusters) %>% dplyr::count()
ggdat <- allMD %>% group_by(DataSetCategory, TissueType, g_seurat_clusters) %>%
    dplyr::count() %>% group_by(g_seurat_clusters) %>%
    mutate(frac = n / sum(n))
ggdat <- ggdat %>% filter(DataSetCategory != "Proliferative")
ggdat$g_seurat_clusters <- factor(ggdat$g_seurat_clusters, levels = c("MAIT-like_c3", "CD8_c12", "CD8_c13", "CD4_c10", "CD4_c9", "CD4_c7", "CD8_c3", "CD4_c6", "CD4_c2", "NKT_c4", "CD8_c2", "NKT_c0", "CD8_c1", "CD8_c7", "CD8_c10", "CD8_c8", "CD8_c6", "CD8_c9", "CD8_c0", "CD8_c11", "CD8_c4", "CD4_c4", "MAIT-like_c1", "MAIT_c2", "CD4_c5", "CD4_c0", "CD8_c5", "CD4_c11", "CD4_c8", "CD4_c3", "CD4_c1"))
ggdat$TissueType <- factor(ggdat$TissueType, levels = c("Healthy donor", "Uninvolved normal tissue", "Primary tumor tissue", "Metastatic tumor tissue"))
g <- ggdat %>% ggplot() +
    geom_bar(aes(x = g_seurat_clusters, y = frac,
                 fill = TissueType), stat = "identity") +
    ## facet_grid(. ~ DataSetCategory, scales = "free", space = "free") +
    theme_classic() +
    theme(
        text = element_text(size = 10),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        legend.direction = "vertical",
        strip.background = element_blank(),
        legend.key.size = unit(2.5, "mm")) +
    guides(fill = guide_legend(ncol = 1, override.aes = list(size = 0.05)))
ggsave(file.path(figurePath, paste0("cluster-TissueType_fraction.pdf")), g,  width = 210, height = 297/4, units = "mm")



#' add module score ##########################################################

CD8_Obj <- readRDS(CD8_path)
Idents(CD8_Obj) <- CD8_Obj$seurat_clusters
CD8FunctionList <- "/rsrch3/home/genomic_med/ychu2/configs/scSeqs/database/Markers/CD8/forFig2"
fileFullPathL = list.files(path = CD8FunctionList, recursive = F)
fileNames = basename(fileFullPathL)
if(length(fileNames) < 1) stop(paste0('Failed to find data in ', CD8FunctionList))
marker.list = list()
for (fi in seq_along(fileNames)) {
    fileFullPath <- file.path(CD8FunctionList, fileNames[fi])
    FunctionName <- tools::file_path_sans_ext(basename(fileFullPath))
    markert <- read_tsv(fileFullPath)
    colnames(markert)[1] <- c("Gene")
    tempGenes <- markert %>% pull(Gene)
    marker.list[[FunctionName]] <- tempGenes
}
marker.list[["DKK3"]] <- c("DKK3")
CD8_Obj <-  AddModuleScore(CD8_Obj,
                           features = marker.list,
                           ctrl = 5,
                           name = "FunctionScore")
for(i in 1:length(marker.list)){
    colnames(CD8_Obj@meta.data)[colnames(CD8_Obj@meta.data) == paste0("FunctionScore", i)] <- names(marker.list)[i]
}
marker.list <- marker.list[names(marker.list) != "DKK3"]
Idents(CD8_Obj) <- CD8_Obj$seurat_clusters
Differentiation <- c("Naive", "Activation:Effector function", "Exhaustion")
Function <- c("TCR Signaling", "Cytotoxicity", "Cytokine:Cytokine receptor",
              "Chemokine:Chemokine receptor", "Senescence", "Anergy",
              "NFKB Signaling", "Stress response", "MAPK Signaling", "Adhesion",
              "IFN Response")
Metabolism <- c("Oxidative phosphorylation", "Glycolysis", "Fatty acid metabolism")
Apoptosis <- c("Pro-apoptosis", "Anti-apoptosis")
MarkerNameVector <- c(Differentiation, Function, Metabolism, Apoptosis)
FunctionScoreMatrix <- matrix(0,
                              ncol = length(unique(CD8_Obj$seurat_clusters)),
                              nrow = length(marker.list))
colnames(FunctionScoreMatrix) <- paste0("CD8_c", 0:13)
rownames(FunctionScoreMatrix) <- MarkerNameVector
for(ci in 1:ncol(FunctionScoreMatrix)){
    for(ri in 1:nrow(FunctionScoreMatrix)){
        FunctionVec <- as_tibble(CD8_Obj@meta.data) %>% pull(MarkerNameVector[ri])
        fv <- mean(FunctionVec[CD8_Obj$seurat_clusters == levels(CD8_Obj$seurat_clusters)[ci]])
        FunctionScoreMatrix[ri, ci] <- fv
    }
}
FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, rescale, to=c(-1, 1)))
targetFunctionName <- c("Naive", "Activation:Effector function", "Cytotoxicity", "Exhaustion", "Stress response", "IFN Response")
CD8FuncScoreMatrix <- FunctionScoreMatrix[targetFunctionName,]
CD4_Obj <- readRDS(CD4_path)
CD4_Obj$seurat_clusters <- droplevels(CD4_Obj$seurat_clusters)
CD4FunctionList <- "/rsrch3/home/genomic_med/ychu2/configs/scSeqs/database/Markers/CD4/forFig3"
fileFullPathL = list.files(path = CD4FunctionList, recursive = F)
fileNames = basename(fileFullPathL)
if(length(fileNames) < 1) stop(paste0('Failed to find data in ', CD4FunctionList))
marker.list = list()
for (fi in seq_along(fileNames)) {
    fileFullPath <- file.path(CD4FunctionList, fileNames[fi])
    FunctionName <- tools::file_path_sans_ext(basename(fileFullPath))
    markert <- read_tsv(fileFullPath)
    colnames(markert)[1] <- c("Gene")
    tempGenes <- markert %>% pull(Gene)
    marker.list[[FunctionName]] <- tempGenes
}
CD4_Obj <-  AddModuleScore(CD4_Obj,
                           features = marker.list,
                           ctrl = 5,
                           name = "FunctionScore")
for(i in 1:length(marker.list)){
    colnames(CD4_Obj@meta.data)[colnames(CD4_Obj@meta.data) == paste0("FunctionScore", i)] <- names(marker.list)[i]
}
Idents(CD4_Obj) <- CD4_Obj$seurat_clusters
Differentiation <- c("Naive", "Activation:Effector function", "Exhaustion")
Function <- c("TCR signaling", "Cytotoxicity", "Cytokine:Cytokine receptor",
              "Chemokine:Chemokine receptor", "Stress response", "Adhesion",
              "IFN response", "Treg signature", "Costimulatory molecules", "Transcription factor")
Metabolism <- c("OXPHOS", "Glycolysis", "Lipid metabolism")
Apoptosis <- c("Pro-apoptosis", "Anti-apoptosis")
MarkerNameVector <- c(Differentiation, Function, Metabolism, Apoptosis)
FunctionScoreMatrix <- matrix(0,
                              ncol = length(unique(CD4_Obj$seurat_clusters)),
                              nrow = length(marker.list))
colnames(FunctionScoreMatrix) <- paste0("CD4_c", levels(CD4_Obj$seurat_clusters))
rownames(FunctionScoreMatrix) <- MarkerNameVector
for(ci in 1:ncol(FunctionScoreMatrix)){
    for(ri in 1:nrow(FunctionScoreMatrix)){
        FunctionVec <- as_tibble(CD4_Obj@meta.data) %>% pull(MarkerNameVector[ri])
        fv <- mean(FunctionVec[CD4_Obj$seurat_clusters == levels(CD4_Obj$seurat_clusters)[ci]])
        FunctionScoreMatrix[ri, ci] <- fv
    }
}
FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, rescale, to=c(-1, 1)))
targetFunctionName <- c("Naive", "Activation:Effector function", "Cytotoxicity", "Exhaustion", "Stress response", "IFN response")
CD4FuncScoreMatrix <- FunctionScoreMatrix[targetFunctionName,]
NK_Obj <- readRDS(NK_path)
NKFunctionList <- "/rsrch3/home/genomic_med/ychu2/configs/scSeqs/database/Markers/CD4/forFig3"
fileFullPathL = list.files(path = NKFunctionList, recursive = F)
fileNames = basename(fileFullPathL)
if(length(fileNames) < 1) stop(paste0('Failed to find data in ', NKFunctionList))
marker.list = list()
for (fi in seq_along(fileNames)) {
    fileFullPath <- file.path(NKFunctionList, fileNames[fi])
    FunctionName <- tools::file_path_sans_ext(basename(fileFullPath))
    markert <- read_tsv(fileFullPath)
    colnames(markert)[1] <- c("Gene")
    tempGenes <- markert %>% pull(Gene)
    marker.list[[FunctionName]] <- tempGenes
}
NK_Obj <-  AddModuleScore(NK_Obj,
                           features = marker.list,
                           ctrl = 5,
                           name = "FunctionScore")
for(i in 1:length(marker.list)){
    colnames(NK_Obj@meta.data)[colnames(NK_Obj@meta.data) == paste0("FunctionScore", i)] <- names(marker.list)[i]
}
Idents(NK_Obj) <- NK_Obj$g_seurat_clusters
Differentiation <- c("Naive", "Activation:Effector function", "Exhaustion")
Function <- c("TCR signaling", "Cytotoxicity", "Cytokine:Cytokine receptor",
              "Chemokine:Chemokine receptor", "Stress response", "Adhesion",
              "IFN response", "Treg signature", "Costimulatory molecules", "Transcription factor")
Metabolism <- c("OXPHOS", "Glycolysis", "Lipid metabolism")
Apoptosis <- c("Pro-apoptosis", "Anti-apoptosis")
MarkerNameVector <- c(Differentiation, Function, Metabolism, Apoptosis)
FunctionScoreMatrix <- matrix(0,
                              ncol = length(unique(NK_Obj$g_seurat_clusters)),
                              nrow = length(marker.list))
colnames(FunctionScoreMatrix) <- unique(NK_Obj$g_seurat_clusters)
rownames(FunctionScoreMatrix) <- MarkerNameVector
for(ci in 1:ncol(FunctionScoreMatrix)){
    for(ri in 1:nrow(FunctionScoreMatrix)){
        FunctionVec <- as_tibble(NK_Obj@meta.data) %>% pull(MarkerNameVector[ri])
        fv <- mean(FunctionVec[NK_Obj$g_seurat_clusters == unique(NK_Obj$g_seurat_clusters)[ci]])
        FunctionScoreMatrix[ri, ci] <- fv
    }
}
FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, rescale, to=c(-1, 1)))
targetFunctionName <- c("Naive", "Activation:Effector function", "Cytotoxicity", "Exhaustion", "Stress response", "IFN response")
NKFuncScoreMatrix <- FunctionScoreMatrix[targetFunctionName,]
rownames(CD8FuncScoreMatrix) <- targetFunctionName
rownames(CD4FuncScoreMatrix) <- targetFunctionName
rownames(NKFuncScoreMatrix) <- targetFunctionName
FuncScoreMatrix <- cbind(CD8FuncScoreMatrix, CD4FuncScoreMatrix)
FuncScoreMatrix <- cbind(FuncScoreMatrix, NKFuncScoreMatrix)
ggplotdata <- reshape2::melt(FuncScoreMatrix)
ggplotdata$Var2 <- factor(ggplotdata$Var2, levels = c("MAIT-like_c3", "CD8_c12", "CD8_c13", "CD4_c10", "CD4_c9", "CD4_c7", "CD8_c3", "CD4_c6", "CD4_c2", "NKT_c4", "CD8_c2", "NKT_c0", "CD8_c1", "CD8_c7", "CD8_c10", "CD8_c8", "CD8_c6", "CD8_c9", "CD8_c0", "CD8_c11", "CD8_c4", "CD4_c4", "MAIT-like_c1", "MAIT_c2", "CD4_c5", "CD4_c0", "CD8_c5", "CD4_c11", "CD4_c8", "CD4_c3", "CD4_c1"))
ggplotdata$Var1 <- factor(ggplotdata$Var1, levels = rev(targetFunctionName))
my.breaks <- c(seq(-1, 0, by=0.1), seq(0.1, 1, by=0.1))
my.colors <- c(
    colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)/2),
    colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks)/2))
g <- ggplotdata %>%
    ggplot() +
    geom_tile(aes(x = Var2, y = Var1, fill = value)) +
    scale_fill_gradient2(high = "red", low="blue", mid = "white") +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    )
ggsave(file.path(figurePath, paste0("temp_function.pdf")), g,  width = 210, height = 297/4, units = "mm")
