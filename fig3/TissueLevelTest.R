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

figurePath <- "/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/scripts/pipelines/step4_figs/FigsForPaper/fig3/out"
if(!dir.exists(figurePath)){
    dir.create(figurePath)
}
colorsForDataType <- c("#6DCCDD", "#EDCAE0", "#F494BE", "#F9B26C", "#A6ADCC", "#C4DA5D")
colorForClass1 <- c("#C4DA5D", "#6DCCDD", "#F494BE", "#EDCAE0")
dataTypeLevel <- c("CD4", "CD8", "NK", "NKT", "MAIT", "Proliferative")

Treg_path <- '/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/validate/Sub_Treg_CD4_V5/nPC_15/UMAP_dist_0.1_nneighbor_20/p1_sub_Treg_CD4_V5_UMAP_dist_0.1_nneighbor_20_CLUSTER_res_0.3/cluster.rds'
Treg_Obj <- readRDS(Treg_path)
TFH_path <- '/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/validate/Sub_TFH_CD4_V6/nPC_10/UMAP_dist_0.01_nneighbor_50/p1_sub_TFH_CD4_V6_UMAP_dist_0.01_nneighbor_50_CLUSTER_res_0.3/cluster.rds'
TFH_Obj <- readRDS(TFH_path)
Treg_Obj <- subset(Treg_Obj, cells = Cells(Treg_Obj)[Cells(Treg_Obj) %in% Cells(CD4_Obj)])
TFH_Obj <- subset(TFH_Obj, cells = Cells(TFH_Obj)[Cells(TFH_Obj) %in% Cells(CD4_Obj)])
for(tempColumn in c("CancerType", "TissueType", "OrganSite", "Disease", "Patient", "Sample")){
    Treg_Obj@meta.data[,tempColumn] <- CD4_Obj@meta.data[match(Cells(Treg_Obj), Cells(CD4_Obj)),tempColumn]
    TFH_Obj@meta.data[,tempColumn] <- CD4_Obj@meta.data[match(Cells(TFH_Obj), Cells(CD4_Obj)),tempColumn]
}


#' filter sample #############################################################
sampleN <- allMD %>% group_by(OrganSite, Sample) %>% count() %>% filter(n > 200)
TypeSampleN <- allMD
TypeSampleN$CancerType[TypeSampleN$TissueType %in% c("Healthy donor", "Uninvolved normal tissue")] = TypeSampleN$OrganSite[TypeSampleN$TissueType %in% c("Healthy donor", "Uninvolved normal tissue")]
TypeSampleN <- TypeSampleN %>% group_by(TissueType, CancerType, Sample) %>% count() %>% filter(Sample %in% sampleN$Sample)
TissueType.labs <- c("H", "U", "P", "M")
names(TissueType.labs) <- c("Healthy donor", "Uninvolved normal tissue", "Primary tumor tissue", "Metastatic tumor tissue")

allObj <- list(CD8 = CD8_Obj,
               CD4 = CD4_Obj,
               NK = NK_Obj,
               Proliferative = Proliferative_Obj,
               Treg = Treg_Obj,
               TFH = TFH_Obj)

for(oi in 1:length(allObj)){
    tempName <- names(allObj[oi])
    tempObj <- allObj[[oi]]
    Idents(tempObj) <- tempObj@meta.data$seurat_clusters
    ggplotdata <- as_tibble(tempObj@meta.data, rownames = NA) %>%
        select(seurat_clusters, Sample, CancerType, TissueType, OrganSite, orig.ident)
    ggplotdata$CancerType[ggplotdata$TissueType %in% c("Healthy donor", "Uninvolved normal tissue")] = ggplotdata$OrganSite[ggplotdata$TissueType %in% c("Healthy donor", "Uninvolved normal tissue")]
    ggplotdata <- ggplotdata %>%
        group_by(seurat_clusters, TissueType, CancerType, Sample) %>%
        dplyr::count() %>%
        filter(Sample %in% sampleN$Sample)
    ggplotdata$frac <- ggplotdata$n / sampleN$n[match(ggplotdata$Sample, sampleN$Sample)]

    totalObj_ggplotdata <- c()
    for(tempSeuratCluster in unique(ggplotdata$seurat_clusters)){
        temp_ggplotdata <- ggplotdata %>%
            filter(seurat_clusters == tempSeuratCluster)
        for(temp_TissueType in unique(TypeSampleN$TissueType)){
            temp_TypeSampleN <- TypeSampleN %>% filter(TissueType == temp_TissueType)
            temp_diff_Samples <- setdiff(unique( temp_TypeSampleN %>% pull(Sample)),
                                         unique(temp_ggplotdata %>% filter(TissueType == temp_TissueType) %>% pull(Sample)))
            if(length(temp_diff_Samples) > 0 ){
                t_rows <- tibble(seurat_clusters = rep(tempSeuratCluster, length(temp_diff_Samples)),
                                 TissueType = rep(temp_TissueType, length(temp_diff_Samples)),
                                 CancerType = temp_TypeSampleN$CancerType[temp_TypeSampleN$Sample %in% temp_diff_Samples],
                                 Sample = temp_diff_Samples,
                                 n = rep(0, length(temp_diff_Samples)),
                                 frac = rep(0, length(temp_diff_Samples)))
                temp_ggplotdata <- bind_rows(temp_ggplotdata, t_rows)
            }
        }
        temp_ggplotdata$TissueType <- factor(temp_ggplotdata$TissueType, levels = c("Healthy donor", "Uninvolved normal tissue", "Primary tumor tissue", "Metastatic tumor tissue"))
        if(tempName == "NK"){
            if(tempSeuratCluster %in% c(2)){
                tempDataSet = "MAIT"
            }else if(tempSeuratCluster %in% c(1,3)){
                tempDataSet = "MAIT-like"
            }else if(tempSeuratCluster %in% c(0, 4)){
                tempDataSet = "NKT"
            }else{
                tempDataSet = "NK"
            }
        }else{
            tempDataSet = tempName
        }
        temp_ggplotdata$g_seurat_clusters <- paste0(tempDataSet, "_c", tempSeuratCluster)
        totalObj_ggplotdata <- bind_rows(totalObj_ggplotdata, temp_ggplotdata)
    }


    ggplotdata <- totalObj_ggplotdata
    ggplotdata$TissueType <- plyr::mapvalues(ggplotdata$TissueType,
                                             c("Healthy donor",
                                               "Uninvolved normal tissue",
                                               "Primary tumor tissue",
                                               "Metastatic tumor tissue"),
                                             c("H", "U", "P", "M"))
    ggplotdata$TissueType <- factor(ggplotdata$TissueType, levels = c("H", "U", "P", "M"))

    g <- ggplotdata %>%
        ggstatsplot::grouped_ggbetweenstats(
                         data = .,
                         x = TissueType,
                         y = frac,
                         grouping.var = g_seurat_clusters,
                         xlab = "Tissue type",
                         ylab = "Sample fraction",
                         pairwise.display = "all", # display only significant pairwise comparisons
                         p.adjust.method = "fdr", # adjust p-values for multiple tests using this method
                         ggtheme = theme_classic(),
                         package = "ggsci",
                         palette = "default_jco",
                         plotgrid.args = list(ncol = 3))
    ggsave(file.path(figurePath, paste0(tempName, "_cluster_all-sample_boxplot_gbs.pdf")),
           g, width = 400, height = 600, units = "mm")

    if(tempName == "NK"){
        ## temp_ggplotdata <- ggplotdata %>%
        temp_ggplotdata <- ggplotdata
        ##     filter(g_seurat_clusters %in% c("NKT_c4", "MAIT-like_c1", "NKT_c0"))
        temp_ggplotdata$g_seurat_clusters <- factor(temp_ggplotdata$g_seurat_clusters,
                                                    levels = c("NKT_c0", "NKT_c4", "MAIT-like_c1", "MAIT_c2", "MAIT-like_c3"))
        g <- temp_ggplotdata %>% ggplot() +
            geom_jitter(aes(x = TissueType,
                            y = frac,
                            fill = TissueType,
                            color = TissueType), size = 0.2) +
            geom_boxplot(aes(x = TissueType,
                             y = frac),
                         fill = NA,
                         outlier.shape = NA,
                         lwd = 0.2) +
            facet_wrap(g_seurat_clusters ~ ., scales ="free", nrow =1) +
            theme_classic() +
            theme(text = element_text(size = 10),
                  legend.position = "none",
                  strip.background = element_rect(colour="white", fill="white"))
        ggsave(file.path(figurePath,
                         paste0(tempName, "_cluster_selected-sample_boxplot.pdf")),
               g, width = 210, height = 297/8, units = "mm")
    }

    if(tempName == "Proliferative"){
        temp_ggplotdata <- ggplotdata %>%
            filter(g_seurat_clusters %in% paste0("Proliferative_c", c(0, 6)))
        temp_ggplotdata$g_seurat_clusters <- plyr::mapvalues(temp_ggplotdata$g_seurat_clusters,
                                                             paste0("Proliferative_c", c(0, 6)),
                                                             paste0("P_c", c(0, 6)))
        temp_ggplotdata$g_seurat_clusters <- factor(temp_ggplotdata$g_seurat_clusters,
                                                    levels = paste0("P_c", c(0, 6)))
        g <- temp_ggplotdata %>% ggplot() +
            geom_jitter(aes(x = TissueType,
                            y = frac,
                            fill = TissueType,
                            color = TissueType), size = 0.2) +
            geom_boxplot(aes(x = TissueType,
                             y = frac),
                         fill = NA,
                         outlier.shape = NA,
                         lwd = 0.2) +
            facet_wrap(g_seurat_clusters ~ ., scales ="free", nrow =1) +
            theme_classic() +
            theme(text = element_text(size = 10),
                  legend.position = "none",
                  strip.background = element_rect(colour="white", fill="white"))
        ggsave(file.path(figurePath,
                         paste0(tempName, "_cluster_selected-sample_boxplot.pdf")),
               g, width = 210 * 2 / 5, height = 297/8, units = "mm")
    }

    if(tempName == "Treg"){
        temp_ggplotdata <- ggplotdata %>%
            filter(g_seurat_clusters %in% paste0("Treg_c", 0:6))
        temp_ggplotdata$g_seurat_clusters <- factor(temp_ggplotdata$g_seurat_clusters,
                                                    levels = paste0("Treg_c", 0:6))
        g <- temp_ggplotdata %>% ggplot() +
            geom_jitter(aes(x = TissueType,
                            y = frac,
                            fill = TissueType,
                            color = TissueType), size = 0.2) +
            geom_boxplot(aes(x = TissueType,
                             y = frac),
                         fill = NA,
                         outlier.shape = NA,
                         lwd = 0.2) +
            facet_wrap(g_seurat_clusters ~ ., scales ="free", nrow =1) +
            theme_classic() +
            theme(text = element_text(size = 10),
                  legend.position = "none",
                  strip.background = element_rect(colour="white", fill="white"))
        ggsave(file.path(figurePath,
                         paste0(tempName, "_cluster_selected-sample_boxplot.pdf")),
               g, width = 210, height = 297/6, units = "mm")
    }
}
