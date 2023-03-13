#--------------------------------------------------------------
# filename : SampleFractionCorrelation.r
# Date : 2022-10-29
# contributor : Yanshuo Chu
# function: SampleFractionCorrelation
# R version: R/4.0.3
#--------------------------------------------------------------

print('<==== SampleFractionCorrelation.r ====>')
rm(list=ls())

library(tidyverse)
library(Seurat )
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
library(gridExtra)
library(grid)
set.seed(123)

CD8_path <- '/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/0_write_sample_info/CD8_2022-10-20.rds'
CD4_path <- '/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/0_write_sample_info/CD4_2022-10-20.rds'
Innate_path <- '/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/0_write_sample_info/Innate_2022-10-20.rds'

## CD8_path <- '/rsrch3/scratch/genomic_med/ychu2/data/Tcellproject/data/ClinicInfoBk/CD8_cluster_2021-07-22.rds'
## CD4_path <- '/rsrch3/scratch/genomic_med/ychu2/data/Tcellproject/data/ClinicInfoBk/CD4_cluster_2021-07-22.rds'
## Innate_path <- '/rsrch3/scratch/genomic_med/ychu2/data/Tcellproject/data/ClinicInfoBk/NKT_cluster_2021-07-22.rds'

CD8_Obj <- readRDS(CD8_path)
CD4_Obj <- readRDS(CD4_path)
Innate_Obj <- readRDS(Innate_path)

allPath <- list(CD8 = CD8_path, CD4 = CD4_path, NK = Innate_path)
allObj <- list(CD8 = CD8_Obj, CD4 = CD4_Obj, NK = Innate_Obj)


common_columns <- c()
for(oi in 1:length(allObj)) {
    tempName <- names(allObj[oi])
    tempObj <- allObj[[oi]]
    tempObj$barcode <- Cells(tempObj)
    if (oi == 1) {
        common_columns <- colnames(tempObj@meta.data)
    } else {
        common_columns <- intersect(common_columns, colnames(tempObj@meta.data))
    }
}

allMD <- c()
for (oi in 1:length(allObj)) {
    tempName <- names(allObj[oi])
    tempObj <- allObj[[oi]]
    tempMD <- as_tibble(tempObj@meta.data, rownames = NA) %>%
        dplyr::select(all_of(common_columns)) %>% mutate(barcode = Cells(tempObj))
    allMD <- bind_rows(allMD, tempMD)
}
colorsForDataType <- c("#6DCCDD", "#EDCAE0", "#F494BE", "#F9B26C", "#A6ADCC", "#C4DA5D")
colorForClass1 <- c("#C4DA5D", "#6DCCDD", "#F494BE", "#EDCAE0")
dataTypeLevel <- c("CD4", "CD8", "MAIT", "Tgd", "NKT", "Proliferative")

figure_path <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/1_main_fig4/")
if (!dir.exists(figure_path)) {
  dir.create(figure_path, recursive = T)
}
setwd(figure_path)


#' sample fraction spearman correlation ######################################
sampleN <- allMD %>% group_by(Sample) %>% count() %>% filter(n > 200)
allMD_filtered <- allMD %>% filter(Sample %in% sampleN$Sample)
allMD_filtered$CancerType[allMD_filtered$TissueType %in% c("Healthy donor", "Uninvolved normal tissue")] <-
    allMD_filtered$OrganSite[allMD_filtered$TissueType %in% c("Healthy donor", "Uninvolved normal tissue")]

#' Step1: initial sample fraction matrix #####################################

## Healthy donor
## Metastatic tumor tissue
## Primary tumor tissue
## Uninvolved normal tissue

tumor_type_groups <- list(
  Primary = "Primary tumor tissue",
  Metastatic = "Metastatic tumor tissue",
  Healthy = "Healthy donor",
  Normal = "Uninvolved normal tissue",
  NormalAndTumor = c(
    "Uninvolved normal tissue",
    "Primary tumor tissue",
    "Metastatic tumor tissue"),
  NormalAndPrimary = c(
    "Uninvolved normal tissue",
    "Primary tumor tissue"),
  NormalAndHealthy = c("Uninvolved normal tissue",
                       "Healthy donor"))


for (tyi in seq_along(tumor_type_groups)) {
  tTy_name <- names(tumor_type_groups[tyi])
  tTy <- tumor_type_groups[[tyi]]

  tMD <- allMD_filtered %>%
    filter(TissueType %in% tTy)

  target_cell_types <- tMD %>%
    group_by(cell.type, Sample) %>%
    count()
  target_cell_types$frac <- target_cell_types$n / sampleN$n[match(target_cell_types$Sample, sampleN$Sample)]
  target_cell_types <- target_cell_types %>%
    ungroup() %>% 
    group_by(cell.type) %>%
    summarise(medianFrac = mean(frac)) %>%
    filter(medianFrac > 0.02) %>%
    pull(cell.type)

  tMD <- tMD %>%
    filter(cell.type %in% target_cell_types)

  matFrac <- matrix(
    rep(0, length(target_cell_types) * length(unique(tMD$Sample))),
    nrow = length(target_cell_types),
    ncol = length(unique(tMD$Sample)))

  rownames(matFrac) <- sort(unique(tMD$cell.type))
  colnames(matFrac) <- sort(unique(tMD$Sample))
  for(tC in unique(target_cell_types)){
    for(tS in unique(tMD$Sample)){
      matFrac[tC, tS] <-
        dim(tMD %>% filter(cell.type == tC) %>% filter(Sample == tS))[1] / sampleN$n[sampleN$Sample == tS]
    }
  }

  mc <- cor(t(matFrac), method = "spearman")
  mc_p_adj <- rstatix::cor_mat(as.tibble(t(matFrac)), method = "spearman") %>%
    rstatix::cor_get_pval() %>%
    rstatix::adjust_pvalue(method = "fdr")
  mc_p_adj <- data.frame(mc_p_adj[,2:dim(mc_p_adj)[2]])
  rownames(mc_p_adj) <- colnames(mc_p_adj)
  for(mcj in 1:dim(mc_p_adj)[2]){
    mc_p_adj[, mcj] <-  gtools::stars.pval(mc_p_adj[, mcj])
  }
  rownames(mc_p_adj) <- rownames(mc)
  colnames(mc_p_adj) <- colnames(mc)

  mmc = mc
  mmc[mmc == 1] = -1
  mc[mc == 1] = max(mmc)
  p_cor <- pheatmap(mc,
                    show_colnames = T,
                    show_rownames = T,
                    cluster_cols = T,
                    cluster_rows = T,
                    border_color = F)
  row_order <- p_cor$tree_row$order
  row_order_name <- rownames(mc)[row_order]
  mc[mmc == -1] <- NA
  mc <- mc[row_order_name, row_order_name]
  pdf(file.path(paste0(tTy_name, "_cor.pdf")))
  pheatmap(mc,
           show_colnames = T,
           show_rownames = T,
           cluster_cols = F,
           cluster_rows = F,
           border_color = F)
  dev.off()


  mc_tibble <- mc %>%
    as.data.frame() %>%
    rownames_to_column(var = "x") %>%
    pivot_longer(cols = -x, names_to = "y") %>%
    mutate(key=paste0(x, "__", y))
  mc_padj_tibble <- mc_p_adj %>%
    rownames_to_column(var = "x") %>%
    pivot_longer(cols = -x, names_to = "y") %>%
    mutate(key=paste0(x, "__", y))
  mc_tibble$star <- mc_padj_tibble$value[match(mc_tibble$key, mc_padj_tibble$key)]
  mc_tibble$star[is.na(mc_tibble$value)] <- ""
  mc_tibble$x <- factor(mc_tibble$x, levels = row_order_name)
  mc_tibble$y <- factor(mc_tibble$y, levels = row_order_name)
  g <- mc_tibble %>%
    ggplot() +
    geom_tile(mapping = aes(x = x, y = y, fill = value)) +
    geom_text(mapping = aes(x = x, y = y, label = star),
              angle = 45,
              vjust = 0.8,
              hjust = 0.5) +
    scale_fill_gradient2(low = "#336699", high = "#FF3333", mid = "white",
                         na.value = "grey",  midpoint = 0) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = -60, vjust = 0.5, hjust=0),
          strip.text.y = element_text(angle = 0),
          strip.background = element_rect(colour=NA, fill=NA)) 
  ggsave(file.path(paste0(tTy_name, "_tile.pdf")), width = 8)
}
