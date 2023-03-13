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

CD8_path <- '/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/0_write_sample_info/CD8_2022-10-20.rds'
CD4_path <- '/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/0_write_sample_info/CD4_2022-10-20.rds'
Innate_path <- '/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/0_write_sample_info/Innate_2022-10-20.rds'
Proliferative_path <- '/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/0_write_sample_info/Proliferative_2022-10-20.rds'

CD8_Obj <- readRDS(CD8_path)
CD4_Obj <- readRDS(CD4_path)
Innate_Obj <- readRDS(Innate_path)
Proliferative_Obj <- readRDS(Proliferative_path)

allPath <- list(CD8 = CD8_path,
                CD4 = CD4_path,
                Innate = Innate_path,
                Proliferative = Proliferative_path)
allObj <- list(CD8 = CD8_Obj,
               CD4 = CD4_Obj,
               Innate = Innate_Obj,
               Proliferative = Proliferative_Obj)
common_columns <- c()
for (oi in 1:length(allObj)) {
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
dataTypeLevel <- c("CD4", "CD8", "Tgd", "NKT", "MAIT", "Proliferative")

figure_path <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/1_main_fig4/")
if (!dir.exists(figure_path)) {
  dir.create(figure_path, recursive = T)
}
setwd(figure_path)


tAllObj <- list(CD8 = CD8_Obj,
                CD4 = CD4_Obj,
                Innate = Innate_Obj)
matrixRoeL <- list()
matrixCellL <- list()
allMD_type <- unique(paste0(allMD$TissueType, "_", allMD$CancerType))
totalMD <- c()
for (oi in 1:length(tAllObj)) {
  tempName <- names(tAllObj[oi])
  tempObj <- tAllObj[[oi]]
  md <- as_tibble(tempObj@meta.data, rownames = NA)
  md <- md %>% filter(TissueType %in% c("Primary tumor tissue",
                                        "Metastatic tumor tissue"))
  md$type <- paste0(md$TissueType, "_", md$CancerType)
  N <- dim(md)[1]
  ctmd <- md %>% group_by(cell.type, type) %>% dplyr::count()
  tmd <- md %>% group_by(type) %>% dplyr::count()
  cmd <- md %>% group_by(cell.type) %>% dplyr::count()
  matRoe <- matrix(rep(0, length(unique(md$cell.type)) *
                          length(allMD_type)),
                   nrow = length(unique(md$cell.type)),
                   ncol = length(allMD_type))
  rownames(matRoe) <- as.character(sort(unique(md$cell.type)))
  colnames(matRoe) <- sort(unique(allMD_type))
  matCell <- matrix(rep(0, length(unique(md$cell.type)) *
                           length(allMD_type)),
                    nrow = length(unique(md$cell.type)),
                    ncol = length(allMD_type))
  rownames(matCell) <- as.character(sort(unique(md$cell.type)))
  colnames(matCell) <- sort(allMD_type)
  for (ci in 1:length(unique(md$cell.type))) {
    cin <- sort(unique(md$cell.type))[ci]
    k <- cmd$n[cmd$cell.type == cin]
    for (ti in 1:length(unique(md %>% pull(type)))) {
      tin <- sort(unique(md %>% pull(type)))[ti]
      ## print(ctmd)
      if (any(ctmd$cell.type == cin & ctmd$type == tin)) {
        n <- ctmd$n[ctmd$cell.type == cin & ctmd$type == tin]
      } else {
        n <- 0
      }
      M <- tmd$n[tmd$type == tin]
      matRoe[as.character(cin), tin] <- (n/M) / (k/N)
      matCell[as.character(cin), tin] <- n
    }
  }
  matRoe[matRoe>10] <- 10

  ## colnames(matRoe) <- paste0(colnames(matRoe), " (", colSums(matCell)[colnames(matRoe)], ")")
  rownames(matRoe) <- paste0(rownames(matRoe), " (", rowSums(matCell)[rownames(matRoe)], ")")

  ## rdx <- rowSums(matCell) >= 200
  ## matRoe <- matRoe[rdx, ]
  ## matCell <- matCell[rdx, ]

  matrixRoeL[[length(matrixRoeL)+1]] <- matRoe
  matrixCellL[[length(matrixCellL)+1]] <- matCell

  tmd <- ctmd %>%
    filter(type %in% colnames(matCell)) %>%
    filter(cell.type %in% rownames(matCell)) %>%
    mutate(DataSet = tempName)
  totalMD <- bind_rows(totalMD, tmd)
}

targetMatRoe <- rbind(matrixRoeL[[1]], matrixRoeL[[2]])
targetMatRoe <- rbind(targetMatRoe, matrixRoeL[[3]])
targetMatCell <- rbind(matrixCellL[[1]], matrixCellL[[2]])
targetMatCell <- rbind(targetMatCell, matrixCellL[[3]])

targetMatRoe <- targetMatRoe[rowSums(targetMatCell) > 100,
                             colSums(targetMatCell) > 100]

group_info <- str_extract(colnames(targetMatRoe), "^.*(?=_)")
## colnames(targetMatRoe) <- str_extract(colnames(targetMatRoe), "(?<=_).*")
col_annotation <- data.frame(TissueType = group_info)
rownames(col_annotation) <- colnames(targetMatRoe)
my.breaks <- seq(0, 1.99, by=0.01)
my.colors <- c(
    colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)/2),
    colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks)/2))
my.breaks.extra <- seq(2, 10, by = (10 - 2)/99)
my.colors.extra <- colorRampPalette(colors = c("#FD9AA0", "#550000"))(length(my.breaks.extra))
my.breaks <- c(my.breaks, my.breaks.extra)
my.colors <- c(my.colors, my.colors.extra)
widthi = ncol(targetMatRoe)/6 + 6
heighti = nrow(targetMatRoe)/10 + 3

pdf(file.path(figure_path, paste0("CD4CD8Innate_combined_roe.pdf")), width = widthi, height = heighti)
print(pheatmap(targetMatRoe,
               color = my.colors,
               breaks = my.breaks,
               ## annotation_col = col_annotation,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = T,
               cluster_rows = T,
               fontsize = 8,
               border_color = F))
dev.off()

###############################################################################
#                                  Figure 4d                                  #
###############################################################################
rm(list=ls())
library(Seurat)
library(tidyverse)
library(ggplot2)
library(stringr)
library(gtools)
library(reshape2)
library(ggdendro)

cmr_T <- read_csv("/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/scripts/pipelines/step4_figs/FigsForPaper/fig6/UKData_2_krupa/correlation_matrix.csv")
cmp_T <- read_csv("/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/scripts/pipelines/step4_figs/FigsForPaper/fig6/UKData_2_krupa/correlation_matrix_updated.csv")
CD4F_T <- read_csv("/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/scripts/pipelines/step4_figs/FigsForPaper/fig6/UKData_2_krupa/cd4_fulldata.csv")
CD8F_T <- read_csv("/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/scripts/pipelines/step4_figs/FigsForPaper/fig6/UKData_2_krupa/cd8_fulldata.csv")

figurePath <- "/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/scripts/pipelines/step4_figs/FigsForPaper/fig6/UKData_2_krupa"
if(!dir.exists(figurePath)){
    dir.create(figurePath)
}
setwd(figurePath)

#' draw cbe p heatmap #########################################################

cmr_T <- cmr_T %>%
    mutate(v1 = variable_1, v2 = variable_2, Coefficient = corr_out) %>%
    select(v1, v2, Coefficient)

cmp_T <- cmp_T %>%
    mutate(v1 = Var1, v2 = variable, P = value) %>%
    select(v1, v2, P)
cmr_T$P <- 1
cmr_T$P <- cmp_T$P

ggplotdata <- cmr_T
ggplotdata$pstars <- stars.pval(ggplotdata$P)
ggplotdata$pstarslevel <- 1
ggplotdata$pstarslevel[ggplotdata$pstars == "."] <- 2
ggplotdata$pstarslevel[ggplotdata$pstars == "*"] <- 3
ggplotdata$pstarslevel[ggplotdata$pstars == "**"] <- 4
ggplotdata$pstarslevel[ggplotdata$pstars == "***"] <- 5
tggplotdata <- ggplotdata
tggplotdata$v1 <- ggplotdata$v2
tggplotdata$v2 <- ggplotdata$v1
ggplotdata <- bind_rows(ggplotdata, tggplotdata) %>% distinct()

CorTable <- ggplotdata %>%
    mutate(ctv = Coefficient) %>%
    select(v1, v2, ctv) %>%
    dcast(v1 ~ v2)

CM <- as.matrix(CorTable[,2:dim(CorTable)[2]])
colnames(CM) <- colnames(CorTable)[2:dim(CorTable)[2]]
rownames(CM) <- CorTable[,1]
CM.dendro <- as.dendrogram(hclust(d = dist(x = CM)))
CM.plot <- ggdendrogram(data = CM.dendro, rotate = FALSE)
pdf(file.path(figurePath, paste0("CM_row_cluster.pdf")))
print(CM.plot)
dev.off()
row.order <- order.dendrogram(CM.dendro)
print(rownames(CM)[row.order])
CM.dendro <- as.dendrogram(hclust(d = dist(x = t(CM))))
CM.plot <- ggdendrogram(data = CM.dendro, rotate = TRUE)
pdf(file.path(figurePath, paste0("CM_col_dendro.pdf")))
print(CM.plot)
dev.off()
col.order <- order.dendrogram(CM.dendro)
print(colnames(CM)[col.order])
ggplotdata$v1 <- factor(x = ggplotdata$v1,
                          levels = rownames(CM)[row.order],
                          ordered = TRUE)
ggplotdata$v2 <- factor(x = ggplotdata$v2,
                              levels = colnames(CM)[col.order],
                              ordered = TRUE)
if(abs(min(ggplotdata$Coefficient)) > abs(max(ggplotdata$Coefficient))){
    colorRadius <- abs(min(ggplotdata$Coefficient))
}else{
    colorRadius <- abs(max(ggplotdata$Coefficient))
}

g <- ggplot(ggplotdata, aes(y = v1, x = v2)) +
    geom_point(aes(size = pstars, fill = Coefficient),
               color = "black", stroke = 0.1, shape = 22) +
    scale_fill_gradientn(
        colours = c("#286513", "#348319", "white", "#833419", "#652713"),
        values =scales::rescale(c(-colorRadius, -1, 0, 1, colorRadius)),
        limits = c(-colorRadius, colorRadius),
        guide = "colourbar") +
    scale_size_manual(values = c(0.1, 1, 2, 3, 4)) +
    theme_classic() +
    theme(text = element_text(size = 10),
          legend.position = "bottom",
          strip.background = element_blank(),
          axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
          legend.box="vertical", legend.margin=margin())
ggsave("correlation.pdf", g, width = 118, height = 140, units = "mm")


