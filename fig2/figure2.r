#--------------------------------------------------------------
# filename : figure2.r
# Date : 2022-10-20
# contributor : Yanshuo Chu
# function: figure2
# R version: R/4.0.3
#--------------------------------------------------------------

print('<==== figure2.r ====>')
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
library(scales)
set.seed(123)

CD8_path <- '/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/0_write_sample_info/CD8_2022-10-20.rds'
CD4_path <- '/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/0_write_sample_info/CD4_2022-10-20.rds'
Innate_path <- '/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/0_write_sample_info/Innate_2022-10-20.rds'
Proliferative_path <- '/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/0_write_sample_info/Proliferative_2022-10-20.rds'
Treg_path <- '/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/0_write_sample_info/Treg_2022-10-20.rds'
TFH_path <- '/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/0_write_sample_info/TFH_2022-10-20.rds'

CD8_Obj <- readRDS(CD8_path)
CD4_Obj <- readRDS(CD4_path)
Innate_Obj <- readRDS(Innate_path)
Proliferative_Obj <- readRDS(Proliferative_path)
Treg_Obj <- readRDS(Treg_path)
TFH_Obj <- readRDS(TFH_path)

allPath <- list(CD8 = CD8_path, CD4 = CD4_path, NK = Innate_path, Proliferative = Proliferative_path)
allObj <- list(CD8 = CD8_Obj, CD4 = CD4_Obj, NK = Innate_Obj, Proliferative = Proliferative_Obj)
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

figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/1_main_fig2/")
if(!dir.exists(figurePath)){
  dir.create(figurePath, recursive = T)
}
setwd(figurePath)


colorsForDataType <- c("#6DCCDD", "#EDCAE0", "#F494BE", "#F9B26C", "#A6ADCC", "#C4DA5D")
colorForClass1 <- c("#C4DA5D", "#6DCCDD", "#F494BE", "#EDCAE0")
dataTypeLevel <- c("CD4", "CD8", "MAIT", "Tgd", "NKT", "Proliferative")


###############################################################################
#'                             Manuscript: figure2A                          '#
###############################################################################
## umapColor <- colorRampPalette(colorsForDataType)(15)
Idents(CD8_Obj) <- CD8_Obj$seurat_clusters
umapColor <- c("#6DCCDD","#9ACBDE","#C8CADF","#EDC6DD","#F092B1","#F27FA5","#F47892", "#F6A395","#F8AD77","#E7B080","#C9AFA2","#ABADC5","#AEB9AC","#B9C984", "#C4DA5D")
png(file.path(figurePath, "figure2A_label.png"), width = 210 / 2, height = 210 / 2, units = "mm", res = 600)
DimPlot(CD8_Obj, label = T) +
    ## scale_color_manual(values = hue_pal( h.start = 0, direction = 1, l = 80)(15)) +
    ## scale_color_manual(values = umapColor) +
    scale_color_manual(values = colorRampPalette(umapColor)(length(unique(CD8_Obj$seurat_clusters)))) +
    scale_x_reverse()+
    theme_void() + theme(text = element_text(size = 6),
                         legend.position = "none")
dev.off()
png(file.path(figurePath, "figure2A.png"), width = 210 / 2, height = 210 / 2, units = "mm", res = 600)
DimPlot(CD8_Obj, label = F) +
    ## scale_color_manual(values = hue_pal( h.start = 0, direction = 1, l = 80)(15)) +
    ## scale_color_manual(values = umapColor) +
    scale_color_manual(values = colorRampPalette(umapColor)(length(unique(CD8_Obj$seurat_clusters)))) +
    scale_x_reverse()+
    theme_void() + theme(text = element_text(size = 6),
                         legend.position = "none")
dev.off()

###############################################################################
#'                             Manuscipt: figure2B                           '#
###############################################################################
Idents(CD8_Obj) <- CD8_Obj$seurat_clusters
markers <- c("GZMK", "GZMB", "PRF1", "CD44", "CD69", "FAS", "FASLG", "PDCD1",
"LAG3", "CTLA4", "FGFBP2", "GZMH", "GNLY", "NR4A1", "BAG3", "HSPA1A", "OSA1",
"IFIT1", "MX1", "DKK3", "CCR4", "EOMES", "CNN2", "LIMD2", "CD27", "LGR6",
"KLRC4", "CD244", "SEMA4A", "ITGA1", "KLRB1", "PRDM1", "CCR7", "SELL",
"TCF7","TRGV5", "TRGV10")
gene <- intersect(markers, rownames(CD8_Obj))
p <- DotPlot(CD8_Obj, features = rev(gene))
data <- p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]
## data$id <- factor(data$id, levels = rev(dataTypeLevel))
plotx <- ggplot(data, aes(x = id, y = features.plot)) +        ## global aes
    geom_point(aes(fill = avg.exp.scaled, size = pct.exp),
               color = 'black',
               shape = 21,
               stroke = 0.01)  +    ## geom_point for circle illusion
                                        #scale_fill_gradientn(colours=rev(color),limits=c(0,max(data$avg.exp)))+       ## color of the corresponding aes
    scale_x_discrete(breaks=0:13, labels=paste0("CD8_c", 0:13)) +
    xlab("") + ylab("") +
    ## scale_fill_gradient2(high = colorsForDataType[3], mid = "white", low = "#6DCCFF")+
    scale_fill_gradientn(
        ## limits = c(-1.5, 1.5),
        colors = c("#5DBCFF", "#6DCCFF", "white", colorsForDataType[3], "#F484AE"))+
    scale_size(range = c(0, 3.2), limits = c(0, 100), breaks = c(0,20,40,60,80,100))+             ## to tune the size of circles
    theme(
        text = element_text(size = 10),
        panel.grid.major = element_line(colour = "grey90", size=0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
        axis.text.y=element_text(angle = -30, vjust = 1, hjust = 1),
        legend.position="bottom",
        legend.title=element_text(size=10)) +
    guides(size = guide_legend(title.position="top",
                               title.hjust = 0.5,
                               ncol = 1,
                               byrow = T,
                               override.aes = list(stroke = 0.4)),
           fill = guide_colourbar(title.position = "top", title.hjust = 0.5))
ggsave(file.path(figurePath, "CD8_bubbleplot.pdf"), plotx, width = 210/3.4,
       height = 297/1.7, units = "mm")




###############################################################################
#'                              Manuscript: Fig2E                            '#
###############################################################################
CD8_Obj <- readRDS(CD8_path)
Idents(CD8_Obj) <- CD8_Obj$seurat_clusters
CD8FunctionList <- "/rsrch3/home/genomic_med/ychu2/configs/public/knowledge/database/Markers/CD8/forFig2"
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
CD8_Obj <-  AddModuleScore(CD8_Obj,
                           features = marker.list,
                           ctrl = 5,
                           name = "FunctionScore")
for(i in 1:length(marker.list)){
    colnames(CD8_Obj@meta.data)[colnames(CD8_Obj@meta.data) == paste0("FunctionScore", i)] <- names(marker.list)[i]
}
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
orderC = c("CD8_c13", "CD8_c3", "CD8_c6", "CD8_c0", "CD8_c11", "CD8_c9", "CD8_c10", "CD8_c12", "CD8_c8", "CD8_c2", "CD8_c7", "CD8_c4", "CD8_c5", "CD8_c1")
FunctionScoreMatrix <- FunctionScoreMatrix[,orderC]
my.breaks <- c(seq(-1, 0, by=0.1), seq(0.1, 1, by=0.1))
my.colors <- c(
    colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)/2),
    colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks)/2))
## cellType_col <- data.frame(cell.type = CD8_Obj_CellType)
## rownames(cellType_col) <- colnames(FunctionScoreMatrix)
signatureType_row <- data.frame(Signature.type = c(
                                    rep("Differentiation", length(Differentiation)),
                                    rep("Function", length(Function)),
                                    rep("Metabolism", length(Metabolism)),
                                    rep("Apoptosis", length(Apoptosis))))
rownames(signatureType_row) <- MarkerNameVector
pheatmap(FunctionScoreMatrix,
         show_colnames = T,
         show_rownames = T,
         ## annotation_col = cellType_col,
         annotation_row = signatureType_row,
         gaps_row = c(3, 14, 17),
         cluster_rows = F,
         cluster_cols = F,
         breaks = my.breaks,
         color = my.colors,
         border_color = "NA",
         fontsize = 8,
         width = 5,
         height = 3.8,
         filename = file.path(figurePath, paste0("CD8_FunctionScore_heatmap.pdf")))

###############################################################################
#'                              Manuscript: Fig2E                           '#
###############################################################################
CD8_Obj <- readRDS(CD8_path)
CD8FunctionList <- "/rsrch3/home/genomic_med/ychu2/configs/public/knowledge/database/Markers/CD8/forFig2"
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
CD8_Obj <-  AddModuleScore(CD8_Obj,
                           features = marker.list,
                           ctrl = 5,
                           name = "FunctionScore")
for(i in 1:length(marker.list)){
    colnames(CD8_Obj@meta.data)[colnames(CD8_Obj@meta.data) == paste0("FunctionScore", i)] <- names(marker.list)[i]
}
coord = Embeddings(object = CD8_Obj, reduction = "umap")
coord = coord[,c(1,2)]
colnames(coord) = c("UMAP1", "UMAP2")
coord = data.frame(ID = rownames(coord), coord)
meta = CD8_Obj@meta.data;
meta = data.frame(ID = rownames(meta), meta,stringsAsFactors = F)
meta = left_join(meta, coord, by = 'ID')
ggplotdata <- as_tibble(meta, rownames = NA) %>% select(UMAP1, UMAP2, Naive, `Activation.Effector.function`, Exhaustion, Cytotoxicity) %>% gather(value = Score, key = Signature, -UMAP1, -UMAP2)
ggplotdata$Signature[ggplotdata$Signature == "Activation.Effector.function"] <- "Activation/Effector"
ggplotdata$Signature <- factor(ggplotdata$Signature, levels = c("Naive", "Activation/Effector", "Cytotoxicity", "Exhaustion"))
a <- min(ggplotdata$Score)
b <- max(ggplotdata$Score)
g <- ggplot() +
    geom_point(data = ggplotdata, aes(x = UMAP1, y = UMAP2, color = Score), size = 0.01) +
    ## geom_point(data = ggplotdata %>% filter(Score <= 0.15, Score >= -0.2), aes(x = UMAP1, y = UMAP2, color = Score), size = 0.01) +
    ## geom_point(data = ggplotdata %>% filter(Score > 0.15), aes(x = UMAP1, y = UMAP2, color = Score), size = 0.01) +
    scale_color_gradientn(colors = c("#0000FF", "#8888FF", "#AAAAFF", "#FFFFFF", "#FF8888", "#FF5555", "#FF0000"), values = c(0, -a/3/(b-a), -a*2/3/(b-a), -a/(b-a), (1-a/(b-a))/3, (1-a/(b-a))*2/3, 1)) +
    ## scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
                          ## values = c(min(unique(ggplotdata$Score)), 0, max(unique(ggplotdata$Score)))) +
    scale_x_reverse()+
    facet_grid(. ~ Signature, scales = "free") +
    theme_void() +
    theme(text = element_text(size = 8),
          strip.background = element_rect(colour=NA, fill=NA),
          ## strip.text = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.2),
          strip.text = element_blank(),
          strip.placement = "none",
          panel.spacing.x=unit(0, "lines"),
          panel.spacing.y=unit(0, "lines"),
          legend.position = "none")
ggsave(file.path(figurePath, paste0("figure2D2.png")), g, device = "png", height = 297/10, width = 210/2, units = "mm")


###############################################################################
#                                  Figure 2F/G                                  #
###############################################################################

library(monocle3)
monocle.cds <- readRDS("/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/validate/CD8_V6/nPC_50/UMAP_dist_0.1_nneighbor_50/p1CD8_V6_UMAP_dist_0.1_nneighbor_50_CLUSTER_res_0.3/data_pseudotime.rds")
monocle.cds <- monocle.cds[, Cells(CD8_Obj)]
png(file.path(figurePath, "figure2trajectory.png"), width = 210 / 2, height = 210 / 2, units = "mm", res = 600)
plot_cells(cds = monocle.cds,
           color_cells_by = "pseudotime",
           show_trajectory_graph = F,
           trajectory_graph_color = "white",
           trajectory_graph_segment_size = 0.5,
           graph_label_size = 2,
           cell_size = 1,
           label_cell_groups = F,
           label_groups_by_cluster = F,
           label_branch_points = F,
           label_roots = F,
           label_leaves = F) +
    scale_color_gradientn(colours =
                          c("#0000DD", "#AAAAFF", "#EEEEEE", "#FFAAAA", "#AA0000"))+
    scale_x_reverse()+
    theme_void() + theme(text = element_text(size = 6),
                         legend.position = "none")
dev.off()



library(monocle3)
CD8_Obj <- readRDS(CD8_path)
monocle.cds <- readRDS(file.path(dirname(CD8_path), "data_pseudotime.rds"))
CD8FunctionList <- "/rsrch3/home/genomic_med/ychu2/configs/public/knowledge/database/Markers/CD8/forFig2"
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
CD8_Obj <-  AddModuleScore(CD8_Obj,
                           features = marker.list,
                           ctrl = 5,
                           name = "FunctionScore")
for(i in 1:length(marker.list)){
    colnames(CD8_Obj@meta.data)[colnames(CD8_Obj@meta.data) == paste0("FunctionScore", i)] <- names(marker.list)[i]
}
CD8_Obj@meta.data$Pseudotime <- monocle.cds@principal_graph_aux@listData$UMAP$pseudotime[match(Cells(CD8_Obj), names(monocle.cds@principal_graph_aux@listData$UMAP$pseudotime))]
CD8_Obj@meta.data$Pseudotime[is.infinite(CD8_Obj@meta.data$pseudotime)] <- sort(unique(CD8_Obj@meta.data$Pseudotime), decreasing = T)[2] + 1
md <- as_tibble(CD8_Obj@meta.data, rownames = NA) %>% select(seurat_clusters, Naive, Cytotoxicity, Exhaustion, `Activation:Effector function`, `Anergy`, `Senescence`, `Chemokine:Chemokine receptor`, `Cytokine:Cytokine receptor`, `IFN Response`, `MAPK Signaling`, `NFKB Signaling`, `TCR Signaling`, `Anti-apoptosis`, `Pro-apoptosis`, `Stress response`, Adhesion, Pseudotime)
path1 <- c(3, 0, 7, 1)
path2 <- c(3, 0, 2)
pathL <- list(path1 = path1, path2 = path2)
signatureL <- list(
    path1 = c("Naive", "Activation:Effector function", "Exhaustion", "Cytotoxicity"),
    path2 = c("Naive", "Activation:Effector function", "Exhaustion", "Cytotoxicity"))
totalMD <- c()
for(li in 1:length(pathL)){
    tempName <- names(pathL[li])
    tempPath <- pathL[[li]]
    tempMD <- md %>% filter(seurat_clusters %in% tempPath) %>% select(-seurat_clusters)
    tempMD <- tempMD %>% dplyr::select(matches(paste(c(signatureL[[li]], "Pseudotime"), collapse = "|")))
    tempMD <- tempMD %>%  gather(value = Value, key = Signature, -Pseudotime) %>% mutate(Path = tempName)
    totalMD <- bind_rows(totalMD, tempMD)
}
md <- as_tibble(CD8_Obj@meta.data, rownames = NA) %>%  select(Naive, Cytotoxicity, Exhaustion, `Activation:Effector function`, `TCR Signaling`, Pseudotime) %>% gather(value = Value, key = Signature, -Pseudotime) %>% mutate(Path = "all")
md <- bind_rows(totalMD, md)
md$Path <- factor(md$Path, levels = c("all", "path1", "path2"))
g <- ggplot(md) +
    stat_smooth(aes(x = Pseudotime, y = Value, color = Signature), method = "lm", formula = y ~ poly(x, 3))+
    facet_wrap(Path ~ .) +
    theme_classic() +
    theme(text = element_text(size = 10),
          legend.position = "right" ,
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.grid = element_line(),
          panel.grid.major = element_line(),   # Major grid lines
          panel.grid.minor = element_line(),   # Minor grid lines
          panel.grid.major.x = element_line(size = 0.5), # Vertical major grid lines
          panel.grid.major.y = element_line(size = 0.5), # Horizontal major grid lines
          panel.grid.minor.x = element_line(size = 0.25), # Vertical minor grid lines
          panel.grid.minor.y = element_line(size = 0.25)  # Vertical major grid lines
          )
ggsave(file.path(figurePath, paste0("figure2F.pdf")), g, width = 210, height = 297/5, units = "mm")

#' Naive, Cytotoxicity, Exhaustion score ######################################

###############################################################################
#                                  Figure 2j                                  #
###############################################################################


mycolors <- c("#99cc99", "#99ccff", "#cc6699", "#ffcc99")
names(mycolors) <- c("Healthy donor", "Uninvolved normal tissue", "Primary tumor tissue", "Metastatic tumor tissue")
md <- as_tibble(CD8_Obj@meta.data)
fracMD <- md %>%
  group_by(cell.type, TissueType) %>%
  count %>%
  ungroup %>%
  group_by(cell.type) %>%
  mutate(frac = n / (sum(n)))
fracMD$TissueType <- factor(fracMD$TissueType,
                            levels = names(mycolors))
g <- fracMD %>%
  ggplot() +
  geom_bar(aes(x = cell.type, y = frac,
               fill = TissueType),
           color = NA,
           position="stack", stat="identity") +
  scale_fill_manual(values = mycolors) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        strip.text.y = element_text(angle = 0),
        strip.background = element_rect(colour=NA, fill=NA))
ggsave(file.path(figurePath, paste0("figure2jup.pdf")), g)


N <- dim(md)[1]
cmd <- md %>% group_by(cell.type) %>% dplyr::count()
ctmd <- md %>% group_by(cell.type, TissueType) %>% dplyr::count()
tmd <- md %>% group_by(TissueType) %>% dplyr::count()
matRoe <- matrix(rep(0, length(unique(md$cell.type)) *
                        length(unique(md$TissueType))),
                  nrow = length(unique(md$TissueType)),
                  ncol = length(unique(md$cell.type)))
rownames(matRoe) <- as.character(names(mycolors))
colnames(matRoe) <- unique(md$cell.type)
matCell <- matrix(rep(0, length(unique(md$cell.type)) *
                        length(unique(md$TissueType))),
                 nrow = length(unique(md$TissueType)),
                 ncol = length(unique(md$cell.type)))
rownames(matCell) <- as.character(names(mycolors))
colnames(matCell) <- unique(md$cell.type)
for(ci in 1:length(unique(md$cell.type))){
  cin <- unique(md$cell.type)[ci]
  k <- cmd$n[cmd$cell.type == cin]
  for(ti in 1:length(unique(md %>% pull(TissueType)))) {
    tin <- unique(md %>% pull(TissueType))[ti]
    if(any(ctmd$cell.type == cin & ctmd$TissueType == tin)) {
      n <- ctmd$n[ctmd$cell.type == cin & ctmd$TissueType == tin]
    } else {
      n <- 0
    }
    M <- tmd$n[tmd$TissueType == tin]
    matRoe[tin, cin] <- (n/M) / (k/N)
    matCell[tin, cin] <- n
  }
}
matRoe[matRoe>3] = 3
my.breaks <- seq(0, 1.99, by=0.01)
my.colors <- c(
  colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)/2),
  colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks)/2))
my.breaks.extra <- seq(2, 3, by = (3 - 2)/99)
my.colors.extra <- colorRampPalette(colors = c("#FD9AA0", "#550000"))(length(my.breaks.extra))
my.breaks <- c(my.breaks, my.breaks.extra)
my.colors <- c(my.colors, my.colors.extra)
widthi = ncol(matRoe)/6 + 4.3
heighti = nrow(matRoe)/10 + 3
pdf(file.path(getwd(), "figure2jdown.pdf"))
print(pheatmap(matRoe,
               color = my.colors,
               breaks = my.breaks,
               cluster_rows = F,
               cluster_cols = F))
dev.off()

###############################################################################
#                                  Figure 2F old version                      #
###############################################################################



cell.type.order <- c( "CD8_c7_p-Tex", "CD8_c1_Tex")
Idents(CD8_Obj) <- CD8_Obj$cell.type
temp_obj <- subset(CD8_Obj, idents = cell.type.order)
markers <- c( "CXCL13", "TIGIT", "GNLY", "PRF1", "GZMK", "GZMH", "GZMB",
  "GZMA", "IRF4", "NR4A2", "NR4A1", "IFNG", "TNF", "CD27", "CD28", "LAYN",
  "CTLA4", "LAG3", "HAVCR2", "PDCD1", "CD38", "CD101", "TNFRSF9", "ENTPD1",
  "BATF", "TOX", "TCF7", "SELL", "CXCR3", "CXCR5", "CCR4", "CX3CR1", "FGFBP2",
  "FCGR3A")
gene <- intersect(markers, rownames(temp_obj))
p <- DotPlot(temp_obj, features = markers)
data <- p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]
data$id <- factor(data$id, levels = rev(cell.type.order))
plotx <- ggplot(data, aes(y = id, x = features.plot)) +
    geom_point(aes(fill = avg.exp.scaled, size = pct.exp),
               color = 'black',
               shape = 21,
               stroke = 0.01) +
    xlab("") + ylab("") +
    scale_fill_gradientn(
        colors = c("#5DBCFF", "#6DCCFF", "white", colorsForDataType[3], "#F484AE")) +
    scale_size(range = c(0, 3.2), limits = c(0, 100), breaks = c(0,20,40,60,80,100)) +
    theme(
        text = element_text(size = 10),
        panel.grid.major = element_line(colour = "grey90", size=0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
        legend.position="bottom",
        legend.title=element_text(size=10)) +
    guides(size = guide_legend(title.position="top",
                               title.hjust = 0.5,
                               byrow = T,
                               override.aes = list(stroke = 0.4)),
           fill = guide_colourbar(title.position = "top", title.hjust = 0.5))
ggsave(file.path(figurePath, "CD8_bubbleplot_3j.pdf"), plotx, width = 210/1.8,
       height = 297/5.5, units = "mm")




###############################################################################
#'                   Manuscript: fig2i fig3 umap and density                 '#
###############################################################################

theme_black <- function(base_size = 12, base_family = "") {
    theme_grey(base_size = base_size, base_family = base_family) %+replace%
        theme(
            axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.background = element_rect(color = NA, fill = "black"),
            legend.key = element_rect(color = "white",  fill = "black"),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = NULL,
            legend.key.width = NULL,
            legend.text = element_text(size = base_size*0.8, color = "white"),
            legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),
            legend.position = "none",
            legend.text.align = NULL,
            legend.title.align = NULL,
            legend.direction = "vertical",
            legend.box = NULL,
            panel.background = element_rect(fill = "black", color  =  NA),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.spacing = unit(0, "lines"),
            strip.background = element_rect(fill = "grey30", color = "grey10"),
            strip.text.x = element_text(size = base_size*0.8, color = "white"),
            strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),
            plot.background = element_rect(color = "black", fill = "black"),
            plot.title = element_text(size = base_size*1.2, color = "white"),
            plot.margin = unit(rep(0, 4), "lines")
        )
}
get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
}
umapColor <- c("#6DCCDD","#9ACBDE","#C8CADF","#EDC6DD","#F092B1","#F27FA5","#F47892",
               "#F6A395","#F8AD77","#E7B080","#C9AFA2","#ABADC5","#AEB9AC","#B9C984",
               "#C4DA5D")
for(oi in 1:length(allObj)){
    tempName <- names(allObj[oi])
    tempObj <- allObj[[oi]]
    Idents(tempObj) <- tempObj$seurat_clusters
    coord = Embeddings(object = tempObj, reduction = "umap")
    coord = coord[,c(1,2)]
    colnames(coord) = c("UMAP_1", "UMAP_2")
    coord = data.frame(ID = rownames(coord), coord)
    meta = tempObj@meta.data
    meta = data.frame(ID = rownames(meta), meta,stringsAsFactors = F)
    meta = left_join(meta, coord, by = 'ID')
    # randomly sample cells of large tissuetype ###############################
    print(tempName)
    print("table(meta$TissueType)")
    print(table(meta$TissueType))
    minNum <- min(table(meta$TissueType))
    print("minNum")
    print(minNum)
    meta <- meta %>%
        group_by(TissueType) %>%
        slice_sample(n = minNum)
    for(ti in 1: length(unique(meta$TissueType))){
        tin <- sort(unique(meta$TissueType))[ti]
        meta_tin <- meta[meta$TissueType == tin,]
        meta_ntin <- meta[meta$TissueType != tin,]
        g <- ggplot(data = meta_tin, mapping = aes(x = UMAP_1, y = UMAP_2)) +
            stat_density_2d(aes(fill = ..density..), geom = "raster", contour = F) +
            geom_point(color = 'white', size = .05) +
            scale_fill_viridis(option="magma") +
            theme_black()
        png(file.path(figurePath, paste0(tempName, "_TissueType-", tin, "_umap_density.png")))
        print(g)
        dev.off()
        g <- ggplot(data = meta_tin, mapping = aes(x = UMAP_1, y = UMAP_2)) +
            geom_point(aes(color = seurat_clusters), size = 0.5) +
            scale_color_manual(values = umapColor) +
            theme_void() +
            theme(text = element_text(size = 6),
                  legend.position = "none")
        png(file.path(figurePath, paste0(tempName, "_TissueType-",tin, "_umap.png")))
        print(g)
        dev.off()
    }
}


###############################################################################
#                       Sample fraction of Fig 2k/j/f/g                       #
###############################################################################

library(ggstatsplot)

# step1 set folder ############################################################
figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/1_main_fig2/")
if(!dir.exists(figurePath)) {
  dir.create(figurePath, recursive = T)
}
setwd(figurePath)

# sample cell number, larger than 200 #########################################
sampleN <- allMD %>%
  group_by(OrganSite, Sample) %>%
  count() %>%
  filter(n > 200)

# TypeSampleN total sample number of CancerType #####################################
TypeSampleN <- allMD
TypeSampleN$CancerType[TypeSampleN$TissueType %in%
                       c("Healthy donor", "Uninvolved normal tissue")] =
  TypeSampleN$OrganSite[TypeSampleN$TissueType %in%
                        c("Healthy donor", "Uninvolved normal tissue")]
TypeSampleN <- TypeSampleN %>%
  group_by(TissueType, CancerType, Sample) %>%
  count() %>%
  filter(Sample %in% sampleN$Sample)
TissueType.labs <- c("H", "U", "P", "M")
names(TissueType.labs) <- c("Healthy donor",
                            "Uninvolved normal tissue",
                            "Primary tumor tissue",
                            "Metastatic tumor tissue")

for(oi in 1:length(allObj)){

  tempName <- names(allObj[oi])
  tempObj <- allObj[[oi]]

  
  # set tempDataSet string ####################################################
  if(tempName == "Innate") {
    if(tempSeuratCluster %in% c(1,2)) {
      tempDataSet = "MAIT"
    }else if(tempSeuratCluster %in% c(3)) {
      tempDataSet = "Tgd"
    }else if(tempSeuratCluster %in% c(0,4)) {
      tempDataSet = "NKT"
    }
  }else {
    tempDataSet = tempName
  }

  # set obj sub folder ######################################################
  tempFigurePath <- file.path(figurePath, tempName)
  if(!dir.exists(tempFigurePath)){
    dir.create(tempFigurePath)
  }

  Idents(tempObj) <- tempObj@meta.data$cell.type

  # get cancer type sample number and fraction of current obj ###############
  ggplotdata <- as_tibble(tempObj@meta.data, rownames = NA) %>%
    select(cell.type,
           Sample,
           CancerType,
           TissueType,
           OrganSite,
           orig.ident)
  ggplotdata$CancerType[ggplotdata$TissueType %in%
                        c("Healthy donor", "Uninvolved normal tissue")] =
    ggplotdata$OrganSite[ggplotdata$TissueType %in%
                         c("Healthy donor", "Uninvolved normal tissue")]
  ggplotdata <- ggplotdata %>%
    group_by(cell.type,
             TissueType,
             CancerType,
             Sample) %>%
    dplyr::count() %>%
    filter(Sample %in% sampleN$Sample)
  ggplotdata$frac <- ggplotdata$n / sampleN$n[match(ggplotdata$Sample,
                                                    sampleN$Sample)]

  totalObj_ggplotdata <- c()

  # add mission value of each sub cluster to existing data #####################################
  for(tempCT in unique(ggplotdata$cell.type)) {
    # get existing data #######################################################
    temp_ggplotdata <- ggplotdata %>% filter(cell.type == tempCT)

    # add missing data ########################################################
    for(temp_TissueType in unique(TypeSampleN$TissueType)) {
      temp_TypeSampleN <- TypeSampleN %>%
        filter(TissueType == temp_TissueType)
      temp_diff_Samples <- setdiff(unique( temp_TypeSampleN %>% pull(Sample)),
                                   unique(temp_ggplotdata %>%
                                          filter(TissueType == temp_TissueType) %>%
                                          pull(Sample) ))
      if(length(temp_diff_Samples) > 0 ) {
        t_rows <- tibble(cell.type = rep(tempCT, length(temp_diff_Samples)),
                         TissueType = rep(temp_TissueType, length(temp_diff_Samples)),
                         CancerType = temp_TypeSampleN$CancerType[temp_TypeSampleN$Sample %in% temp_diff_Samples],
                         Sample = temp_diff_Samples,
                         n = rep(0, length(temp_diff_Samples)),
                         frac = rep(0, length(temp_diff_Samples)))
        temp_ggplotdata <- bind_rows(temp_ggplotdata, t_rows)
      }
    }
    temp_ggplotdata$TissueType <- factor(temp_ggplotdata$TissueType,
                                         levels = c("Healthy donor",
                                                    "Uninvolved normal tissue",
                                                    "Primary tumor tissue",
                                                    "Metastatic tumor tissue"))

    # add to totalObj_ggplotdata ##############################################
    totalObj_ggplotdata <- bind_rows(totalObj_ggplotdata, temp_ggplotdata)
  }

  if(tempName == "CD8") {
    # generat cell fraction plot of each sub CD8 clusters #####################

    # add x key ###############################################################
    temp_ggplotdata <- totalObj_ggplotdata %>%
      filter(cell.type %in% c("CD8_c3_Tn", "CD8_c0_t-Teff", "CD8_c1_Tex")) %>% 
      mutate(key = paste0(TissueType, "_", CancerType))
    temp_ggplotdata$cell.type <- factor(temp_ggplotdata$cell.type, levels = c("CD8_c3_Tn",
                                                                              "CD8_c0_t-Teff",
                                                                              "CD8_c1_Tex"))

    # sample number larger than 2 of each x key ###############################
    selectedT <- temp_ggplotdata %>%
      group_by(TissueType, CancerType, cell.type) %>%
      dplyr::count() %>%
      filter(n >= 3) %>%
      mutate(key = paste0(TissueType, "_", CancerType))
    temp_ggplotdata <- temp_ggplotdata %>% filter(key %in% selectedT$key)

    temp_ggplotdata$names <- paste0(temp_ggplotdata$TissueType, "_", temp_ggplotdata$CancerType)


    temp_ggplotdata <- temp_ggplotdata %>%
      mutate(TissueType = fct_recode(TissueType,
                                     H = "Healthy donor",
                                     U = "Uninvolved normal tissue",
                                     P = "Primary tumor tissue",
                                     M = "Metastatic tumor tissue"))

    temp_ggplotdata$CancerType <- factor(temp_ggplotdata$CancerType,
                                         levels = c(
                                           "BM",
                                           "Blood",
                                           "LN",
                                           "Breast",
                                           "Lung",
                                           "colon",
                                           "Tongue/Tonsil",
                                           "LGG",
                                           "AML",
                                           "HCC",
                                           "FL",
                                           "GBM",
                                           "STAD",
                                           "NSCLC",
                                           "BCC",
                                           "BRCA",
                                           "CRC",
                                           "DLBCL",
                                           "HNSC",
                                           "SKCM"))

    g <- ggplot(data = temp_ggplotdata) +
      geom_boxplot(aes(x = CancerType, y = frac, color = CancerType),
                   size = 0.3, outlier.shape = NA, width = 0.618, lwd = 0.1) +
      geom_jitter(shape = 16, position = position_jitter(0.2),
                  size = 0.5,  alpha = 0.9,
                  aes(x = CancerType, y = frac, fill = CancerType,
                      color = CancerType)) +
      scale_color_manual(values = colorRampPalette(colorsForDataType)(length(unique(temp_ggplotdata$CancerType)))) +
      facet_grid(cell.type ~ TissueType, scales = "free", space = "free_x") +
      xlab("") + ylab("") +
      theme_classic() +
      theme( text = element_text(size = 10),
            strip.background = element_blank(),
            axis.text.x = element_text(angle = 60 , vjust = 1, hjust=1),
            legend.position = "none")
    ggsave(file.path(tempFigurePath, paste0("Figure2j.pdf")),
           g, width = 122, height = 297/4, units = "mm")

    
    totalObj_ggplotdata$TissueTypeShort <- plyr::mapvalues(totalObj_ggplotdata$TissueType,
                                                      c("Healthy donor",
                                                        "Uninvolved normal tissue",
                                                        "Primary tumor tissue",
                                                        "Metastatic tumor tissue"),
                                                      c("H", "U", "P", "M"))
    g <- totalObj_ggplotdata %>%
      ggplot() +
      geom_boxplot(aes(x = TissueTypeShort,
                       y = frac,
                       color = TissueType),
                   size = 0.2, outlier.shape = NA) +
      geom_jitter(shape=16,
                  position=position_jitter(0.2),
                  size = 0.2,
                  aes(x = TissueTypeShort,
                      y = frac,
                      fill = CancerType,
                      color = CancerType)) +
      facet_wrap(cell.type ~ ., scales = "free") +
      xlab("") + ylab("") +
      scale_color_manual(values = mycolors) +
      theme_classic() +
      theme(text = element_text(size = 10),
            strip.background = element_blank(),
            legend.position = "none")
    ggsave(file.path(tempFigurePath,
                     paste0("figure2i.pdf")),
           g, width = 210, height = 200, units = "mm")

  }

    ## if(tempName == "CD4"){

    ##     temp_ggplotdata <- totalObj_ggplotdata %>%
    ##         filter(g_seurat_clusters %in% c("CD4_c1", "CD4_c3", "CD4_c4", "CD4_c8")) %>% mutate(key = paste0(TissueType, "_", CancerType))
    ##     selectedT <- temp_ggplotdata %>% group_by(TissueType, CancerType, g_seurat_clusters) %>% dplyr::count() %>% filter(n >= 3) %>% mutate(key = paste0(TissueType, "_", CancerType))
    ##     temp_ggplotdata <- temp_ggplotdata %>% filter(key %in% selectedT$key)
    ##     temp_ggplotdata$names <- paste(temp_ggplotdata$TissueType, temp_ggplotdata$CancerType)
    ##     temp_ggplotdata$CancerType <- factor(temp_ggplotdata$CancerType,
    ##                                          levels = c("Lymph node",
    ##                                                     "Breast",
    ##                                                     "Head and neck",
    ##                                                     "Lung",
    ##                                                     "Blood",
    ##                                                     "Brain",
    ##                                                     "BM",
    ##                                                     "colon",
    ##                                                     "LUSC",
    ##                                                     "HNSC",
    ##                                                     "BCC",
    ##                                                     "LUAD",
    ##                                                     "FL",
    ##                                                     "NSCLC",
    ##                                                     "BRCA",
    ##                                                     "PDAC",
    ##                                                     "SKCM",
    ##                                                     "OV",
    ##                                                     "STAD",
    ##                                                     "HCC",
    ##                                                     "GBM",
    ##                                                     "MCL",
    ##                                                     "AML",
    ##                                                     "LBCL",
    ##                                                     "CRC"))
    ##     temp_ggplotdata$g_seurat_clusters <- factor(temp_ggplotdata$g_seurat_clusters, levels = c("CD4_c1", "CD4_c3", "CD4_c4", "CD4_c8"))
    ##     temp_ggplotdata <- temp_ggplotdata %>% mutate(TissueType = fct_recode(TissueType, H = "Healthy donor", U = "Uninvolved normal tissue", P = "Primary tumor tissue", M = "Metastatic tumor tissue"))
    ##     g <- ggplot(data = temp_ggplotdata) +
    ##         geom_boxplot(aes(x = CancerType, y = frac, color = CancerType),
    ##                      size = 0.3, outlier.shape = NA, width = 0.618, lwd = 0.1) +
    ##         ## geom_violin(aes(x = CancerType, y = frac, color = CancerType)) +
    ##         geom_jitter(shape = 16, position = position_jitter(0.2),
    ##                     size = 0.5,  alpha = 0.9,
    ##                     aes(x = CancerType, y = frac, fill = CancerType,
    ##                         color = CancerType)) +
    ##         scale_color_manual(values = colorRampPalette(colorsForDataType)(length(unique(temp_ggplotdata$CancerType)))) +
    ##         facet_grid(g_seurat_clusters ~ TissueType,
    ##                    scales = "free", space = "free_x") +
    ##         xlab("") + ylab("") +
    ##         ## scale_y_log10() +
    ##         theme_classic() +
    ##         theme( text = element_text(size = 10),
    ##               strip.background = element_blank(),
    ##               axis.text.x = element_text(angle = 60 , vjust = 1, hjust=1),
    ##               legend.position = "none")
    ##     ## scale_x_discrete(breaks = temp_ggplotdata$names,
    ##     ##                  labels = temp_ggplotdata$CancerType)
    ##     ggsave(file.path(tempFigurePath, paste0("CD4_figure2_sample_boxplot.pdf")), g, width = 122, height = 297/4, units = "mm")

    ##     for(temp_TissueType in unique(temp_ggplotdata$TissueType)){
    ##         for(temp_CancerType in unique(temp_ggplotdata %>%
    ##                                       filter(TissueType == temp_TissueType) %>%
    ##                                       pull(CancerType))){
    ##             if(all(temp_ggplotdata %>% filter(TissueType == temp_TissueType) %>%
    ##                    filter(CancerType == temp_CancerType) %>% pull(frac) == 0)){
    ##                 next
    ##             }
    ##             g <- temp_ggplotdata %>% filter(TissueType == temp_TissueType) %>%
    ##                 filter(CancerType == temp_CancerType) %>%
    ##                 ggstatsplot::ggbetweenstats(
    ##                                  data = .,
    ##                                  x = g_seurat_clusters,
    ##                                  y = frac,
    ##                                  ## grouping.var = TissueType,
    ##                                  xlab = "Tissue type",
    ##                                  ylab = "Sample fraction",
    ##                                  ## pairwise.display = "all", # display only significant pairwise comparisons
    ##                                  p.adjust.method = "fdr", # adjust p-values for multiple tests using this method
    ##                                  ggtheme = theme_classic()
    ##                                  ## package = "ggsci",
    ##                                  ## palette = "default_jco",
    ##                                  ## plotgrid.args = list(nrow = 4)
    ##                              )
    ##             ggsave(file.path(tempFigurePath, paste0(temp_TissueType, "-", temp_CancerType, "_cancertype_boxplot_gbs.pdf")), g, width = 100, height = 200, units = "mm")
    ##         }
    ##     }


    ##     g <- totalObj_ggplotdata %>%
    ##         ggplot(data = .) +
    ##         geom_boxplot(aes(x = CancerType,
    ##                          y = frac,
    ##                          color = CancerType),
    ##                      size = 0.2, outlier.shape = NA) +
    ##         geom_jitter(shape=16, position=position_jitter(0.2),
    ##                     size = 0.2, aes(x = CancerType, y = frac,
    ##                                     fill = CancerType,
    ##                                     color = CancerType)) +
    ##         facet_grid(g_seurat_clusters ~ TissueType,
    ##                    scales = "free", space = "free_x") +
    ##         xlab("") + ylab("") +
    ##         scale_color_manual(values = colorRampPalette(colorsForDataType)(length(unique(totalObj_ggplotdata$CancerType)))) +
    ##         theme_classic() +
    ##         theme( text = element_text(size = 10),
    ##               strip.background = element_blank(),
    ##               axis.text.x = element_text(angle = 60 , vjust = 1, hjust=1),
    ##               legend.position = "none")
    ##     ## scale_x_discrete(breaks = temp_ggplotdata$names,
    ##     ##                  labels = temp_ggplotdata$CancerType)
    ##     ggsave(file.path(tempFigurePath, paste0("CD8_figure2_cancertype_sample_boxplot.pdf")), g, width = 400, height = 600, units = "mm")


    ##     for(temp_TissueType in unique(totalObj_ggplotdata$TissueType)){
    ##         temp_ggplotdata <- totalObj_ggplotdata %>% filter(TissueType == temp_TissueType) %>% mutate(key = paste0(TissueType, "_", CancerType))
    ##         temp_ggplotdata$CancerType[temp_ggplotdata$CancerType == "LUAD"] <- "NSCLC"
    ##         temp_ggplotdata$CancerType[temp_ggplotdata$CancerType == "LUSC"] <- "NSCLC"
    ##         selectedT <- temp_ggplotdata %>% group_by(TissueType, CancerType, g_seurat_clusters) %>% dplyr::count() %>% filter(n >= 3) %>% mutate(key = paste0(TissueType, "_", CancerType))
    ##         temp_ggplotdata <- temp_ggplotdata %>% filter(key %in% selectedT$key)
    ##         temp_mean_ggplotdata <- temp_ggplotdata %>% group_by(CancerType, g_seurat_clusters) %>% summarise(meanFrac = mean(frac))
    ##         CancerType_names <- unique(temp_ggplotdata$CancerType)
    ##         psimMatrix <- matrix(rep(NA, length(unique(temp_ggplotdata$CancerType)) ^ 2),
    ##                              nrow = length(unique(temp_ggplotdata$CancerType)),
    ##                              ncol = length(unique(temp_ggplotdata$CancerType)))
    ##         rownames(psimMatrix) <- sort(unique(temp_ggplotdata$CancerType))
    ##         colnames(psimMatrix) <- sort(unique(temp_ggplotdata$CancerType))
    ##             for(i in 1: (length(CancerType_names) - 1)){
    ##                 for(j in (i + 1):length(CancerType_names)){
    ##                     i_CancerType_name <- CancerType_names[i]
    ##                     j_CancerType_name <- CancerType_names[j]
    ##                     i_fracs <- temp_mean_ggplotdata %>% filter(CancerType == i_CancerType_name) %>% pull(meanFrac)
    ##                     j_fracs <- temp_mean_ggplotdata %>% filter(CancerType == j_CancerType_name) %>% pull(meanFrac)
    ##                     t.res <- t.test(i_fracs, j_fracs, paired = TRUE, alternative = "two.sided")
    ##                     psimMatrix[i_CancerType_name, j_CancerType_name] <- t.res[[3]]
    ##                     psimMatrix[j_CancerType_name, i_CancerType_name] <- t.res[[3]]
    ##                 }
    ##             }

    ##         psimMatrix <- psimMatrix[rownames(psimMatrix) != "MCL", colnames(psimMatrix) != "MCL"]
    ##         my.breaks0 <- seq(min(psimMatrix[!is.na(psimMatrix)]), 0.05, by=0.002)
    ##         my.breaks1 <- seq(0.05, max(psimMatrix[!is.na(psimMatrix)]), by=0.002)
    ##         my.colors <- c(
    ##             colorRampPalette(colors = c("red", "white"))(length(my.breaks0)),
    ##             colorRampPalette(colors = c("white", "blue"))(length(my.breaks1)))
    ##         my.breaks = c(my.breaks0, my.breaks1)
    ##         pdf(file.path(tempFigurePath, paste0(temp_TissueType, "_CancerType_pSimMatrix.pdf")))
    ##         print(pheatmap(psimMatrix,
    ##                        color = my.colors,
    ##                        breaks = my.breaks,
    ##                        ## annotation_row = cellType_col,
    ##                        show_colnames = T,
    ##                        show_rownames = T,
    ##                        cluster_cols = T,
    ##                        cluster_rows = T,
    ##                        border_color = F))
    ##         dev.off()

    ##         CancerType_cell_number_frac <- allMD %>% filter(DataSet == "CD8") %>% mutate(g_seurat_clusters = paste0(DataSet, "_c", seurat_clusters)) %>% filter(TissueType == temp_TissueType) %>% filter(CancerType %in% unique(temp_ggplotdata$CancerType)) %>% select(g_seurat_clusters, CancerType, barcode, DataSet) %>% group_by(CancerType, g_seurat_clusters) %>% dplyr::count() %>% group_by(CancerType) %>% mutate(frac = n / sum(n))
    ##         g <- ggplot(CancerType_cell_number_frac) +
    ##             geom_bar(aes(x = CancerType, y = frac,
    ##                          fill = g_seurat_clusters),
    ##                      color = "white",
    ##                      position="stack", stat="identity") +
    ##             theme_classic()
    ##         ggsave(file.path(tempFigurePath, paste0(temp_TissueType, "_stackfrac_cancertype.pdf")), g)


    ##         fracMatrix <- matrix(rep(0, length(unique(CancerType_cell_number_frac$CancerType)) *
    ##                                     length(unique(CancerType_cell_number_frac$g_seurat_clusters))),
    ##                              ncol = length(unique(CancerType_cell_number_frac$CancerType)),
    ##                              nrow = length(unique(CancerType_cell_number_frac$g_seurat_clusters)))
    ##         rownames(fracMatrix) <- sort(unique(CancerType_cell_number_frac$g_seurat_clusters))
    ##         colnames(fracMatrix) <- sort(unique(CancerType_cell_number_frac$CancerType))
    ##         for(temp_rn in sort(unique(CancerType_cell_number_frac$g_seurat_clusters))){
    ##             for(temp_cl in sort(unique(CancerType_cell_number_frac$CancerType))){
    ##                 if(!any(CancerType_cell_number_frac$CancerType == temp_cl & CancerType_cell_number_frac$g_seurat_clusters == temp_rn)){
    ##                     next
    ##                 }
    ##                 fracMatrix[temp_rn, temp_cl] <- CancerType_cell_number_frac$frac[CancerType_cell_number_frac$CancerType == temp_cl & CancerType_cell_number_frac$g_seurat_clusters == temp_rn]
    ##             }
    ##         }
    ##         my.breaks <- seq(min(fracMatrix[!is.na(fracMatrix)]), max(fracMatrix[!is.na(fracMatrix)]), by=0.002)
    ##         my.colors <- c(
    ##             colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)/2),
    ##             colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks)/2))
    ##         pdf(file.path(tempFigurePath, paste0(temp_TissueType, "_CancerType_cluster_fraction_heatmap.pdf")))
    ##         print(pheatmap(fracMatrix,
    ##                        color = my.colors,
    ##                        breaks = my.breaks,
    ##                        ## annotation_row = cellType_col,
    ##                        show_colnames = T,
    ##                        show_rownames = T,
    ##                        cluster_cols = F,
    ##                        cluster_rows = F,
    ##                        border_color = F))
    ##         dev.off()
    ##     }


    ##     #' plus skmc  #########################################################
    ##     temp_ggplotdata <- totalObj_ggplotdata %>% filter(TissueType == "Primary tumor tissue") %>% mutate(key = paste0(TissueType, "_", CancerType))
    ##     temp_ggplotdata$CancerType[temp_ggplotdata$CancerType == "LUAD"] <- "NSCLC"
    ##     temp_ggplotdata$CancerType[temp_ggplotdata$CancerType == "LUSC"] <- "NSCLC"
    ##     SKCM_temp_ggplotdata <- totalObj_ggplotdata %>% filter(CancerType == "SKCM") %>% mutate(key = paste0(TissueType, "_", CancerType))
    ##     temp_ggplotdata <- bind_rows(temp_ggplotdata, SKCM_temp_ggplotdata)
    ##     selectedT <- temp_ggplotdata %>% group_by(TissueType, CancerType, g_seurat_clusters) %>% dplyr::count() %>% filter(n >= 3) %>% mutate(key = paste0(TissueType, "_", CancerType))
    ##     temp_ggplotdata <- temp_ggplotdata %>% filter(key %in% selectedT$key)
    ##     temp_mean_ggplotdata <- temp_ggplotdata %>% group_by(CancerType, g_seurat_clusters) %>% summarise(meanFrac = mean(frac))
    ##     CancerType_names <- unique(temp_ggplotdata$CancerType)
    ##     psimMatrix <- matrix(rep(NA, length(unique(temp_ggplotdata$CancerType)) ^ 2),
    ##                          nrow = length(unique(temp_ggplotdata$CancerType)),
    ##                          ncol = length(unique(temp_ggplotdata$CancerType)))
    ##     rownames(psimMatrix) <- sort(unique(temp_ggplotdata$CancerType))
    ##     colnames(psimMatrix) <- sort(unique(temp_ggplotdata$CancerType))
    ##     for(i in 1: (length(CancerType_names) - 1)){
    ##         for(j in (i + 1):length(CancerType_names)){
    ##             i_CancerType_name <- CancerType_names[i]
    ##             j_CancerType_name <- CancerType_names[j]
    ##             i_fracs <- temp_mean_ggplotdata %>% filter(CancerType == i_CancerType_name) %>% pull(meanFrac)
    ##             j_fracs <- temp_mean_ggplotdata %>% filter(CancerType == j_CancerType_name) %>% pull(meanFrac)
    ##             t.res <- t.test(i_fracs, j_fracs, paired = TRUE, alternative = "two.sided")
    ##             psimMatrix[i_CancerType_name, j_CancerType_name] <- t.res[[3]]
    ##             psimMatrix[j_CancerType_name, i_CancerType_name] <- t.res[[3]]
    ##         }
    ##     }

    ##     psimMatrix <- psimMatrix[rownames(psimMatrix) != "MCL", colnames(psimMatrix) != "MCL"]
    ##     my.breaks <- seq(min(psimMatrix[!is.na(psimMatrix)]), max(psimMatrix[!is.na(psimMatrix)]), by=0.002)
    ##     my.colors <- c(
    ##         colorRampPalette(colors = c("#FD9AA0", "white"))(length(my.breaks)/2),
    ##         colorRampPalette(colors = c("white", "#6DCCFD"))(length(my.breaks)/2))
    ##     pdf(file.path(tempFigurePath, paste0("CancerType_pSimMatrix.pdf")))
    ##     print(pheatmap(psimMatrix,
    ##                    color = my.colors,
    ##                    breaks = my.breaks,
    ##                    ## annotation_row = cellType_col,
    ##                    show_colnames = T,
    ##                    show_rownames = T,
    ##                    cluster_cols = T,
    ##                    cluster_rows = T,
    ##                    border_color = F))
    ##     dev.off()
    ## }

}

