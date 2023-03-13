#--------------------------------------------------------------
# filename : figure3.r
# Date : 2022-10-24
# contributor : Yanshuo Chu
# function: figure3
# R version: R/4.0.3
#--------------------------------------------------------------

print('<==== figure3.r ====>')
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
for(oi in 1:length(allObj)) {
    tempName <- names(allObj[oi])
    tempObj <- allObj[[oi]]
    tempObj$barcode <- Cells(tempObj)
    if(oi == 1) {
        common_columns <- colnames(tempObj@meta.data)
    }else {
        common_columns <- intersect(common_columns, colnames(tempObj@meta.data))
    }
}
allMD <- c()
for(oi in 1:length(allObj)) {
    tempName <- names(allObj[oi])
    tempObj <- allObj[[oi]]
    tempMD <- as_tibble(tempObj@meta.data, rownames = NA) %>%
        dplyr::select(all_of(common_columns)) %>% mutate(barcode = Cells(tempObj))
    allMD <- bind_rows(allMD, tempMD)
}

figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/1_main_fig3/")
if(!dir.exists(figurePath)) {
  dir.create(figurePath, recursive = T)
}
setwd(figurePath)

colorsForDataType <- c("#6DCCDD", "#EDCAE0", "#F494BE", "#F9B26C", "#A6ADCC", "#C4DA5D")
colorForClass1 <- c("#C4DA5D", "#6DCCDD", "#F494BE", "#EDCAE0")
dataTypeLevel <- c("CD4", "CD8", "MAIT", "Tgd", "NKT", "Proliferative")

###############################################################################
#                                  Figure 3a                                  #
###############################################################################

umapColor <- colorRampPalette(colorsForDataType)(12)
Idents(CD4_Obj) <- CD4_Obj$cell.type
png(file.path(figurePath, "figure3A.png"), width = 210 / 2, height = 210 / 2, units = "mm", res = 600)
DimPlot(CD4_Obj, label = T) +
    scale_color_manual(values = umapColor) +
    theme_void() + theme(text = element_text(size = 6),
                         legend.position = "none")
dev.off()

###############################################################################
#                                  Figure 3b                                  #
###############################################################################

Idents(CD4_Obj) <- CD4_Obj$seurat_clusters
markers <- c(
"IL7R",
"GPR183",
"CD69",
"FOXP3",
"IL2RA",
"CTLA4",
"TNFRSF4",
"RPL31",
"RPL21",
"CXCL13",
"PDCD1",
"TOX",
"ICOS",
"BCL6",
"NR4A1",
"BAG3",
"FOS",
"JUN",
"IFNG",
"GZMA",
"GZMH",
"GZMB",
"GZMK",
"NKG7",
"PRF1",
"FHIT",
"CCR7",
"LEF1",
"SELL",
"TCF7",
"TCEA3",
"CDC25B",
"IL17A",
"IL17F",
"RORA",
"KLRB1",
"CCR6",
"SLC40A1",
"ANKRD55",
"ISG15",
"IFI44L",
"IFIT1",
"EOMES")
gene <- intersect(markers, rownames(CD4_Obj))
p <- DotPlot(CD4_Obj, features = rev(gene))
data <- p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]
## data$id <- factor(data$id, levels = rev(dataTypeLevel))
plotx <- ggplot(data, aes(y = id, x = features.plot)) +        ## global aes
    geom_point(aes(fill = avg.exp.scaled, size = pct.exp),
               color = 'black',
               shape = 21,
               stroke = 0.005) + ## geom_point for circle illusion
                                  # scale_fill_gradientn(colours=rev(color),limits=c(0,max(data$avg.exp))) +
    ## color of the corresponding aes
    scale_y_discrete(breaks=0:11, labels=paste0("CD4_c", 0:11), limits = rev) +
    scale_x_discrete(limits = rev) +
    xlab("") + ylab("") +
    ## scale_fill_gradient2(high = colorsForDataType[3], mid = "white", low = "#6DCCFF")+
    scale_fill_gradientn(
        ## limits = c(-1.5, 1.5),
        colors = c("#5DBCFF", "#6DCCFF", "white", colorsForDataType[3], "#F484AE"))+
    scale_size(range = c(0, 3.5), limits = c(0, 100), breaks = c(0,20,40,60,80,100))+             ## to tune the size of circles
    theme(
        text = element_text(size = 10),
        panel.grid.major = element_line(colour = "grey90", size=0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
        legend.position="bottom",
        legend.title=element_text(size=8)) +
    guides(size = guide_legend(title.position="top",
                               title.hjust = 0.5,
                               override.aes = list(stroke = 0.4)),
           fill = guide_colourbar(title.position = "top", title.hjust = 0.5))
ggsave(file.path(figurePath, "CD4_bubbleplot.pdf"), plotx, height = 210/2.6, width = 140, units = "mm")

###############################################################################
#                                  Figure 3c                                  #
###############################################################################

library(scales)
CD4_Obj <- readRDS(CD4_path)
CD4FunctionList <- "/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/knowledge/public/Markers/CD4/forFig3"
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
for(i in 1:length(marker.list)) {
    colnames(CD4_Obj@meta.data)[colnames(CD4_Obj@meta.data) == paste0("FunctionScore", i)] <- names(marker.list)[i]
}
Idents(CD4_Obj) <- CD4_Obj$seurat_clusters
Differentiation <- c("Naive", "Activation:Effector function", "Exhaustion")
Function <- c("TCR signaling", "Cytotoxicity", "Cytokine:Cytokine receptor",
              "Chemokine:Chemokine receptor", "Stress response", "Adhesion",
              "IFN response", "Treg signature", "Costimulatory molecules")
Metabolism <- c("OXPHOS", "Glycolysis", "Lipid metabolism")
Apoptosis <- c("Pro-apoptosis", "Anti-apoptosis")
MarkerNameVector <- c(Differentiation, Function, Metabolism, Apoptosis)
FunctionScoreMatrix <- matrix(0,
                              ncol = length(unique(CD4_Obj$seurat_clusters)),
                              nrow = length(MarkerNameVector))
colnames(FunctionScoreMatrix) <- paste0("CD4_c", unique(CD4_Obj$seurat_clusters))
rownames(FunctionScoreMatrix) <- MarkerNameVector
for(ci in 1:ncol(FunctionScoreMatrix)) {
    for(ri in 1:nrow(FunctionScoreMatrix)) {
        FunctionVec <- as_tibble(CD4_Obj@meta.data) %>% pull(MarkerNameVector[ri])
        fv <- mean(FunctionVec[CD4_Obj$seurat_clusters == unique(CD4_Obj$seurat_clusters)[ci]])
        FunctionScoreMatrix[ri, ci] <- fv
    }
}
FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, rescale, to=c(-1, 1)))
my.breaks <- c(seq(-1, 0, by=0.1), seq(0.1, 1, by=0.1))
my.colors <- c(
    colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)/2),
    colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks)/2))
signatureType_row <- data.frame(Signature.type = c(
                                    rep("Differentiation", length(Differentiation)),
                                    rep("Function", length(Function)),
                                    rep("Metabolism", length(Metabolism)),
                                    rep("Apoptosis", length(Apoptosis))))
rownames(signatureType_row) <- MarkerNameVector
col_name_v <- paste0("CD4_c", c(1,3,8,7,9,10,2,6,4,11,0,5))
FunctionScoreMatrix <- FunctionScoreMatrix[, col_name_v]
pheatmap(FunctionScoreMatrix,
         show_colnames = T,
         show_rownames = T,
         annotation_row = signatureType_row,
         gaps_row = c(3, 12, 15),
         cluster_rows = F,
         cluster_cols = F,
         breaks = my.breaks,
         color = my.colors,
         border_color = "NA",
         fontsize = 8,
         width = 5,
         height = 3.8,
         filename = file.path(figurePath, paste0("CD4_FunctionScore_heatmap.pdf")))



## library(ROGUE)
## totalRogue <- c()
## for(ct in as.vector(unique(CD4_Obj@meta.data$seurat_clusters))) {
##     targetMatrix <- CD4_Obj@assays$RNA@data[, rownames(CD4_Obj@meta.data[CD4_Obj@meta.data$seurat_clusters == ct,])]
##     if(dim(targetMatrix)[2] > 10000) {
##         ## for(i in 1:100) {
##             targetMatrix <- targetMatrix[, sample(colnames(targetMatrix), 10000, replace = F)]
##             targetMatrix <- matr.filter(as.matrix(targetMatrix), min.cells = 10, min.genes = 10)
##             targetSE <- SE_fun(targetMatrix)
##             ct.rogue.value <- CalculateRogue(targetSE, platform = "UMI")
##             tempTable <- tibble(seurat_clusters = paste0("CD4_c", ct), rogue = ct.rogue.value)
##             totalRogue <- bind_rows(totalRogue, tempTable)
##         ## }
##     }
##     else{
##         targetMatrix <- matr.filter(as.matrix(targetMatrix), min.cells = 10, min.genes = 10)
##         targetSE <- SE_fun(targetMatrix)
##         ct.rogue.value <- CalculateRogue(targetSE, platform = "UMI")
##         tempTable <- tibble(seurat_clusters = paste0("CD4_c", ct), rogue = ct.rogue.value)
##         totalRogue <- bind_rows(totalRogue, tempTable)
##     }
## }
## ## write_tsv(totalRogue, file.path(figurePath, paste0('totalRogue', "_", Sys.Date(), '.tsv')))
## pdf(file.path(figurePath, "rogue_CD4_linepoint.pdf"))
## ggplot(totalRogue, aes(x = seurat_clusters, y = rogue)) +
##     geom_point()
## dev.off()



md <- as_tibble(CD4_Obj@meta.data, rownames = NA) %>% mutate(barcode = Cells(CD4_Obj)) %>%
    select(barcode, seurat_clusters)
lmd <- md %>% filter(seurat_clusters %in% 0:5) %>% group_by(seurat_clusters) %>%  sample_n(10000)
md <- md %>% filter(! seurat_clusters %in% 0:5)
md <- bind_rows(lmd, md)

rogue.res <- rogue(CD4_Obj@assays$RNA@counts[,md$barcode], labels = CD4_Obj@meta.data$seurat_clusters[Cells(CD4_Obj) %in% md$barcode], samples = CD4_Obj@meta.data$batch[Cells(CD4_Obj) %in% md$barcode], platform = "UMI", span = 0.6)

pdf(file.path(figurePath, "rogue_CD4_boxplot.pdf"))
rogue.boxplot(rogue.res)
dev.off()

###############################################################################
#'                   Figure 3d                                               '#
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
for(oi in 1:length(allObj)) {
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
    for(ti in 1: length(unique(meta$TissueType))) {
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
#                       Figure 3F and 3G                                      #
###############################################################################

library(ggstatsplot)

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

for(oi in 1:length(allObj)) {

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
  if(!dir.exists(tempFigurePath)) {
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

  if(tempName == "CD4") {
    # generat cell fraction plot of each sub CD8 clusters #####################

    # add x key ###############################################################
    temp_ggplotdata <- totalObj_ggplotdata %>%
      filter(cell.type %in% c("CD4_c1_Treg", "CD4_c3_Tfh", "CD4_c0_Tcm")) %>% 
      mutate(key = paste0(TissueType, "_", CancerType))
    temp_ggplotdata$cell.type <- factor(temp_ggplotdata$cell.type, levels = c("CD4_c1_Treg",
                                                                              "CD4_c3_Tfh",
                                                                              "CD4_c0_Tcm"))

    # sample number larger than 2 of each x key ###############################
    selectedT <- temp_ggplotdata %>%
      group_by(TissueType, CancerType, cell.type) %>%
      dplyr::count() %>%
      filter(n >= 3) %>%
      mutate(key = paste0(TissueType, "_", CancerType))
    temp_ggplotdata <- temp_ggplotdata %>% filter(key %in% union(selectedT$key, "Primary tumor tissue_STAD"))

    temp_ggplotdata$names <- paste0(temp_ggplotdata$TissueType, "_", temp_ggplotdata$CancerType)


    temp_ggplotdata <- temp_ggplotdata %>%
      mutate(TissueType = fct_recode(TissueType,
                                     H = "Healthy donor",
                                     U = "Uninvolved normal tissue",
                                     P = "Primary tumor tissue",
                                     M = "Metastatic tumor tissue"))

    temp_ggplotdata$CancerType <- factor(temp_ggplotdata$CancerType,
                                         levels = c(
                                           "LN",
                                           "Breast",
                                           "Tongue/Tonsil",
                                           "Lung",
                                           "Blood",
                                           "BM",
                                           "colon",
                                           "Brain",
                                           "HNSC",
                                           "BCC",
                                           "NSCLC",
                                           "SKCM",
                                           "FL",
                                           "BRCA",
                                           "GBM",
                                           "HCC",
                                           "AML",
                                           "DLBCL",
                                           "LGG",
                                           "STAD",
                                           "CRC"))

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
    ggsave(file.path(tempFigurePath, paste0("Figure3g.pdf")),
           g, width = 122, height = 297/4, units = "mm")

    

    mycolors <- c("#99cc99", "#99ccff", "#cc6699", "#ffcc99")
    names(mycolors) <- c("Healthy donor", "Uninvolved normal tissue", "Primary tumor tissue", "Metastatic tumor tissue")
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
                     paste0("figure3f.pdf")),
           g, width = 210, height = 200, units = "mm")

  }
}




###############################################################################
#'                                    Treg                                   '#
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
library(scales)
#' load seurat object #########################################################
Treg_path <- '/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/0_write_sample_info/Treg_2022-10-20.rds'
Treg_Obj <- readRDS(Treg_path)
TFH_path <- '/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/0_write_sample_info/TFH_2022-10-20.rds'
TFH_Obj <- readRDS(TFH_path)
CD4_path <- '/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/0_write_sample_info/CD4_2022-10-20.rds'
CD4_Obj <- readRDS(CD4_path)
CD8_path <- '/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/0_write_sample_info/CD8_2022-10-20.rds'
CD8_Obj <- readRDS(CD8_path)
Innate_path <- '/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/0_write_sample_info/Innate_2022-10-20.rds'
Innate_Obj <- readRDS(Innate_path)
P_path <- '/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/0_write_sample_info/Proliferative_2022-10-20.rds'
P_Obj <- readRDS(P_path)


figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/1_main_fig3/")
if(!dir.exists(figurePath)){
  dir.create(figurePath, recursive = T)
}
setwd(figurePath)

colorsForDataType <- c("#6DCCDD", "#EDCAE0", "#F494BE", "#F9B26C", "#A6ADCC", "#C4DA5D")
colorForClass1 <- c("#C4DA5D", "#6DCCDD", "#F494BE", "#EDCAE0")
dataTypeLevel <- c("CD4", "CD8", "NK", "NKT", "MAIT", "Proliferative")

Treg_Obj <- subset(Treg_Obj, cells = Cells(Treg_Obj)[Cells(Treg_Obj) %in% Cells(CD4_Obj)])
TFH_Obj <- subset(TFH_Obj, cells = Cells(TFH_Obj)[Cells(TFH_Obj) %in% Cells(CD4_Obj)])
for(tempColumn in c("CancerType", "TissueType", "OrganSite", "Disease", "Patient", "Sample")) {
    Treg_Obj@meta.data[,tempColumn] <- CD4_Obj@meta.data[match(Cells(Treg_Obj), Cells(CD4_Obj)),tempColumn]
    TFH_Obj@meta.data[,tempColumn] <- CD4_Obj@meta.data[match(Cells(TFH_Obj), Cells(CD4_Obj)),tempColumn]
}

allObj <- list(Treg = Treg_Obj,
               TFH = TFH_Obj,
               Innate = Innate_Obj)

###############################################################################
#                                  Figure 3h                                  #
###############################################################################
umapColor <- colorRampPalette(colorsForDataType)(length(unique(Treg_Obj$seurat_clusters)))
png(file.path(figurePath, "Treg_UMAP.png"), width = 210/3.5, height = 210/3.5, units = "mm", res = 600)
DimPlot(Treg_Obj, label = T) +
    scale_color_manual(values = umapColor) +
    theme_void() + theme(text = element_text(size = 5),
                         legend.position = "none")
dev.off()

###############################################################################
#                                  Figure 3i                                  #
###############################################################################
umapColor <- colorRampPalette(colorsForDataType)(length(unique(TFH_Obj$seurat_clusters)))
Idents(TFH_Obj) <- TFH_Obj$seurat_clusters
png(file.path(figurePath, "TFH_UMAP.png"), width = 210/3.5,
    height = 210/3.5, units = "mm", res = 600)
DimPlot(TFH_Obj, label = T) +
    scale_color_manual(values = umapColor) +
    theme_void() + theme(text = element_text(size = 5),
                         legend.position = "none")
dev.off()

###############################################################################
#                                  Figure 3j                                  #
###############################################################################
Idents(Treg_Obj) <- Treg_Obj$seurat_clusters
markers <- c(
    "KLF2",
    "LEF1",
    "SELL",
    "IL1R2",
    "IL1R1",
    "TNFRSF4",
    "TNFRSF18",
    "TNFRSF9",
    "IL2RA",
    "BATF",
    "IL21R",
    "FOXP3",
    "FANK1",
    "GATA3",
    "TNFRSF13B",
    "GPX1",
    "IL2RA",
    "HSPA1A",
    "HSPA1B",
    "BAG3",
    "HSPE1",
    "HSPB1",
    "HSPD1",
    "NR4A1",
    "IL7R",
    "LRRC32",
    "CCR7",
    "ICA1",
    "CXCL13",
    "TCF7",
    "KLRB1")
gene <- intersect(markers, rownames(Treg_Obj))
p <- DotPlot(Treg_Obj, features = rev(gene))
data <- p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]
## data$id <- factor(data$id, levels = rev(dataTypeLevel))
plotx <- ggplot(data, aes(y = id, x = features.plot)) +        ## global aes
    geom_point(aes(fill = avg.exp.scaled, size = pct.exp),
               color = 'black',
               shape = 21,
               stroke = 0.005) +
    scale_y_discrete(breaks=0:11, labels=paste0("Treg_c", 0:11), limits = rev) +
    scale_x_discrete(limits = rev) +
    xlab("") + ylab("") +
    scale_fill_gradientn(
        colors = c("#5DBCFF", "#6DCCFF", "white", colorsForDataType[3], "#F484AE"))+
    scale_size(range = c(0, 3.5), limits = c(0, 100), breaks = c(0,20,40,60,80,100))+             ## to tune the size of circles
    ## scale_x_reverse() + scale_y_reverse() +
    theme(
        text = element_text(size = 10),
        panel.grid.major = element_line(colour = "grey90", size=0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
        legend.position="bottom",
        legend.title=element_text(size=8)) +
    guides(size = guide_legend(title.position="top",
                               title.hjust = 0.5,
                               override.aes = list(stroke = 0.4)),
           fill = guide_colourbar(title.position = "top", title.hjust = 0.5))
ggsave(file.path(figurePath, "Treg_bubbleplot.pdf"), plotx, width = 210/2.15,
       height = 67, units = "mm")


###############################################################################
#                                  Figure 3i                                  #
###############################################################################
Idents(TFH_Obj) <- TFH_Obj$seurat_clusters
umapColor <- colorRampPalette(colorsForDataType)(length(unique(TFH_Obj$seurat_clusters)))
umapColor <- c("#6DCCDD","#9ACBDE","#C8CADF","#EDC6DD","#F092B1","#F27FA5","#F47892", "#F6A395","#F8AD77","#E7B080","#C9AFA2","#ABADC5","#AEB9AC","#B9C984", "#C4DA5D")
png(file.path(figurePath, "TFH_UMAP_figure3D.png"), width = 210 / 2, height = 210 / 2, units = "mm", res = 600)
DimPlot(TFH_Obj, label = T) +
    scale_color_manual(values = colorRampPalette(umapColor)(length(unique(TFH_Obj$seurat_clusters)))) +
    theme_void() + theme(text = element_text(size = 6),
                         legend.position = "none")
dev.off()


###############################################################################
#                                  Figure 3k                                  #
###############################################################################
Idents(TFH_Obj) <- TFH_Obj$seurat_clusters
markers <- c(
    "NMB",
    "IL7R",
    "CXCL13",
    "CD200",
    "CCR7",
    "PDCD1",
    "TOX2",
    "SH2D1A",
    "IL6R",
    "CD40LG",
    "TIGIT",
    "CXCR5",
    "ASCL2",
    "ICOS",
    "ASCL2",
    "KLRB1",
    "GPR183",
    "CD69",
    "KLF2",
    "CCR5",
    "TGFB1",
    "STAT5",
    "PRF1",
    "PRDM1",
    "NKG7",
    "MAF",
    "LAG3",
    "KLRB1",
    "IL2RG",
    "IL2RB",
    "IL10RA",
    "IFNG",
    "HAVCR2",
    "GZMB",
    "GZMA",
    "GNLY",
    "ENTPD1",
    "IRF4",
    "BATF",
    "TNFRSF4",
    "IL2RA",
    "NRN1"
)
gene <- intersect(markers, rownames(TFH_Obj))
p <- DotPlot(TFH_Obj, features = rev(gene))
data <- p$data[ ,c('id','features.plot','pct.exp','avg.exp.scaled')]
## data$id <- factor(data$id, levels = rev(dataTypeLevel))
plotx <- ggplot(data, aes(y = id, x = features.plot)) +        ## global aes
    geom_point(aes(fill = avg.exp.scaled, size = pct.exp),
               color = 'black',
               shape = 21,
               stroke = 0.005)  +    ## geom_point for circle illusion
                                        #scale_fill_gradientn(colours=rev(color),limits=c(0,max(data$avg.exp)))+       ## color of the corresponding aes
    scale_y_discrete(breaks=0:11, labels=paste0("TFH_c", 0:11), limits = rev) +
    scale_x_discrete(limits = rev) +
    xlab("") + ylab("") +
    ## scale_fill_gradient2(high = colorsForDataType[3], mid = "white", low = "#6DCCFF")+
    scale_fill_gradientn(
        ## limits = c(-1.5, 1.5),
        colors = c("#5DBCFF", "#6DCCFF", "white", colorsForDataType[3], "#F484AE"))+
    scale_size(range = c(0, 3.5), limits = c(0, 100), breaks = c(0,20,40,60,80,100))+             ## to tune the size of circles
    theme(
        text = element_text(size = 10),
        panel.grid.major = element_line(colour = "grey90", size=0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
        legend.position="bottom",
        legend.title=element_text(size=8)) +
    guides(size = guide_legend(title.position="top",
                               title.hjust = 0.5,
                               override.aes = list(stroke = 0.4)),
           fill = guide_colourbar(title.position = "bottom", title.hjust = 0.5))
ggsave(file.path(figurePath, "TFH_bubbleplot.pdf"), plotx, width = 210/1.8,
       height = 297/5, units = "mm")


###############################################################################
#'                                NK MAIT GDT                                '#
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
library(scales)
#' load seurat object #########################################################
NK_path <- '/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/validate/NKGDT_V5/nPC_10/UMAP_dist_0.1_nneighbor_35/p1NKGDT_V4_10_UMAP_dist_0.1_nneighbor_35_CLUSTER_res_0.3/cluster.rds'
NK_Obj <- readRDS(NK_path)
figurePath <- "/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/figures2"
colorsForDataType <- c("#6DCCDD", "#EDCAE0", "#F494BE", "#F9B26C", "#A6ADCC", "#C4DA5D")
colorForClass1 <- c("#C4DA5D", "#6DCCDD", "#F494BE", "#EDCAE0")
dataTypeLevel <- c("CD4", "CD8", "NK", "NKT", "MAIT", "Proliferative")


umapColor <- colorRampPalette(colorsForDataType)(length(unique(NK_Obj$seurat_clusters)))
## umapColor <- c("#6DCCDD","#9ACBDE","#C8CADF","#EDC6DD","#F092B1","#F27FA5","#F47892", "#F6A395","#F8AD77","#E7B080","#C9AFA2","#ABADC5","#AEB9AC","#B9C984", "#C4DA5D")
png(file.path(figurePath, "figure4A.png"), width = 210 / 2, height = 210 / 2, units = "mm", res = 600)
DimPlot(NK_Obj, label = T) +
    ## scale_color_manual(values = hue_pal( h.start = 0, direction = 1, l = 80)(15)) +
    scale_color_manual(values = umapColor) +
    ## scale_x_reverse()+
    theme_void() + theme(text = element_text(size = 6),
                         legend.position = "none")
dev.off()


Idents(NK_Obj) <- NK_Obj$seurat_clusters
markers <- c(
"CD3D",
"CD3G",
"KLRF1",
"LTB",
"KLRB1",
"IL7R",
"CCR6",
"RORC",
"KLRC1",
"XCL1",
"XCL2",
"FCER1G",
"IL2RB",
"CD7",
"ITGB1",
"PRSS23",
"TRAV1-2",
"SLC4A10",
"IL23R",
"IL18R1",
"KLRG1",
"CXCR6",
"ZBTB16",
"ITGAM",
"KLRD1",
"FGFBP2",
"FCGR3A",
"SELL",
"NCAM1",
"TBX21",
"ZNF683",
"CX3CR1",
"CD8B",
"KLRF1",
"LAG3",
"KIR2DL1",
"KIR3DL1",
"KIR2DL3",
"KIR3DL2",
"TIGIT",
"EOMES"
)
gene <- intersect(markers, rownames(NK_Obj))
p <- DotPlot(NK_Obj, features = rev(gene))
data <- p$data[ ,c('id','features.plot','pct.exp','avg.exp.scaled')]
## data$id <- factor(data$id, levels = rev(dataTypeLevel))
plotx <- ggplot(data, aes(x = id, y = features.plot)) +        ## global aes
    geom_point(aes(fill = avg.exp.scaled, size = pct.exp),
               color = 'black',
               shape = 21,
               stroke = 0.005)  +    ## geom_point for circle illusion
                                        #scale_fill_gradientn(colours=rev(color),limits=c(0,max(data$avg.exp)))+       ## color of the corresponding aes
    scale_x_discrete(breaks=0:11, labels=paste0("NK_c", 0:11)) +
    xlab("") + ylab("") +
    ## scale_fill_gradient2(high = colorsForDataType[3], mid = "white", low = "#6DCCFF")+
    scale_fill_gradientn(
        ## limits = c(-1.5, 1.5),
        colors = c("#5DBCFF", "#6DCCFF", "white", colorsForDataType[3], "#F484AE"))+
    scale_size(range = c(0, 3.5), limits = c(0, 100), breaks = c(0,20,40,60,80,100))+             ## to tune the size of circles
    theme(
        text = element_text(size = 10),
        panel.grid.major = element_line(colour = "grey90", size=0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 1),
        ## legend.position="bottom",
        legend.title=element_text(size=8)) +
    guides(size = guide_legend(title.position="top",
                               title.hjust = 0.5,
                               override.aes = list(stroke = 0.4)),
           fill = guide_colourbar(title.position = "top", title.hjust = 0.5))
ggsave(file.path(figurePath, "NK_bubbleplot.pdf"), plotx, width = 210/2.5,
       height = 297/2, units = "mm")


Idents(NK_Obj) <- NK_Obj$seurat_clusters
markers <- c(
    "TNF",
    "IFNG",
    "GZMH",
    "GZMK",
    "GNLY",
    "PRF1",
    "GZMB",
    "GZMA",
    "GZMM"
)
gene <- intersect(markers, rownames(NK_Obj))
p <- DotPlot(NK_Obj, features = rev(gene))
data <- p$data[ ,c('id','features.plot','pct.exp','avg.exp.scaled')]
## data$id <- factor(data$id, levels = rev(dataTypeLevel))
plotx <- ggplot(data, aes(x = id, y = features.plot)) +        ## global aes
    geom_point(aes(fill = avg.exp.scaled, size = pct.exp),
               color = 'black',
               shape = 21,
               stroke = 0.005)  +    ## geom_point for circle illusion
                                        #scale_fill_gradientn(colours=rev(color),limits=c(0,max(data$avg.exp)))+       ## color of the corresponding aes
    scale_x_discrete(breaks=0:11, labels=paste0("NK_c", 0:11)) +
    xlab("") + ylab("") +
    ## scale_fill_gradient2(high = colorsForDataType[3], mid = "white", low = "#6DCCFF")+
    scale_fill_gradientn(
        ## limits = c(-1.5, 1.5),
        colors = c("#5DBCFF", "#6DCCFF", "white", colorsForDataType[3], "#F484AE"))+
    scale_size(range = c(0, 3.5), limits = c(0, 100), breaks = c(0,20,40,60,80,100))+             ## to tune the size of circles
    theme(
        text = element_text(size = 10),
        panel.grid.major = element_line(colour = "grey90", size=0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 1),
        ## legend.position="bottom",
        legend.title=element_text(size=8)) +
    guides(size = guide_legend(title.position="top",
                               title.hjust = 0.5,
                               override.aes = list(stroke = 0.4)),
           fill = guide_colourbar(title.position = "top", title.hjust = 0.5))
ggsave(file.path(figurePath, "NK_cytotoxicity_bubbleplot.pdf"), plotx, width = 210/2.5,
       height = 297/4, units = "mm")




###############################################################################
#'                                  figure 5                                 '#
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
CD8_path <- '/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/validate/CD8_V6/nPC_50/UMAP_dist_0.1_nneighbor_50/p1CD8_V6_UMAP_dist_0.1_nneighbor_50_CLUSTER_res_0.3/cluster.rds'
CD4_path <- '/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/validate/CD4_V7/nPC_50/UMAP_dist_0.1_nneighbor_50/p1CD4_V7_UMAP_dist_0.1_nneighbor_50_CLUSTER_res_0.3/cluster.rds'
NK_path <- '/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/validate/NKGDT_V5/nPC_10/UMAP_dist_0.1_nneighbor_35/p1NKGDT_V4_10_UMAP_dist_0.1_nneighbor_35_CLUSTER_res_0.3/cluster.rds'
Proliferative_path <- '/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/validate/Proliferating_V3/nPC_30/UMAP_dist_0.1_nneighbor_50/p1Proliferating_V3_UMAP_dist_0.1_nneighbor_50_CLUSTER_res_0.8/cluster.rds'
CD8_Obj <- readRDS(CD8_path)
CD4_Obj <- readRDS(CD4_path)
NK_Obj <- readRDS(NK_path)
Proliferative_Obj <- readRDS(Proliferative_path)
CD4_Obj@meta.data$DataSet <- "CD4"
CD8_Obj@meta.data$DataSet <- "CD8"
NK_Obj@meta.data$DataSet <- "NK"
NK_Obj@meta.data$DataSet[NK_Obj$seurat_clusters %in% c(0, 2, 5)] <- "NK"
NK_Obj@meta.data$DataSet[NK_Obj$seurat_clusters %in% c(4)] <- "MAIT"
NK_Obj@meta.data$DataSet[NK_Obj$seurat_clusters %in% c(1, 6)] <- "MAIT-like"
NK_Obj@meta.data$DataSet[NK_Obj$seurat_clusters %in% c(3, 7, 8)] <- "NKT"
Proliferative_Obj@meta.data$DataSet <- "Proliferative"
allPath <- list(CD8 = CD8_path, CD4 = CD4_path, NK = NK_path, Proliferative = Proliferative_path)
allObj <- list(CD8 = CD8_Obj, CD4 = CD4_Obj, NK = NK_Obj, Proliferative = Proliferative_Obj)
common_columns <- c()
for(oi in 1:length(allObj)) {
    tempName <- names(allObj[oi])
    tempObj <- allObj[[oi]]
    if(oi == 1) {
        common_columns <- colnames(tempObj@meta.data)
    }else{
        common_columns <- intersect(common_columns, colnames(tempObj@meta.data))
    }
}
allMD <- c()
for(oi in 1:length(allObj)) {
    tempName <- names(allObj[oi])
    tempObj <- allObj[[oi]]
    tempMD <- as_tibble(tempObj@meta.data, rownames = NA) %>% select(common_columns) %>% mutate(DataSetType = tempName)
    allMD <- bind_rows(allMD, tempMD)
}
isHealthyDonor <- allMD$Disease == "Healthy donor"
isPrimaryTumor <- allMD$TissueType == "Primary tumor tissue"
isUninvolvedNormal <- allMD$`Tissue-Type-1` == "Uninvolved normal tissue"
isMetastaticTumor <- allMD$Disease %in% c("HNSC_Metastatic", "STAD_Metastatic", "LUAD_Metastatic", "Melanoma_Metastatic", "OV_Metastatic")
allMD$Class_Level_1 <- ""
allMD$Class_Level_1[isHealthyDonor] <- "Healthy donor"
allMD$Class_Level_1[isPrimaryTumor] <- "Primary tumor tissue"
allMD$Class_Level_1[isUninvolvedNormal] <- "Uninvolved normal tissue"
allMD$Class_Level_1[isMetastaticTumor] <- "Metastatic tumor tissue"
allMD$Class_Level_2 <- ""
allMD$Class_Level_2[allMD$Class_Level_1 == "Healthy donor" &
                    allMD$`Tissue-type-2` == "BM_Healthy normal"] <- "BM"
allMD$Class_Level_2[allMD$Class_Level_1 == "Healthy donor" &
                    allMD$`Tissue-type-2` == "LN_Healthy"] <- "LN"
allMD$Class_Level_2[allMD$Class_Level_1 == "Healthy donor" &
                    allMD$`Tissue-type-2` == "PBMC_Healthy normal"] <- "PBMC"
allMD$Class_Level_2[allMD$Class_Level_1 == "Uninvolved normal tissue"] <-  allMD$OrganSite[allMD$Class_Level_1 == "Uninvolved normal tissue"]
allMD$Class_Level_2[allMD$Class_Level_1 == "Primary tumor tissue" &
                    allMD$batch == "AML.T.rds"] <- "AML"
allMD$Class_Level_2[allMD$Class_Level_1 == "Primary tumor tissue" &
                    allMD$batch == "BCC_Yost.rds"] <- "BCC"
allMD$Class_Level_2[allMD$Class_Level_1 == "Primary tumor tissue" &
                    allMD$batch %in% c("Breast_MDACC.rds", "Breast_GSE114725.rds")] <- "BRCA"
allMD$Class_Level_2[allMD$Class_Level_1 == "Primary tumor tissue" &
                    allMD$batch == "Gastric_JQin"] <- "STAD"
allMD$Class_Level_2[allMD$Class_Level_1 == "Primary tumor tissue" &
                    allMD$batch == "Glioma.Tcell.scs.rds"] <- "GBM"
allMD$Class_Level_2[allMD$Class_Level_1 == "Primary tumor tissue" &
                    allMD$batch %in% c("HeadNeck_Puram.rds", "HN.Tcell.scs.rds")] <- "HNSC"
allMD$Class_Level_2[allMD$Class_Level_1 == "Primary tumor tissue" &
                    allMD$batch %in% c("Liver_GSE125449.rds")] <- "HCC"
allMD$Class_Level_2[allMD$Class_Level_1 == "Primary tumor tissue" &
                    allMD$Disease %in% c("LUAD_Primary")] <- "LUAD"
allMD$Class_Level_2[allMD$Class_Level_1 == "Primary tumor tissue" &
                    allMD$Disease %in% c("LUSC_Primary")] <- "LUSC"
allMD$Class_Level_2[allMD$Class_Level_1 == "Primary tumor tissue" &
                    allMD$Disease %in% c("SCLC_Primary")] <- "SCLC"
allMD$Class_Level_2[allMD$Class_Level_1 == "Primary tumor tissue" &
                    allMD$Disease %in% c("NSCLC_Primary")] <- "NSCLC"
allMD$Class_Level_2[allMD$Class_Level_1 == "Primary tumor tissue" &
                    allMD$Disease %in% c("DHL", "DLBCL", "THL")] <- "LBCL"
allMD$Class_Level_2[allMD$Class_Level_1 == "Primary tumor tissue" &
                    allMD$Disease %in% c("FL")] <- "FL"
allMD$Class_Level_2[allMD$Class_Level_1 == "Primary tumor tissue" &
                    allMD$batch %in% c("MCL.T.Shaojun.rds", "MCL.T.David.rds")] <- "MCL"
allMD$Class_Level_2[allMD$Class_Level_1 == "Primary tumor tissue" &
                    allMD$Disease %in% c("Melanoma_Primary")] <- "SKCM"
allMD$Class_Level_2[allMD$Class_Level_1 == "Primary tumor tissue" &
                    allMD$Disease %in% c("PDAC_Primary")] <- "PDAC"
allMD$Class_Level_2[allMD$Class_Level_1 == "Metastatic tumor tissue" &
                    allMD$Disease %in% c("STAD_Metastatic")] <- "STAD"
allMD$Class_Level_2[allMD$Class_Level_1 == "Metastatic tumor tissue" &
                    allMD$Disease %in% c("HNSC_Metastatic")] <- "HNSC"
allMD$Class_Level_2[allMD$Class_Level_1 == "Metastatic tumor tissue" &
                    allMD$Disease %in% c("LUAD_Metastatic")] <- "LUAD"
allMD$Class_Level_2[allMD$Class_Level_1 == "Metastatic tumor tissue" &
                    allMD$Disease %in% c("Melanoma_Metastatic")] <- "SKCM"
allMD$Class_Level_2[allMD$Class_Level_1 == "Metastatic tumor tissue" &
                    allMD$Disease %in% c("OV_Metastatic")] <- "OV"
colorsForDataType <- c("#6DCCDD", "#EDCAE0", "#F494BE", "#F9B26C", "#A6ADCC", "#C4DA5D")
colorForClass1 <- c("#C4DA5D", "#6DCCDD", "#F494BE", "#EDCAE0")
dataTypeLevel <- c("CD4", "CD8", "NK", "NKT", "MAIT", "Proliferative")
figurePath <- "/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/figures4"
if(!dir.exists(figurePath)) {
    dir.create(figurePath)
}


for(oi in 1:length(allObj)) {
    tempName <- names(allObj[oi])
    tempObj <- allObj[[oi]]
    md <- as_tibble(tempObj@meta.data, rownames = NA)
    N <- dim(md)[1]
    for(typei in c("Disease", "OrganSite", "TissueType")) {
        ctmd <- md %>% group_by_("seurat_clusters", typei) %>% dplyr::count()
        tmd <- md %>% group_by_(typei) %>% dplyr::count()
        cmd <- md %>% group_by(seurat_clusters) %>% dplyr::count()
        matRoe <- matrix(rep(0, length(unique(md$seurat_clusters)) *
                                length(unique(md %>% pull(typei)))),
                         nrow = length(unique(md$seurat_clusters)),
                         ncol = length(unique(md %>% pull(typei))))
        for(ci in 1:length(unique(md$seurat_clusters))) {
            cin <- sort(unique(md$seurat_clusters))[ci]
            k <- cmd$n[cmd$seurat_clusters == cin]
            for(ti in 1:length(unique(md %>% pull(typei)))) {
                tin <- sort(unique(md %>% pull(typei)))[ti]
                if(any(ctmd$seurat_clusters == cin & ctmd[,typei] == tin)) {
                    n <- ctmd$n[ctmd$seurat_clusters == cin & ctmd[,typei] == tin]
                }else{
                    n <- 0
                }
                M <- tmd$n[tmd[,typei] == tin]
                matRoe[ci, ti] <- (n/M) / (k/N)
            }
        }
        rownames(matRoe) <- paste0(tempName, "_", sort(unique(md$seurat_clusters)))
        colnames(matRoe) <- sort(unique(md %>% pull(typei)))
        ## my.breaks <- c(seq(0, 5, by=0.01), seq(5.0001, 30, by=0.3))
        ## my.colors <- c(colorRampPalette(colors = c("blue", "white", "red"))(length(my.breaks)/2),
        ##                colorRampPalette(colors = c("red", "darkred"))(length(my.breaks)/2))
        matRoe[matRoe > 3] <- 3
        my.breaks <- seq(0, 3, by=0.02)
        my.colors <- c(
            colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)/2),
            colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks)/2))
        widthi = ncol(matRoe)/6 + 4.3
        heighti = nrow(matRoe)/10 + 3
        pdf(file.path(figurePath, paste0(typei, "_", tempName, "_roe.pdf")), width = widthi, height = heighti)
        print(pheatmap(matRoe,
                       color = my.colors,
                       breaks = my.breaks,
                       ## annotation_row = cellType_col,
                       show_colnames = T,
                       show_rownames = T,
                       cluster_cols = T,
                       cluster_rows = T,
                       border_color = F))
        dev.off()
    }
}

#' mark each cluster name #####################################################
allMD$g_seurat_clusters[allMD$DataSet != "Proliferative"] <- paste0(allMD$DataSet[allMD$DataSet != "Proliferative"], "_c", allMD$seurat_clusters[allMD$DataSet != "Proliferative"])
allMD$g_seurat_clusters[allMD$DataSet == "Proliferative"] <- allMD$DataSet[allMD$DataSet == "Proliferative"]
ClusterT <- allMD %>% group_by(g_seurat_clusters) %>% dplyr::count()
N <- dim(allMD)[1]
CL1OrganSiteT <- allMD %>% group_by(Class_Level_1, OrganSite) %>% dplyr::count()
CL1OrganSiteT$g_organsite <- paste0(CL1OrganSiteT$Class_Level_1, "_", CL1OrganSiteT$OrganSite)
matRoe <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) *
                        dim(CL1OrganSiteT)[1]),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = dim(CL1OrganSiteT)[1])
rownames(matRoe) <- sort(unique(allMD$g_seurat_clusters))
colnames(matRoe) <- sort(unique(CL1OrganSiteT$g_organsite))
for(ci in unique(allMD$Class_Level_1)) {
    OrganSiteT <- allMD %>% filter(Class_Level_1 == ci) %>% group_by(OrganSite) %>% dplyr::count()
    for(osi in unique(OrganSiteT$OrganSite)) {
        OrganSiteClusterT <- allMD %>% filter(Class_Level_1 == ci) %>% group_by(OrganSite, g_seurat_clusters) %>% dplyr::count()
        for(gsci in unique(allMD$g_seurat_clusters)) {
            k <- ClusterT$n[ClusterT$g_seurat_clusters == gsci]
            if(gsci %in% OrganSiteClusterT$g_seurat_clusters[OrganSiteClusterT$OrganSite == osi]) {
                n <- OrganSiteClusterT$n[OrganSiteClusterT$OrganSite == osi & OrganSiteClusterT$g_seurat_clusters == gsci]
            } else{
                n <- 0
            }
            M <- OrganSiteT$n[OrganSiteT$OrganSite == osi]
            tempRowName <- gsci
            tempColName <- paste0(ci, "_", osi)
            matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
        }
    }
}
matRoe[matRoe > 2] <- 2
my.breaks <- seq(0, 2, by=0.01)
my.colors <- c(
    colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)/2),
    colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks)/2))
widthi = ncol(matRoe)/11 + 4.3
heighti = nrow(matRoe)/6 + 3
pdf(file.path(figurePath, paste0("all_clusterVStissue_roe.pdf")), width = widthi, height = heighti)
print(pheatmap(matRoe,
               color = my.colors,
               breaks = my.breaks,
               ## annotation_row = cellType_col,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = T,
               cluster_rows = T,
               border_color = F))
dev.off()


ggplotdata <- allMD %>% group_by(Class_Level_1, OrganSite, Class_Level_2) %>% dplyr::count() %>% group_by(Class_Level_1, OrganSite) %>% arrange(desc(Class_Level_2)) %>% mutate(frac = n / sum(n), pos = cumsum(frac) - frac / 2)
g <- ggplot(ggplotdata) +
    geom_bar(aes(x = OrganSite, y = frac, fill = Class_Level_2), stat = "identity", color = "white") +
    geom_text(aes(OrganSite, pos, group=Class_Level_2, label=Class_Level_2),
              size=2, hjust=0.5,vjust=0.5,angle=90) +
    facet_grid(. ~ Class_Level_1, scales = "free") +
    theme_classic() +
    theme(text = element_text(size = 8),
          legend.title=element_text(size=8),
          legend.position = "bottom",
          axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(file.path(figurePath, "bar.pdf"), g, width = 210, height= 150, units = "mm")




#' mark each cluster name #####################################################
allMD$g_seurat_clusters[allMD$DataSet != "Proliferative"] <- paste0(allMD$DataSet[allMD$DataSet != "Proliferative"], "_c", allMD$seurat_clusters[allMD$DataSet != "Proliferative"])
allMD$g_seurat_clusters[allMD$DataSet == "Proliferative"] <- allMD$DataSet[allMD$DataSet == "Proliferative"]
ClusterT <- allMD %>% group_by(g_seurat_clusters) %>% dplyr::count()
N <- dim(allMD)[1]
matRoe <- matrix(rep(0, length(unique(allMD$g_seurat_clusters)) *
                        length(unique(allMD$Class_Level_2))),
                 nrow = length(unique(allMD$g_seurat_clusters)),
                 ncol = length(unique(allMD$Class_Level_2)))
rownames(matRoe) <- sort(unique(allMD$g_seurat_clusters))
colnames(matRoe) <- sort(unique(allMD$Class_Level_2))
for(ci in unique(allMD$Class_Level_2)) {
    CancerT <- allMD %>% filter(Class_Level_2 == ci) %>% dplyr::count()
    for(gsci in unique(allMD$g_seurat_clusters)) {
        k <- ClusterT$n[ClusterT$g_seurat_clusters == gsci]
        if(gsci %in% allMD$g_seurat_clusters[allMD$Class_Level_2 == ci]) {
            CancerClusterT <- allMD %>% filter(Class_Level_2 == ci, g_seurat_clusters == gsci) %>% dplyr::count()
            n <- CancerClusterT$n
        }else{
            n <- 0
        }
        M <- CancerT$n
        tempRowName <- gsci
        tempColName <- ci
        matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
    }
}
matRoe[matRoe > 2] <- 2
my.breaks <- seq(0, 2, by=0.01)
my.colors <- c(
    colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)/2),
    colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks)/2))
widthi = ncol(matRoe)/11 + 4.3
heighti = nrow(matRoe)/6 + 3
pdf(file.path(figurePath, paste0("all_clusterVScancer_roe.pdf")), width = widthi, height = heighti)
print(pheatmap(matRoe,
               color = my.colors,
               breaks = my.breaks,
               ## annotation_row = cellType_col,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = T,
               cluster_rows = T,
               border_color = F))
dev.off()







#' draw boxplot for clusters in Cancer type ################################################
allMD$g_seurat_clusters[allMD$DataSet != "Proliferative"] <- paste0(allMD$DataSet[allMD$DataSet != "Proliferative"], "_c", allMD$seurat_clusters[allMD$DataSet != "Proliferative"])
allMD$g_seurat_clusters[allMD$DataSet == "Proliferative"] <- allMD$DataSet[allMD$DataSet == "Proliferative"]
ggplotdata <- c()
for(ci in unique(allMD$Class_Level_1)) {
    ggplotdataTemp <- allMD %>% filter(Class_Level_1 == ci) %>% group_by(Class_Level_2, g_seurat_clusters, OrganSite) %>% dplyr::count() %>% group_by(g_seurat_clusters) %>% mutate(frac = n / sum(n)) %>% select(Class_Level_2, g_seurat_clusters, OrganSite, frac) %>% mutate(Type = ci)
    ggplotdata <- bind_rows(ggplotdata, ggplotdataTemp)
}
ggplotdata$g_tissue <- paste0(ggplotdata$Type, "_", ggplotdata$Class_Level_2)
g <- ggplot(ggplotdata) +
    geom_boxplot(aes(x = g_tissue, y = frac, color = g_tissue)) +
    facet_grid(g_seurat_clusters ~ ., scales = "free") +
    theme_classic()+
    theme(text = element_text(size = 8),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
          )
ggsave(file.path(figurePath, "all_clusterVScancer_BoxPlot.pdf"), g, width = 210, height = 497, units = "mm")

ggplotdata <- allMD %>% group_by(Class_Level_2, g_seurat_clusters, OrganSite) %>% dplyr::count() %>% group_by(g_seurat_clusters) %>% mutate(frac = n / sum(n)) %>% select(Class_Level_2, g_seurat_clusters, OrganSite, frac)
g <- ggplot(ggplotdata) +
    geom_boxplot(aes(x = OrganSite, y = frac, color = OrganSite)) +
    facet_grid(g_seurat_clusters ~ ., scales = "free") +
    theme_classic()+
    theme(text = element_text(size = 8),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file.path(figurePath, "all_clusterVScancer_BoxPlot_onlyOrganSite.pdf"), g, width = 210, height = 497, units = "mm")
