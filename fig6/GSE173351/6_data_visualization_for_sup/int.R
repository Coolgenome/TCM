#--------------------------------------------------------------
# filename : int.R
# Date : 2022-10-27
# contributor : Yanshuo Chu
# function: int
# R version: R/4.0.3
#--------------------------------------------------------------

print('<==== int.R ====>')
rm(list=ls())

library(optparse)
library(tidyverse)
library(Seurat)

figure_path <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE173351/6_data_visualization_for_sup/")
if (!dir.exists(figure_path)) {
  dir.create(figure_path, recursive = T)
}
setwd(figure_path)

seurat_obj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE173351/4_harmony/nPC_30/UMAP_dist_0.1_nneighbor_50/GSE173351_UMAP_dist_0.1_nneighbor_50_CLUSTER_res_0.3/cluster.rds")

# step1 remove bad clusters  ##################################################
## bad_clusters <- 14:19
## seurat_obj <- subset(seurat_obj, idents = bad_clusters, invert = T)
## saveRDS(seurat_obj, file.path(getwd(), paste0('clusters_remove_bad_clusters.rds')))

# step2 draw umap #############################################################

colorsForDataType <- c("#6DCCDD", "#EDCAE0", "#F494BE", "#F9B26C", "#A6ADCC", "#C4DA5D")

umapColor <- c("#6DCCDD","#9ACBDE","#C8CADF","#EDC6DD","#F092B1","#F27FA5","#F47892", "#F6A395","#F8AD77","#E7B080","#C9AFA2","#ABADC5","#AEB9AC","#B9C984", "#C4DA5D")
png(file.path(getwd(), "UMAP.png"), width = 210 / 2, height = 210 / 2, units = "mm", res = 600)
DimPlot(seurat_obj, label = T) +
  scale_color_manual(values = colorRampPalette(umapColor)(length(unique(seurat_obj$seurat_clusters)))) +
  theme_void() + theme(text = element_text(size = 6),
                       legend.position = "none")
dev.off()


markers <- c("CD3D",
             "CD4",
             "CD40LG",
             "CD8A",
             "CD8B",
             "GZMK",
             "SLC4A10",
             "TRAV1-2",
             "GZMB",
             "GZMH",
             "FCGR3A",
             "KLRD1",
             "NCAM1",
             "TRDV2",
             "TRGV9",
             "TRGV10",
             "MKI67",
             "TRGV2",
             "CD19",
             "MS4A1",
             "NKG7",
             "S100A8",
             "SELL",
             "TNFRSF4",
             "SESN3",
             "CCR7",
             "MAGEH1",
             "CCR8",
             "BATF")

gene <- intersect(markers, rownames(seurat_obj))
p <- DotPlot(seurat_obj, features = markers)
data <- p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]
data$id <- factor(data$id, levels = rev(cell.type.order))
plotx <- ggplot(data, aes(y = id, x = features.plot)) +
  geom_point(aes(fill = avg.exp.scaled, size = pct.exp),
             color = 'black',
             shape = 21,
             stroke = 0.01) +
  xlab("") + ylab("") +
  scale_fill_gradientn(colors = c("white", "#9370DB", "#000000")) +
  scale_size(range = c(0, 3.2), limits = c(0, 100),
             breaks = c(0, 20, 40, 60, 80, 100)) +
  theme(
    text = element_text(size = 10),
    panel.grid.major = element_line(colour = "grey90", size=0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
    ## legend.position="bottom",
    legend.title=element_text(size=10))
  ## guides(size = guide_legend(title.position="top",
  ##                            title.hjust = 0.5,
  ##                            byrow = T,
  ##                            override.aes = list(stroke = 0.4)),
         ## fill = guide_colourbar(title.position = "top", title.hjust = 0.5))

ggsave(file.path(figure_path, "bubbleplot.pdf"),
       plotx,
       width = 210,
       height = 297/4,
       units = "mm")

## 10 Proliferative T

