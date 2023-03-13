#--------------------------------------------------------------
# filename : merge.R
# Date : 2022-11-05
# contributor : Yanshuo Chu
# function: merge
# R version: R/4.0.3
#--------------------------------------------------------------

print('<==== merge.R ====>')
rm(list=ls())

tenx.data = Read10X("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/data/SCP1288/expression/60c76a18771a5b0ba10ea91b")
seurat_obj = CreateSeuratObject(counts = tenx.data,
                           min.cells = 3,
                           min.features = 200,
                           project = "SCP1288")

metaInfo <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/data/SCP1288/metadata/Final_SCP_Metadata.txt")
metaInfo <- metaInfo[2:dim(metaInfo)[1],]
clusterInfo <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/data/SCP1288/cluster/Final_SCP_ClusterFile.txt")
clusterInfo <- clusterInfo[2:dim(clusterInfo)[1],]


seurat_obj$cell.type <- clusterInfo$FinalCellType[match(Cells(seurat_obj), clusterInfo$NAME)]

for (cln in colnames(metaInfo)[2:length(colnames(metaInfo))]) {
  seurat_obj@meta.data[,cln] <- ""
  seurat_obj@meta.data[,cln] <- metaInfo[match(Cells(seurat_obj), metaInfo$NAME), cln]
}

CD8_clusters <- c(
  "41BB-Hi CD8+ T cell",
  "41BB-Lo CD8+ T cell",
  "Cycling CD8+ T cell",
  "MitoHigh CD8+ T cell",
  "MX1-Hi CD8+ T cell")


CD4_clusters <- c(
  "Effector T-Helper",
  "Memory T-Helper",
  "MitoHigh T-Helper",
  "T-Reg" )

figure_path <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/SCP1288/0_merge/")
if (!dir.exists(figure_path)) {
  dir.create(figure_path, recursive = T)
}
setwd(figure_path)

Idents(seurat_obj) <- seurat_obj$FinalCellType
CD8_obj <- subset(seurat_obj, idents = CD8_clusters)
saveRDS(CD8_obj, file.path(getwd(), paste0('CD8.rds')))

CD4_obj <- subset(seurat_obj, idents = CD4_clusters)
saveRDS(CD4_obj, file.path(getwd(), paste0('CD4.rds')))
