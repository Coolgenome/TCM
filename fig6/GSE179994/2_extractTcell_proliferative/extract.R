#--------------------------------------------------------------
# filename : extract.R
# Date : 2022-02-16
# contributor : Yanshuo Chu
# function: extract
#--------------------------------------------------------------

print('<==== extract.R ====>')
rm(list=ls())

library(data.table)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(harmony)
library(ggstatsplot)

figure_path <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE179994/2_extractTcell_proliferative/")
if (!dir.exists(figure_path)) {
  dir.create(figure_path, recursive = T)
}
setwd(figure_path)

seuratObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE179994/merged/merged.obj")

md <- as_tibble(seuratObj@meta.data)
md %>%
  group_by(celltype, cluster) %>%
  count %>%
  as.data.frame

CD4_clusters <- c("CD4_C1-Naive", "CD4_C2-Tcm", "CD4_C3-Tem", "CD4_C4-CD69", "CD4_C5-ISG15", "CD4_C6-RPL", "CD4_C7-Th1-like", "CD4_C8-Treg", "XCL1")
CD4_prolif_clusters <- c("CD4_C9-Prolif.")
CD8_clusters <- c("Non-exhausted", "Tex")
CD8_prolif_clusters <- c("Prolif.")

Idents(seuratObj) <- seuratObj$cluster

CD4SeuratObj <- subset(seuratObj, idents = CD4_clusters)
saveRDS(CD4SeuratObj, paste0('CD4SeuratObj', "_", Sys.Date(), '.rds'))
CD4ProlifSeuratObj <- subset(seuratObj, idents = CD4_prolif_clusters)
saveRDS(CD4ProlifSeuratObj, paste0('CD4ProlifSeuratObj', "_", Sys.Date(), '.rds'))

CD8SeuratObj <- subset(seuratObj, idents = CD8_clusters)
saveRDS(CD8SeuratObj, paste0('CD8SeuratObj', "_", Sys.Date(), '.rds'))
CD8ProlifSeuratObj <- subset(seuratObj, idents = CD8_prolif_clusters)
saveRDS(CD8ProlifSeuratObj, paste0('CD8ProlifSeuratObj', "_", Sys.Date(), '.rds'))



clinicT <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/data/GSE179994/ClinicData.txt") %>%
    filter(`Treatment Hx` %in% c("On treatment"))
noResponseSamples <- clinicT %>%
    filter(Response == "No") %>%
    pull(`Sample Name`)
ResponseSamples <- clinicT %>%
    filter(Response == "Yes") %>%
    pull(`Sample Name`)

md <- seuratObj@meta.data
TotalSampleCellNum <- md %>%
    group_by(sample) %>%
    count()

totalT <- c()
for(tempCluster in unique(md$cluster)){
    TNR <- md %>%
        filter(cluster == tempCluster) %>%
        group_by(sample) %>%
        count() %>%
        filter(sample %in% c(noResponseSamples, ResponseSamples))

    TNR$Frac <- 0.0
    TNR$Frac <- TNR$n / TotalSampleCellNum$n[match(TNR$sample, TotalSampleCellNum$sample)]

    TNR$isResponse <- "NO"
    TNR$isResponse[TNR$sample %in% ResponseSamples] <- "YES"

    TNR$cluster <- tempCluster

    totalT <- bind_rows(totalT, TNR)
}


g <- totalT %>%
    ggstatsplot::grouped_ggbetweenstats(
                     data = .,
                     x = isResponse,
                     y = Frac,
                     grouping.var = cluster,
                     xlab = "",
                     ylab = "Sample fraction",
                     ## pairwise.display = "aiwl", # display only significant pairwise comparisons
                     p.adjust.method = "fdr", # adjust p-values for multiple tests using this method
                     ggtheme = theme_classic(),
                     package = "ggsci",
                     palette = "default_jco",
                     plotgrid.args = list(ncol = 1))


ggsave(file.path(paste0("response_bar.pdf")), g, width = 200, height = 1200, units = "mm")
