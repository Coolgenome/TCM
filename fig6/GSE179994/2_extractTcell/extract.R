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

seuratObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE179994/merged/merged.obj")

Idents(seuratObj) <- seuratObj$celltype
CD4SeuratObj <- subset(seuratObj, idents = "CD4")
saveRDS(CD4SeuratObj, file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE179994/", paste0('CD4SeuratObj', "_", Sys.Date(), '.rds')))
CD8SeuratObj <- subset(seuratObj, idents = "CD8")
saveRDS(CD8SeuratObj, file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE179994/", paste0('CD8SeuratObj', "_", Sys.Date(), '.rds')))

## md <- as_tibble(seuratObj@meta.data)
## md %>%
##   group_by(celltype, cluster) %>%
##   count %>%
##   as.data.frame


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

figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE179994/2_extractTcell/outs")
if(!dir.exists(figurePath)){
    dir.create(figurePath, recursive = T)
}
setwd(figurePath)

ggsave(file.path(paste0("response_bar.pdf")), g, width = 200, height = 1200, units = "mm")
