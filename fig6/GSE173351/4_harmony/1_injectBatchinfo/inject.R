library(optparse)
library(tidyverse)
library(Seurat)
library(GEOquery)


## gse <- getGEO(filename = "/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/data/GSE173351/GSE176021_series_matrix.txt.gz",
##               GSEMatrix = TRUE,
##               getGPL = FALSE)
## metaInfo <- gse@phenoData@data %>%
##     select(title, `response status:ch1`, source_name_ch1)
## write_tsv(metaInfo, file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/data/GSE173351", paste0('metaInfo', "_", Sys.Date(), '.tsv')))


metaInfo <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/data/GSE173351/metaInfo_2022-05-13.tsv")

seuratObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE173351/3_pca/pca.rds")

seuratObj$response <- metaInfo$response[match(seuratObj$orig.ident, metaInfo$orig.ident)]
seuratObj$tissue <- metaInfo$tissue[match(seuratObj$orig.ident, metaInfo$orig.ident)]
seuratObj$batch <- metaInfo$batch[match(seuratObj$orig.ident, metaInfo$orig.ident)]

outDir <- "/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE173351/4_harmony"
if(!dir.exists(outDir)){
    dir.create(outDir, recursive = T)
}
saveRDS(seuratObj, file.path(outDir, paste0('harmony_input.rds')))



