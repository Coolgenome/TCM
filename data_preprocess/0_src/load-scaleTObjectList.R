#' filename : load-scaleTObjectList.R
#' Date : 2020-07-07
#' contributor : Yanshuo Chu
#' function: load-scaleTObjectList


suppressMessages({library(optparse)
library(tidyverse)
library(Seurat)
library(future)
library(future.apply)
})
print('----loading cellranger outputs----')
##CLI parsing
option_list = list(
    make_option(c("-d","--data"),
                type = 'character',
                help = 'dataFolder',
                metavar = 'character'),
    make_option(c("-o",'--out'),
                type = 'character',
                default = 'cellranger.rds',
                help = 'result file name [default = %default]',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

options(future.globals.maxSize = 36 * 1000 * 1024^2)

print(paste0("Loading data from ", opt$data))
samples = list.files(path = opt$data,recursive = F)
samples = basename(samples)
if(length(samples) < 1) stop(paste0('Failed to find data in ',opt$data))

obj.list = list()
for (ids in seq_along(samples)) {
    ds = samples[ids]
    print(ds)
    seuratObj <- readRDS(file.path(opt$data, ds))

    print(colnames(seuratObj@meta.data))
    if("cell.types" %in% colnames(seuratObj@meta.data)){
        seuratObj@meta.data <- seuratObj@meta.data[,c("orig.ident", "nCount_RNA", "nFeature_RNA", "cell.types")]
    }else{
        seuratObj@meta.data <- seuratObj@meta.data[,c("orig.ident", "nCount_RNA", "nFeature_RNA")]
    }
    seuratObj@meta.data$batch <- ds
    if(str_detect(ds, "^Gastric_JQin*")){
        seuratObj@meta.data$batch <- "Gastric_JQin"
    }
    if(str_detect(ds, "^Lung_Le*")){
        seuratObj@meta.data$batch <- "Lung_Le"
    }

    ##Human: MT; Mouse: mt
    mito.features = grep(pattern = '^MT-|^mt-', x = rownames(x = seuratObj), value = T)
    percent.mito <- Matrix::colSums(x = GetAssayData(object = seuratObj, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seuratObj, slot = 'counts'))
    seuratObj[['percent.mito']] = percent.mito

    seuratObj <- subset(x = seuratObj,
                        subset = nFeature_RNA > 200 &
                            nFeature_RNA < 8000 &
                            percent.mito < 0.15)

    obj.list[[length(obj.list)+1]] = seuratObj
}


plan("multiprocess", workers = 20)
obj.list <- lapply(X = obj.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = obj.list)

obj.list <- lapply(X = obj.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = obj.list, reference=c(1,2), reduction = "rpca", dims = 1:50)
seuratObj <- IntegrateData(anchorset = anchors, dims = 1:50)
seuratObj <- ScaleData(seuratObj, verbose = FALSE)

plan("sequential")

seuratObj <- FindVariableFeatures(object = seuratObj, selection.method = 'vst', nfeatures = 3000)
hvg = VariableFeatures(object = seuratObj)
gene.pattern <- c("MALAT1", "^MT-", "^RPL", "^RPS", "^LOC(0-9)", "^TR(A|B|G|D)V", "^MTRNR")
hvg <- hvg[!hvg %in% grep(paste0(gene.pattern, collapse = "|"), hvg, value = T)]

seuratObj <- RunPCA(object = seuratObj, features= hvg, npcs=150, verbose = FALSE)
VariableFeatures(seuratObj) <- hvg

##save data
print("saving output")
saveRDS(seuratObj, file = opt$out)
print('----end----')
