#' filename : determinePC.R
#' Date : 2020-07-07
#' contributor : Yanshuo Chu
#' function: determinePC

print('----determinePC----')

##libraries
suppressMessages({
    library(optparse)
    library(future)
    library(readr)
    library(rjson)
    library(Seurat)
})

##CLI parsing
option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                default = NULL,
                help = "r data file input(after normalization",
                metavar = 'character'),
    make_option(c("-o",'--out'),
                type = 'character',
                default = 'seuratObj.rds',
                help = 'output of seuratObj to determinePC',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$data)) {
    print_help(opt_parser)
    stop("Input data must be provided", call. = F)
}

##Load data
seuratObj <- readRDS(opt$data)

options(future.globals.maxSize = 10 * 1000 * 1024^2)

plan("multiprocess", workers = 12)
seuratObj <- NormalizeData(object = seuratObj,
                           normalization.method = "LogNormalize",
                           scale.factor = 1e4)
plan("sequential")

seuratObj <- FindVariableFeatures(object = seuratObj, selection.method = 'vst', nfeatures = 2000)

hvg = VariableFeatures(object = seuratObj)

plan("multiprocess", workers = 12)
seuratObj <- ScaleData(object = seuratObj,
                       features = hvg,
                       vars.to.regress = c("nCount_RNA", "percent.mito"))
plan("sequential")

seuratObj <- RunPCA(object = seuratObj, features = hvg, npcs=150, verbose = FALSE)

##generate PCA loadings plot
pdf(paste0(opt$out, ".elbowplot.pdf"))
ElbowPlot(object = seuratObj, ndims = 150, reduction = 'pca')
dev.off()

saveRDS(seuratObj, opt$out)
print('----end----')
