#' filename : visualize_batch.R
#' Date : 2020-07-07
#' contributor : Yanshuo Chu
#' function: visualize_batch

##libraries
suppressMessages({library(optparse)
library(readr)
library(rjson)
library(Seurat)
library(dplyr)
library(rlist)
library(ggpubr)
library(ggplot2)
})
print('---visualize embeding---')
##CLI parsing
option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                help = "r data file input(after runtsne/runumap)",
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

extraInWidth <- round((length(unique(seuratObj@meta.data$batch)) - 10) / 10)
isPDF <- length(Cells(seuratObj)) < 50000
if(isPDF){
    print(file.path(dirname(opt$data),'batch.pdf'))
    pdf(file.path(dirname(opt$data),'batch.pdf'), height = 9, width = (9 + extraInWidth))
    print(DimPlot(object=seuratObj, reduction="umap", group.by='batch', label=TRUE))
    dev.off()
}else{
    print(file.path(dirname(opt$data),'batch.png'))
    png(file.path(dirname(opt$data),'batch.png'), height = 9, width = (9 + extraInWidth), units = "in", res = 300)
    print(DimPlot(object=seuratObj, reduction="umap", group.by='batch', label=TRUE))
    dev.off()
}

print('---end---')
