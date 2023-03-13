#' filename : umap-harmony.R
#' Date : 2020-04-23
#' contributor : Yanshuo Chu
#' function: run umap for harmony data

##libraries
suppressMessages({library(optparse)
library(readr)
library(rjson)
library(SeuratData)
library(harmony)
library(Seurat)})
print('---snn clustering---')
##CLI parsing
option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                default = NULL,
                help = "r data file input(after normalization",
                metavar = 'character'),
    make_option(c("-o",'--out'),
                type = 'character',
                default = 'snn-harmony.rds',
                help = 'output file name for the r data file [default = %default]',
                metavar = 'character'),
    make_option(c("-r",'--reduction'),
                type = 'character',
                default = 'harmony',
                help = 'reduction method harmony',
                metavar = 'character'),
    make_option(c("-n",'--npc'),
                type = 'integer',
                default = 40,
                help = 'npc default 4 for dims',
                metavar = 'integer'),
    make_option(c("-e",'--resolution'),
                type = 'double',
                default = 0.2,
                help = 'resolution default 0.4',
                metavar = 'double')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$data)) {
    print_help(opt_parser)
    stop("Input data must be provided", call. = F)
}

##Load data
seuratObj <- readRDS(opt$data)


seuratObj <- FindNeighbors(object = seuratObj,
                           reduction=opt$reduction,
                           dims = 1:opt$npc)

seuratObj <- FindClusters(object = seuratObj,
                          resolution = opt$resolution)


DefaultAssay(seuratObj) <- "RNA"

saveRDS(seuratObj, file = opt$out)
print('---end---')
