#' filename : RunCCA.R
#' Date : 2020-07-07
#' contributor : Yanshuo Chu
#' function: RunCCA


suppressMessages({library(optparse)
library(ggplot2)
library(readr)
library(stringr)
library(rjson)
library(tidyverse)
library(Seurat)})
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

print(paste0("Loading data from ", opt$data))
samples = list.files(path = opt$data, recursive = F)
samples = basename(samples)
if(length(samples) < 1) stop(paste0('Failed to find data in ', opt$data))

obj.list = list()
genes.use = list()
for (ids in seq_along(samples)) {
    ds = samples[ids]
    print(ds)
    seuratObj <- readRDS(file.path(opt$data, ds))
    genes.use <- c(genes.use, head(rownames(seuratObj@hvg.info), 2000))

    obj.list[[length(obj.list)+1]] = seuratObj
}

for (i in 1:length(obj.list)) {
    genes.use <- genes.use[genes.use %in% rownames(obj.list[[i]]@scale.data)]
}

seuratObj <- RunMultiCCA(obj.list, genes.use = genes.use)

##save data
print("saving output")
saveRDS(combined.data, file = opt$out)
print('----end----')
