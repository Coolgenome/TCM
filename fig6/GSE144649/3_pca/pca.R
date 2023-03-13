#'--------------------------------------------------------------
#' filename : pca.R
#' Date : 2022-05-04
#' contributor : Yanshuo Chu
#' function: pca
#'--------------------------------------------------------------

print('<==== pca.R ====>')

suppressMessages({
    library(optparse)
    library(Seurat)
    library(tidyverse)
    library(ggplot2)
})

option_list = list(
    make_option(c("-d","--data"),
                type = 'character',
                help = 'data.rds',
                metavar = 'character'),
    make_option(c("-o","--out"),
                type = 'character',
                help = 'folder',
                metavar = 'character'),
    make_option(c("-n",'--npc'),
                type = 'integer',
                default = 100,
                help = 'dims = 1:npc',
                metavar = 'integer')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

seuratObj <- readRDS(opt$data) %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(npcs = opt$npc, verbose = FALSE)

saveRDS(seuratObj, file.path(opt$out, "pca.rds"))


