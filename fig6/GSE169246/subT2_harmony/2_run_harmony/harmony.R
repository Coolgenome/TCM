#--------------------------------------------------------------
# filename : harmony.R
# Date : 2022-09-01
# contributor : Yanshuo Chu
# function: harmony
#--------------------------------------------------------------

print('<==== harmony.R ====>')

suppressMessages({
    library(optparse)
    library(tidyverse)
    library(harmony)
    library(Seurat)})

option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                default = NULL,
                help = "r data file input(after normalization",
                metavar = 'character'),
    make_option(c("-o",'--out'),
                type = 'character',
                default = 'harmony.pdf',
                help = 'output file name for the r data file [default = %default]',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$data)) {
    print_help(opt_parser)
    stop("Input data must be provided", call. = F)
}

##Load data
seuratObj <- readRDS(opt$data) %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = 100, verbose = FALSE)

## By default, use all pc
seuratObj <- RunHarmony(seuratObj, "batch")

saveRDS(seuratObj, opt$out)

print('---end---')
