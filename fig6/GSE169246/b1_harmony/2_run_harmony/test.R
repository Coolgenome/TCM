##libraries
library(Seurat)
library(tidyverse)
library(harmony)

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

## seuratObj <- readRDS(opt$data)
seuratObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/b1_harmony/1_injectBatchinfo/outs/harmony_input.rds")

## By default, use all pc
seuratObj <- RunHarmony(seuratObj, "batch", assay.use = "RNA", reduction = "pca", verbose = T)
saveRDS(seuratObj, "/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/b1_harmony/2_run_harmony/harmony_out.obj")

print('---end---')
