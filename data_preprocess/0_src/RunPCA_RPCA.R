#'-----------------------------------
#' filename : RunPCA_RPCA.R
#' Date : 2020-07-26
#' contributor : Yanshuo Chu
#' function: RunPCA_RPCA
#'-----------------------------------

print('<== Run PCA ==>')

##libraries
suppressMessages({library(optparse)
    library(tidyverse)
    library(Seurat)})



##CLI parsing
option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                default = NULL,
                help = "r data file input(after normalization",
                metavar = 'character'),
    ## make_option(c("-t", "--table"),
    ##             type = "character",
    ##             default = NULL,
    ##             help = "bad gene table",
    ##             metavar = 'character'),
    make_option(c("-o",'--out'),
                type = 'character',
                default = 'snn-harmony.rds',
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
seuratObj <- readRDS(opt$data)
DefaultAssay(seuratObj) <- "integrated"


cellCycleGeneT1 <- read_tsv("/rsrch3/home/genomic_med/ychu2/configs/scSeqs/database/general/cell-cycle-gene-list.txt")
cellCycleGeneT2 <- read_tsv("/rsrch3/home/genomic_med/ychu2/configs/scSeqs/database/general/regev_lab_cell_cycle_genes.txt")
seuratObj <- FindVariableFeatures(object = seuratObj, selection.method = 'vst', nfeatures = 3000)
hvg = VariableFeatures(object = seuratObj)
gene.pattern <- c("MALAT1", "^MT-", "^RPL", "^RPS", "^LOC(0-9)", "^TR(A|B|G|D)V", "^MTRNR")
hvg <- hvg[!hvg %in% grep(paste0(gene.pattern, collapse = "|"), hvg, value = T)]
## hvg <- setdiff(hvg, cellCycleGeneT1$marker)
## hvg <- setdiff(hvg, cellCycleGeneT2$marker)
## print("move out proliferative markers")
seuratObj <- RunPCA(object = seuratObj, features= hvg, npcs=150, verbose = FALSE)
VariableFeatures(seuratObj) <- hvg


## seuratObj <- FindVariableFeatures(object = seuratObj, selection.method = 'vst', nfeatures = 2200)
## hvg = VariableFeatures(object = seuratObj)
## gene.pattern <- c("MALAT1", "^MT-", "^RPL", "^RPS", "^LOC(0-9)", "^TR(A|B|G|D)V", "^MTRNR")
## hvg <- hvg[!hvg %in% grep(paste0(gene.pattern, collapse = "|"), hvg, value = T)]
## print("length(hvg)")
## print(length(hvg))
## seuratObj <- RunPCA(object = seuratObj, features= hvg, npcs=150, verbose = FALSE)
## VariableFeatures(seuratObj) <- hvg

saveRDS(seuratObj, file = opt$out)
print('---end---')

