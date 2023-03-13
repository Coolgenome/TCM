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
    make_option(c("-i",'--dist'),
                type = 'double',
                default = 0.4,
                help = 'dist default 0.4',
                metavar = 'double'),
    make_option(c("-e",'--nneighbors'),
                type = 'integer',
                default = 30,
                help = 'n neighbors, default 30',
                metavar = 'integer')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$data)) {
    print_help(opt_parser)
    stop("Input data must be provided", call. = F)
}

##Load data
seuratObj <- readRDS(opt$data)

## ## ###############################################################################
## ## #'                  scale regress out proliferative markers                  '#
## ## ###############################################################################
## cellCycleGeneT1 <- read_tsv("/rsrch3/home/genomic_med/ychu2/configs/scSeqs/database/general/cell-cycle-gene-list.txt")
## cellCycleGeneT2 <- read_tsv("/rsrch3/home/genomic_med/ychu2/configs/scSeqs/database/general/regev_lab_cell_cycle_genes.txt")
## seuratObj <- CellCycleScoring(seuratObj,
##                               s.features = s.genes,
##                               g2m.features = g2m.genes,
##                               set.ident = TRUE)
## seuratObj <- ScaleData(seuratObj,
##                        vars.to.regress =c("S.Score", "G2M.Score"),
##                        features = rownames(seuratObj))
## seuratObj <- FindVariableFeatures(object = seuratObj, selection.method = 'vst', nfeatures = 3000)
## hvg = VariableFeatures(object = seuratObj)
## gene.pattern <- c("MALAT1", "^MT-", "^RPL", "^RPS", "^LOC(0-9)", "^TR(A|B|G|D)V", "^MTRNR")
## hvg <- hvg[!hvg %in% grep(paste0(gene.pattern, collapse = "|"), hvg, value = T)]
## hvg <- setdiff(hvg, cellCycleGeneT1)
## hvg <- setdiff(hvg, cellCycleGeneT2)
## seuratObj <- RunPCA(object = seuratObj, features= hvg, npcs=150, verbose = FALSE)
## VariableFeatures(seuratObj) <- hvg
## ## ###############################################################################

##run snn clustering
seuratObj <- RunUMAP(object = seuratObj,
                     reduction = opt$reduction,
                     dims = 1:opt$npc,
                     min.dist = opt$dist,
                     n.neighbors = opt$nneighbors)

saveRDS(seuratObj, file = opt$out)
print('---end---')
