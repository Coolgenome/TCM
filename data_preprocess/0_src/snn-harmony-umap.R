##libraries
suppressMessages({library(optparse)
library(readr)
library(rjson)
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
                default = 'snn.rds',
                help = 'output file name for the r data file [default = %default]',
                metavar = 'character'),
    make_option(c("-c","--param"),
                type = 'character',
                help = 'json file name contain function parameters',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$data)) {
    print_help(opt_parser)
    stop("Input data must be provided", call. = F)
}
if(is.null(opt$param)) {
    print_help(opt_parser)
    stop("json file name (containing parameters) must be provided", call. = F)
}

##load param
param <- fromJSON(file = opt$param)

##Load data
norm.data <- readRDS(opt$data)
norm.data <- FindVariableFeatures(object = norm.data, selection.method = 'vst',
                                  nfeatures = 2000)
length(x = VariableFeatures(object = norm.data))
hvg = VariableFeatures(object = norm.data)
norm.data <- ScaleData(object = norm.data, features = hvg,
                       vars.to.regress = c("nCount_RNA", "percent.mito"))
norm.data <- RunPCA(object = norm.data, features = hvg, verbose = FALSE)
norm.data <- RunHarmony(norm.data, "batch")

npc = param$npc
##run snn clustering
norm.data <- RunUMAP(object = norm.data,
                     reduction = "harmony",
                     dims = 1:npc,
                     min.dist = param$dist,
                     n.neighbors = param$nneigh)

norm.data <- FindNeighbors(object = norm.data, reduction="harmony",
                           dims = 1:npc, k.param = param$k)
snn.obj <- FindClusters(object = norm.data, resolution = param$res)

saveRDS(snn.obj,file = opt$out)
print('---end---')
