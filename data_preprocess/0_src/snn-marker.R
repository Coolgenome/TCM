#' cluster-markers.r
#' 
#' 19-10-07 09:53:55
#' 
#' contributor: guangchun
#'
#' cluster markers analysis
#' 

##libraries
suppressMessages({
    library(optparse)
    library(readr)
    library(rjson)
    library(Seurat)
})
print('---calculate snn cluster makers---')
##CLI parsing
option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                default = NULL,
                help = "r data file input(after snn clustering)",
                metavar = 'character'),
    make_option(c("-o",'--out'),
                type = 'character',
                default = 'markers.txt',
                help = 'output file name for the markers [default = %default]',
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
seuratObj <- readRDS(opt$data)
Idents(seuratObj) <- seuratObj@meta.data$seurat_clusters

## hvg <- VariableFeatures(seuratObj)
## gene.pattern <- c("MALAT1", "^MT-", "^RPL", "^RPS", "^LOC(0-9)", "^TR(A|B|G|D)V", "^MTRNR")
## hvg <- hvg[!hvg %in% grep(paste0(gene.pattern, collapse = "|"), hvg, value = T)]

markers <- FindAllMarkers(object = seuratObj)
## markers <- FindAllMarkers(object = seuratObj,
##                           only.pos = TRUE,
##                           min.pct = 0.25,
##                           logfc.threshold = 0.25)

write_tsv(data.frame(markers), opt$out)
print('---end---')
