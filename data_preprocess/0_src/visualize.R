#' filename : visualize.R
#' Date : 2020-07-07
#' contributor : Yanshuo Chu
#' function: visualize

##libraries


suppressMessages({library(optparse)
library(readr)
library(rjson)
library(Seurat)
library(dplyr)
library(rlist)
library(ggpubr)
library(ggplot2)})
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

coord = Embeddings(object = seuratObj, reduction = "umap")
coord = coord[,c(1,2)]
colnames(coord) = c("dim1", "dim2")
coord = data.frame(ID = rownames(coord), coord)
meta = seuratObj@meta.data;
meta = data.frame(ID = rownames(meta), meta,stringsAsFactors = F)
meta = left_join(meta, coord, by = 'ID')
write_tsv(meta, file.path(dirname(opt$data),'visualization_coordinates.tsv'))

isPDF <- length(Cells(seuratObj)) < 50000
extraInWidth <- round((length(unique(seuratObj@meta.data$seurat_clusters)) - 10) / 10)
if(isPDF){
    print(file.path(dirname(opt$data),'umap.pdf'))
    pdf(file.path(dirname(opt$data),'umap.pdf'), height = 9, width = (9 + extraInWidth))
    print(DimPlot(object=seuratObj, reduction="umap", group.by='seurat_clusters', label=TRUE))
    dev.off()
}else{
    print(file.path(dirname(opt$data),'umap.png'))
    png(file.path(dirname(opt$data),'umap.png'), height = 9, width = (9 + extraInWidth), units = "in", res=600)
    print(DimPlot(object=seuratObj, reduction="umap", group.by='seurat_clusters', label=TRUE))
    dev.off()
}

print('---end---')
