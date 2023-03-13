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
                default = 'umaps.pdf',
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
norm.data <- readRDS(opt$data)

##run snn clustering

for(dist in c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001)){
    norm.umap.data <- RunUMAP(object = norm.data,
                          reduction = "harmony",
                          dims = 1:40,
                          min.dist = dist)

    norm.umap.data <- FindNeighbors(object = norm.umap.data, reduction="harmony", dims = 1:40)
    snn.obj <- FindClusters(object = norm.umap.data, resolution = 0.2)

    pdf(paste0(opt$out, 'batch: dist(', dist, ').pdf'))
    print(DimPlot(object=snn.obj, reduction="umap", group.by='batch', label=TRUE, plot.title=paste0('batch: dist(', dist)))
    print(DimPlot(object=snn.obj, reduction="umap", group.by='control', label=TRUE, plot.title=paste0('batch: dist(', dist)))
    print(DimPlot(object=snn.obj, reduction="umap", group.by='seurat_clusters', label=TRUE, plot.title=paste0('batch: dist(', dist)))
    dev.off()
}

norm.umap.data <- RunUMAP(object = norm.data,
                          reduction = "harmony",
                          dims = 1:40,
                          min.dist = 0.1)
norm.umap.data <- FindNeighbors(object = norm.umap.data, reduction="harmony", dims = 1:40)

for(res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)){
    snn.obj <- FindClusters(object = norm.umap.data, resolution = res)
    print(paste0("plot:",opt$out, 'batch: res(', res, ').pdf'))
    pdf(paste0(opt$out, 'batch: res(', res, ').pdf'))
    print(DimPlot(object=snn.obj, reduction="umap", group.by='batch', label=TRUE, plot.title=paste0('res', res)))
    print(DimPlot(object=snn.obj, reduction="umap", group.by='control', label=TRUE, plot.title=paste0('res', res)))
    print(DimPlot(object=snn.obj, reduction="umap", group.by='seurat_clusters', label=TRUE, plot.title=paste0('res', res)))
    dev.off()
}

# saveRDS(snn.obj,file = opt$out)
print('---end---')
