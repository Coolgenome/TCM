#' monocle3.r
#' 
#' 19-10-21 10:53:54
#' 
#' contributor: guangchun
#'
#' run monocle3 for trajectory analysis
#' 

suppressMessages({library(optparse)
library(readr)
library(rjson)
library(Seurat)
library(monocle3)})
print('---monocle3---')
##CLI parsing
option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                help = "seurat r data file after normalization",
                metavar = 'character'),
    make_option(c('--out-prefix'),
                type = 'character',
                default = 'infercnv',
                help = 'prefix of output files [default = %default]',
                metavar = 'character'),
    make_option(c('--ncore'),
                type = 'integer',
                default = 1,
                help = 'number of cores to use [default = %default]',
                metavar = 'character'),
    make_option(c('-c','--param'),
                type = 'character',
                help = 'json file contain parameters',
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
    stop("json file name (containing user defined genes) must be provided", call. = F)
}
##load param
param <- fromJSON(file = opt$param)

print('perform monocle3')
obj = readRDS(opt$d)
data <- Seurat::GetAssayData(obj, assay="RNA", slot = 'counts')
pd <- data.frame(data = obj@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- fData
##Construct monocle3 cds
monocle_cds <- new_cell_data_set(data,
                              cell_metadata = pd,
                              gene_metadata = fd)
cds = monocle_cds
print('mc3: preproc')
cds = preprocess_cds(cds, num_dim = 100)
pdf(paste0(opt$`out-prefix`,'monocle3-PCA-diagnosis-',Sys.Date(),'.pdf'))
plot_pc_variance_explained(cds)
dev.off()

print('mc3: reduce dimension')
cds <- reduce_dimension(cds,reduction_method = 'tSNE') ##tsne
pdf(paste0(opt$`out-prefix`,'monocle3-tsne-',Sys.Date(),'.pdf'))
plot_cells(cds, reduction_method="tSNE", color_cells_by=paste0('data.',
                                                               param$color_by),
           show_trajectory_graph=FALSE)
dev.off()

cds <- reduce_dimension(cds) ##umap
pdf(paste0(opt$`out-prefix`,'monocle3-umap-',Sys.Date(),'.pdf'))
plot_cells(cds, reduction_method="UMAP", color_cells_by=paste0('data.',
                                                               param$color_by),
           show_trajectory_graph = F)
dev.off()

print('mc3: cluster cells')
cds <- cluster_cells(cds)
pdf(paste0(opt$`out-prefix`,'monocle3-groupcells-diagnosis-',Sys.Date(),'.pdf'))
plot_cells(cds,show_trajectory_graph = F)
plot_cells(cds, color_cells_by="partition", group_cells_by="partition",
           show_trajectory_graph = F)
dev.off()

##find marker genes for each partition
## print('mc3: find marker genes')
## marker_test_res <- top_markers(cds, group_cells_by="cluster", 
##                                reference_cells=param$reference_cells,
##                                cores=opt$ncore)

## top_specific_markers <- marker_test_res %>%
##                             filter(fraction_expressing >= 0.10) %>%
##                             group_by(cell_group) %>%
##                             top_n(param$top_n, pseudo_R2)

## top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

## pdf(paste0(opt$`out-prefix`,
##            'monocle3-topmarker-diagnosis-',param$top_n,'-',Sys.Date(),'.pdf'))
## plot_genes_by_group(cds,
##                     top_specific_marker_ids,
##                     group_cells_by="cluster",
##                     ordering_type="maximal_on_diag",
##                     max.size=3)
## dev.off()

##learn trajectory path
print('mc3: learn trajectory path')
cds <- learn_graph(cds)

##plot trajectory
pdf(paste0(opt$`out-prefix`,'monocle3-trajectory-',Sys.Date(),'.pdf'))
plot_cells(cds,
           color_cells_by = "cluster",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           show_trajectory_graph = F,
           graph_label_size=1.5)
plot_cells(cds,
           color_cells_by = "data.orig.ident",
           label_cell_groups=FALSE,
           show_trajectory_graph = F,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

dev.off()
cds <- order_cells(cds)
pdf(paste0(opt$`out-prefix`,'monocle3-trajectory-pseudotime-',Sys.Date(),'.pdf'))

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

dev.off()
saveRDS(cds,file= paste0(opt$`out-prefix`,'.monocle3Obj.rds'))
print('---end---')
