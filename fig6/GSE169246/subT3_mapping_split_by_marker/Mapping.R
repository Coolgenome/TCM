#'--------------------------------------------------------------
#' filename : Mapping.R
#' Date : 2021-08-08
#' contributor : Yanshuo Chu
#' function: Mapping
#'--------------------------------------------------------------

print('<==== Mapping.R ====>')

suppressMessages({
    library(optparse)
    library(tidyverse)
    library(Seurat)
    library(SeuratObject)
    library(cowplot)
})

option_list = list(
    make_option(c("-r","--referenceData"),
                type = 'character',
                help = 'data.rds',
                metavar = 'character'),
    make_option(c("-q","--queryData"),
                type = 'character',
                help = 'data.rds',
                metavar = 'character'),
    make_option(c("-o","--out"),
                type = 'character',
                help = 'out',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

refSeuratObj <- readRDS(opt$referenceData)
querySeuratObj <- readRDS(opt$queryData)

DefaultAssay(refSeuratObj) <- "RNA"
DefaultAssay(querySeuratObj) <- "RNA"

refSeuratObj <- refSeuratObj %>%
    NormalizeData(verbose = T) %>%
    FindVariableFeatures(selection.method = "vst")
hvgR = VariableFeatures(object = refSeuratObj)
refSeuratObj <- refSeuratObj %>%
    ScaleData(verbose = T) %>%
    RunPCA(verbose = T, features = hvgR)

querySeuratObj <- querySeuratObj %>%
NormalizeData(verbose = T) %>%
FindVariableFeatures(selection.method = "vst")
hvgQ = VariableFeatures(object = querySeuratObj)
querySeuratObj <- querySeuratObj %>%
    ScaleData(verbose = T) %>%
    RunPCA(verbose = T, features = hvgQ)

refSeuratObj$reference.cell.type <- refSeuratObj$seurat_clusters
temp.anchors <- FindTransferAnchors(reference = refSeuratObj,
                                    query = querySeuratObj,
                                    reduction = "pcaproject",
                                    reference.reduction = "pca",
                                    reference.assay = "RNA",
                                    query.assay = "RNA",
                                    features = intersect(hvgR, hvgQ))

#' create a new umap model, exactly the same as the existing one #############
refSeuratObj[["umap.new"]] <- CreateDimReducObject(embeddings = refSeuratObj[["umap"]]@cell.embeddings, key = "UMAPnew_")
umap.new.model <- list()
umap.new.model$n_epochs <- 500
umap.new.model$alpha <-1
umap.new.model$method <- "umap"
umap.new.model$negative_sample_rate <- 5
umap.new.model$gamma <- 1
umap.new.model$approx_pow <- 0
umap.new.model$metric$cosine <- list()
umap.new.model$embedding <- refSeuratObj[["umap.new"]]@cell.embeddings
ab_param <- uwot:::find_ab_params(spread = 1, min_dist = 0.3)
umap.new.model$a <- ab_param["a"]
umap.new.model$b <- ab_param["b"]
refSeuratObj[["umap.new"]]@misc$model <- umap.new.model

querySeuratObj <- MapQuery(anchorset = temp.anchors,
                          reference = refSeuratObj,
                          query = querySeuratObj,
                          refdata = list(celltype = "reference.cell.type"),
                          reference.reduction = "pca",
                          reduction.model = "umap.new")

Idents(refSeuratObj) <- refSeuratObj$seurat_clusters
Idents(querySeuratObj) <- querySeuratObj$predicted.celltype
colorsForDataType <- c("#6DCCDD", "#EDCAE0", "#F494BE", "#F9B26C", "#A6ADCC", "#C4DA5D")
umapColor <- colorRampPalette(colorsForDataType)(length(unique(refSeuratObj$seurat_clusters)))

p1 <- DimPlot(refSeuratObj,
              reduction = "umap",
              group.by = "seurat_clusters",
              label = F,
              label.size = 3,
              pt.size = 0.1,
              repel = F) +
    scale_color_manual(values = umapColor) +
    theme_void() +
    theme(text = element_text(size = 8),
          legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank())

querySeuratObj$predicted.celltype <-
    factor(querySeuratObj$predicted.celltype,
           levels = levels(refSeuratObj$seurat_clusters))

p2 <- DimPlot(querySeuratObj,
              reduction = "ref.umap",
              group.by = "predicted.celltype",
              label = F,
              pt.size = 0.1) +
    scale_color_manual(values = umapColor) +
    theme_void() +
    theme(text = element_text(size = 8),
          legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank())

coord = Embeddings(
    object = querySeuratObj,
    reduction = "ref.umap")
coord = coord[,c(1,2)]
colnames(coord) = c("dim1", "dim2")
coord = data.frame(ID = rownames(coord), coord)
meta = querySeuratObj@meta.data
meta = data.frame(ID = rownames(meta),
                  meta,
                  stringsAsFactors = F)
meta = left_join(meta, coord, by = 'ID')

if(str_detect(opt$referenceData, "CD8")){
    p1 <- p1 + scale_x_reverse()
    p2 <- p2 + scale_x_reverse()
}
png(file.path(opt$out, "reference.png"))
print(p1)
dev.off()
png(file.path(opt$out, "query_mapped.png"))
print(p2)
dev.off()

## for(i in 1:10){
##     g <- ggplot(meta) +
##         geom_point(
##             aes(x = dim1,
##                 y = dim2,
##                 color = predicted.celltype),
##             size = 0.01 * i) +
##         scale_color_manual(values = umapColor) +
##         theme_void() +
##         theme(text = element_text(size = 8),
##               legend.position = "none",
##               axis.text.x = element_blank(),
##               axis.text.y = element_blank(),
##               axis.ticks = element_blank())
##     if(str_detect(opt$referenceData, "CD8")){
##         g  <- g  + scale_x_reverse()
##     }
##     png(file.path(opt$out, paste0("query2_mapped_", i, ".png")))
##     print(g)
##     dev.off()
## }

write_tsv(meta, file.path(opt$out, paste0('meta', "_", Sys.Date(), '.tsv')))

## saveRDS(querySeuratObj, file.path(opt$out, paste0('querySeuratObj', "_", Sys.Date(), '.rds')))
