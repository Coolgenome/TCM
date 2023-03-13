#'--------------------------------------------------------------
#' filename : monocle2.R
#' Date : 2021-03-18
#' contributor : Yanshuo Chu
#' function: monocle2
#'--------------------------------------------------------------

print('<==== monocle2 ====>')

suppressMessages({
    library(optparse)
    library(readr)
    library(rjson)
    library(Seurat)
    library(monocle)
    library(tidyverse)
})
print('---monocle---')
##CLI parsing
option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                help = "seurat r data file after normalization",
                metavar = 'character'),
    make_option(c("-f", "--frac"),
                type = 'double',
                default = 0.2,
                help = 'fraction of cells',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$data)) {
    print_help(opt_parser)
    stop("Input data must be provided", call. = F)
}


print('perform monocle')
outputFolder <- file.path(dirname(opt$d), "monocle2")
if(!dir.exists(outputFolder)){
    dir.create(outputFolder)
}

## seuratObj = readRDS("/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/validate/Sub_Treg_CD4_V5/nPC_15/UMAP_dist_0.1_nneighbor_20/p1_sub_Treg_CD4_V5_UMAP_dist_0.1_nneighbor_20_CLUSTER_res_0.3/cluster.rds")

seuratObj = readRDS(opt$d)

## markerList <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/validate/Sub_Treg_CD4_V5/nPC_15/UMAP_dist_0.1_nneighbor_20/p1_sub_Treg_CD4_V5_UMAP_dist_0.1_nneighbor_20_CLUSTER_res_0.3/snn-single-markers.tsv")

if(file.exists(file.path(dirname(opt$d), "snn-single-markers.tsv"))){
    markerList <- read_tsv(file.path(dirname(opt$d), "snn-single-markers.tsv"))
}else if(file.exists(file.path(dirname(opt$d), "snn-markers.tsv"))){
    markerList <- read_tsv(file.path(dirname(opt$d), "snn-markers.tsv"))
}


gene.pattern <- c("MALAT1", "^HSPA", "^MT-", "^MT\\.", "^RPL", "^RPS", "^LOC(0-9)", "^TR(A|B|G|D)V", "^MTRNR")
noiseGenesByPattern <- grep(paste0(gene.pattern, collapse = "|"), markerList$gene, value = T)
noiseGenes <- union(c("ZNF14", "VIPAS39", "ZNF408", "TADA1", "ARMCX5",
                      "DCUN1D4", "ATG4C", "ZNF559", "TRIB3", "ZNF337", "IFIT5",
                      "C1orf216", "DHX58", "S100B", "EGR2", "ZKSCAN8", "PDIK1L",
                      "ZBTB6", "SLC40A1", "KBTBD7"), noiseGenesByPattern)
## avg 0.7 pct.ratio 1.5
markerList <- markerList %>% mutate(pct.ratio = pct.1 / pct.2) %>%
    filter(! gene %in% noiseGenes & avg_logFC >= 0.5 & pct.ratio >= 1.5 & pct.1 >= 0.2) %>%
    arrange(desc(cluster), desc(pct.ratio), avg_logFC) %>% group_by(cluster)

md <- as_tibble(seuratObj@meta.data, rownames = NA)
md$barcode <- rownames(md)

md <- md %>% group_by(seurat_clusters) %>% sample_frac(opt$f)
## md <- md %>% group_by(seurat_clusters) %>% sample_frac(0.01)

seuratObj <- subset(seuratObj, cells = md$barcode)
Idents(seuratObj) <- seuratObj$seurat_clusters

data <- as(seuratObj@assays$RNA@counts, 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = seuratObj@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

##Construct monocle cds
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

cds = monocle_cds
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)
cds <- setOrderingFilter(cds, unique(markerList$gene))
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree', auto_param_selection = F)
cds <- orderCells(cds)

plots <- list(
    Pseudotime = plot_cell_trajectory(
        cds,
        color_by = "Pseudotime",
        cell_size=0.1),
    seurat_clusters = plot_cell_trajectory(
        cds, color_by = "seurat_clusters",
        cell_size=0.1)+
        geom_text(aes(label = seurat_clusters, color = seurat_clusters)),
    seurat_clusters_withoutlabel = plot_cell_trajectory(
        cds, color_by = "seurat_clusters",
        cell_size=0.1),
    seurat_clusters_withoutlabel_split = plot_cell_trajectory(
        cds, color_by = "seurat_clusters",
        show_state_number = F,
        show_branch_points = F,
        cell_size=1) + facet_wrap(seurat_clusters ~ .)
)


for(i in 1:length(plots))
{
    pdf(file.path(outputFolder, paste0("f_", opt$f, names(plots[i]), 'monocle-plots.pdf')))
    print(plots[[i]])
    dev.off()
}

saveRDS(cds,file= file.path(outputFolder, paste0(opt$f, '_monocleobj.rds')))
print('---end---')
