#--------------------------------------------------------------
# filename : Mapping.R
# Date : 2022-05-02
# contributor : Yanshuo Chu
# function: Mapping
#--------------------------------------------------------------

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


## CD8_ReferenceDataPath="/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/validate/CD8_V6/nPC_50/UMAP_dist_0.1_nneighbor_50/p1CD8_V6_UMAP_dist_0.1_nneighbor_50_CLUSTER_res_0.3/cluster.rds"
## CD4_ReferenceDataPath="/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/validate/CD4_V7/nPC_50/UMAP_dist_0.1_nneighbor_50/p1CD4_V7_UMAP_dist_0.1_nneighbor_50_CLUSTER_res_0.3/cluster.rds"
## QUERYDATAFOLDER="/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE186144/1_split/A"

refSeuratObj <- readRDS(opt$referenceData)
querySeuratObj <- readRDS(opt$queryData)

## refSeuratObj <- readRDS(CD8_ReferenceDataPath)
## querySeuratObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE144649/5_extractT/CD8.rds")

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
hvgR = VariableFeatures(object = querySeuratObj)
querySeuratObj <- querySeuratObj %>%
    ScaleData(verbose = T) %>%
    RunPCA(verbose = T, features = hvgR)


temp.anchors <- FindTransferAnchors(reference = refSeuratObj,
                                    query = querySeuratObj,
                                    reference.reduction = "pca",
                                    k.filter = NA,
                                    dims = 1:20,
                                    features = intersect(rownames(refSeuratObj), rownames(querySeuratObj)))

querySeuratObj <- MapQuery(anchorset = temp.anchors,
                          reference = refSeuratObj,
                          query = querySeuratObj,
                          refdata = refSeuratObj$seurat_clusters)

querySeuratObj$predicted.id <-
    factor(querySeuratObj$predicted.id,
           levels = levels(refSeuratObj$seurat_clusters))

saveRDS(querySeuratObj, file.path(opt$out, paste0('querySeuratObj', "_", Sys.Date(), '.rds')))
