#--------------------------------------------------------------
# filename : Mapping.R
# Date : 2022-04-28
# contributor : Yanshuo Chu
# function: MultiMap
#--------------------------------------------------------------

print('<==== Mapping.R ====>')

suppressMessages({
    library(optparse)
    library(tidyverse)
    library(Seurat)
    library(SeuratObject)
    library(cowplot)
    library(MultiMap)
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

## queryPath <- "/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE179994/TCells/ForMapping/CD8.rds"

refSeuratObj <- readRDS(opt$referenceData)
querySeuratObj <- readRDS(opt$queryData)
## refSeuratObj <- readRDS(CD8_ReferenceDataPath)
## querySeuratObj <- readRDS(queryPath)

DefaultAssay(refSeuratObj) <- "RNA"
DefaultAssay(querySeuratObj) <- "RNA"

A.L <- MultiMap::FindTransferAnchors_MB(reference = refSeuratObj,
                                        query = querySeuratObj,
                                        features = intersect(rownames(refSeuratObj), rownames(querySeuratObj)),
                                        reduction = "pcaproject")
gamma.L <- MultiMap::GetRefBatchWeight_MB(A.L, refSeuratObj@meta.data$seurat_clusters)

querySeuratObj <- querySeuratObj %>%
    NormalizeData(verbose = T) %>%
    FindVariableFeatures(selection.method = "vst") %>%
    ScaleData(verbose = T) %>%
    RunPCA(verbose = T, npcs = 80)

sortedBatch <- gamma.L %>%
    slice_max(Gamma, prop = 0.75) %>%
    pull(batch)

Idents(refSeuratObj) <- refSeuratObj$batch
predicted.L <- list()
for(i in 1:length(sortedBatch)){
    RefBatchName <- sortedBatch[i]
    tempRefObj <- subset(refSeuratObj, idents = RefBatchName)

    tempRefObj <- tempRefObj %>%
        NormalizeData(verbose = T) %>%
        FindVariableFeatures(selection.method = "vst") %>%
        ScaleData(verbose = T) %>%
        RunPCA(verbose = T, npc = 80)

    temp.anchors <- FindTransferAnchors(reference = tempRefObj,
                                        query = querySeuratObj,
                                        features = intersect(rownames(tempRefObj), rownames(querySeuratObj)),
                                        dims = 1:30,
                                        reference.reduction = "pca")

    tempQueryObj <- MapQuery(anchorset = temp.anchors,
                             reference = tempRefObj,
                             query = querySeuratObj,
                             refdata = tempRefObj$seurat_clusters)

    predicted.L[[i]] <- tempQueryObj@meta.data$predicted.id
    names(predicted.L)[i] <- RefBatchName
}

predictedT <- bind_cols(predicted.L)
predictedT <- predictedT %>% mutate(MostFreq = apply(predictedT[, 1:dim(predictedT)[2]],
                            MARGIN =1,
                            function(x){return(names(which.max(table(x))))}))

saveRDS(predictedT, file.path(opt$out, paste0('predictedT', "_", Sys.Date(), '.rds')))
