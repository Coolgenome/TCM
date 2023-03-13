#'--------------------------------------------------------------
#' filename : SeuratMapping.R
#' Date : 2021-07-21
#' contributor : Yanshuo Chu
#' function: SeuratMapping
#' R: 4.0.3
#' seurat: 4
#'--------------------------------------------------------------

print('<==== SeuratMapping.R ====>')

suppressMessages({
    library(optparse)
    library(tidyverse)
    library(Seurat)
    library(cowplot)
})

option_list = list(
    make_option(c("-d","--data"),
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

seuratObj <- readRDS(opt$data)

## NKT_PATH="/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/validate/NKTMAIT_V6/nPC_5/UMAP_dist_0.1_nneighbor_35/p1_NKTMAIT_v6_UMAP_dist_0.1_nneighbor_35_CLUSTER_res_0.3/cluster.rds"
## seuratObj <- readRDS(NKT_PATH)

DEGsT <- read_tsv(file.path(dirname(opt$data), "snn-single-markers.tsv"))
DEGs <- unique(DEGsT %>% pull(gene))

Idents(seuratObj) <- seuratObj$batch
sortedBatch <- sort(table(seuratObj$batch), decreasing = T)
totalT <- c()

###############################################################################
#               Here Map to Each Batch, find out best reference               #
###############################################################################
# Or based on a score method to map query to reference ########################
# Compare the score method and together method, check the difference###########


for(i in 1:length(sortedBatch)){
    tempRefObj <- subset(seuratObj, idents = names(sortedBatch)[i])
    for(j in 1:length(sortedBatch)){

        if(i == j){next}
        ## if(sortedBatch[i] < 500){next}
        tryCatch({
            tempQueryObj <- subset(seuratObj,
                                   idents = names(sortedBatch)[j])

            DefaultAssay(tempRefObj) <- "RNA"
            DefaultAssay(tempQueryObj) <- "RNA"

            tempRefObj$reference.cell.type <- tempRefObj$seurat_clusters

            temp.anchors <- FindTransferAnchors(reference = tempRefObj,
                                                query = tempQueryObj,
                                                features = DEGs,
                                                k.filter = NA,
                                                reduction = "pcaproject",
                                                reference.reduction = "pca",
                                                reference.assay = "RNA",
                                                query.assay = "RNA")

            ## Error: No features to use in finding transfer anchors. To troubleshoot, try explicitly providing features to the features \
            ## 2 parameter and ensure that they are present in both reference and query assays.
            ## 1 Execution halted

            #' create a new umap model, exactly the same as the existing one #############

            tempRefObj[["umap.new"]] <- CreateDimReducObject(
                embeddings = tempRefObj[["umap"]]@cell.embeddings, key = "UMAPnew_")
            umap.new.model <- list()
            umap.new.model$n_epochs <- 500
            umap.new.model$alpha <-1
            umap.new.model$method <- "umap"
            umap.new.model$negative_sample_rate <- 5
            umap.new.model$gamma <- 1
            umap.new.model$approx_pow <- 0
            umap.new.model$metric$cosine <- list()
            umap.new.model$embedding <- tempRefObj[["umap.new"]]@cell.embeddings
            ab_param <- uwot:::find_ab_params(spread = 1, min_dist = 0.3)
            umap.new.model$a <- ab_param["a"]
            umap.new.model$b <- ab_param["b"]
            tempRefObj[["umap.new"]]@misc$model <- umap.new.model


            tempQueryObj <- MapQuery(anchorset = temp.anchors,
                                     reference = tempRefObj,
                                     query = tempQueryObj,
                                     refdata = list(celltype = "reference.cell.type"),
                                     reference.reduction = "pca",
                                     reduction.model = "umap.new")

            correctPosition <- as.numeric(
                tempQueryObj$predicted.celltype == tempQueryObj$seurat_clusters)
            wrongPosition <- as.numeric(
                tempQueryObj$predicted.celltype != tempQueryObj$seurat_clusters)

            tempTibble <- tibble(
                QueryDataSet = names(sortedBatch)[j],
                RefDataSet = names(sortedBatch)[i],
                MatchedCellNum = sum(correctPosition),
                QueryDataCellNum = length(correctPosition),
                RefDataCellNum = length(Cells(tempRefObj)),
                ACC = sum(correctPosition)/length(correctPosition))

            totalT <- bind_rows(totalT, tempTibble)

            p1 <- DimPlot(tempRefObj,
                          reduction = "umap",
                          group.by = "reference.cell.type",
                          label = TRUE,
                          label.size = 3,
                          repel = TRUE) +
                NoLegend() +
                ggtitle("Reference annotations")

            p2 <- DimPlot(tempQueryObj,
                          reduction = "ref.umap",
                          group.by = "predicted.celltype", label = TRUE,
                          label.size = 3, repel = TRUE) +
                NoLegend() +
                ggtitle("Query transferred labels")

            p3 <- DimPlot(subset(tempQueryObj, cells = Cells(tempQueryObj)[as.logical(correctPosition)]),
                          reduction = "umap",
                          group.by = "predicted.celltype", label = TRUE,
                          label.size = 3, repel = TRUE) +
                NoLegend() +
                ggtitle(paste0("Query correct, ACC=", round(sum(correctPosition)/length(correctPosition), 2)))

            p4 <- DimPlot(subset(tempQueryObj, cells = Cells(tempQueryObj)[as.logical(wrongPosition)]),
                          reduction = "umap",
                          group.by = "seurat_clusters", label = TRUE,
                          label.size = 3, repel = TRUE) +
                NoLegend() +
                ggtitle("Query wrong")

            if(str_detect(opt$data, "CD8")){
                p1 <- p1 + scale_x_reverse()
                p2 <- p2 + scale_x_reverse()
                p3 <- p3 + scale_x_reverse()
                p4 <- p4 + scale_x_reverse()
            }

            g <- plot_grid(p1,p2,p3,p4, nrow = 1, axis = "bt", align = 'h')
            pdf(file.path(opt$out, paste0("Query-", names(sortedBatch)[j],
                                          "_Ref-", names(sortedBatch)[i], "_umap.pdf")),
                width = 12, height=3)
            print(g)
            dev.off()
        },
        error = function(e){print(e)}
        )
    }
}

g <- totalT %>% ggplot() +
    geom_point(aes(x = RefDataSet,
                   y = QueryDataSet,
                   size = ACC,
                   color = RefDataCellNum,
                   fill = RefDataCellNum), shape = 22) +
    xlab("RefDataSet") + ylab("QueryDataSet") +
    theme_classic() +
    theme(text = element_text(size = 10),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pdf(file.path(opt$out, "allQuerysACC_point.pdf"))
print(g)
dev.off()

write_tsv(totalT, file.path(opt$out, paste0('totalT', "_", Sys.Date(), '.tsv')))
