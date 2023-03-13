#'--------------------------------------------------------------
#' filename : DEG-bubble-plot.R
#' Date : 2020-08-25
#' contributor : Yanshuo Chu
#' function: DEG-bubble-plot
#'--------------------------------------------------------------

print('<==== DEG-bubble-plot ====>')

suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(ggplot2)
    library(optparse)
    library(ggplot2)
    library(RColorBrewer)
    library(ggpubr)
    library(viridis)
    library(rlist)
})


option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                help = "r data file input(after runtsne/runumap)",
                metavar = 'character'),
    make_option(c("-g","--DEGs"),
                type = 'character',
                help = 'DEG table file (by findallmarkers)',
                metavar = 'character'),
    make_option(c("-n","--number"),
                type = 'integer',
                default = 5,
                help = 'number of DEGs',
                metavar = 'integer')
);


opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$data)) {
    print_help(opt_parser)
    stop("Input data must be provided", call. = F)
}

if(is.null(opt$DEGs)) {
    print_help(opt_parser)
    stop("DEGs must be provided", call. = F)
}

bubblePlotFolder <- file.path(dirname(opt$data), "DEG_BubblePlot")
if(!dir.exists(bubblePlotFolder)){
    dir.create(bubblePlotFolder)
}

featurePlotFolder <- file.path(dirname(opt$data), "DEG_FeaturePlot")
if(!dir.exists(featurePlotFolder)){
    dir.create(featurePlotFolder)
}

seuratObj <- readRDS(opt$data)

## seuratObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/CD4_V4/nPC_50/UMAP_dist_0.1_nneighbor_50/p1CD4_V4_50_UMAP_dist_0.1_nneighbor_50_CLUSTER_res_0.3/cluster.rds")

Idents(seuratObj) <- seuratObj@meta.data$seurat_clusters

#' filter out noise genes #####################################################


markerList <- read_tsv(opt$DEGs)

## markerList <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/CD4_V4/nPC_50/UMAP_dist_0.1_nneighbor_50/p1CD4_V4_50_UMAP_dist_0.1_nneighbor_50_CLUSTER_res_0.3/snn-single-markers.tsv")

gene.pattern <- c("MALAT1", "^HSPA", "^MT-", "^MT\\.", "^RPL", "^RPS", "^LOC(0-9)", "^TR(A|B|G|D)V", "^MTRNR")
noiseGenesByPattern <- grep(paste0(gene.pattern, collapse = "|"), markerList$gene, value = T)
noiseGenes <- union(c("ZNF14", "VIPAS39", "ZNF408", "TADA1", "ARMCX5", "DCUN1D4", "ATG4C", "ZNF559", "TRIB3", "ZNF337", "IFIT5", "C1orf216", "DHX58",
                "S100B", "EGR2", "ZKSCAN8", "PDIK1L", "ZBTB6", "SLC40A1", "KBTBD7"), noiseGenesByPattern)

## avg 0.7 pct.ratio 1.5
markerList <- markerList %>% mutate(pct.ratio = pct.1 / pct.2) %>%
    filter(! gene %in% noiseGenes & avg_logFC >= 0.5 & pct.ratio >= 1.5 & pct.1 >= 0.2) %>%
    arrange(desc(cluster), desc(pct.ratio), avg_logFC) %>% group_by(cluster) %>%
    slice_max(order_by = pct.1, n = opt$number)

write_tsv(markerList, file.path(dirname(opt$DEGs), paste0('top20-DEGs', "_", Sys.Date(), '.tsv')))

## markerList <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/CD4_V4/nPC_50/UMAP_dist_0.1_nneighbor_50/p1CD4_V4_50_UMAP_dist_0.1_nneighbor_50_CLUSTER_res_0.3/snn-single-markers.tsv") %>% filter(!gene %in% noiseGenes) %>% arrange(desc(cluster), avg_logFC) %>% group_by(cluster) %>% slice_max(order_by = avg_logFC, n = 5)

#' generate list with name ####################################################
genes.list <- list()
for(cci in unique(markerList$cluster)){
    genes.list[[as.character(cci)]] <- unique(markerList$gene[markerList$cluster == cci])
}

#' add total gene to list to generate a total bubble plot #####################
## DEGFileName <- tools::file_path_sans_ext(basename(opt$DEGs))
## genes.list[[DEGFileName]] <- markerList$gene


for(i in 1:length(genes.list)){
    gene.list.name <- names(genes.list[i])
    gene.list <- genes.list[[i]]

    gene<-intersect(gene.list,rownames(seuratObj))
    p<-DotPlot(seuratObj, features = gene)
    data<-p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]
    #data[data$avg.exp>7,'avg.exp']<-7
    plotx<-ggplot(data, aes(x = features.plot,y = id)) +        ## global aes
      geom_point(aes(fill = avg.exp.scaled,size =pct.exp),color='black',shape=21)  +    ## geom_point for circle illusion
      #scale_fill_gradientn(colours=rev(color),limits=c(0,max(data$avg.exp)))+       ## color of the corresponding aes
      scale_fill_viridis()+
      scale_size(range = c(0,6), limits = c(0, 100), breaks = c(0,20,40,60,80,100))+             ## to tune the size of circles
      theme(panel.grid.major = element_line(colour = "grey90",size=0.2), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.text.x=element_text(angle = 90,vjust = 1,hjust = 1) )
    ht <- round(8/31 * length(unique(data$id)))
    wh <- round(4/10 * length(gene.list))
    if(ht < 3){ht = 3;}
    if(wh < 5){wh = 5;}
    ggsave(
        file.path(bubblePlotFolder,
                  paste0(tools::file_path_sans_ext(basename(opt$DEGs)),
                         '.', gene.list.name, '.marker.bubble.pdf')),
           plotx, height=ht, width=wh, limitsize=F)

    gp.list=list()
    for(genei in 1: length(gene)){
        gp.list = list.append(gp.list,
                              annotate_figure(p = FeaturePlot(seuratObj,
                                                              reduction = "umap",
                                                              features = as.vector(gene[genei]),
                                                              cols=c("lightgray", "blue", "black")),
                                              top = ggpubr::text_grob(label = gene.list.name, face="bold", size = 20, color="red")
                                              )
                              )
    }

    combined.gp <- do.call(ggarrange, c(gp.list, ncol = 2, nrow = 2))
    pdf(
        file.path(featurePlotFolder,
                  paste0(tools::file_path_sans_ext(basename(opt$DEGs)),
                         '.', gene.list.name, '.feature.bubble.pdf'))
    )
    print(combined.gp)
    dev.off()
}
