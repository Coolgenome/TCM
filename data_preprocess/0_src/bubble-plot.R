#' filename : bubble-plot.R
#' Date : 2020-06-22
#' contributor : Yanshuo Chu
#' function: bubble-plot

suppressMessages({library(optparse)
library(rjson)
library(rlist)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(viridis)
library(stringr)
})
print('---visualize embeding---')
##CLI parsing
option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                help = "r data file input(after runtsne/runumap)",
                metavar = 'character'),
    make_option(c("-o",'--out'),
                type = 'character',
                default = 'visualization.pdf',
                help = 'output file name for the plot [default = %default]',
                metavar = 'character'),
    make_option(c("-m","--markers"),
                type = 'character',
                help = 'markers table file (by findallmarkers)',
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
## seuratObj <- readRDS("")

Idents(seuratObj) <- seuratObj@meta.data$seurat_clusters

markerList <- read_tsv(opt$markers)

## markerList <- read_tsv("/rsrch3/home/genomic_med/ychu2/configs/scSeqs/database/Markers/CD4/markers/CD4_naive_clusters_comparison.txt")
colnames(markerList)[1:2] <- c("marker", "celltype")

if(any(c("Cd3d", "Cd79a", "Cd79b", "S100a8") %in% rownames(seuratObj))){
    firstup <- function(x) {
        substr(x, 1, 1) <- toupper(substr(x, 1, 1))
        substr(x, 2, nchar(x)) <- tolower(substr(x, 2, nchar(x)))
        x
    }
    markerList$marker <- firstup(markerList$marker)
}

targetMarkers <- intersect(markerList$marker, rownames(seuratObj))
markerList <- markerList %>% filter(marker %in% targetMarkers)

genes.list <- list()

if(length(unique(markerList$celltype)) > 1) {
    for(ct in unique(markerList$celltype)){
        genes.list[[str_replace(ct, " ", "-")]] <- markerList$marker[markerList$celltype == ct]
    }
}

totalName <- tools::file_path_sans_ext(basename(opt$markers))
## totalName <- tools::file_path_sans_ext("/rsrch3/home/genomic_med/ychu2/configs/scSeqs/database/CD8Anne.txt")

genes.list[[totalName]] <- markerList$marker

for(i in 1:length(genes.list)){
    gene.list.name <- names(genes.list[i])
    gene.list <- genes.list[[i]]

    if(length(gene.list) >= 50) {
        partN <- round(length(gene.list) / 25)
        cellN <- ceiling(length(gene.list) / partN)
        for(pni in 1:partN){
            startIdx <- (pni - 1) * cellN + 1
            endIdx <- pni * cellN
            if(endIdx > length(gene.list)){
                endIdx <- length(gene.list)
            }

            print(paste0("length(gene.list): ", length(gene.list)))
            print("startIdx")
            print(startIdx)
            print("endIdx")
            print(endIdx)

            temp.gene.list <- gene.list[startIdx:endIdx]

            gene <- intersect(temp.gene.list,rownames(seuratObj))
            p <- DotPlot(seuratObj, features = gene)
            data <- p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]
                                        #data[data$avg.exp>7,'avg.exp']<-7
            plotx <- ggplot(data, aes(x = features.plot,y = id)) +        ## global aes
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
                file.path(opt$out,
                          paste0(tools::file_path_sans_ext(basename(opt$markers)), '.',
                                 gene.list.name, '.part', pni, '.marker.bubble.pdf')),
                plotx, height=ht, width=wh, limitsize = F)
        }
    }

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
        file.path(opt$out,
                  paste0(tools::file_path_sans_ext(basename(opt$markers)), '.', gene.list.name, '.marker.bubble.pdf')),
        plotx, height=ht, width=wh, limitsize = F)
}

