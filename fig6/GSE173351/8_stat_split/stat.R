#--------------------------------------------------------------
# filename : stat.R
# Date : 2022-03-04
# contributor : Yanshuo Chu
# function: stat
#--------------------------------------------------------------

print('<==== stat.R ====>')
rm(list=ls())

suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(ggplot2)
    library(ggpubr)
    library(stringr)
    library(pheatmap)
})

figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE173351/8_stat_split/outs")
if(!dir.exists(figurePath)){
    dir.create(figurePath, recursive = T)
}
setwd(figurePath)

CD4PredictedObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE173351/7_mapping/CD4/querySeuratObj_2022-05-14.rds")
CD8PredictedObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE173351/7_mapping/CD8/querySeuratObj_2022-05-14.rds")



CD4_meta <- CD4PredictedObj@meta.data
CD4_meta$Patient <- stringr::str_extract(CD4_meta$orig.ident, "^.*-[^_]+(?=_.*)")
CD4_meta$BarcodeInOrigIdent <- stringr::str_extract(pattern = "[A-Z]+-[0-9]+", rownames(CD4_meta))
CD4_meta$isResponse = "NR"
CD4_meta$isResponse[CD4_meta$response %in% "MPR"] = "R"
CD4_meta$type = ""
CD4_meta$cdr3 = ""

CD8_meta <- CD8PredictedObj@meta.data
CD8_meta$Patient <- stringr::str_extract(CD8_meta$orig.ident, "^.*-[^_]+(?=_.*)")
CD8_meta$BarcodeInOrigIdent <- stringr::str_extract(pattern = "[A-Z]+-[0-9]+", rownames(CD8_meta))
CD8_meta$isResponse = "NR"
CD8_meta$isResponse[CD8_meta$response %in% "MPR"] = "R"
CD8_meta$type = ""
CD8_meta$cdr3 = ""

###############################################################################
#                             load TCR annotation                             #
###############################################################################
TCR_annotation_T <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/data/GSE173351/TCR_MANA_VIRAL.txt")

TCR_Contig_Path <- "/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/data/GSE173351/GSE176021/TCR"
print(paste0("Loading data from ", TCR_Contig_Path))
samplesDir = list.dirs(path = TCR_Contig_Path, recursive = F)
samplesName = basename(samplesDir)
for (ids in seq_along(samplesName)) {
    ds = samplesName[ids]
    print(ds)
    tempPatient <- stringr::str_extract(ds, "^.*-[^_]+(?=_.*)")
    tempOrigIdent <- ds
    fcaT = read_csv(file.path(TCR_Contig_Path, ds, "filtered_contig_annotations.csv")) %>%
        select(barcode, cdr3) %>%
        mutate(barcodeID = paste0(tempOrigIdent, "_", barcode))
    print(paste(ds, ":", dim(fcaT)[1]))
    CD4_meta$cdr3[rownames(CD4_meta) %in% fcaT$barcodeID] <- fcaT$cdr3[match(rownames(CD4_meta)[rownames(CD4_meta) %in% fcaT$barcodeID], fcaT$barcodeID)]
    CD8_meta$cdr3[rownames(CD8_meta) %in% fcaT$barcodeID] <- fcaT$cdr3[match(rownames(CD8_meta)[rownames(CD8_meta) %in% fcaT$barcodeID], fcaT$barcodeID)]
}

for(tempPatient in unique(TCR_annotation_T$Patient)){
    print(tempPatient)
    tat <- TCR_annotation_T %>%
        filter(Patient == tempPatient)
    for(i in 1:dim(tat)[1]){
        CD4_meta$type[CD4_meta$Patient == tempPatient &
                      CD4_meta$cdr3 == tat$TCR[i]] <- tat$Type[i]
        CD8_meta$type[CD8_meta$Patient == tempPatient &
                      CD8_meta$cdr3 == tat$TCR[i]] <- tat$Type[i]
    }
}

CD8_meta$type <- plyr::mapvalues(CD8_meta$type,
                                        c("MANA", "Viral (CMV, EBV, Influenza A)", "Viral (Influenza)"),
                                        c("MANA", "Viral", "Viral"))
CD8_meta <- CD8_meta %>%
    mutate(group = paste0(type, "-", isResponse))

allMD <- list(CD8 = CD8_meta)
allT <- bind_rows(CD4_meta, CD8_meta)

allTSampleCN <- allT %>%
    group_by(orig.ident) %>%
    count

## targetTreatments = c("lpilimumab+Nivolumab", "lpilimumab")
targetTissues = unique(allT$tissue)
targetCellSet = c("AllT", "Split")

for(tCS in targetCellSet){
    for(tT in c("tumor")){
        for(ami in 1:length(allMD)){
            tempName <- names(allMD[ami])
            tempMD <- allMD[[ami]] %>%
                filter(tissue == tT) %>%
                filter(type != "") %>%
                mutate(group = paste0(type, "-", isResponse))
            if(tCS == "Split"){
                TotalSampleCellNum <- tempMD %>%
                    group_by(orig.ident) %>%
                    count()
            }else{
                TotalSampleCellNum <- allTSampleCN %>%
                    filter(orig.ident %in% tempMD$orig.ident)
            }
            ClusterT <- tempMD %>%
                group_by(predicted.id) %>%
                dplyr::count()
            N <- dim(tempMD)[1]
            matRoe <- matrix(rep(0, length(unique(tempMD$predicted.id)) * 4),
                             nrow = length(unique(tempMD$predicted.id)),
                             ncol = 4)
            matCell <- matrix(rep(0, length(unique(tempMD$predicted.id)) * 4),
                              nrow = length(unique(tempMD$predicted.id)),
                              ncol = 4)
            rownames(matRoe) <- sort(unique(tempMD$predicted.id))
            colnames(matRoe) <- unique(tempMD$group)
            rownames(matCell) <- sort(unique(tempMD$predicted.id))
            colnames(matCell) <- unique(tempMD$group)
            ggplotdata <- c()
            for(tempGroup in unique(tempMD$group)){
                targetTcellMD <- tempMD %>% filter(group == tempGroup)
                M <- dim(targetTcellMD)[1]
                for(gsci in unique(tempMD$predicted.id)){
                    k <- ClusterT$n[ClusterT$predicted.id == gsci]
                    if(gsci %in% targetTcellMD$predicted.id){
                        n <- dim(targetTcellMD %>% filter(predicted.id == gsci))[1]
                    } else{
                        n <- 0
                    }
                    tempData <- tibble(group = tempGroup,
                                       Cluster = gsci,
                                       Roe = (n/M) / (k/N))
                    ggplotdata <- bind_rows(ggplotdata, tempData)
                    tempRowName <- as.character(gsci)
                    tempColName <- tempGroup
                    matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
                    matCell[tempRowName, tempColName] <- n
                }
            }
            matCell <- matCell[,c("MANA-NR",
                                  "MANA-R",
                                  "Viral-NR",
                                  "Viral-R")]
            matRoe <- matRoe[,c("MANA-NR",
                                "MANA-R",
                                "Viral-NR",
                                "Viral-R")]
            matRoe = matRoe[rowSums(matCell) > 100,]
            matRoe[matRoe>10] = 10

            ## my.breaks <- seq(0, 0.99, by=0.01)
            ## my.colors <- c(colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)))
            ## my.breaks.extra <- seq(1, 1.99, by = (1.99 - 1)/99)
            ## my.colors.extra <- colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks.extra))
            ## my.breaks.extra2 <- seq(2, 10, by = (10 - 2)/99)
            ## my.colors.extra2 <- colorRampPalette(colors = c("#FD9AA0", "red"))(length(my.breaks.extra2))
            ## my.breaks <- c(my.breaks, my.breaks.extra, my.breaks.extra2)
            ## my.colors <- c(my.colors, my.colors.extra, my.colors.extra2)


            my.breaks <- seq(0, 1.99, by=0.01)
            my.colors <- c(
                colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)/2),
                colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks)/2))
            my.breaks.extra <- seq(2, 10, by = (10 - 2)/99)
            my.colors.extra <- colorRampPalette(colors = c("#FD9AA0", "#550000"))(length(my.breaks.extra))
            my.breaks <- c(my.breaks, my.breaks.extra)
            my.colors <- c(my.colors, my.colors.extra)

            widthi = ncol(matRoe)/12 + 1.8
            heighti = nrow(matRoe)/7 + 2
            print(rownames(matRoe))
            pdf(file.path(figurePath, paste0(tempName, "-", tCS, "_Tissue-", tT, "_response_roe_heatmap.pdf")), width = widthi, height = heighti)
            print(pheatmap(matRoe,
                           color = my.colors,
                           breaks = my.breaks,
                           show_colnames = T,
                           show_rownames = T,
                           cluster_cols = F,
                           cluster_rows = F,
                           border_color = F))
            dev.off()
            pdf(file.path(figurePath, paste0(tempName, "-", tCS, "_Tissue-", tT, "_response_cell_heatmap.pdf")), width = widthi * 1.5, height = heighti)
            print(pheatmap(matCell,
                           display_numbers = T,
                           show_colnames = T,
                           show_rownames = T,
                           cluster_cols = F,
                           cluster_rows = F,
                           number_format = "%.0f",
                           border_color = F))
            dev.off()
        }
    }
}

CD8_meta <- CD8_meta %>%
    mutate(group = paste0(type, "-", isResponse))

CD8PredictedObj$group <- CD8_meta$group[match(Cells(CD8PredictedObj), rownames(CD8_meta))]
CD8PredictedObj$type <- CD8_meta$type[match(Cells(CD8PredictedObj), rownames(CD8_meta))]

groups <- c(
    "MANA-NR",
    "MANA-R",
    "Viral-NR",
    "Viral-R")
CD8PredictedObj$group <- factor(CD8PredictedObj$group, levels = groups)

types <- c(
    "MANA",
    "Viral")
CD8PredictedObj$type <- factor(CD8PredictedObj$type, levels = types)


groups <- c("MANA-NR", "MANA-R", "Viral-NR", "Viral-R")
CD8PredictedObj$group <- factor(CD8PredictedObj$group, levels = groups)

types <- c("MANA", "Viral")
CD8PredictedObj$type <- factor(CD8PredictedObj$type, levels = types)

hs_genes <- c("HSPA1A", "HSPA1B")

Idents(CD8PredictedObj) <- CD8PredictedObj$group
tempObj <- subset(CD8PredictedObj, idents = groups)

## pdf(file.path(getwd(), "CD8_HSPA1A_violinplot.pdf"), width = 5, height = 4)
## VlnPlot(tempObj, features = "HSPA1A", pt.size = 0) +
##     geom_boxplot(color = "black", fill = NA, outlier.shape = NA, width = 0.2)
## dev.off()
## pdf(file.path(getwd(), "CD8_HSPA1B_violinplot.pdf"), width = 5, height = 4)
## VlnPlot(tempObj, features = "HSPA1B", pt.size = 0) +
##     geom_boxplot(color = "black", fill = NA, outlier.shape = NA, width = 0.2)
## dev.off()


Idents(tempObj) <- tempObj$type
pdf(file.path(getwd(), "CD8_HSPA1A_type_violinplot.pdf"), width = 5, height = 4)
VlnPlot(tempObj, features = "HSPA1A", pt.size = 0) +
    geom_boxplot(color = "black", fill = NA, outlier.shape = NA, width = 0.2) +
    stat_compare_means()
dev.off()
pdf(file.path(getwd(), "CD8_HSPA1B_type_violinplot.pdf"), width = 5, height = 4)
VlnPlot(tempObj, features = "HSPA1B", pt.size = 0) +
    geom_boxplot(color = "black", fill = NA, outlier.shape = NA, width = 0.2) +
    stat_compare_means()
dev.off()



pdf(file.path(getwd(), "CD8_groups_ridgeplot.pdf"), width = 8, height = 4)
RidgePlot(tempObj, features = hs_genes)
dev.off()

targetMD <- tempObj@meta.data
targetMD$HSPA1A <- 0
targetMD$HSPA1A <- tempObj@assays$RNA@data["HSPA1A",]
targetMD$HSPA1B <- 0
targetMD$HSPA1B <- tempObj@assays$RNA@data["HSPA1B",]

targetMD <- as_tibble(targetMD) %>%
    select(group, HSPA1A, HSPA1B) %>%
    gather(key = "Marker", value = "Expression", -group)

## g <- ggstatsplot::ggbetweenstats(
##                       data = targetMD %>% filter(Marker == "HSPA1A"),
##                       pairwise.display = "all",
##                       p.adjust.method = "fdr",
##                       x = group,
##                       y = Expression)
## ggsave(file.path(figurePath, paste0("CD8_HSPA1A_group_ggbetweenstats.pdf")), g, width = 210, height = 210, units = "mm")


g <- targetMD %>%
    ggplot() +
    geom_violin(aes(x = group, y = Expression, fill = group), lwd = 0.1) +
    geom_boxplot(aes(x = group, y = Expression),
                 color = "black",
                 fill = NA,
                 width = 0.1,
                 lwd = 0.5,
                 outlier.shape = NA) +
    facet_grid(. ~ Marker) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file.path(getwd(), "CD8_groups_violinplot.pdf"), g, width = 6, height = 4)




g <- ggstatsplot::ggbetweenstats(
                      data = targetMD %>% filter(Marker == "HSPA1B"),
                      pairwise.display = "all",
                      p.adjust.method = "fdr",
                      x = group,
                      y = Expression)
ggsave(file.path(figurePath, paste0("CD8_HSPA1B_group_ggbetweenstats.pdf")), g, width = 210, height = 210, units = "mm")

## tempMeta <- as_tibble(tempObj@meta.data) %>%
##     group_by(group, predicted.id) %>%
##     count %>%
##     ungroup %>%
##     group_by(group) %>%
##     mutate(frac = n / sum(n)) %>%
##     filter(predicted.id == 4)

Idents(CD8PredictedObj) <- CD8PredictedObj$type
tempObj <- subset(CD8PredictedObj, idents = types)

pdf(file.path(getwd(), "CD8_types_violinplot.pdf"), width = 8, height = 4)
VlnPlot(tempObj, features = hs_genes, pt.size = 0)
dev.off()
