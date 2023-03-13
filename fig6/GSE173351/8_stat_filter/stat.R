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
    library(stringr)
    library(pheatmap)
})
figure_path <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE173351/8_stat_filter/")
if (!dir.exists(figure_path)) {
  dir.create(figure_path, recursive = T)
}
setwd(figure_path)


CD4PredictedObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE173351/7_mapping_filter/CD4/querySeuratObj_2022-05-25.rds")
CD8PredictedObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE173351/7_mapping_filter/CD8/querySeuratObj_2022-05-25.rds")

Idents(CD8PredictedObj) <- CD8PredictedObj$predicted.id
pdf(file.path(getwd(), "bubbleplot_prof.pdf"))
DotPlot(CD8PredictedObj, features = c("MKI67", "TOP2A"))
dev.off()




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



allMD <- list(CD8 = CD8_meta, CD4 = CD4_meta)
allT <- bind_rows(CD4_meta, CD8_meta)

allTSampleCN <- allT %>%
    group_by(orig.ident) %>%
    count

## targetTreatments = c("lpilimumab+Nivolumab", "lpilimumab")
targetTissues = unique(allT$tissue)
targetCellSet = c("AllT", "Split")

for(tCS in targetCellSet){
    for(tT in targetTissues){
        for(ami in 1:length(allMD)){
            tempName <- names(allMD[ami])
            tempMD <- allMD[[ami]] %>% filter(tissue == tT)
            if(tCS == "Split"){
                TotalSampleCellNum <- tempMD %>%
                    group_by(orig.ident) %>%
                    count()
            }else{
                TotalSampleCellNum <- allTSampleCN %>%
                    filter(orig.ident %in% tempMD$orig.ident)
            }
            totalT <- c()
            for(tempCluster in as.numeric(unique(tempMD$predicted.id))){
                TNR <- tempMD %>%
                    filter(predicted.id == tempCluster) %>%
                    group_by(orig.ident, isResponse) %>%
                    count()
                TNR$Frac <- 0.0
                TNR$Frac <- TNR$n / TotalSampleCellNum$n[match(TNR$orig.ident, TotalSampleCellNum$orig.ident)]

                TNR$cluster <- tempCluster
                totalT <- bind_rows(totalT, TNR)
                missedSamples <- setdiff(
                    unique(TotalSampleCellNum$orig.ident),
                    unique(TNR$orig.ident))

                if(length(missedSamples) >0){
                    isResponse = rep("No", length(missedSamples))
                    isResponse[missedSamples %in% c("MPR")] = "Yes"
                    tempTNR <- tibble(
                        orig.ident = missedSamples,
                        isResponse = isResponse,
                        Frac = 0,
                        cluster = tempCluster
                    )
                    totalT <- bind_rows(totalT, tempTNR)
                }
            }
            g <- totalT %>%
                ggstatsplot::grouped_ggbetweenstats(
                                 data = .,
                                 x = isResponse,
                                 y = Frac,
                                 grouping.var = cluster,
                                 xlab = "",
                                 ylab = "Sample fraction",
                                 ## pairwise.display = "aiwl", # display only significant pairwise comparisons
                                 p.adjust.method = "fdr", # adjust p-values for multiple tests using this method
                                 ggtheme = theme_classic(),
                                 package = "ggsci",
                                 palette = "default_jco",
                                 plotgrid.args = list(ncol = 1))
            figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE173351/7_mapping_stat_filter",
                                    tempName)
            if(!dir.exists(figurePath)){
                dir.create(figurePath, recursive = T)
            }
            ggsave(file.path(figurePath, paste0(tempName, "-", tCS, "_Tissue-", tT,"_response_bar.pdf")), g, width = 200, height = 1200, units = "mm")
            ClusterT <- tempMD %>%
                group_by(predicted.id) %>%
                dplyr::count()
            N <- dim(tempMD)[1]
            YesMD <- tempMD %>%
                filter(isResponse == "Yes") %>%
                rename(Group = isResponse)
            NoMD <- tempMD %>%
                filter(isResponse == "No") %>%
                rename(Group = isResponse)
            totalTargetTcellMD <- bind_rows(YesMD, NoMD)
            matRoe <- matrix(rep(0, length(unique(tempMD$predicted.id)) * 2),
                             nrow = length(unique(tempMD$predicted.id)),
                             ncol = 2)
            matCell <- matrix(rep(0, length(unique(tempMD$predicted.id)) * 2),
                              nrow = length(unique(tempMD$predicted.id)),
                              ncol = 2)
            rownames(matRoe) <- sort(unique(tempMD$predicted.id))
            colnames(matRoe) <- c("Yes", "No")
            rownames(matCell) <- sort(unique(tempMD$predicted.id))
            colnames(matCell) <- c("Yes", "No")
            ggplotdata <- c()
            for(tempGroup in c("Yes", "No")){
                targetTcellMD <- totalTargetTcellMD %>% filter(Group == tempGroup)
                M <- dim(targetTcellMD)[1]
                for(gsci in unique(tempMD$predicted.id)){
                    k <- ClusterT$n[ClusterT$predicted.id == gsci]
                    if(gsci %in% targetTcellMD$predicted.id){
                        n <- dim(targetTcellMD %>% filter(predicted.id == gsci))[1]
                    } else{
                        n <- 0
                    }
                    tempData <- tibble(isResponse = tempGroup,
                                       Cluster = gsci,
                                       Roe = (n/M) / (k/N))
                    ggplotdata <- bind_rows(ggplotdata, tempData)
                    tempRowName <- as.character(gsci)
                    tempColName <- tempGroup
                    matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
                    matCell[tempRowName, tempColName] <- n
                }
            }
            matRoe[matRoe>10] = 10
            ## colnames(matRoe) <- paste0(colnames(matRoe), " (", colSums(matCell)[colnames(matRoe)], ")")
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
                           cluster_cols = T,
                           cluster_rows = F,
                           number_format = "%.0f",
                           border_color = F))
            dev.off()
        }
    }
}
