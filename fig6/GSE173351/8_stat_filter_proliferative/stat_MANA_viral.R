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
    library(pheatmap)})

midRed <- 1.2
threshold <- 5

figure_path <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE173351/8_stat_filter_proliferative/")
if (!dir.exists(figure_path)) {
  dir.create(figure_path, recursive = T)
}
setwd(figure_path)

CD8PredictedObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE173351/7_mapping_filter_proliferative/CD8/querySeuratObj_2022-11-01.rds")
CD4PredictedObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE173351/7_mapping_filter_proliferative/CD4/querySeuratObj_2022-11-01.rds")
PPredictedObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE173351/7_mapping_filter_proliferative/P/querySeuratObj_2022-11-01.rds")

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

P_meta <- PPredictedObj@meta.data
P_meta$Patient <- stringr::str_extract(P_meta$orig.ident, "^.*-[^_]+(?=_.*)")
P_meta$BarcodeInOrigIdent <- stringr::str_extract(pattern = "[A-Z]+-[0-9]+", rownames(P_meta))
P_meta$isResponse = "NR"
P_meta$isResponse[P_meta$response %in% "MPR"] = "R"
P_meta$type = ""
P_meta$cdr3 = ""

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
    P_meta$cdr3[rownames(P_meta) %in% fcaT$barcodeID] <- fcaT$cdr3[match(rownames(P_meta)[rownames(P_meta) %in% fcaT$barcodeID], fcaT$barcodeID)]
}

for (tempPatient in unique(TCR_annotation_T$Patient)) {
    print(tempPatient)
    tat <- TCR_annotation_T %>%
        filter(Patient == tempPatient)
    for (i in 1:dim(tat)[1]) {
        CD4_meta$type[CD4_meta$Patient == tempPatient &
                      CD4_meta$cdr3 == tat$TCR[i]] <- tat$Type[i]
        CD8_meta$type[CD8_meta$Patient == tempPatient &
                      CD8_meta$cdr3 == tat$TCR[i]] <- tat$Type[i]
        P_meta$type[P_meta$Patient == tempPatient &
                      P_meta$cdr3 == tat$TCR[i]] <- tat$Type[i]
    }
}

CD8_meta$type <- plyr::mapvalues(CD8_meta$type,
                                        c("MANA", "Viral (CMV, EBV, Influenza A)", "Viral (Influenza)"),
                                        c("MANA", "Viral", "Viral"))
CD8_CellType <- c("CD8_c0_t-Teff",
                  "CD8_c1_Tex",
                  "CD8_c2_Teff",
                  "CD8_c3_Tn",
                  "CD8_c4_Tstr",
                  "CD8_c5_Tisg",
                  "CD8_c6_Tcm",
                  "CD8_c7_p-Tex",
                  "CD8_c8_Teff_KLRG1",
                  "CD8_c9_Tsen",
                  "CD8_c10_Teff_CD244",
                  "CD8_c11_Teff_SEMA4A",
                  "CD8_c12_Trm",
                  "CD8_c13_Tn_TCF7")
CD8_meta$predicted.id <- plyr::mapvalues(CD8_meta$predicted.id,
                                         0:(length(CD8_CellType)-1),
                                         CD8_CellType)
CD8_meta$predicted.id <- droplevels(CD8_meta$predicted.id)

P_meta$type <- plyr::mapvalues(P_meta$type,
                                 c("MANA", "Viral (CMV, EBV, Influenza A)", "Viral (Influenza)"),
                                 c("MANA", "Viral", "Viral"))
P_meta$predicted.id <- plyr::mapvalues(P_meta$predicted.id,
                                       0:7,
                                       c("P_CD8",
                                         "P_CD4",
                                         "P_Else",
                                         "P_Else",
                                         "P_CD8",
                                         "P_CD8",
                                         "P_Treg",
                                         "P_CD8"))

CD8_meta <- CD8_meta %>%
    mutate(group = paste0(type, "-", isResponse))
P_meta <- P_meta %>%
  mutate(group = paste0(type, "-", isResponse))

allMD <- list(CD8 = CD8_meta, Proliferative = P_meta)
allT <- bind_rows(CD4_meta, CD8_meta)
allT <- bind_rows(allT, P_meta)

allTSampleCN <- allT %>%
    group_by(orig.ident) %>%
    count


###############################################################################
#                               get sample info                               #
###############################################################################

totalMD <- bind_rows(CD8_meta, P_meta)

sampleMD <- totalMD %>%
        filter(tissue == "tumor") %>%
        filter(type != "") %>%
  mutate(group = paste0(type, "-", isResponse)) %>%
  group_by(group,
           orig.ident) %>%
  dplyr::count() %>%
  ungroup() %>%
  group_by(group) %>%
  dplyr::count()
write_tsv(sampleMD, file.path(figurePath,
                              paste0('NSCLC_caushi_sampleNum',
                                     "_", Sys.Date(), '.tsv')))

groupCellNum <- totalMD %>%
        filter(tissue == "tumor") %>%
        filter(type != "") %>%
  mutate(group = paste0(type, "-", isResponse)) %>%
  group_by(group) %>%
  dplyr::count()
write_tsv(groupCellNum, file.path(figurePath,
                                  paste0('NSCLC_caushi_cellNum',
                                         "_", Sys.Date(), '.tsv')))

###############################################################################
#                            finish get sample info                           #
###############################################################################


## targetTreatments = c("lpilimumab+Nivolumab", "lpilimumab")

rns <- unique(totalMD$predicted.id)
matRoe <- matrix(rep(0, length(rns) * 4), nrow = length(rns), ncol = 4)
matCell <- matrix(rep(0, length(rns) * 4), nrow = length(rns), ncol = 4)
rownames(matRoe) <- sort(rns)
colnames(matRoe) <- c("MANA-NR", "MANA-R", "Viral-NR", "Viral-R")
rownames(matCell) <- sort(rns)
colnames(matCell) <- c("MANA-NR", "MANA-R", "Viral-NR", "Viral-R")


for(ami in 1:length(allMD)) {
    tempName <- names(allMD[ami])
    tempMD <- allMD[[ami]] %>%
        filter(tissue == "tumor") %>%
        filter(type != "") %>%
        mutate(group = paste0(type, "-", isResponse))
    TotalSampleCellNum <- allTSampleCN %>%
        filter(orig.ident %in% tempMD$orig.ident)
    ClusterT <- tempMD %>%
        group_by(predicted.id) %>%
        dplyr::count()
    N <- dim(tempMD)[1]

    for(tempGroup in unique(tempMD$group)) {
        targetTcellMD <- tempMD %>% filter(group == tempGroup)
        M <- dim(targetTcellMD)[1]
        for(gsci in unique(tempMD$predicted.id)) {
            k <- ClusterT$n[ClusterT$predicted.id == gsci]
            if(gsci %in% targetTcellMD$predicted.id) {
                n <- dim(targetTcellMD %>% filter(predicted.id == gsci))[1]
            } else {
                n <- 0
            }
            tempData <- tibble(group = tempGroup,
                               Cluster = gsci,
                               Roe = (n/M) / (k/N))
            tempRowName <- gsci
            tempColName <- tempGroup
            matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
            matCell[tempRowName, tempColName] <- n
        }
    }
}


tMR <- matRoe
tMC <- matCell
tMR[tMR > threshold] = threshold
my.breaks.low <- seq(0, 1, by=0.01)
my.colors.low <- colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks.low))
my.breaks.mid <- seq(1.01, midRed, by=0.01)
my.colors.mid <- colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks.mid))
my.breaks.up <- seq(midRed + 0.01, threshold, by = (threshold - midRed)/99)
my.colors.up <- colorRampPalette(colors = c("#FD9AA0", "#550000"))(length(my.breaks.up))
my.breaks <- c(my.breaks.low, my.breaks.mid, my.breaks.up)
my.colors <- c(my.colors.low, my.colors.mid, my.colors.up)
widthi = ncol(tMR)/12 + 4
heighti = nrow(tMR)/7 + 3
pdf(file.path(figure_path, paste0("figure6k_Roe.pdf")), width = widthi, height = heighti)
print(pheatmap(tMR,
               color = my.colors,
               breaks = my.breaks,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figure_path, paste0("figure6k_Cell.pdf")), width = widthi, height = heighti)
print(pheatmap(tMC,
               show_colnames = T,
               show_rownames = T,
               display_numbers = T,
               number_format = "%.0f",
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()


tMR <- matRoe
tMC <- matCell
tMR <- tMR[rowSums(tMC) > 30,]
tMC <- tMC[rowSums(tMC) > 30,]
tMR[tMR > threshold] = threshold
tMR <- tMR[c("CD8_c4_Tstr", "CD8_c1_Tex", "CD8_c5_Tisg", "CD8_c2_Teff", "CD8_c3_Tn", "CD8_c0_t-Teff", "P_CD8"),]
tMC <- tMC[c("CD8_c4_Tstr", "CD8_c1_Tex", "CD8_c5_Tisg", "CD8_c2_Teff", "CD8_c3_Tn", "CD8_c0_t-Teff", "P_CD8"),]
my.breaks.low <- seq(0, 1, by=0.01)
my.colors.low <- colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks.low))
my.breaks.mid <- seq(1.01, midRed, by=0.01)
my.colors.mid <- colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks.mid))
my.breaks.up <- seq(midRed + 0.01, threshold, by = (threshold - midRed)/99)
my.colors.up <- colorRampPalette(colors = c("#FD9AA0", "#550000"))(length(my.breaks.up))
my.breaks <- c(my.breaks.low, my.breaks.mid, my.breaks.up)
my.colors <- c(my.colors.low, my.colors.mid, my.colors.up)
widthi = ncol(tMR)/12 + 4
heighti = nrow(tMR)/7 + 3
pdf(file.path(figure_path, paste0("figure6k_Roe.pdf")), width = widthi, height = heighti)
print(pheatmap(tMR,
               color = my.colors,
               breaks = my.breaks,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figure_path, paste0("figure6k_Cell.pdf")), width = widthi, height = heighti)
print(pheatmap(tMC,
               show_colnames = T,
               show_rownames = T,
               display_numbers = T,
               number_format = "%.0f",
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
