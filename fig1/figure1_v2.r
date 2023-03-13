#--------------------------------------------------------------
# filename : figure1_v2.r
# Date : 2022-10-20
# contributor : Yanshuo Chu
# function: figure1_v2
# R version: R/4.0.3
#--------------------------------------------------------------

print('<==== figure1_v2.r ====>')

rm(list=ls())

library(tidyverse)
library(Seurat)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(wesanderson)
library(viridis)
library(ggsci)
library(tidytext)
library(ggpubr)
library(cowplot)
library(facetscales)
library(latex2exp)
library(ggstatsplot)
set.seed(123)

CD8_path <- '/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/0_write_sample_info/CD8_2022-10-20.rds'
CD4_path <- '/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/0_write_sample_info/CD4_2022-10-20.rds'
Innate_path <- '/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/0_write_sample_info/Innate_2022-10-20.rds'
Proliferative_path <- '/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/0_write_sample_info/Proliferative_2022-10-20.rds'

CD8_Obj <- readRDS(CD8_path)
CD4_Obj <- readRDS(CD4_path)
Innate_Obj <- readRDS(Innate_path)
Proliferative_Obj <- readRDS(Proliferative_path)
allPath <- list(CD8 = CD8_path, CD4 = CD4_path, NK = Innate_path, Proliferative = Proliferative_path)
allObj <- list(CD8 = CD8_Obj, CD4 = CD4_Obj, NK = Innate_Obj, Proliferative = Proliferative_Obj)
common_columns <- c()
for(oi in 1:length(allObj)){
    tempName <- names(allObj[oi])
    tempObj <- allObj[[oi]]
    tempObj$barcode <- Cells(tempObj)
    if(oi == 1){
        common_columns <- colnames(tempObj@meta.data)
    }else{
        common_columns <- intersect(common_columns, colnames(tempObj@meta.data))
    }
}
allMD <- c()
for(oi in 1:length(allObj)){
    tempName <- names(allObj[oi])
    tempObj <- allObj[[oi]]
    tempMD <- as_tibble(tempObj@meta.data, rownames = NA) %>%
        dplyr::select(all_of(common_columns)) %>% mutate(barcode = Cells(tempObj))
    allMD <- bind_rows(allMD, tempMD)
}
colorsForDataType <- c("#6DCCDD", "#EDCAE0", "#F494BE", "#F9B26C", "#A6ADCC", "#C4DA5D")
colorForClass1 <- c("#C4DA5D", "#6DCCDD", "#F494BE", "#EDCAE0")
dataTypeLevel <- c("CD4", "CD8", "MAIT", "Tgd", "NKT", "Proliferative")

figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/1_main_fig1/")
if(!dir.exists(figurePath)){
  dir.create(figurePath, recursive = T)
}
setwd(figurePath)


## fig1a: statistics

allMD %>% dim()
allMD %>% group_by(Sample) %>%
    dplyr::count() %>% dim()
unique(allMD$Sample)
allMD %>% group_by(Patient) %>%
    dplyr::count() %>% dim()

t <- allMD %>% group_by(OrganSite) %>% dplyr::count()

t <- allMD %>%
    filter(TissueType %in% c("Primary tumor tissue", "Metastatic tumor tissue")) %>%
    group_by(CancerType) %>% dplyr::count()

table(allMD$cell.type)

###############################################################################
#                                sup figure 1b                                #
###############################################################################
allMD$DataSet <- plyr::mapvalues(allMD$cell.type, c( "CD8_c0_t-Teff",
"CD8_c1_Tex", "CD8_c2_Teff", "CD8_c3_Tn", "CD8_c4_Tstr", "CD8_c5_Tisg",
"CD8_c6_Tcm", "CD8_c7_p-Tex", "CD8_c8_Teff_KLRG1", "CD8_c9_Tsen",
"CD8_c10_Teff_CD244", "CD8_c11_Teff_SEMA4A", "CD8_c12_Trm", "CD8_c13_Tn_TCF7",
"CD4_c0_Tcm", "CD4_c1_Treg", "CD4_c2_Tn", "CD4_c3_Tfh", "CD4_c4_Tstr",
"CD4_c5_Tctl", "CD4_c6_Tn_FHIT", "CD4_c7_Tn_TCEA3", "CD4_c8_Th17",
"CD4_c9_Tn_TCF7_SLC40A1", "CD4_c10_Tn_LEF1_ANKRD55", "CD4_c11_Tisg",
"NKT_c0_FCGR3A_GZMB", "MAIT-like_c1", "MAIT_c2_KLRB1", "Tgd_c3_CX3CR1",
"NKT_c4_KIR_TIGIT_CXCR6", "P_c0_CD8_KLRD1_GZMB_CCL4/5", "P_c1_CD4_CD40LG_IL7R",
"P_c2_DNT", "P_c3_DNT_GZMK+", "P_c4_CD8_C1QBP_MT1G/X/E",
"P_c5_CD8_CXCL13_CCL3/4/5_PD-1", "P_c6_Treg", "P_c7_CD8_GZMK_NKG7_LAG3_PRF1"),
c( "CD8", "CD8", "CD8", "CD8", "CD8", "CD8", "CD8", "CD8", "CD8", "CD8", "CD8",
"CD8", "CD8", "CD8", "CD4", "CD4", "CD4", "CD4", "CD4", "CD4", "CD4", "CD4",
"CD4", "CD4", "CD4", "CD4", "NKT", "MAIT", "MAIT", "Tgd", "NKT", "P", "P", "P",
"Proliferative", "Proliferative", "Proliferative", "Proliferative", "Proliferative"))


ggdat <- allMD %>%
  group_by(DataSet,TissueType, CancerType) %>%
    dplyr::count()
ggdat$TissueType <- factor(ggdat$TissueType, levels = c("Healthy donor", "Uninvolved normal tissue", "Primary tumor tissue", "Metastatic tumor tissue"))
g <- ggdat %>% ggplot() +
  geom_bar(aes(x = CancerType, y = n, fill = DataSet), stat = "identity", position = "dodge") +
  scale_fill_manual(values=c(
                      "#ffcccc", #CD8
                      "#66cccc", #CD4
                      "#9999cc", #NKT
                      "#ff99cc", #MAIT
                      "#ff9966", #Tgd
                      "#cccc66" #P
                    )) +
  scale_y_continuous(trans='log10') +
  facet_grid(. ~ TissueType, scales = "free", space = "free") +
  ylab("Cell number") +
  theme_classic() +
  theme(
    text = element_text(size = 10),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
    legend.direction = "vertical",
    strip.background = element_blank(),
    legend.key.size = unit(2.5, "mm"))
ggsave(file.path(paste0("figureS1b.pdf")), g,  width = 210, height = 297/4, units = "mm")

###############################################################################
#'                                    manuscript: fig 1b                     '#
###############################################################################

tempAllMD <- allMD
tempAllMD$OrganSite[tempAllMD$OrganSite %in% c("Kidney", "Spleen", "Muscular", "Pelvis")] <- "Others"
ggplotdata_subject <- tempAllMD %>%  group_by(OrganSite, Patient) %>% dplyr::count() %>% group_by(OrganSite) %>% dplyr::count() %>% rename(Subject = n)
ggplotdata_sample <- tempAllMD %>%  group_by(OrganSite, Sample) %>% dplyr::count() %>% group_by(OrganSite) %>% dplyr::count() %>% rename(Sample = n)
ggplotdata_cell <- tempAllMD %>%  group_by(OrganSite) %>% dplyr::count() %>% rename(Cell = n) %>% arrange(Cell)
ggplotdata <- full_join(ggplotdata_subject, ggplotdata_sample, by = c("OrganSite"))
ggplotdata <- full_join(ggplotdata, ggplotdata_cell, by = c("OrganSite"))
ggplotdata <- ggplotdata %>% gather(key = "Type", value = "Number", -OrganSite )
ggplotdata$Type <- factor(ggplotdata$Type, levels = c("Cell", "Sample","Subject"))
ggplotdata$OrganSite <- factor(ggplotdata$OrganSite, levels = ggplotdata_cell %>% pull(OrganSite))
scales_x <- list(
    `Cell` = scale_x_log10(limits = c(100, 100000), breaks = c(100, 1000, 10000, 100000), label = c(TeX("10^2"), TeX("10^3"), TeX("10^4"), TeX("10^5"))),
    `Sample` = scale_x_log10(limits = c(1, 100), breaks = c(1, 10, 100), label = c(TeX("1"), TeX("10"), TeX("10^2"))),
    `Subject` = scale_x_log10(limits = c(1, 100), breaks = c(1, 10, 100), label = c(TeX("1"), TeX("10"), TeX("10^2"))))
g_Organ_figure1B <- ggplot(data = ggplotdata) +
    geom_bar(aes(x = Number, y = OrganSite),fill = "#66F0FF", stat = "identity") +
    xlab("") + ylab("") +
    facet_grid_sc(cols = vars(Type), scales = list(x = scales_x)) +
    theme_classic() +
    theme(
        text = element_text(size = 10),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
write_tsv(ggplotdata, file.path(figurePath, paste0('figure1B', "_", Sys.Date(), '.tsv')))
pdf(file.path(figurePath, paste0("Organ_figure1B.pdf")), height = 3, width = 5)
print(g_Organ_figure1B)
dev.off()
ggplotdata <- allMD %>% filter(!OrganSite %in% c("Spleen", "Muscular")) %>% group_by(OrganSite, TissueType) %>% dplyr::count() %>% group_by(OrganSite) %>% mutate(frac = n / sum(n))
ggplotdata$ForOrder <- ggplotdata$frac
for(oi in unique(ggplotdata$OrganSite)){
    for(ci in unique(ggplotdata$TissueType[ggplotdata$OrganSite == oi])){
        if("Primary tumor tissue" %in% unique(ggplotdata$TissueType[ggplotdata$OrganSite == oi])){
            ggplotdata$ForOrder[ggplotdata$TissueType == ci & ggplotdata$OrganSite == oi] <- ggplotdata$frac[ggplotdata$OrganSite == oi & ggplotdata$TissueType == "Primary tumor tissue"]
        }else{
            ggplotdata$ForOrder[ggplotdata$TissueType == ci & ggplotdata$OrganSite == oi] <- 0
        }
    }
}
ggplotdata$TissueType <- factor(ggplotdata$TissueType, levels = c("Healthy donor", "Uninvolved normal tissue", "Metastatic tumor tissue", "Primary tumor tissue"))
ggplotdata <- ggplotdata %>% rename(Type = TissueType)
g_Organ_Class1_figure1B <- ggplot(data =
                                      ## ggplotdata %>% filter(!(
                                      ##     Type == "Uninvolved normal tissue" & OrganSite == "Pancreas"))
                                      ggplotdata) +
    geom_bar(aes(x = reorder(OrganSite, ForOrder), y = frac, fill = Type), stat = "identity") +
    scale_fill_manual(values=colorForClass1[c(1,2,4,3)]) +
    coord_flip() +
    ylab("Cell fraction") + xlab("") +
    theme_classic() +
    theme(
        text = element_text(size = 10),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        legend.direction = "vertical",
        legend.key.size = unit(2.5, "mm")
    ) +
    guides(fill = guide_legend(ncol = 1, override.aes = list(size = 0.05)))
## pdf(file.path(figurePath, paste0("Organ_Class1_figure1B.pdf")), height = 5, width = 3)
## print(g)
## dev.off()
gn <- ggplotdata %>% select(OrganSite, Type, n) %>% spread(Type, n) %>% replace(is.na(.), 0)
write_tsv(gn, file.path(figurePath, paste0('figure1C_n', "_", Sys.Date(), '.tsv')))
gf <- ggplotdata %>% select(OrganSite, Type, frac) %>% spread(Type, frac) %>% replace(is.na(.), 0)
write_tsv(gf, file.path(figurePath, paste0('figure1C_frac', "_", Sys.Date(), '.tsv')))
g <- plot_grid(
    g_Organ_figure1B + theme(legend.position = "none"),
    g_Organ_Class1_figure1B,
    nrow = 1,
    axis = "bt",
    align = 'h',
    rel_widths = c(1, 1.7))
ggsave(file.path(figurePath, paste0("figure1B.pdf")), g,  width = 210/7 *5, height = 297/3.2, units = "mm")

###############################################################################
#'                        Manuscript:figure1D                                '#
###############################################################################

seuratObj<- readRDS("/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/whole_seuratObj_2021-03-03.rds")
Idents(seuratObj) <- seuratObj$DataSet

markers <- c("CD3D",
             "CD4",
             "CD40LG",
             "CD8A",
             "CD8B",
             "GZMK",
             "GZMB",
             "GZMH",
             "FCGR3A",
             "KLRD1",
             "KLRF1",
             "NCAM1",
             "TRAV1-2",
             "SLC4A10",
             "MKI67")
gene <- intersect(markers, rownames(seuratObj))
p <- DotPlot(seuratObj, features = gene)
data <- p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]
data$id <- factor(data$id, levels = rev(dataTypeLevel))
plotx <- ggplot(data, aes(x = features.plot, y = id)) +        ## global aes
    geom_point(aes(fill = avg.exp.scaled, size = pct.exp),
               color = 'black',
               shape = 21,
               stroke = 0.01) +
    scale_fill_gradient2(high = colorsForDataType[3], mid = "#F0D0D0", low = colorsForDataType[1])+
    scale_size(range = c(0,6), limits = c(0, 100), breaks = c(0,20,40,60,80,100)) +
    theme(
        text = element_text(size = 10),
        panel.grid.major = element_line(colour = "grey90", size=0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
        legend.position="bottom",
        legend.title=element_text(size=10)) +
    guides(size = guide_legend(title.position="top",
                               title.hjust = 0.5,
                               override.aes = list(stroke = 0.4)),
           fill = guide_colourbar(title.position="top", title.hjust = 0.5))
ggsave(file.path(figurePath, "whole_bubbleplot.pdf"), plotx, width = 210/2, height = 297/4.2, units = "mm")

###############################################################################
#'                     Manuscript:   figure1E
###############################################################################

tempAllMD <- allMD
tempAllMD$DataSet[tempAllMD$cell.type %in% c("MAIT-like_c1",
                                             "MAIT_c2_KLRB1")] <- "MAIT"
tempAllMD$DataSet[tempAllMD$cell.type %in% c("Tgd_c3_CX3CR1")] <- "Tgd"
tempAllMD$DataSet[tempAllMD$cell.type %in% c("NKT_c0_FCGR3A_GZMB", "NKT_c4_KIR_TIGIT_CXCR6")] <- "NKT"
ggplotdata <- tempAllMD %>% group_by(DataSet) %>% count() %>% mutate(frac = n / sum(ggplotdata$n))
write_tsv(ggplotdata, file.path(getwd(), paste0('figure1E', "_", Sys.Date(), '.tsv')))
colnTable <- as.vector(ggplotdata$n)
names(colnTable) <- c("CD4", "CD8", "MAIT", "NKT", "Proliferative", "Tgd")
## colnTable <- colnTable[dataTypeLevel]
cellNumber <- colnTable
colnTable <- colnTable / sum(colnTable) * 100
colnTable <- round(colnTable, digits = 2)
lbls <- paste(names(colnTable), "\n", cellNumber, ", ", colnTable, "%", sep="")
pdf(file.path(figurePath, paste0("DataSet_pieplot.pdf")), width = 5, height = 5)
pie(colnTable, labels = lbls, col = colorsForDataType, main=paste0("Total cell ", sum(cellNumber)), border = NA, clockwise = T, init.angle = 180)
dev.off()


###############################################################################
#'                      Manuscript: figure 1F                                '#
###############################################################################

tempAllMD <- allMD
tempAllMD$DataSet[tempAllMD$cell.type %in% c("MAIT-like_c1",
                                             "MAIT_c2_KLRB1")] <- "MAIT"
tempAllMD$DataSet[tempAllMD$cell.type %in% c("Tgd_c3_CX3CR1")] <- "Tgd"
tempAllMD$DataSet[tempAllMD$cell.type %in% c("NKT_c0_FCGR3A_GZMB", "NKT_c4_KIR_TIGIT_CXCR6")] <- "NKT"
ggplotdata <- c()
for(ci in unique(tempAllMD$TissueType)){
    if(ci %in% c("Healthy donor", "Uninvolved normal tissue")){
        ggplotdataTemp <- tempAllMD %>% filter(TissueType == ci) %>%
            mutate(DetailOrganSite = paste0(OrganSite)) %>%
            group_by(DetailOrganSite, DataSet) %>% dplyr::count() %>%
            group_by(DetailOrganSite) %>% mutate(frac = n / sum(n)) %>%
            select(DetailOrganSite, DataSet, frac, n) %>% mutate(Type = ci) %>%
            group_by(DetailOrganSite, Type) %>% mutate(ColSum = sum(n))
    }else{
        ggplotdataTemp <- tempAllMD %>% filter(TissueType == ci) %>%
            mutate(DetailOrganSite = paste0(CancerType, "_", OrganSite)) %>%
            group_by(DetailOrganSite, DataSet) %>% dplyr::count() %>%
            group_by(DetailOrganSite) %>% mutate(frac = n / sum(n)) %>%
            select(DetailOrganSite, DataSet, frac, n) %>% mutate(Type = ci) %>%
            group_by(DetailOrganSite, Type) %>% mutate(ColSum = sum(n))
    }
    ggplotdata <- bind_rows(ggplotdata, ggplotdataTemp)
}
origin.ggplotdata <- ggplotdata
ggplotdata <- ggplotdata %>% filter(ColSum >= 200)
ggplotdata$ForOrder <- ggplotdata$frac
for(fi in unique(ggplotdata$Type)){
    for(oi in unique(ggplotdata$DetailOrganSite[ggplotdata$Type == fi])){
        ggplotdata$ForOrder[ggplotdata$Type == fi & ggplotdata$DetailOrganSite == oi] <- ggplotdata$frac[ggplotdata$DetailOrganSite == oi & ggplotdata$Type == fi & ggplotdata$DataSet == "CD4"]
    }
}
ggplotdata$Type <- factor(ggplotdata$Type, levels = c("Healthy donor", "Uninvolved normal tissue", "Primary tumor tissue", "Metastatic tumor tissue"))
ggplotdata$DataSet <- factor(ggplotdata$DataSet, levels = dataTypeLevel)
g_OrganSite_figure1F <- ggplot(data = ggplotdata) +
    geom_bar(aes(x = reorder_within(DetailOrganSite, -ForOrder, Type), y = frac, fill = DataSet), stat = "identity") +
    ylab("Cell fraction") + xlab("") +
    scale_fill_manual(values = colorsForDataType) +
    scale_x_reordered() +
    facet_grid(. ~ Type , scales = "free_x", space = "free_x")+
    theme_classic() +
    theme(
        text = element_text(size = 10),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        strip.background = element_blank(),
          strip.text.x = element_text(size = 10))
pdf(file.path(figurePath, paste0("DetailOrganSite_figure1F.pdf")), width = 8, height = 3)
print(g_OrganSite_figure1F)
dev.off()
ggsave(file.path(figurePath, paste0("DetailOrganSite_figure1F.pdf")), g_OrganSite_figure1F,  width = 210 * 1.1, height = 497/7, units = "mm")

write_tsv(origin.ggplotdata, file.path(figurePath, paste0('figure1F_c_detail', "_", Sys.Date(), '.tsv')))

for(tti in unique(origin.ggplotdata$Type)){
    tempT <- origin.ggplotdata %>% filter(Type == tti)

    tempTn <- tempT %>% select(DetailOrganSite, DataSet, n) %>%
        spread(DataSet, n) %>% replace(is.na(.), 0)
    write_tsv(tempTn, file.path(figurePath, paste0(tti, '_n', "_", Sys.Date(), '.tsv')))

    tempTfrac <- tempT %>% select(DetailOrganSite, DataSet, frac) %>%
        spread(DataSet, frac) %>% replace(is.na(.), 0)
    write_tsv(tempTfrac, file.path(figurePath, paste0(tti, '_frac', "_", Sys.Date(), '.tsv')))
}

#' statistics ################################################################

## set.seed(123)
## library(ggstatsplot)
## figurePath <- "/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/figures4"
## figurePath <- file.path(figurePath, "tissue_sample_fraction_all_boxplot")
## if(!dir.exists(figurePath)){
##     dir.create(figurePath)
## }
## sampleN <- allMD %>% group_by(OrganSite, Sample) %>% count() %>% filter(n > 200)
## TypeSampleN <- allMD
## TypeSampleN$CancerType[TypeSampleN$TissueType %in% c("Healthy donor", "Uninvolved normal tissue")] = TypeSampleN$OrganSite[TypeSampleN$TissueType %in% c("Healthy donor", "Uninvolved normal tissue")]
## TypeSampleN <- TypeSampleN %>% group_by(TissueType, CancerType, Sample) %>% count() %>% filter(Sample %in% sampleN$Sample)
## TissueType.labs <- c("H", "U", "P", "M")
## names(TissueType.labs) <- c("Healthy donor", "Uninvolved normal tissue", "Primary tumor tissue", "Metastatic tumor tissue")
## for(oi in 1:length(allObj)){
##     tempName <- names(allObj[oi])
##     tempObj <- allObj[[oi]]
##     tempFigurePath <- file.path(figurePath, tempName)
##     if(!dir.exists(tempFigurePath)){
##         dir.create(tempFigurePath)
##     }
##     Idents(tempObj) <- tempObj@meta.data$seurat_clusters
##     ggplotdata <- as_tibble(tempObj@meta.data, rownames = NA) %>% select(seurat_clusters, Sample, CancerType, TissueType, OrganSite, orig.ident)
##     ggplotdata$CancerType[ggplotdata$TissueType %in% c("Healthy donor", "Uninvolved normal tissue")] = ggplotdata$OrganSite[ggplotdata$TissueType %in% c("Healthy donor", "Uninvolved normal tissue")]
##     ggplotdata <- ggplotdata %>% group_by(seurat_clusters, TissueType, CancerType, Sample) %>% dplyr::count() %>% filter(Sample %in% sampleN$Sample)
##     ggplotdata$frac <- ggplotdata$n / sampleN$n[match(ggplotdata$Sample, sampleN$Sample)]

##     totalObj_ggplotdata <- c()
##     for(tempSeuratCluster in unique(ggplotdata$seurat_clusters)){
##         temp_ggplotdata <- ggplotdata %>% filter(seurat_clusters == tempSeuratCluster)
##         for(temp_TissueType in unique(TypeSampleN$TissueType)){
##             temp_TypeSampleN <- TypeSampleN %>% filter(TissueType == temp_TissueType)
##             temp_diff_Samples <- setdiff(unique( temp_TypeSampleN %>% pull(Sample)),
##                                          unique(temp_ggplotdata %>% filter(TissueType == temp_TissueType) %>% pull(Sample)))
##             if(length(temp_diff_Samples) > 0 ){
##                 t_rows <- tibble(seurat_clusters = rep(tempSeuratCluster, length(temp_diff_Samples)),
##                                  TissueType = rep(temp_TissueType, length(temp_diff_Samples)),
##                                  CancerType = temp_TypeSampleN$CancerType[temp_TypeSampleN$Sample %in% temp_diff_Samples],
##                                  Sample = temp_diff_Samples,
##                                  n = rep(0, length(temp_diff_Samples)),
##                                  frac = rep(0, length(temp_diff_Samples)))
##                 temp_ggplotdata <- bind_rows(temp_ggplotdata, t_rows)
##             }
##         }
##         temp_ggplotdata$TissueType <- factor(temp_ggplotdata$TissueType, levels = c("Healthy donor", "Uninvolved normal tissue", "Primary tumor tissue", "Metastatic tumor tissue"))
##         ## g <- ggplot(data = temp_ggplotdata) +
##         ##     geom_boxplot(aes(x = TissueType, y = frac, color = CancerType), size = 0.2, outlier.shape = NA) +
##         ##     geom_jitter(shape=16, position=position_jitter(0.2), size = 0.2,
##         ##                 aes(x = TissueType, y = frac, fill = CancerType, color = CancerType)) +
##         ##     xlab("") + ylab("") +
##         ##     ## scale_y_log10() +
##         ##     theme_classic() +
##         ##     theme(text = element_text(size = 8),
##         ##           strip.background = element_blank(),
##         ##           axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
##         ##           legend.position = "none")
##         if(tempName == "NK"){
##             if(tempSeuratCluster %in% c(4)){
##                 tempDataSet = "MAIT"
##             }else if(tempSeuratCluster %in% c(1,6)){
##                 tempDataSet = "MAIT-like"
##             }else if(tempSeuratCluster %in% c(3,7,8)){
##                 tempDataSet = "NKT"
##             }else{
##                 tempDataSet = "NK"
##             }
##         }else{
##             tempDataSet = tempName
##         }
##         ## ggsave(file.path(tempFigurePath, paste0(tempDataSet, "_c", tempSeuratCluster, "_sample_boxplot.pdf")), g, width = 210/8, height = 297/8, units = "mm")
##         temp_ggplotdata$g_seurat_clusters <- paste0(tempDataSet, "_c", tempSeuratCluster)
##         totalObj_ggplotdata <- bind_rows(totalObj_ggplotdata, temp_ggplotdata)
##         ## g <- ggstatsplot::ggbetweenstats(
##         ##                       data = temp_ggplotdata,
##         ##                       x = TissueType,
##         ##                       y = frac,
##         ##                       xlab = "Tissue type",
##         ##                       ylab = "Sample fraction",
##         ##                       ggtheme = theme_classic(),
##         ##                       plot.type = "box") +
##         ##     theme(text = element_text(size = 8))
##         ## ggsave(file.path(tempFigurePath, paste0(tempDataSet, "_c", tempSeuratCluster, "_sample_boxplot_gbs.pdf")), g, width = 210/1.5, height = 297/2, units = "mm")
##     }
##     if(tempName == "CD8"){
##         g <- totalObj_ggplotdata %>%
##             ggstatsplot::grouped_ggbetweenstats(
##                              data = .,
##                              x = TissueType,
##                              y = frac,
##                              grouping.var = g_seurat_clusters,
##                              xlab = "Tissue type",
##                              ylab = "Sample fraction",
##                              pairwise.display = "all", # display only significant pairwise comparisons
##                              p.adjust.method = "fdr", # adjust p-values for multiple tests using this method
##                              ggtheme = theme_classic(),
##                              package = "ggsci",
##                              palette = "default_jco",
##                              plotgrid.args = list(ncol = 3))
##         ggsave(file.path(tempFigurePath, paste0("cluster_all-sample_boxplot_gbs.pdf")), g, width = 400, height = 500, units = "mm")
        
##         ggplotdata <- totalObj_ggplotdata
##         ggplotdata$TissueType <- plyr::mapvalues(ggplotdata$TissueType,
##                                                  c("Healthy donor",
##                                                    "Uninvolved normal tissue",
##                                                    "Primary tumor tissue",
##                                                    "Metastatic tumor tissue"),
##                                                  c("H", "U", "P", "M"))
##         ggplotdata$TissueType <- factor(ggplotdata$TissueType, levels = c("H", "U", "P", "M"))
##         temp_ggplotdata <- ggplotdata %>%
##             filter(g_seurat_clusters %in% paste0("CD8_c", c(5, 9, 10, 11)))
##         levels(temp_ggplotdata$g_seurat_clusters) <-
##             paste0("CD8_c", c(5, 9, 10, 11))
##         temp_ggplotdata$g_seurat_clusters <- factor(temp_ggplotdata$g_seurat_clusters,
##                                                     levels = paste0("CD8_c", c(5, 9, 10, 11)))
##         g <- temp_ggplotdata %>% ggplot() +
##             geom_jitter(aes(x = TissueType,
##                             y = frac,
##                             fill = TissueType,
##                             color = TissueType), size = 0.2) +
##             geom_boxplot(aes(x = TissueType,
##                              y = frac),
##                          fill = NA,
##                          outlier.shape = NA,
##                          lwd = 0.2) +
##             facet_wrap(g_seurat_clusters ~ ., scales ="free", nrow =2) +
##             theme_classic() +
##             theme(text = element_text(size = 10),
##                   legend.position = "none",
##                   strip.background = element_rect(colour="white", fill="white"))
##         ggsave(file.path(figurePath,
##                          paste0(tempName, "_cluster_selected-sample_boxplot.pdf")),
##                g, width = 210 * 2/5, height = 290/4, units = "mm")

##         ## totalObj_ggplotdata$TissueType <- plyr::mapvalues(totalObj_ggplotdata$TissueType,
##         ##                                                   c("Healthy donor",
##         ##                                                     "Uninvolved normal tissue",
##         ##                                                     "Primary tumor tissue",
##         ##                                                     "Metastatic tumor tissue"),
##         ##                                                   c("H", "U", "P", "M"))
##         ## totalObj_ggplotdata$g_seurat_clusters <- factor(totalObj_ggplotdata$g_seurat_clusters,
##         ##                                                 levels = paste0("CD8_c", 0:13))
##         ## g <- totalObj_ggplotdata %>%
##         ##     ggplot() +
##         ##     geom_jitter(aes(x = TissueType,
##         ##                     y = frac,
##         ##                     fill = TissueType,
##         ##                     color = TissueType)) +
##         ##     geom_boxplot(aes(x = TissueType,
##         ##                      y = frac),
##         ##                  fill = NA,
##         ##                  outlier.shape = NA) +
##         ##     facet_wrap(g_seurat_clusters ~ ., scales ="free") +
##         ##     theme_classic() +
##         ##     theme(text = element_text(size = 10),
##         ##         strip.background = element_rect(colour="white", fill="white"))
##         ## ggsave(file.path(tempFigurePath, paste0("cluster_all-sample_boxplot.pdf")), g, width = 210, height = 290/2, units = "mm")


##         temp_ggplotdata <- totalObj_ggplotdata %>% filter(g_seurat_clusters %in% c("CD8_c1", "CD8_c2", "CD8_c3", "CD8_c0"))
##         temp_ggplotdata$g_seurat_clusters <- factor(temp_ggplotdata$g_seurat_clusters, levels = c("CD8_c3", "CD8_c0", "CD8_c2", "CD8_c1"))
##         g <- temp_ggplotdata %>%
##             ggstatsplot::grouped_ggbetweenstats(
##                              data = .,
##                              x = TissueType,
##                              y = frac,
##                              grouping.var = g_seurat_clusters,
##                              xlab = "Tissue type",
##                              ylab = "Sample fraction",
##                              ## pairwise.display = "aiwl", # display only significant pairwise comparisons
##                              p.adjust.method = "fdr", # adjust p-values for multiple tests using this method
##                              ggtheme = theme_classic(),
##                              package = "ggsci",
##                              palette = "default_jco",
##                              plotgrid.args = list(nrow = 1))
##         ggsave(file.path(tempFigurePath, paste0("cluster3_0_2_1-sample_boxplot_gbs.pdf")), g, width = 400, height = 292/2, units = "mm")


##         testForFig2J <- temp_ggplotdata %>%
##             filter(g_seurat_clusters == "CD8_c1") 

##         tg1 <- testForFig2J %>% filter((TissueType == "Primary tumor tissue" & CancerType == "NSCLC") |
##                                        (TissueType == "Uninvolved normal tissue" & CancerType == "Lung"))
##         tg2 <- testForFig2J %>% filter((TissueType == "Primary tumor tissue" & CancerType == "HNSC") |
##                                        (TissueType == "Uninvolved normal tissue" & CancerType == "Tongue/Tonsil"))
##         tg3 <- testForFig2J %>% filter((CancerType == "FL") |
##                                        (TissueType == "Healthy donor" & CancerType == "LN"))
##         tg4 <- testForFig2J %>% filter((CancerType == "DLBCL") |
##                                        (TissueType == "Healthy donor" & CancerType == "LN"))
##         tg5 <- testForFig2J %>% filter((CancerType == "DLBCL" | CancerType == "FL") |
##                                        (TissueType == "Healthy donor" & CancerType == "LN"))
##         tg5$CancerType[tg5$CancerType %in% c("DLBCL", "FL")] <- "DLBCL+FL"
##         tgL <- list(tg1 = tg1, tg2 = tg2, tg3 = tg3, tg4 = tg4, tg5 = tg5)

##         for(tgli in 1:length(tgL)){
##             tempName <- names(tgL[tgli])
##             tempTg <- tgL[[tgli]]

            
##             g <- ggstatsplot::ggbetweenstats(
##                                   data = tempTg,
##                                   x = CancerType,
##                                   y = frac,
##                                   xlab = "Tissue type",
##                                   ylab = "Sample fraction",
##                                   ggtheme = theme_classic(),
##                                   plot.type = "box") +
##                 theme(text = element_text(size = 8))
##             ggsave(file.path("/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/scripts/pipelines/step4_figs/FigsForPaper/supplement/test/outs", paste0(tempName, "_sample_boxplot_gbs.pdf")), g, width = 210/1.5, height = 297/2, units = "mm")
##         }




        

##         temp_ggplotdata %>% select(g_seurat_clusters, TissueType, frac) %>% group_by(g_seurat_clusters, TissueType) %>% mutate(sd = sd(frac)) %>% select(g_seurat_clusters, TissueType, sd) %>% group_by(g_seurat_clusters, TissueType, sd) %>% dplyr::count()
##         g <- totalObj_ggplotdata %>%
##             ggstatsplot::grouped_ggbetweenstats(
##                              data = .,
##                              x = g_seurat_clusters,
##                              y = frac,
##                              grouping.var = TissueType,
##                              xlab = "Tissue type",
##                              ylab = "Sample fraction",
##                              ## pairwise.display = "all", # display only significant pairwise comparisons
##                              p.adjust.method = "fdr", # adjust p-values for multiple tests using this method
##                              ggtheme = theme_classic(),
##                              package = "ggsci",
##                              palette = "default_jco",
##                              plotgrid.args = list(nrow = 4))
##         ggsave(file.path(tempFigurePath, paste0("Tissue_all-sample_boxplot_gbs.pdf")), g, width = 200, height = 800, units = "mm")

##         allMD <- allMD %>% mutate(g_seurat_clusters = paste0(DataSet, "_c", seurat_clusters))
##         temp_ggplotdata <- allMD %>% filter(g_seurat_clusters %in% c("CD8_c3", "CD8_c0", "CD8_c2", "CD8_c1")) %>% group_by(g_seurat_clusters, TissueType) %>% dplyr::count() %>% group_by(g_seurat_clusters) %>% mutate(frac = n / sum(n))

##         temp_ggplotdata$g_seurat_clusters <- factor(temp_ggplotdata$g_seurat_clusters, levels = c("CD8_c3", "CD8_c0", "CD8_c2", "CD8_c1"))

##         g <- ggplot(temp_ggplotdata, aes(x="", y=frac, fill=TissueType)) +
##             geom_bar(stat="identity", width=1, color="white") +
##             coord_polar("y", start=0) +
##             facet_grid(g_seurat_clusters ~ .)
##             theme_classic() # remove background, grid, numeric labels
##         ggsave(file.path(tempFigurePath, paste0("Tissue_sample3-0-2-1_pieplot_gbs.pdf")), g, width = 200, height = 200, units = "mm")
##     }
##     if(tempName == "CD4"){

##         ggplotdata <- totalObj_ggplotdata
##         ggplotdata$TissueType <- plyr::mapvalues(ggplotdata$TissueType,
##                                                  c("Healthy donor",
##                                                    "Uninvolved normal tissue",
##                                                    "Primary tumor tissue",
##                                                    "Metastatic tumor tissue"),
##                                                  c("H", "U", "P", "M"))
##         ggplotdata$TissueType <- factor(ggplotdata$TissueType, levels = c("H", "U", "P", "M"))
##         temp_ggplotdata <- ggplotdata %>%
##             filter(g_seurat_clusters %in% paste0("CD4_c", c(0, 5, 6, 7, 9, 10)))
##         temp_ggplotdata$g_seurat_clusters <- factor(temp_ggplotdata$g_seurat_clusters,
##                                                     levels = paste0("CD4_c", c(0, 5, 6, 7, 9, 10)))
##         g <- temp_ggplotdata %>%
##             ggplot() +
##             geom_jitter(aes(x = TissueType,
##                             y = frac,
##                             fill = TissueType,
##                             color = TissueType), size = 0.2) +
##             geom_boxplot(aes(x = TissueType,
##                              y = frac),
##                          fill = NA,
##                          outlier.shape = NA,
##                          lwd = 0.2) +
##             facet_wrap(g_seurat_clusters ~ ., scales ="free", nrow =1) +
##             theme_classic() +
##             theme(text = element_text(size = 10),
##                   legend.position = "none",
##                   strip.background = element_rect(colour="white", fill="white"))
##         ggsave(file.path(figurePath,
##                          paste0(tempName, "_cluster_selected-sample_boxplot.pdf")),
##                g, width = 210, height = 290/5, units = "mm")
        
##         g <- totalObj_ggplotdata %>%
##             ggstatsplot::grouped_ggbetweenstats(
##                              data = .,
##                              x = TissueType,
##                              y = frac,
##                              grouping.var = g_seurat_clusters,
##                              xlab = "Tissue type",
##                              ylab = "Sample fraction",
##                              pairwise.display = "all", # display only significant pairwise comparisons
##                              p.adjust.method = "fdr", # adjust p-values for multiple tests using this method
##                              ggtheme = theme_classic(),
##                              package = "ggsci",
##                              palette = "default_jco",
##                              plotgrid.args = list(ncol = 3))
##         ggsave(file.path(tempFigurePath, paste0("cluster_all-sample_boxplot_gbs.pdf")), g, width = 400, height = 500, units = "mm")

##         totalObj_ggplotdata$TissueType <- plyr::mapvalues(totalObj_ggplotdata$TissueType,
##                                                           c("Healthy donor",
##                                                             "Uninvolved normal tissue",
##                                                             "Primary tumor tissue",
##                                                             "Metastatic tumor tissue"),
##                                                           c("H", "U", "P", "M"))
##         totalObj_ggplotdata$g_seurat_clusters <- factor(totalObj_ggplotdata$g_seurat_clusters,
##                                                         levels = paste0("CD4_c", 0:11))
##         g <- totalObj_ggplotdata %>%
##             ggplot() +
##             geom_jitter(aes(x = TissueType,
##                             y = frac,
##                             fill = TissueType,
##                             color = TissueType)) +
##             geom_boxplot(aes(x = TissueType,
##                              y = frac),
##                          fill = NA,
##                          outlier.shape = NA) +
##             facet_wrap(g_seurat_clusters ~ ., scales ="free") +
##             theme_classic() +
##             theme(text = element_text(size = 10),
##                 strip.background = element_rect(colour="white", fill="white"))
##         ggsave(file.path(tempFigurePath, paste0("cluster_all-sample_boxplot.pdf")), g, width = 210, height = 290/2, units = "mm")

        

##         temp_ggplotdata <- totalObj_ggplotdata %>% filter(g_seurat_clusters %in% c("CD4_c1", "CD4_c3", "CD4_c4", "CD4_c8"))
##         temp_ggplotdata$g_seurat_clusters <- factor(temp_ggplotdata$g_seurat_clusters, levels = c("CD4_c1", "CD4_c3", "CD4_c4", "CD4_c8"))

##         g <- temp_ggplotdata %>%
##             ggstatsplot::grouped_ggbetweenstats(
##                              data = .,
##                              x = TissueType,
##                              y = frac,
##                              grouping.var = g_seurat_clusters,
##                              xlab = "Tissue type",
##                              ylab = "Sample fraction",
##                              pairwise.display = "all", # display only significant pairwise comparisons
##                              p.adjust.method = "fdr", # adjust p-values for multiple tests using this method
##                              ggtheme = theme_classic(),
##                              package = "ggsci",
##                              palette = "default_jco",
##                              plotgrid.args = list(nrow = 1))
##         ggsave(file.path(tempFigurePath, paste0("cluster1_3_4_8-sample_boxplot_gbs.pdf")), g, width = 400, height = 292/2, units = "mm")

##         g <- temp_ggplotdata %>%
##             ggstatsplot::grouped_ggbetweenstats(
##                              data = .,
##                              x = g_seurat_clusters,
##                              y = frac,
##                              grouping.var = TissueType,
##                              xlab = "Cluster",
##                              ylab = "Sample fraction",
##                              ## pairwise.display = "aiwl", # display only significant pairwise comparisons
##                              p.adjust.method = "fdr", # adjust p-values for multiple tests using this method
##                              ggtheme = theme_classic(),
##                              package = "ggsci",
##                              palette = "default_jco",
##                              plotgrid.args = list(nrow = 1))
##         ggsave(file.path(tempFigurePath, paste0("cluster1_3_4_8-sample_boxplot_TissueType_gbs.pdf")), g, width = 400, height = 292/2, units = "mm")


##         temp_ggplotdata %>% select(g_seurat_clusters, TissueType, frac) %>% group_by(g_seurat_clusters, TissueType) %>% mutate(sd = sd(frac)) %>% select(g_seurat_clusters, TissueType, sd) %>% group_by(g_seurat_clusters, TissueType, sd) %>% dplyr::count()

##         temp_ggplotdata %>% group_by(g_seurat_clusters, TissueType) %>% count()


##         g <- totalObj_ggplotdata %>%
##             ggstatsplot::grouped_ggbetweenstats(
##                              data = .,
##                              x = g_seurat_clusters,
##                              y = frac,
##                              grouping.var = TissueType,
##                              xlab = "Tissue type",
##                              ylab = "Sample fraction",
##                              ## pairwise.display = "all", # display only significant pairwise comparisons
##                              p.adjust.method = "fdr", # adjust p-values for multiple tests using this method
##                              ggtheme = theme_classic(),
##                              package = "ggsci",
##                              palette = "default_jco",
##                              plotgrid.args = list(nrow = 4))
##         ggsave(file.path(tempFigurePath, paste0("Tissue_all-sample_boxplot_gbs.pdf")), g, width = 200, height = 800, units = "mm")

##         library(ggforce)
##         allMD <- allMD %>% mutate(g_seurat_clusters = paste0(DataSet, "_c", seurat_clusters))
##         temp_ggplotdata <- allMD %>% filter(g_seurat_clusters %in% c("CD4_c1", "CD4_c3", "CD4_c0", "CD4_c8")) %>% group_by(g_seurat_clusters, TissueType) %>% dplyr::count() %>% group_by(g_seurat_clusters) %>% mutate(frac = n / sum(n))
##         temp_ggplotdata$TissueType <- factor(temp_ggplotdata$TissueType, levels = c("Healthy donor", "Uninvolved normal tissue", "Primary tumor tissue", "Metastatic tumor tissue"))
##         g <- ggplot(temp_ggplotdata) +
##             geom_bar(aes(x = n, y = TissueType, fill = TissueType), stat = "identity")  +
##             facet_grid(g_seurat_clusters ~ .) +
##             theme_classic() +
##             scale_x_continuous(trans=trans_reverser('log10'))
##         ggsave(file.path(tempFigurePath, paste0("Tissue_sample1-3-0-8_barplot_gbs.pdf")), g, width = 210, height = 297/4, units = "mm")


##         temp_ggplotdata$g_seurat_clusters <- factor(temp_ggplotdata$g_seurat_clusters, levels = c("CD4_c1", "CD4_c3", "CD4_c0", "CD4_c8"))
##         g <- ggplot(temp_ggplotdata, aes(x="", y=frac, fill=TissueType)) +
##             geom_bar(stat="identity", width=1, color="white") +
##             coord_polar("y", start=0) +
##             facet_grid(g_seurat_clusters ~ .)
##         theme_classic() # remove background, grid, numeric labels
##         ggsave(file.path(tempFigurePath, paste0("Tissue_sample1-3-0-8_pieplot_gbs.pdf")), g, width = 200, height = 200, units = "mm")
##     }
##     ## g <- ggplot(data = totalObj_ggplotdata) +
##     ##     geom_boxplot(aes(x = reorder(CancerType, frac, FUN=median),
##     ##                      y = frac,
##     ##                      color = CancerType),
##     ##                  size = 0.2,
##     ##                  outlier.shape = NA) +
##     ##     geom_jitter(shape=16,
##     ##                 position=position_jitter(0.2),
##     ##                 size = 0.2,
##     ##                 aes(x = reorder(CancerType, frac, FUN=median),
##     ##                     y = frac,
##     ##                     fill = CancerType,
##     ##                     color = CancerType)) +
##     ##     facet_grid(g_seurat_clusters ~ TissueType,
##     ##                        scales = "free_x", space = "free",
##     ##                        labeller = labeller(TissueType = TissueType.labs)) +
##     ##     xlab("") + ylab("") +
##     ##     ## scale_y_log10() +
##     ##     theme_classic() +
##     ##     theme( text = element_text(size = 8),
##     ##           strip.background = element_blank(),
##     ##           axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
##     ##           legend.position = "none")
##     ## ggsave(file.path(tempFigurePath, paste0("all_boxplot.pdf")), g, width = 210/2, height = 297, units = "mm")
##     ## g <- ggplot(data = totalObj_ggplotdata) +
##     ##     geom_boxplot(aes(x = reorder(CancerType, frac, FUN=median),
##     ##                      y = frac,
##     ##                      color = g_seurat_clusters),
##     ##                  size = 0.2,
##     ##                  outlier.shape = NA) +
##     ##     geom_jitter(shape=16,
##     ##                 position=position_jitter(0.2),
##     ##                 size = 0.2,
##     ##                 aes(x = reorder(CancerType, frac, FUN=median),
##     ##                     y = frac,
##     ##                     fill = g_seurat_clusters,
##     ##                     color = g_seurat_clusters)) +
##     ##     facet_grid(TissueType ~ .,
##     ##                        scales = "free_x", space = "free",
##     ##                        labeller = labeller(TissueType = TissueType.labs)) +
##     ##     xlab("") + ylab("") +
##     ##     ## scale_y_log10() +
##     ##     theme_classic() +
##     ##     theme( text = element_text(size = 8),
##     ##           strip.background = element_blank(),
##     ##           axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
##     ##           legend.position = "none")
##     ## ggsave(file.path(tempFigurePath, paste0("all_boxplot_layout2.pdf")), g, width = 297*4, height = 210/4, units = "mm")
##     ## for(tempTissueType in unique(totalObj_ggplotdata$TissueType)){
##     ##     temp_ggplotdata <- totalObj_ggplotdata %>% filter(TissueType == tempTissueType)
##     ##     g <- ggplot(data = temp_ggplotdata) +
##     ##         geom_boxplot(aes(x = g_seurat_clusters,
##     ##                          y = frac,
##     ##                          color = CancerType),
##     ##                      size = 0.2,
##     ##                      outlier.shape = NA) +
##     ##         geom_jitter(shape=16,
##     ##                     position=position_jitter(0.2),
##     ##                     size = 0.2,
##     ##                     aes(x = g_seurat_clusters,
##     ##                         y = frac,
##     ##                         fill = CancerType,
##     ##                         color = CancerType)) +
##     ##         xlab("") + ylab("") +
##     ##         facet_grid(CancerType ~ .) +
##     ##         ## scale_y_log10() +
##     ##         theme_classic() +
##     ##         theme( text = element_text(size = 8),
##     ##               strip.background = element_blank(),
##     ##               axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
##     ##               legend.position = "bottom")
##     ##     ggsave(file.path(tempFigurePath, paste0(tempTissueType, "_all_boxplot_layout2.pdf")), g, width = 210/2, height = 297, units = "mm")
##     ## }
## }






## figurePath <- "/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/figures4"
## figurePath <- file.path(figurePath, "cluster_sample_fraction_all_boxplot")
## if(!dir.exists(figurePath)){
##     dir.create(figurePath)
## }
## sampleN <- allMD %>% group_by(OrganSite, Sample) %>% count() %>% filter(n > 200)
## TypeSampleN <- allMD
## TypeSampleN$CancerType[TypeSampleN$TissueType %in% c("Healthy donor", "Uninvolved normal tissue")] = TypeSampleN$OrganSite[TypeSampleN$TissueType %in% c("Healthy donor", "Uninvolved normal tissue")]
## TypeSampleN <- TypeSampleN %>% group_by(TissueType, CancerType, Sample) %>% count() %>% filter(Sample %in% sampleN$Sample)
## TissueType.labs <- c("H", "U", "P", "M")
## names(TissueType.labs) <- c("Healthy donor", "Uninvolved normal tissue", "Primary tumor tissue", "Metastatic tumor tissue")
## for(oi in 1:length(allObj)){
##     tempName <- names(allObj[oi])
##     tempObj <- allObj[[oi]]
##     tempFigurePath <- file.path(figurePath, tempName)
##     if(!dir.exists(tempFigurePath)){
##         dir.create(tempFigurePath)
##     }
##     Idents(tempObj) <- tempObj@meta.data$seurat_clusters
##     ggplotdata <- as_tibble(tempObj@meta.data, rownames = NA) %>% select(seurat_clusters, Sample, CancerType, TissueType, OrganSite, orig.ident)
##     ggplotdata$CancerType[ggplotdata$TissueType %in% c("Healthy donor", "Uninvolved normal tissue")] = ggplotdata$OrganSite[ggplotdata$TissueType %in% c("Healthy donor", "Uninvolved normal tissue")]
##     ggplotdata <- ggplotdata %>% group_by(seurat_clusters, TissueType, CancerType, Sample) %>% dplyr::count() %>% filter(Sample %in% sampleN$Sample)
##     ggplotdata$frac <- ggplotdata$n / sampleN$n[match(ggplotdata$Sample, sampleN$Sample)]
##     totalObj_ggplotdata <- c()
##     for(tempSeuratCluster in unique(ggplotdata$seurat_clusters)){
##         temp_ggplotdata <- ggplotdata %>% filter(seurat_clusters == tempSeuratCluster)
##         for(temp_TissueType in unique(TypeSampleN$TissueType))
##         {
##             temp_TypeSampleN <- TypeSampleN %>% filter(TissueType == temp_TissueType)
##             temp_diff_Samples <- setdiff(unique( temp_TypeSampleN %>% pull(Sample)),
##                                          unique(temp_ggplotdata %>% filter(TissueType == temp_TissueType) %>% pull(Sample) ))
##             if(length(temp_diff_Samples) > 0 ){
##                 t_rows <- tibble(seurat_clusters = rep(tempSeuratCluster, length(temp_diff_Samples)),
##                                  TissueType = rep(temp_TissueType, length(temp_diff_Samples)),
##                                  CancerType = temp_TypeSampleN$CancerType[temp_TypeSampleN$Sample %in% temp_diff_Samples],
##                                  Sample = temp_diff_Samples,
##                                  n = rep(0, length(temp_diff_Samples)),
##                                  frac = rep(0, length(temp_diff_Samples)))
##                 temp_ggplotdata <- bind_rows(temp_ggplotdata, t_rows)
##             }
##         }
##         temp_ggplotdata$TissueType <- factor(temp_ggplotdata$TissueType, levels = c("Healthy donor", "Uninvolved normal tissue", "Primary tumor tissue", "Metastatic tumor tissue"))
##         ## g <- ggplot(data = temp_ggplotdata) +
##         ##     geom_boxplot(aes(x = reorder(CancerType, -frac, FUN=median), y = frac, color = CancerType), size = 0.2, outlier.shape = NA) +
##         ##     geom_jitter(shape=16, position=position_jitter(0.2), size = 0.2, aes(x = reorder(CancerType, frac, FUN=median), y = frac, fill = CancerType, color = CancerType)) +
##         ##     ggforce::facet_row(TissueType ~ ., scales = "free_x", space = "free",
##         ##                        labeller = labeller(TissueType = TissueType.labs)) +
##         ##     xlab("") + ylab("") +
##         ##     ## scale_y_log10() +
##         ##     theme_classic() +
##         ##     theme( text = element_text(size = 8),
##         ##           strip.background = element_blank(),
##         ##         axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
##         ##           legend.position = "none")
##         if(tempName == "NK"){
##             if(tempSeuratCluster %in% c(4)){
##                 tempDataSet = "MAIT"
##             }else if(tempSeuratCluster %in% c(1,6)){
##                 tempDataSet = "MAIT-like"
##             }else if(tempSeuratCluster %in% c(3,7,8)){
##                 tempDataSet = "NKT"
##             }else{
##                 tempDataSet = "NK"
##             }
##         }else{
##             tempDataSet = tempName
##         }
##         ## ggsave(file.path(tempFigurePath, paste0(tempDataSet, "_c", tempSeuratCluster, "_sample_boxplot.pdf")), g, width = 210/2, height = 297/8, units = "mm")
##         temp_ggplotdata$g_seurat_clusters <- paste0(tempDataSet, "_c", tempSeuratCluster)
##         totalObj_ggplotdata <- bind_rows(totalObj_ggplotdata, temp_ggplotdata)
##         ## temp_ggplotdata <- temp_ggplotdata %>% mutate(key = paste0(TissueType, "_", CancerType))
##         ## selectedT <- temp_ggplotdata %>% group_by(TissueType, CancerType, g_seurat_clusters) %>% dplyr::count() %>% filter(n >= 3) %>% mutate(key = paste0(TissueType, "_", CancerType))
##         ## g <- temp_ggplotdata %>% filter(key %in% selectedT$key) %>% 
##         ##     ggstatsplot::grouped_ggbetweenstats(
##         ##                      data = .,
##         ##                      ## x = TissueType,
##         ##                      x = CancerType,
##         ##                      y = frac,
##         ##                      grouping.var = TissueType,
##         ##                      xlab = "Cancer type",
##         ##                      ylab = "Sample fraction",
##         ##                      ## pairwise.display = "all", # display only significant pairwise comparisons
##         ##                      p.adjust.method = "fdr", # adjust p-values for multiple tests using this method
##         ##                      ggtheme = theme_classic(),
##         ##                      package = "ggsci",
##         ##                      palette = "default_jco",
##         ##                      plotgrid.args = list(nrow = 4))
##         ## ggsave(file.path(tempFigurePath, paste0(tempDataSet, "_c", tempSeuratCluster, "_cancer_type_sample_boxplot.pdf")), g, width = 210, height = 600, units = "mm")
##     }

##     if(tempName == "CD8"){
##         temp_ggplotdata <- totalObj_ggplotdata %>% filter(g_seurat_clusters %in% c("CD8_c3", "CD8_c0", "CD8_c2", "CD8_c1")) %>% mutate(key = paste0(TissueType, "_", CancerType))
##         selectedT <- temp_ggplotdata %>% group_by(TissueType, CancerType, g_seurat_clusters) %>% dplyr::count() %>% filter(n >= 3) %>% mutate(key = paste0(TissueType, "_", CancerType))
##         temp_ggplotdata <- temp_ggplotdata %>% filter(key %in% selectedT$key)
##         temp_ggplotdata$names <- paste0(temp_ggplotdata$TissueType, "_", temp_ggplotdata$CancerType)
##         temp_ggplotdata$names <- factor(temp_ggplotdata$names,
##                                         levels = c("Healthy donor_BM",
##                                                    "Healthy donor Lymph_node",
##                                                    "Healthy donor_Blood",
##                                                    "Uninvolved normal tissue_Brain",
##                                                    "Uninvolved normal tissue_Breast",
##                                                    "Uninvolved normal tissue_Lung",
##                                                    "Uninvolved normal tissue_Blood",
##                                                    "Uninvolved normal tissue_colon",
##                                                    "Uninvolved normal tissue_Head and neck",
##                                                    "Primary tumor tissue_NSCLC",
##                                                    "Primary tumor tissue_HNSC",
##                                                    "Primary tumor tissue_LBCL",
##                                                    "Primary tumor tissue_LUSC",
##                                                    "Primary tumor tissue_BRCA",
##                                                    "Primary tumor tissue_FL",
##                                                    "Primary tumor tissue_LUAD",
##                                                    "Primary tumor tissue_STAD",
##                                                    "Primary tumor tissue_BCC",
##                                                    "Primary tumor tissue_HCC",
##                                                    "Primary tumor tissue_GBM",
##                                                    "Primary tumor tissue_MCL",
##                                                    "Primary tumor tissue_AML",
##                                                    "Primary tumor tissue_CRC",
##                                                    "Primary tumor tissue_PDAC",
##                                                    "Metastatic tumor tissue_SKCM",
##                                                    "Metastatic tumor tissue_OV",
##                                                    "Metastatic tumor tissue_LUAD",
##                                                    "Metastatic tumor tissue_STAD"))
##         temp_ggplotdata$CancerType <- factor(temp_ggplotdata$CancerType,
##                                         levels = c("BM",
##                                                    "Lymph node",
##                                                    "Brain",
##                                                    "Breast",
##                                                    "Lung",
##                                                    "Blood",
##                                                    "colon",
##                                                    "Head and neck",
##                                                    "SKCM",
##                                                    "OV",
##                                                    "NSCLC",
##                                                    "HNSC",
##                                                    "LBCL",
##                                                    "LUSC",
##                                                    "BRCA",
##                                                    "FL",
##                                                    "LUAD",
##                                                    "STAD",
##                                                    "BCC",
##                                                    "HCC",
##                                                    "GBM",
##                                                    "MCL",
##                                                    "AML",
##                                                    "CRC",
##                                                    "PDAC"))
##         temp_ggplotdata$g_seurat_clusters <- factor(temp_ggplotdata$g_seurat_clusters, levels = c("CD8_c3", "CD8_c0", "CD8_c2", "CD8_c1"))
##         temp_ggplotdata <- temp_ggplotdata %>% mutate(TissueType = fct_recode(TissueType, H = "Healthy donor", U = "Uninvolved normal tissue", P = "Primary tumor tissue", M = "Metastatic tumor tissue"))
##         g <- ggplot(data = temp_ggplotdata) +
##             geom_boxplot(aes(x = CancerType, y = frac, color = CancerType),
##                          size = 0.3, outlier.shape = NA, width = 0.618, lwd = 0.1) +
##             ## geom_violin(aes(x = CancerType, y = frac, color = CancerType)) +
##             geom_jitter(shape = 16, position = position_jitter(0.2),
##                         size = 0.5,  alpha = 0.9,
##                         aes(x = CancerType, y = frac, fill = CancerType,
##                             color = CancerType)) +
##             scale_color_manual(values = colorRampPalette(colorsForDataType)(length(unique(temp_ggplotdata$CancerType)))) +
##             facet_grid(g_seurat_clusters ~ TissueType,
##                        scales = "free", space = "free_x") +
##             xlab("") + ylab("") +
##             ## scale_y_log10() +
##             theme_classic() +
##             theme( text = element_text(size = 10),
##                   strip.background = element_blank(),
##                   axis.text.x = element_text(angle = 60 , vjust = 1, hjust=1),
##                   legend.position = "none")
##             ## scale_x_discrete(breaks = temp_ggplotdata$names,
##             ##                  labels = temp_ggplotdata$CancerType)
##         ggsave(file.path(tempFigurePath, paste0("CD8_figure2_sample_boxplot.pdf")), g, width = 122, height = 297/4, units = "mm")
##         ## for(temp_TissueType in unique(temp_ggplotdata$TissueType)){
##         ##     for(temp_CancerType in unique(temp_ggplotdata %>%
##         ##                                   filter(TissueType == temp_TissueType) %>%
##         ##                                   pull(CancerType))){
##         ##         g <- temp_ggplotdata %>% filter(TissueType == temp_TissueType) %>%
##         ##             filter(CancerType == temp_CancerType) %>% 
##         ##             ggstatsplot::ggbetweenstats(
##         ##                              data = .,
##         ##                              x = g_seurat_clusters,
##         ##                              y = frac,
##         ##                              ## grouping.var = TissueType,
##         ##                              xlab = "Tissue type",
##         ##                              ylab = "Sample fraction",
##         ##                              ## pairwise.display = "all", # display only significant pairwise comparisons
##         ##                              p.adjust.method = "fdr", # adjust p-values for multiple tests using this method
##         ##                              ggtheme = theme_classic()
##         ##                              ## package = "ggsci",
##         ##                              ## palette = "default_jco",
##         ##                              ## plotgrid.args = list(nrow = 4)
##         ##                          )
##         ##         ggsave(file.path(tempFigurePath, paste0(temp_TissueType, "-", temp_CancerType, "_cancertype_boxplot_gbs.pdf")), g, width = 100, height = 200, units = "mm")
##         ##     }
##         ## }

##         g <- totalObj_ggplotdata %>%
##             ggplot(data = .) +
##             geom_boxplot(aes(x = CancerType,
##                              y = frac,
##                              color = CancerType),
##                          size = 0.2, outlier.shape = NA) +
##             geom_jitter(shape=16, position=position_jitter(0.2),
##                         size = 0.2, aes(x = CancerType, y = frac,
##                                         fill = CancerType,
##                                         color = CancerType)) +
##             facet_grid(g_seurat_clusters ~ TissueType,
##                        scales = "free", space = "free_x") +
##             xlab("") + ylab("") +
##             scale_color_manual(values = colorRampPalette(colorsForDataType)(length(unique(totalObj_ggplotdata$CancerType)))) +
##             theme_classic() +
##             theme( text = element_text(size = 10),
##                   strip.background = element_blank(),
##                   axis.text.x = element_text(angle = 60 , vjust = 1, hjust=1),
##                   legend.position = "none")
##         ## scale_x_discrete(breaks = temp_ggplotdata$names,
##         ##                  labels = temp_ggplotdata$CancerType)
##         ggsave(file.path(tempFigurePath, paste0("CD8_figure2_cancertype_sample_boxplot.pdf")), g, width = 400, height = 600, units = "mm")



##         for(temp_TissueType in unique(totalObj_ggplotdata$TissueType)){

##             temp_ggplotdata <- totalObj_ggplotdata %>% filter(TissueType == temp_TissueType) %>% mutate(key = paste0(TissueType, "_", CancerType))
##             temp_ggplotdata$CancerType[temp_ggplotdata$CancerType == "LUAD"] <- "NSCLC"
##             temp_ggplotdata$CancerType[temp_ggplotdata$CancerType == "LUSC"] <- "NSCLC"
##             selectedT <- temp_ggplotdata %>% group_by(TissueType, CancerType, g_seurat_clusters) %>% dplyr::count() %>% filter(n >= 3) %>% mutate(key = paste0(TissueType, "_", CancerType))
##             temp_ggplotdata <- temp_ggplotdata %>% filter(key %in% selectedT$key)
##             temp_mean_ggplotdata <- temp_ggplotdata %>% group_by(CancerType, g_seurat_clusters) %>% summarise(meanFrac = mean(frac))
##             CancerType_names <- unique(temp_ggplotdata$CancerType)

##             psimMatrix <- matrix(rep(NA, length(unique(temp_ggplotdata$CancerType)) ^ 2),
##                                  nrow = length(unique(temp_ggplotdata$CancerType)),
##                                  ncol = length(unique(temp_ggplotdata$CancerType)))
##             rownames(psimMatrix) <- sort(unique(temp_ggplotdata$CancerType))
##             colnames(psimMatrix) <- sort(unique(temp_ggplotdata$CancerType))


##                 for(i in 1: (length(CancerType_names) - 1)){
##                     for(j in (i + 1):length(CancerType_names)){

##                         i_CancerType_name <- CancerType_names[i]
##                         j_CancerType_name <- CancerType_names[j]

##                         i_fracs <- temp_mean_ggplotdata %>% filter(CancerType == i_CancerType_name) %>% pull(meanFrac)
##                         j_fracs <- temp_mean_ggplotdata %>% filter(CancerType == j_CancerType_name) %>% pull(meanFrac)
##                         ## t.res <- t.test(i_fracs, j_fracs, paired = T, alternative = "two.sided")
##                         ## psimMatrix[i_CancerType_name, j_CancerType_name] <- t.res[[3]]
##                         ## psimMatrix[j_CancerType_name, i_CancerType_name] <- t.res[[3]]

##                         res <- lsa::cosine(i_fracs, j_fracs)
##                         psimMatrix[i_CancerType_name, j_CancerType_name] <- res
##                         psimMatrix[j_CancerType_name, i_CancerType_name] <- res
##                     }
##                 }

##             psimMatrix <- psimMatrix[!rownames(psimMatrix) %in% c("MCL", "CRC"), !colnames(psimMatrix) %in% c("MCL", "CRC")]



##             pdf(file.path(tempFigurePath, paste0(temp_TissueType, "_CancerType_pSimMatrix.pdf")))
##             print(pheatmap(psimMatrix,
##                            na_col = "#EEEEEE",
##                            show_colnames = T,
##                            show_rownames = T,
##                            cluster_cols = T,
##                            cluster_rows = T,
##                            border_color = F)
##             dev.off()


##             my.breaks0 <- seq(min(psimMatrix[!is.na(psimMatrix)]), 0.05, by=0.002)
##             my.breaks1 <- seq(0.05, max(psimMatrix[!is.na(psimMatrix)]), by=0.002)
##             my.colors <- c(
##                 colorRampPalette(colors = c(colorsForDataType[3], "white"))(length(my.breaks0)),
##                 colorRampPalette(colors = c("white", "#6DCCFF"))(length(my.breaks1)))
##             my.breaks <- c(my.breaks0, my.breaks1) ^ (1/10000)
##             pdf(file.path(tempFigurePath, paste0(temp_TissueType, "_CancerType_pSimMatrix.pdf")))
##             print(pheatmap(psimMatrix^ (1/10000),
##                            color = my.colors,
##                            breaks = my.breaks,
##                            ## annotation_row = cellType_col,
##                            legend_breaks = c(0.032, 0.05, 0.1, 0.3, 0.96) ^ (1 /10000),
##                            legend_labels = c(0.032, 0.05, 0.1, 0.3, 0.96),
##                            na_col = "#EEEEEE",
##                            show_colnames = T,
##                            show_rownames = T,
##                            cluster_cols = T,
##                            cluster_rows = T,
##                            border_color = F) +
##                   scale_fill_gradient(trans = "log",
##                                   breaks = my.breaks,
##                                   labels = my.breaks))
##             dev.off()



##             fracMatrix <- reshape2::dcast(temp_mean_ggplotdata, CancerType ~ g_seurat_clusters)
##             rownames(fracMatrix) <- fracMatrix$CancerType
##             fracMatrix <- as.matrix(fracMatrix[, 2:dim(fracMatrix)[2]])

##             fracMatrix <- fracMatrix[!rownames(fracMatrix) %in% c("MCL", "CRC"), !colnames(fracMatrix) %in% c("MCL", "CRC")]

##             pdf(file.path(tempFigurePath, paste0(temp_TissueType, "frac_heatmap.pdf")))
##             print(pheatmap(fracMatrix,
##                            na_col = "#EEEEEE",
##                            show_colnames = T,
##                            show_rownames = T,
##                            cluster_cols = T,
##                            cluster_rows = T,
##                            border_color = F))
##             dev.off()

##             library(reshape2)
##             library(gtools)
##             library(ggdendro)
##             ssimMatrix <- melt(psimMatrix)
##             ssimMatrix$pstars <- stars.pval(ssimMatrix$value)
##             ssimMatrix$pstarslevel <- 1
##             ssimMatrix$pstarslevel[ssimMatrix$pstars == "."] <- 2
##             ssimMatrix$pstarslevel[ssimMatrix$pstars == "*"] <- 3
##             ssimMatrix$pstarslevel[ssimMatrix$pstars == "**"] <- 4
##             ssimMatrix$pstarslevel[ssimMatrix$pstars == "***"] <- 5
##             ssimMatrix$pstarslevel[ssimMatrix$pstars == ""] <- NA
##             CM <- psimMatrix
##             CM.dendro <- as.dendrogram(hclust(d = dist(x = CM)))
##             CM.plot <- ggdendrogram(data = CM.dendro, rotate = FALSE)
##             pdf(file.path(tempFigurePath, paste0(temp_TissueType, "_CancerType_pSimMatrix_hc.pdf")))
##             print(CM.plot)
##             dev.off()
##             row.order <- order.dendrogram(CM.dendro)

##             ssimMatrix$Var1 <- factor(x = ssimMatrix$Var1,
##                                       levels = rev(rownames(psimMatrix)[row.order]),
##                                       ordered = TRUE)
##             ssimMatrix$Var2 <- factor(x = ssimMatrix$Var2,
##                                       levels = colnames(psimMatrix)[row.order],
##                                       ordered = TRUE)

##             if(abs(min(ssimMatrix$value[!is.na(ssimMatrix$value)])) >
##                abs(max(ssimMatrix$value[!is.na(ssimMatrix$value)]))){
##                 colorRadius <- abs(min(ssimMatrix$value))
##             }else{
##                 colorRadius <- abs(max(ssimMatrix$value))
##             }
##             g <- ggplot(ssimMatrix, aes(y = Var1, x = Var2)) +
##                 geom_tile(aes(fill = value),
##                            color = "black", stroke = 0.1, shape = 22) +
##                 scale_fill_gradientn(name = "count", trans = "log",
##                                     breaks = c(0.01, 0.05, 0.1, 0.3, 1),
##                                     labels = c(0.01, 0.05, 0.1, 0.3, 1)) +
##                 theme_void() +
##                 theme(text = element_text(size = 10),
##                       legend.position = "bottom",
##                       strip.background = element_blank(),
##                       axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
##                       legend.box="vertical", legend.margin=margin())
##             ggsave(file.path(tempFigurePath,paste0(temp_TissueType, "_CancerType_pSimMatrix_point.pdf")), g, width = 210, height = 297/1.52, units = "mm")



##             CancerType_cell_number_frac <- allMD %>% filter(DataSet == "CD8") %>% mutate(g_seurat_clusters = paste0(DataSet, "_c", seurat_clusters)) %>% filter(TissueType == temp_TissueType) %>% filter(CancerType %in% unique(temp_ggplotdata$CancerType)) %>% select(g_seurat_clusters, CancerType, barcode, DataSet) %>% group_by(CancerType, g_seurat_clusters) %>% dplyr::count() %>% group_by(CancerType) %>% mutate(frac = n / sum(n))
##             g <- ggplot(CancerType_cell_number_frac) +
##                 geom_bar(aes(x = CancerType, y = frac,
##                              fill = g_seurat_clusters),
##                          color = "white",
##                          position="stack", stat="identity") +
##                 theme_classic()
##             ggsave(file.path(tempFigurePath, paste0(temp_TissueType, "_stackfrac_cancertype.pdf")), g)


##             fracMatrix <- matrix(rep(0, length(unique(CancerType_cell_number_frac$CancerType)) *
##                                         length(unique(CancerType_cell_number_frac$g_seurat_clusters))),
##                                  ncol = length(unique(CancerType_cell_number_frac$CancerType)),
##                                  nrow = length(unique(CancerType_cell_number_frac$g_seurat_clusters)))
##             rownames(fracMatrix) <- sort(unique(CancerType_cell_number_frac$g_seurat_clusters))
##             colnames(fracMatrix) <- sort(unique(CancerType_cell_number_frac$CancerType))
##             for(temp_rn in sort(unique(CancerType_cell_number_frac$g_seurat_clusters))){
##                 for(temp_cl in sort(unique(CancerType_cell_number_frac$CancerType))){
##                     if(!any(CancerType_cell_number_frac$CancerType == temp_cl & CancerType_cell_number_frac$g_seurat_clusters == temp_rn)){
##                         next
##                     }
##                     fracMatrix[temp_rn, temp_cl] <- CancerType_cell_number_frac$frac[CancerType_cell_number_frac$CancerType == temp_cl & CancerType_cell_number_frac$g_seurat_clusters == temp_rn]
##                 }
##             }
##             my.breaks <- seq(min(fracMatrix[!is.na(fracMatrix)]), max(fracMatrix[!is.na(fracMatrix)]), by=0.002)
##             my.colors <- c(
##                 colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)/2),
##                 colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks)/2))
##             pdf(file.path(tempFigurePath, paste0(temp_TissueType, "_CancerType_cluster_fraction_heatmap.pdf")))
##             print(pheatmap(fracMatrix,
##                            color = my.colors,
##                            breaks = my.breaks,
##                            ## annotation_row = cellType_col,
##                            show_colnames = T,
##                            show_rownames = T,
##                            cluster_cols = F,
##                            cluster_rows = F,
##                            border_color = F))
##             dev.off()
##         }


##         #' plus skmc  #########################################################
##         temp_ggplotdata <- totalObj_ggplotdata %>% filter(TissueType == "Primary tumor tissue") %>% mutate(key = paste0(TissueType, "_", CancerType))
##         temp_ggplotdata$CancerType[temp_ggplotdata$CancerType == "LUAD"] <- "NSCLC"
##         temp_ggplotdata$CancerType[temp_ggplotdata$CancerType == "LUSC"] <- "NSCLC"
##         SKCM_temp_ggplotdata <- totalObj_ggplotdata %>% filter(CancerType == "SKCM") %>% mutate(key = paste0(TissueType, "_", CancerType))
##         temp_ggplotdata <- bind_rows(temp_ggplotdata, SKCM_temp_ggplotdata)
##         selectedT <- temp_ggplotdata %>% group_by(TissueType, CancerType, g_seurat_clusters) %>% dplyr::count() %>% filter(n >= 3) %>% mutate(key = paste0(TissueType, "_", CancerType))
##         temp_ggplotdata <- temp_ggplotdata %>% filter(key %in% selectedT$key)
##         temp_mean_ggplotdata <- temp_ggplotdata %>% group_by(CancerType, g_seurat_clusters) %>% summarise(meanFrac = mean(frac))
##         CancerType_names <- unique(temp_ggplotdata$CancerType)
##         psimMatrix <- matrix(rep(NA, length(unique(temp_ggplotdata$CancerType)) ^ 2),
##                              nrow = length(unique(temp_ggplotdata$CancerType)),
##                              ncol = length(unique(temp_ggplotdata$CancerType)))
##         rownames(psimMatrix) <- sort(unique(temp_ggplotdata$CancerType))
##         colnames(psimMatrix) <- sort(unique(temp_ggplotdata$CancerType))
##         for(i in 1: (length(CancerType_names) - 1)){
##             for(j in (i + 1):length(CancerType_names)){
##                 i_CancerType_name <- CancerType_names[i]
##                 j_CancerType_name <- CancerType_names[j]
##                 i_fracs <- temp_mean_ggplotdata %>% filter(CancerType == i_CancerType_name) %>% pull(meanFrac)
##                 j_fracs <- temp_mean_ggplotdata %>% filter(CancerType == j_CancerType_name) %>% pull(meanFrac)
##                 t.res <- t.test(i_fracs, j_fracs, paired = TRUE, alternative = "two.sided")
##                 psimMatrix[i_CancerType_name, j_CancerType_name] <- t.res[[3]]
##                 psimMatrix[j_CancerType_name, i_CancerType_name] <- t.res[[3]]
##             }
##         }

##         psimMatrix <- psimMatrix[rownames(psimMatrix) != "MCL", colnames(psimMatrix) != "MCL"]
##         my.breaks <- seq(min(psimMatrix[!is.na(psimMatrix)]), max(psimMatrix[!is.na(psimMatrix)]), by=0.002)
##         my.colors <- c(
##             colorRampPalette(colors = c("#FD9AA0", "white"))(length(my.breaks)/2),
##             colorRampPalette(colors = c("white", "#6DCCFD"))(length(my.breaks)/2))
##         pdf(file.path(tempFigurePath, paste0("CancerType_pSimMatrix.pdf")))
##         print(pheatmap(psimMatrix,
##                        color = my.colors,
##                        breaks = my.breaks,
##                        ## annotation_row = cellType_col,
##                        show_colnames = T,
##                        show_rownames = T,
##                        cluster_cols = T,
##                        cluster_rows = T,
##                        border_color = F))
##         dev.off()

##     }

##     if(tempName == "CD4"){

##         temp_ggplotdata <- totalObj_ggplotdata %>%
##             filter(g_seurat_clusters %in% c("CD4_c1", "CD4_c3", "CD4_c4", "CD4_c8")) %>% mutate(key = paste0(TissueType, "_", CancerType))
##         selectedT <- temp_ggplotdata %>% group_by(TissueType, CancerType, g_seurat_clusters) %>% dplyr::count() %>% filter(n >= 3) %>% mutate(key = paste0(TissueType, "_", CancerType))
##         temp_ggplotdata <- temp_ggplotdata %>% filter(key %in% selectedT$key)
##         temp_ggplotdata$names <- paste(temp_ggplotdata$TissueType, temp_ggplotdata$CancerType)
##         temp_ggplotdata$CancerType <- factor(temp_ggplotdata$CancerType,
##                                              levels = c("Lymph node",
##                                                         "Breast",
##                                                         "Head and neck",
##                                                         "Lung",
##                                                         "Blood",
##                                                         "Brain",
##                                                         "BM",
##                                                         "colon",
##                                                         "LUSC",
##                                                         "HNSC",
##                                                         "BCC",
##                                                         "LUAD",
##                                                         "FL",
##                                                         "NSCLC",
##                                                         "BRCA",
##                                                         "PDAC",
##                                                         "SKCM",
##                                                         "OV",
##                                                         "STAD",
##                                                         "HCC",
##                                                         "GBM",
##                                                         "MCL",
##                                                         "AML",
##                                                         "LBCL",
##                                                         "CRC"))
##         temp_ggplotdata$g_seurat_clusters <- factor(temp_ggplotdata$g_seurat_clusters, levels = c("CD4_c1", "CD4_c3", "CD4_c4", "CD4_c8"))
##         temp_ggplotdata <- temp_ggplotdata %>% mutate(TissueType = fct_recode(TissueType, H = "Healthy donor", U = "Uninvolved normal tissue", P = "Primary tumor tissue", M = "Metastatic tumor tissue"))
##         g <- ggplot(data = temp_ggplotdata) +
##             geom_boxplot(aes(x = CancerType, y = frac, color = CancerType),
##                          size = 0.3, outlier.shape = NA, width = 0.618, lwd = 0.1) +
##             ## geom_violin(aes(x = CancerType, y = frac, color = CancerType)) +
##             geom_jitter(shape = 16, position = position_jitter(0.2),
##                         size = 0.5,  alpha = 0.9,
##                         aes(x = CancerType, y = frac, fill = CancerType,
##                             color = CancerType)) +
##             scale_color_manual(values = colorRampPalette(colorsForDataType)(length(unique(temp_ggplotdata$CancerType)))) +
##             facet_grid(g_seurat_clusters ~ TissueType,
##                        scales = "free", space = "free_x") +
##             xlab("") + ylab("") +
##             ## scale_y_log10() +
##             theme_classic() +
##             theme( text = element_text(size = 10),
##                   strip.background = element_blank(),
##                   axis.text.x = element_text(angle = 60 , vjust = 1, hjust=1),
##                   legend.position = "none")
##         ## scale_x_discrete(breaks = temp_ggplotdata$names,
##         ##                  labels = temp_ggplotdata$CancerType)
##         ggsave(file.path(tempFigurePath, paste0("CD4_figure2_sample_boxplot.pdf")), g, width = 122, height = 297/4, units = "mm")

##         for(temp_TissueType in unique(temp_ggplotdata$TissueType)){
##             for(temp_CancerType in unique(temp_ggplotdata %>%
##                                           filter(TissueType == temp_TissueType) %>%
##                                           pull(CancerType))){
##                 if(all(temp_ggplotdata %>% filter(TissueType == temp_TissueType) %>%
##                        filter(CancerType == temp_CancerType) %>% pull(frac) == 0)){
##                     next
##                 }
##                 g <- temp_ggplotdata %>% filter(TissueType == temp_TissueType) %>%
##                     filter(CancerType == temp_CancerType) %>%
##                     ggstatsplot::ggbetweenstats(
##                                      data = .,
##                                      x = g_seurat_clusters,
##                                      y = frac,
##                                      ## grouping.var = TissueType,
##                                      xlab = "Tissue type",
##                                      ylab = "Sample fraction",
##                                      ## pairwise.display = "all", # display only significant pairwise comparisons
##                                      p.adjust.method = "fdr", # adjust p-values for multiple tests using this method
##                                      ggtheme = theme_classic()
##                                      ## package = "ggsci",
##                                      ## palette = "default_jco",
##                                      ## plotgrid.args = list(nrow = 4)
##                                  )
##                 ggsave(file.path(tempFigurePath, paste0(temp_TissueType, "-", temp_CancerType, "_cancertype_boxplot_gbs.pdf")), g, width = 100, height = 200, units = "mm")
##             }
##         }


##         g <- totalObj_ggplotdata %>%
##             ggplot(data = .) +
##             geom_boxplot(aes(x = CancerType,
##                              y = frac,
##                              color = CancerType),
##                          size = 0.2, outlier.shape = NA) +
##             geom_jitter(shape=16, position=position_jitter(0.2),
##                         size = 0.2, aes(x = CancerType, y = frac,
##                                         fill = CancerType,
##                                         color = CancerType)) +
##             facet_grid(g_seurat_clusters ~ TissueType,
##                        scales = "free", space = "free_x") +
##             xlab("") + ylab("") +
##             scale_color_manual(values = colorRampPalette(colorsForDataType)(length(unique(totalObj_ggplotdata$CancerType)))) +
##             theme_classic() +
##             theme( text = element_text(size = 10),
##                   strip.background = element_blank(),
##                   axis.text.x = element_text(angle = 60 , vjust = 1, hjust=1),
##                   legend.position = "none")
##         ## scale_x_discrete(breaks = temp_ggplotdata$names,
##         ##                  labels = temp_ggplotdata$CancerType)
##         ggsave(file.path(tempFigurePath, paste0("CD8_figure2_cancertype_sample_boxplot.pdf")), g, width = 400, height = 600, units = "mm")


##         for(temp_TissueType in unique(totalObj_ggplotdata$TissueType)){
##             temp_ggplotdata <- totalObj_ggplotdata %>% filter(TissueType == temp_TissueType) %>% mutate(key = paste0(TissueType, "_", CancerType))
##             temp_ggplotdata$CancerType[temp_ggplotdata$CancerType == "LUAD"] <- "NSCLC"
##             temp_ggplotdata$CancerType[temp_ggplotdata$CancerType == "LUSC"] <- "NSCLC"
##             selectedT <- temp_ggplotdata %>% group_by(TissueType, CancerType, g_seurat_clusters) %>% dplyr::count() %>% filter(n >= 3) %>% mutate(key = paste0(TissueType, "_", CancerType))
##             temp_ggplotdata <- temp_ggplotdata %>% filter(key %in% selectedT$key)
##             temp_mean_ggplotdata <- temp_ggplotdata %>% group_by(CancerType, g_seurat_clusters) %>% summarise(meanFrac = mean(frac))
##             CancerType_names <- unique(temp_ggplotdata$CancerType)
##             psimMatrix <- matrix(rep(NA, length(unique(temp_ggplotdata$CancerType)) ^ 2),
##                                  nrow = length(unique(temp_ggplotdata$CancerType)),
##                                  ncol = length(unique(temp_ggplotdata$CancerType)))
##             rownames(psimMatrix) <- sort(unique(temp_ggplotdata$CancerType))
##             colnames(psimMatrix) <- sort(unique(temp_ggplotdata$CancerType))
##                 for(i in 1: (length(CancerType_names) - 1)){
##                     for(j in (i + 1):length(CancerType_names)){
##                         i_CancerType_name <- CancerType_names[i]
##                         j_CancerType_name <- CancerType_names[j]
##                         i_fracs <- temp_mean_ggplotdata %>% filter(CancerType == i_CancerType_name) %>% pull(meanFrac)
##                         j_fracs <- temp_mean_ggplotdata %>% filter(CancerType == j_CancerType_name) %>% pull(meanFrac)
##                         t.res <- t.test(i_fracs, j_fracs, paired = TRUE, alternative = "two.sided")
##                         psimMatrix[i_CancerType_name, j_CancerType_name] <- t.res[[3]]
##                         psimMatrix[j_CancerType_name, i_CancerType_name] <- t.res[[3]]
##                     }
##                 }

##             psimMatrix <- psimMatrix[rownames(psimMatrix) != "MCL", colnames(psimMatrix) != "MCL"]
##             my.breaks0 <- seq(min(psimMatrix[!is.na(psimMatrix)]), 0.05, by=0.002)
##             my.breaks1 <- seq(0.05, max(psimMatrix[!is.na(psimMatrix)]), by=0.002)
##             my.colors <- c(
##                 colorRampPalette(colors = c("red", "white"))(length(my.breaks0)),
##                 colorRampPalette(colors = c("white", "blue"))(length(my.breaks1)))
##             my.breaks = c(my.breaks0, my.breaks1)
##             pdf(file.path(tempFigurePath, paste0(temp_TissueType, "_CancerType_pSimMatrix.pdf")))
##             print(pheatmap(psimMatrix,
##                            color = my.colors,
##                            breaks = my.breaks,
##                            ## annotation_row = cellType_col,
##                            show_colnames = T,
##                            show_rownames = T,
##                            cluster_cols = T,
##                            cluster_rows = T,
##                            border_color = F))
##             dev.off()

##             CancerType_cell_number_frac <- allMD %>% filter(DataSet == "CD8") %>% mutate(g_seurat_clusters = paste0(DataSet, "_c", seurat_clusters)) %>% filter(TissueType == temp_TissueType) %>% filter(CancerType %in% unique(temp_ggplotdata$CancerType)) %>% select(g_seurat_clusters, CancerType, barcode, DataSet) %>% group_by(CancerType, g_seurat_clusters) %>% dplyr::count() %>% group_by(CancerType) %>% mutate(frac = n / sum(n))
##             g <- ggplot(CancerType_cell_number_frac) +
##                 geom_bar(aes(x = CancerType, y = frac,
##                              fill = g_seurat_clusters),
##                          color = "white",
##                          position="stack", stat="identity") +
##                 theme_classic()
##             ggsave(file.path(tempFigurePath, paste0(temp_TissueType, "_stackfrac_cancertype.pdf")), g)


##             fracMatrix <- matrix(rep(0, length(unique(CancerType_cell_number_frac$CancerType)) *
##                                         length(unique(CancerType_cell_number_frac$g_seurat_clusters))),
##                                  ncol = length(unique(CancerType_cell_number_frac$CancerType)),
##                                  nrow = length(unique(CancerType_cell_number_frac$g_seurat_clusters)))
##             rownames(fracMatrix) <- sort(unique(CancerType_cell_number_frac$g_seurat_clusters))
##             colnames(fracMatrix) <- sort(unique(CancerType_cell_number_frac$CancerType))
##             for(temp_rn in sort(unique(CancerType_cell_number_frac$g_seurat_clusters))){
##                 for(temp_cl in sort(unique(CancerType_cell_number_frac$CancerType))){
##                     if(!any(CancerType_cell_number_frac$CancerType == temp_cl & CancerType_cell_number_frac$g_seurat_clusters == temp_rn)){
##                         next
##                     }
##                     fracMatrix[temp_rn, temp_cl] <- CancerType_cell_number_frac$frac[CancerType_cell_number_frac$CancerType == temp_cl & CancerType_cell_number_frac$g_seurat_clusters == temp_rn]
##                 }
##             }
##             my.breaks <- seq(min(fracMatrix[!is.na(fracMatrix)]), max(fracMatrix[!is.na(fracMatrix)]), by=0.002)
##             my.colors <- c(
##                 colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)/2),
##                 colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks)/2))
##             pdf(file.path(tempFigurePath, paste0(temp_TissueType, "_CancerType_cluster_fraction_heatmap.pdf")))
##             print(pheatmap(fracMatrix,
##                            color = my.colors,
##                            breaks = my.breaks,
##                            ## annotation_row = cellType_col,
##                            show_colnames = T,
##                            show_rownames = T,
##                            cluster_cols = F,
##                            cluster_rows = F,
##                            border_color = F))
##             dev.off()
##         }


##         #' plus skmc  #########################################################
##         temp_ggplotdata <- totalObj_ggplotdata %>% filter(TissueType == "Primary tumor tissue") %>% mutate(key = paste0(TissueType, "_", CancerType))
##         temp_ggplotdata$CancerType[temp_ggplotdata$CancerType == "LUAD"] <- "NSCLC"
##         temp_ggplotdata$CancerType[temp_ggplotdata$CancerType == "LUSC"] <- "NSCLC"
##         SKCM_temp_ggplotdata <- totalObj_ggplotdata %>% filter(CancerType == "SKCM") %>% mutate(key = paste0(TissueType, "_", CancerType))
##         temp_ggplotdata <- bind_rows(temp_ggplotdata, SKCM_temp_ggplotdata)
##         selectedT <- temp_ggplotdata %>% group_by(TissueType, CancerType, g_seurat_clusters) %>% dplyr::count() %>% filter(n >= 3) %>% mutate(key = paste0(TissueType, "_", CancerType))
##         temp_ggplotdata <- temp_ggplotdata %>% filter(key %in% selectedT$key)
##         temp_mean_ggplotdata <- temp_ggplotdata %>% group_by(CancerType, g_seurat_clusters) %>% summarise(meanFrac = mean(frac))
##         CancerType_names <- unique(temp_ggplotdata$CancerType)
##         psimMatrix <- matrix(rep(NA, length(unique(temp_ggplotdata$CancerType)) ^ 2),
##                              nrow = length(unique(temp_ggplotdata$CancerType)),
##                              ncol = length(unique(temp_ggplotdata$CancerType)))
##         rownames(psimMatrix) <- sort(unique(temp_ggplotdata$CancerType))
##         colnames(psimMatrix) <- sort(unique(temp_ggplotdata$CancerType))
##         for(i in 1: (length(CancerType_names) - 1)){
##             for(j in (i + 1):length(CancerType_names)){
##                 i_CancerType_name <- CancerType_names[i]
##                 j_CancerType_name <- CancerType_names[j]
##                 i_fracs <- temp_mean_ggplotdata %>% filter(CancerType == i_CancerType_name) %>% pull(meanFrac)
##                 j_fracs <- temp_mean_ggplotdata %>% filter(CancerType == j_CancerType_name) %>% pull(meanFrac)
##                 t.res <- t.test(i_fracs, j_fracs, paired = TRUE, alternative = "two.sided")
##                 psimMatrix[i_CancerType_name, j_CancerType_name] <- t.res[[3]]
##                 psimMatrix[j_CancerType_name, i_CancerType_name] <- t.res[[3]]
##             }
##         }

##         psimMatrix <- psimMatrix[rownames(psimMatrix) != "MCL", colnames(psimMatrix) != "MCL"]
##         my.breaks <- seq(min(psimMatrix[!is.na(psimMatrix)]), max(psimMatrix[!is.na(psimMatrix)]), by=0.002)
##         my.colors <- c(
##             colorRampPalette(colors = c("#FD9AA0", "white"))(length(my.breaks)/2),
##             colorRampPalette(colors = c("white", "#6DCCFD"))(length(my.breaks)/2))
##         pdf(file.path(tempFigurePath, paste0("CancerType_pSimMatrix.pdf")))
##         print(pheatmap(psimMatrix,
##                        color = my.colors,
##                        breaks = my.breaks,
##                        ## annotation_row = cellType_col,
##                        show_colnames = T,
##                        show_rownames = T,
##                        cluster_cols = T,
##                        cluster_rows = T,
##                        border_color = F))
##         dev.off()
##     }
##     ## if(tempName == "NK"){
##     ##     TissueType.list <- list(NK = c("Innate_c0", "Innate_c2", "Innate_c5"),
##     ##                             NKT = c("NKT_c3", "NKT_c7", "NKT_c8"),
##     ##                             MAIT = c("MAIT_c4", "MAIT-like_c1", "MAIT-like_c6"))
##     ##     for(TTli in 1:length(TissueType.list)){
##     ##         tempTissueTypeName <- names(TissueType.list[TTli])
##     ##         tempTissueTypes <- TissueType.list[[TTli]]
##     ##         temp_ggplotdata <- totalObj_ggplotdata %>% filter(g_seurat_clusters %in% tempTissueTypes) %>% mutate(key = paste0(TissueType, "_", CancerType))
##     ##         selectedT <- temp_ggplotdata %>% group_by(TissueType, CancerType, g_seurat_clusters) %>% dplyr::count() %>% filter(n >= 3) %>% mutate(key = paste0(TissueType, "_", CancerType))
##     ##         temp_ggplotdata <- temp_ggplotdata %>% filter(key %in% selectedT$key)
##     ##         temp_ggplotdata$names <- paste0(temp_ggplotdata$TissueType, "_", temp_ggplotdata$CancerType)
##     ##         if(tempTissueTypeName == "NK"){
##     ##             temp_ggplotdata$CancerType <- factor(temp_ggplotdata$CancerType,
##     ##                                                  levels = c(
##     ##                                                      "Lung",
##     ##                                                      "Breast",
##     ##                                                      "Head and neck",
##     ##                                                      "Blood",
##     ##                                                      "BM",
##     ##                                                      "Lymph node",
##     ##                                                      "colon",
##     ##                                                      "BRCA",
##     ##                                                      "HNSC",
##     ##                                                      "HCC",
##     ##                                                      "NSCLC",
##     ##                                                      "OV",
##     ##                                                      "STAD",
##     ##                                                      "SKCM",
##     ##                                                      "GBM",
##     ##                                                      "AML",
##     ##                                                      "PDAC",
##     ##                                                      "LBCL",
##     ##                                                      "LUSC",
##     ##                                                      "FL",
##     ##                                                      "BCC",
##     ##                                                      "CRC",
##     ##                                                      "LUAD",
##     ##                                                      "MCL"
##     ##                                                  ))
##     ##             temp_ggplotdata$g_seurat_clusters <- factor(temp_ggplotdata$g_seurat_clusters, levels = c("Innate_c0", "Innate_c2", "Innate_c5"))
##     ##             temp_ggplotdata <- temp_ggplotdata %>%
##     ##                 filter(names %in% c("Healthy donor_BM",
##     ##                                       "Uninvolved normal tissue_Lung",
##     ##                                       "Primary tumor tissue_BRCA",
##     ##                                       "Primary tumor tissue_HNSC",
##     ##                                       "Primary tumor tissue_HCC",
##     ##                                       "Primary tumor tissue_NSCLC",
##     ##                                       "Primary tumor tissue_GBM",
##     ##                                       "Primary tumor tissue_AML",
##     ##                                       "Primary tumor tissue_MCL",
##     ##                                       "Metastatic tumor tissue_SKCM"))
##     ##         }
##     ##         if(tempTissueTypeName == "NKT"){
##     ##             temp_ggplotdata <- temp_ggplotdata %>%
##     ##                 filter(names %in% c(
##     ##                                       "Healthy donor_BM",
##     ##                                       "Uninvolved normal tissue_Lung",
##     ##                                       "Primary tumor tissue_NSCLC",
##     ##                                       "Metastatic tumor tissue_SKCM"))
##     ##         }
##     ##         if(tempTissueTypeName == "MAIT"){
##     ##             temp_ggplotdata <- temp_ggplotdata %>%
##     ##                 filter(names %in% c(
##     ##                                       "Healthy donor_BM",
##     ##                                       "Uninvolved normal tissue_Lung",
##     ##                                       "Primary tumor tissue_GBM",
##     ##                                       "Primary tumor tissue_HNSC",
##     ##                                       "Primary tumor tissue_NSCLC",
##     ##                                       "Primary tumor tissue_LUAD",
##     ##                                       "Primary tumor tissue_MCL"))
##     ##         }

##     ##         ## temp_ggplotdata$CancerType <- factor(temp_ggplotdata$CancerType,
##     ##         ##                                 levels = c("BM",
##     ##         ##                                            "Lymph node",
##     ##         ##                                            "Brain",
##     ##         ##                                            "Breast",
##     ##         ##                                            "Lung",
##     ##         ##                                            "Blood",
##     ##         ##                                            "colon",
##     ##         ##                                            "Head and neck",
##     ##         ##                                            "SKCM",
##     ##         ##                                            "OV",
##     ##         ##                                            "NSCLC",
##     ##         ##                                            "HNSC",
##     ##         ##                                            "LBCL",
##     ##         ##                                            "LUSC",
##     ##         ##                                            "BRCA",
##     ##         ##                                            "FL",
##     ##         ##                                            "LUAD",
##     ##         ##                                            "STAD",
##     ##         ##                                            "BCC",
##     ##         ##                                            "HCC",
##     ##         ##                                            "GBM",
##     ##         ##                                            "MCL",
##     ##         ##                                            "AML",
##     ##         ##                                            "CRC",
##     ##         ##                                            "PDAC"))
##     ##         ## temp_ggplotdata$g_seurat_clusters <- factor(temp_ggplotdata$g_seurat_clusters, levels = c("CD8_c3", "CD8_c0", "CD8_c2", "CD8_c1"))

##     ##         temp_ggplotdata <- temp_ggplotdata %>% mutate(TissueType = fct_recode(TissueType, H = "Healthy donor", U = "Uninvolved normal tissue", P = "Primary tumor tissue", M = "Metastatic tumor tissue"))

##     ##         if(tempTissueTypeName == "NK"){
##     ##             g <- ggplot(data = temp_ggplotdata) +
##     ##                 geom_boxplot(aes(x = CancerType,
##     ##                                  y = frac, color = CancerType),
##     ##                              size = 0.3, outlier.shape = NA, width = 0.618, lwd = 0.1)
##     ##         }else{
##     ##             g <- ggplot(data = temp_ggplotdata) +
##     ##                 geom_boxplot(aes(x = reorder(CancerType, -frac, FUN=median),
##     ##                                  y = frac, color = CancerType),
##     ##                              size = 0.3, outlier.shape = NA, width = 0.618, lwd = 0.1)
##     ##         }
##     ##         g <- g +
##     ##             ## geom_violin(aes(x = CancerType, y = frac, color = CancerType)) +
##     ##             geom_jitter(shape = 16, position = position_jitter(0.2),
##     ##                         size = 0.5,  alpha = 0.9,
##     ##                         aes(x = CancerType, y = frac, fill = CancerType,
##     ##                             color = CancerType)) +
##     ##             scale_color_manual(values = colorRampPalette(colorsForDataType)(length(unique(temp_ggplotdata$CancerType)))) +
##     ##             facet_grid(g_seurat_clusters ~ TissueType,
##     ##                        scales = "free", space = "free_x") +
##     ##             xlab("") + ylab("") +
##     ##             ## scale_y_log10() +
##     ##             theme_classic() +
##     ##             theme( text = element_text(size = 10),
##     ##                   strip.background = element_blank(),
##     ##                   axis.text.x = element_text(angle = 60 , vjust = 1, hjust=1),
##     ##                   legend.position = "none")

##     ##         ## scale_x_discrete(breaks = temp_ggplotdata$names,
##     ##         ##                  labels = temp_ggplotdata$CancerType)
##     ##         ggsave(file.path(tempFigurePath, paste0(tempTissueTypeName, "_innate_figure2_sample_boxplot.pdf")), g, width = (100 - 30) * (length(unique(temp_ggplotdata$names))/21) + 30, height = 297/4, units = "mm")
##     ##         for(temp_TissueType in unique(temp_ggplotdata$TissueType)){
##     ##             for(temp_CancerType in unique(temp_ggplotdata %>%
##     ##                                           filter(TissueType == temp_TissueType) %>%
##     ##                                           pull(CancerType))){
##     ##                 if(all(temp_ggplotdata %>% filter(TissueType == temp_TissueType) %>%
##     ##                        filter(CancerType == temp_CancerType) %>% pull(frac) == 0)){next}
##     ##                 g <- temp_ggplotdata %>% filter(TissueType == temp_TissueType) %>%
##     ##                     filter(CancerType == temp_CancerType) %>%
##     ##                     ggstatsplot::ggbetweenstats(
##     ##                                      data = .,
##     ##                                      x = g_seurat_clusters,
##     ##                                      y = frac,
##     ##                                      ## grouping.var = TissueType,
##     ##                                      xlab = "Tissue type",
##     ##                                      ylab = "Sample fraction",
##     ##                                      ## pairwise.display = "all", # display only significant pairwise comparisons
##     ##                                      p.adjust.method = "fdr", # adjust p-values for multiple tests using this method
##     ##                                      ggtheme = theme_classic()
##     ##                                      ## package = "ggsci",
##     ##                                      ## palette = "default_jco",
##     ##                                      ## plotgrid.args = list(nrow = 4)
##     ##                                  )
##     ##                 ggsave(file.path(tempFigurePath, paste0(tempTissueTypeName, "_", temp_TissueType, "-", temp_CancerType, "_cancertype_boxplot_gbs.pdf")), g, width = 100, height = 200, units = "mm")
##     ##             }
##     ##         }
##     ##     }



##     ##     for(temp_TissueType in unique(totalObj_ggplotdata$TissueType)){
##     ##         temp_ggplotdata <- totalObj_ggplotdata %>% filter(TissueType == temp_TissueType) %>% mutate(key = paste0(TissueType, "_", CancerType))
##     ##         selectedT <- temp_ggplotdata %>% group_by(TissueType, CancerType, g_seurat_clusters) %>% dplyr::count() %>% filter(n >= 3) %>% mutate(key = paste0(TissueType, "_", CancerType))
##     ##         temp_ggplotdata <- temp_ggplotdata %>% filter(key %in% selectedT$key)
##     ##         temp_mean_ggplotdata <- temp_ggplotdata %>% group_by(CancerType, g_seurat_clusters) %>% summarise(meanFrac = mean(frac))
##     ##         CancerType_names <- unique(temp_ggplotdata$CancerType)
##     ##         psimMatrix <- matrix(rep(NA, length(unique(temp_ggplotdata$CancerType)) ^ 2),
##     ##                              nrow = length(unique(temp_ggplotdata$CancerType)),
##     ##                              ncol = length(unique(temp_ggplotdata$CancerType)))
##     ##         rownames(psimMatrix) <- sort(unique(temp_ggplotdata$CancerType))
##     ##         colnames(psimMatrix) <- sort(unique(temp_ggplotdata$CancerType))
##     ##             for(i in 1: (length(CancerType_names) - 1)){
##     ##                 for(j in (i + 1):length(CancerType_names)){
##     ##                     i_CancerType_name <- CancerType_names[i]
##     ##                     j_CancerType_name <- CancerType_names[j]
##     ##                     i_fracs <- temp_mean_ggplotdata %>% filter(CancerType == i_CancerType_name) %>% pull(meanFrac)
##     ##                     j_fracs <- temp_mean_ggplotdata %>% filter(CancerType == j_CancerType_name) %>% pull(meanFrac)
##     ##                     t.res <- t.test(i_fracs, j_fracs, paired = TRUE, alternative = "two.sided")
##     ##                     psimMatrix[i_CancerType_name, j_CancerType_name] <- t.res[[3]]
##     ##                     psimMatrix[j_CancerType_name, i_CancerType_name] <- t.res[[3]]
##     ##                 }
##     ##             }
##     ##         my.breaks <- seq(min(psimMatrix[!is.na(psimMatrix)]), max(psimMatrix[!is.na(psimMatrix)]), by=0.002)
##     ##         my.colors <- c(
##     ##             colorRampPalette(colors = c("#FD9AA0", "white"))(length(my.breaks)/2),
##     ##             colorRampPalette(colors = c("white", "#6DCCFD"))(length(my.breaks)/2))
##     ##         pdf(file.path(tempFigurePath, paste0(temp_TissueType, "_CancerType_pSimMatrix.pdf")))
##     ##         print(pheatmap(psimMatrix,
##     ##                        color = my.colors,
##     ##                        breaks = my.breaks,
##     ##                        ## annotation_row = cellType_col,
##     ##                        show_colnames = T,
##     ##                        show_rownames = T,
##     ##                        cluster_cols = T,
##     ##                        cluster_rows = T,
##     ##                        border_color = F))
##     ##         dev.off()

##     ##         CancerType_cell_number_frac <- allMD %>% filter(DataSet == "NK") %>% mutate(g_seurat_clusters = paste0(DataSet, "_c", seurat_clusters)) %>% filter(TissueType == temp_TissueType) %>% filter(CancerType %in% unique(temp_ggplotdata$CancerType)) %>% select(g_seurat_clusters, CancerType, barcode, DataSet) %>% group_by(CancerType, g_seurat_clusters) %>% dplyr::count() %>% group_by(CancerType) %>% mutate(frac = n / sum(n))
##     ##         g <- ggplot(CancerType_cell_number_frac) +
##     ##             geom_bar(aes(x = CancerType, y = frac,
##     ##                          fill = g_seurat_clusters),
##     ##                      color = "white",
##     ##                      position="stack", stat="identity") +
##     ##             theme_classic()
##     ##         ggsave(file.path(tempFigurePath, paste0(temp_TissueType, "_stackfrac_cancertype.pdf")), g)


##     ##         fracMatrix <- matrix(rep(0, length(unique(CancerType_cell_number_frac$CancerType)) *
##     ##                                     length(unique(CancerType_cell_number_frac$g_seurat_clusters))),
##     ##                              ncol = length(unique(CancerType_cell_number_frac$CancerType)),
##     ##                              nrow = length(unique(CancerType_cell_number_frac$g_seurat_clusters)))
##     ##         rownames(fracMatrix) <- sort(unique(CancerType_cell_number_frac$g_seurat_clusters))
##     ##         colnames(fracMatrix) <- sort(unique(CancerType_cell_number_frac$CancerType))
##     ##         for(temp_rn in sort(unique(CancerType_cell_number_frac$g_seurat_clusters))){
##     ##             for(temp_cl in sort(unique(CancerType_cell_number_frac$CancerType))){
##     ##                 if(!any(CancerType_cell_number_frac$CancerType == temp_cl & CancerType_cell_number_frac$g_seurat_clusters == temp_rn)){
##     ##                     next
##     ##                 }
##     ##                 fracMatrix[temp_rn, temp_cl] <- CancerType_cell_number_frac$frac[CancerType_cell_number_frac$CancerType == temp_cl & CancerType_cell_number_frac$g_seurat_clusters == temp_rn]
##     ##             }
##     ##         }
##     ##         my.breaks <- seq(min(fracMatrix[!is.na(fracMatrix)]), max(fracMatrix[!is.na(fracMatrix)]), by=0.002)
##     ##         my.colors <- c(
##     ##             colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)/2),
##     ##             colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks)/2))
##     ##         pdf(file.path(tempFigurePath, paste0(temp_TissueType, "_CancerType_cluster_fraction_heatmap.pdf")))
##     ##         print(pheatmap(fracMatrix,
##     ##                        color = my.colors,
##     ##                        breaks = my.breaks,
##     ##                        ## annotation_row = cellType_col,
##     ##                        show_colnames = T,
##     ##                        show_rownames = T,
##     ##                        cluster_cols = F,
##     ##                        cluster_rows = F,
##     ##                        border_color = F))
##     ##         dev.off()
##     ##     }
##     ## }
##     ## g <- ggplot(data = totalObj_ggplotdata) +
##     ##     geom_boxplot(aes(x = reorder(CancerType, frac, FUN=median),
##     ##                      y = frac,
##     ##                      color = CancerType),
##     ##                  size = 0.2,
##     ##                  outlier.shape = NA) +
##     ##     geom_jitter(shape=16,
##     ##                 position=position_jitter(0.2),
##     ##                 size = 0.2,
##     ##                 aes(x = reorder(CancerType, frac, FUN=median),
##     ##                     y = frac,
##     ##                     fill = CancerType,
##     ##                     color = CancerType)) +
##     ##     facet_grid(g_seurat_clusters ~ TissueType,
##     ##                        scales = "free_x", space = "free",
##     ##                        labeller = labeller(TissueType = TissueType.labs)) +
##     ##     xlab("") + ylab("") +
##     ##     ## scale_y_log10() +
##     ##     theme_classic() +
##     ##     theme( text = element_text(size = 8),
##     ##           strip.background = element_blank(),
##     ##           axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
##     ##           legend.position = "none")
##     ## ggsave(file.path(tempFigurePath, paste0("all_boxplot.pdf")), g, width = 210/2, height = 297, units = "mm")
##     ## g <- ggplot(data = totalObj_ggplotdata) +
##     ##     geom_boxplot(aes(x = reorder(CancerType, frac, FUN=median),
##     ##                      y = frac,
##     ##                      color = g_seurat_clusters),
##     ##                  size = 0.2,
##     ##                  outlier.shape = NA) +
##     ##     geom_jitter(shape=16,
##     ##                 position=position_jitter(0.2),
##     ##                 size = 0.2,
##     ##                 aes(x = reorder(CancerType, frac, FUN=median),
##     ##                     y = frac,
##     ##                     fill = g_seurat_clusters,
##     ##                     color = g_seurat_clusters)) +
##     ##     facet_grid(TissueType ~ .,
##     ##                        scales = "free_x", space = "free",
##     ##                        labeller = labeller(TissueType = TissueType.labs)) +
##     ##     xlab("") + ylab("") +
##     ##     ## scale_y_log10() +
##     ##     theme_classic() +
##     ##     theme( text = element_text(size = 8),
##     ##           strip.background = element_blank(),
##     ##           axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
##     ##           legend.position = "none")
##     ## ggsave(file.path(tempFigurePath, paste0("all_boxplot_layout2.pdf")), g, width = 297*4, height = 210/4, units = "mm")
##     ## for(tempTissueType in unique(totalObj_ggplotdata$TissueType)){
##     ##     temp_ggplotdata <- totalObj_ggplotdata %>% filter(TissueType == tempTissueType)
##     ##     g <- ggplot(data = temp_ggplotdata) +
##     ##         geom_boxplot(aes(x = g_seurat_clusters,
##     ##                          y = frac,
##     ##                          color = CancerType),
##     ##                      size = 0.2,
##     ##                      outlier.shape = NA) +
##     ##         geom_jitter(shape=16,
##     ##                     position=position_jitter(0.2),
##     ##                     size = 0.2,
##     ##                     aes(x = g_seurat_clusters,
##     ##                         y = frac,
##     ##                         fill = CancerType,
##     ##                         color = CancerType)) +
##     ##         xlab("") + ylab("") +
##     ##         facet_grid(CancerType ~ .) +
##     ##         ## scale_y_log10() +
##     ##         theme_classic() +
##     ##         theme( text = element_text(size = 8),
##     ##               strip.background = element_blank(),
##     ##               axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
##     ##               legend.position = "bottom")
##     ##     ggsave(file.path(tempFigurePath, paste0(tempTissueType, "_all_boxplot_layout2.pdf")), g, width = 210/2, height = 297, units = "mm")
##     ## }
## }
