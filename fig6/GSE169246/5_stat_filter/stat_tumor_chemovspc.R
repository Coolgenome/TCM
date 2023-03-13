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

threshold <- 5
midRed <- 1.2

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

CD4_CellType <- c("CD4_c0_Tcm",
                  "CD4_c1_Treg",
                  "CD4_c2_Tn",
                  "CD4_c3_Tfh",
                  "CD4_c4_Tstr",
                  "CD4_c5_Tctl",
                  "CD4_c6_Tn_FHIT",
                  "CD4_c7_Tn_TCEA3",
                  "CD4_c8_Th17",
                  "CD4_c9_Tn_TCF7_SLC40A1",
                  "CD4_c10_Tn_LEF1_ANKRD55",
                  "CD4_c11_Tisg")

Innate_CellType <- c("NKT_c0_FCGR3A_GZMB",
                     "MAIT-like_c1",
                     "MAIT_c2_KLRB1",
                     "Tgd_c3_CX3CR1",
                     "NKT_c4_KIR_TIGIT_CXCR6")

Proliferative_CellType <- c("P_c0_CD8_KLRD1_GZMB_CCL4/5",
                            "P_c1_CD4_CD40LG_IL7R",
                            "P_c2_DNT",
                            "P_c3_DNT_GZMK+",
                            "P_c4_CD8_C1QBP_MT1G/X/E",
                            "P_c5_CD8_CXCL13_CCL3/4/5_PD-1",
                            "P_c6_Treg",
                            "P_c7_CD8_GZMK_NKG7_LAG3_PRF1")


TotalAntiPDL1ChemoPatients <- c("P019",
                                "P010",
                                "P012",
                                "P007",
                                "P017",
                                "P001",
                                "P002",
                                "P014",
                                "P004",
                                "P005",
                                "P016")

TotalChemoPatients <- c("P022",
                        "P011",
                        "P020",
                        "P008",
                        "P013",
                        "P025",
                        "P018",
                        "P023",
                        "P024",
                        "P003",
                        "P028")

AllResponsePatients <- c("P019",
                         "P010",
                         "P012",
                         "P007",
                         "P022",
                         "P011",
                         "P020",
                         "P008",
                         "P013")

figure_path <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/5_stat_filter/")
if (!dir.exists(figure_path)) {
  dir.create(figure_path, recursive = T)
}
setwd(figure_path)


CD8_meta <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/MappingResult_filter/CD8/meta_2022-11-03.tsv")
CD8_meta$Sample <- stringr::str_extract(CD8_meta$ID, "(?<=^.{10,20}\\.).+")
CD8_meta$Patient <- stringr::str_extract(CD8_meta$Sample, "P\\d+")
CD8_meta$Tissue <- stringr::str_extract(CD8_meta$Sample, "\\w$")
CD8_meta$TumorTreatment <- stringr::str_extract(CD8_meta$Sample, "^[a-zA-Z]+")
CD8_meta$isResponse <- "NR"
CD8_meta$isResponse[CD8_meta$Patient %in% AllResponsePatients] <- "R"
CD8_meta$TreatmentType <- "PDL1+Chemo"
CD8_meta$TreatmentType[CD8_meta$Patient %in% TotalChemoPatients] <- "Chemo"
CD8_meta$group <- paste0(CD8_meta$Tissue, "_", CD8_meta$TreatmentType, "_", CD8_meta$TumorTreatment, "-", CD8_meta$isResponse)
CD8_meta$predicted.celltype <- plyr::mapvalues(CD8_meta$predicted.celltype,
                                               0:(length(CD8_CellType) - 1),
                                               CD8_CellType)

CD8_obj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/TCells/forMapping/CD8.rds")
CD8_obj$Sample <- CD8_meta$Sample[match(Cells(CD8_obj), CD8_meta$ID)]
CD8_obj$Patient <- CD8_meta$Patient[match(Cells(CD8_obj), CD8_meta$ID)]
CD8_obj$Tissue <- CD8_meta$Tissue[match(Cells(CD8_obj), CD8_meta$ID)]
CD8_obj$TumorTreatment <- CD8_meta$TumorTreatment[match(Cells(CD8_obj), CD8_meta$ID)]
CD8_obj$isResponse <- CD8_meta$isResponse[match(Cells(CD8_obj), CD8_meta$ID)]
CD8_obj$TreatmentType <- CD8_meta$TreatmentType[match(Cells(CD8_obj), CD8_meta$ID)]
CD8_obj$group <- paste0(CD8_obj$Tissue, "_", CD8_obj$TreatmentType, "_", CD8_obj$TumorTreatment, "-", CD8_obj$isResponse)

CD4_meta <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/MappingResult_filter/CD4/meta_2022-11-04.tsv")
CD4_meta$Sample <- stringr::str_extract(CD4_meta$ID, "(?<=^.{10,20}\\.).+")
CD4_meta$Patient <- stringr::str_extract(CD4_meta$Sample, "P\\d+")
CD4_meta$Tissue <- stringr::str_extract(CD4_meta$Sample, "\\w$")
CD4_meta$TumorTreatment <- stringr::str_extract(CD4_meta$Sample, "^[a-zA-Z]+")
CD4_meta$isResponse <- "NR"
CD4_meta$isResponse[CD4_meta$Patient %in% AllResponsePatients] <- "R"
CD4_meta$TreatmentType <- "PDL1+Chemo"
CD4_meta$TreatmentType[CD4_meta$Patient %in% TotalChemoPatients] <- "Chemo"
CD4_meta$group <- paste0(CD4_meta$Tissue, "_", CD4_meta$TreatmentType, "_", CD4_meta$TumorTreatment, "-", CD4_meta$isResponse)
CD4_meta$predicted.celltype <- plyr::mapvalues(CD4_meta$predicted.celltype,
                                               0:(length(CD4_CellType) - 1),
                                               CD4_CellType)

CD4_obj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/TCells/forMapping/CD4.rds")
CD4_obj$Sample <- CD4_meta$Sample[match(Cells(CD4_obj), CD4_meta$ID)]
CD4_obj$Patient <- CD4_meta$Patient[match(Cells(CD4_obj), CD4_meta$ID)]
CD4_obj$Tissue <- CD4_meta$Tissue[match(Cells(CD4_obj), CD4_meta$ID)]
CD4_obj$TumorTreatment <- CD4_meta$TumorTreatment[match(Cells(CD4_obj), CD4_meta$ID)]
CD4_obj$isResponse <- CD4_meta$isResponse[match(Cells(CD4_obj), CD4_meta$ID)]
CD4_obj$TreatmentType <- CD4_meta$TreatmentType[match(Cells(CD4_obj), CD4_meta$ID)]
CD4_obj$group <- paste0(CD4_obj$Tissue, "_", CD4_obj$TreatmentType, "_", CD4_obj$TumorTreatment, "-", CD4_obj$isResponse)

P_meta <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/MappingResult_filter/P/meta_2022-11-03.tsv")
P_meta$Sample <- stringr::str_extract(P_meta$ID, "(?<=^.{10,20}\\.).+")
P_meta$Patient <- stringr::str_extract(P_meta$Sample, "P\\d+")
P_meta$Tissue <- stringr::str_extract(P_meta$Sample, "\\w$")
P_meta$TumorTreatment <- stringr::str_extract(P_meta$Sample, "^[a-zA-Z]+")
P_meta$isResponse <- "NR"
P_meta$isResponse[P_meta$Patient %in% AllResponsePatients] <- "R"
P_meta$TreatmentType <- "PDL1+Chemo"
P_meta$TreatmentType[P_meta$Patient %in% TotalChemoPatients] <- "Chemo"
P_meta$group <- paste0(P_meta$Tissue, "_", P_meta$TreatmentType, "_", P_meta$TumorTreatment, "-", P_meta$isResponse)
P_meta$predicted.celltype <- plyr::mapvalues(P_meta$predicted.celltype,
                                             0:(length(Proliferative_CellType) - 1),
                                             c("P_CD8",
                                               "P_CD4",
                                               "P_Else",
                                               "P_Else",
                                               "P_CD8",
                                               "P_CD8",
                                               "P_Treg",
                                               "P_CD8"))

P_obj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/TCells/forMapping/P.rds")
P_obj@meta.data$Sample <- stringr::str_extract(Cells(P_obj), "(?<=^.{10,20}\\.).+")
P_obj@meta.data$Patient <- stringr::str_extract(P_obj@meta.data$Sample, "P\\d+")
P_obj@meta.data$Tissue <- stringr::str_extract(P_obj@meta.data$Sample, "\\w$")
P_obj@meta.data$TumorTreatment <- stringr::str_extract(P_obj@meta.data$Sample, "^[a-zA-Z]+")
P_obj@meta.data$isResponse <- "NR"
P_obj@meta.data$isResponse[P_obj@meta.data$Patient %in% AllResponsePatients] <- "R"
P_obj@meta.data$TreatmentType <- "PDL1+Chemo"
P_obj@meta.data$TreatmentType[P_obj@meta.data$Patient %in% TotalChemoPatients] <- "Chemo"
P_obj$group <- paste0(P_obj$Tissue, "_", P_obj$TreatmentType, "_", P_obj$TumorTreatment, "-", P_obj$isResponse)

all_obj_L <- list(CD4 = CD4_obj, CD8 = CD8_obj, P = P_obj)
for (aoli in seq_along(all_obj_L)) {
  obj_name <- names(all_obj_L[aoli])
  obj <- all_obj_L[[aoli]]
  
  Idents(obj) <- obj$group
  pdf(file.path(getwd(), paste0(obj_name, "_bubbleplot.pdf")))
  g <- DotPlot(obj, features = c("CD3D", "CD3E", "CD4", "CD8A", "CD8B", "HSPA1A", "HSPA1B", "MKI67")) +
   theme_classic() +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
          strip.text.y = element_text(angle = 0),
          strip.background = element_rect(colour=NA, fill=NA)) 
  print(g)
  dev.off()
}

allT <- bind_rows(CD4_meta, CD8_meta)
allT <- bind_rows(allT, P_meta)

N <- dim(allT)[1]
ClusterT <- allT %>%
  group_by(predicted.celltype) %>%
  count()
allTSampleCN <- allT %>%
    group_by(Patient, Sample) %>%
    count

allMD <- allT %>% 
  filter(TumorTreatment != "Prog")

matRoe <- matrix(rep(0, length(unique(allMD$predicted.celltype)) * length(unique(allMD$group))),
                 nrow = length(unique(allMD$predicted.celltype)),
                 ncol = length(unique(allMD$group)))
matCell <- matrix(rep(0, length(unique(allMD$predicted.celltype)) * length(unique(allMD$group))),
                  nrow = length(unique(allMD$predicted.celltype)),
                  ncol = length(unique(allMD$group)))
rownames(matRoe) <- sort(unique(allMD$predicted.celltype))
colnames(matRoe) <- unique(allMD$group)
rownames(matCell) <- sort(unique(allMD$predicted.celltype))
colnames(matCell) <- unique(allMD$group)

for(tempGroup in unique(allMD$group)){
  targetTcellMD <- allMD %>% filter(group == tempGroup)
  M <- dim(targetTcellMD)[1]
  for(gsci in unique(allMD$predicted.celltype)){
    k <- ClusterT$n[ClusterT$predicted.celltype == gsci]
    if(gsci %in% targetTcellMD$predicted.celltype){
      n <- dim(targetTcellMD %>% filter(predicted.celltype == gsci))[1]
    } else{
      n <- 0
    }
    tempRowName <- as.character(gsci)
    tempColName <- tempGroup
    matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
    matCell[tempRowName, tempColName] <- n
  }
}

## cns <- c("Treatment-naive-Primary", "Treatment-naive-LN-Met", "Primary-Pre-R",
##          "Metastatic-Pre-R", "Metastatic-On-R", "Metastatic-On-NR")
pre_rns <- c(
  "CD4_c4_Tstr",
  "CD8_c4_Tstr",
  "CD8_c1_Tex",
  "CD4_c3_TFH",
  "CD4_c1_Treg",
  "CD8_c5_Tisg",
  "CD8_c2_Teff",
  "CD8_c3_Tn",
  "CD4_c2_Tn",
  "CD8_c0_t-Teff",
  "CD4_c0_Tcm")

rns <- c(pre_rns, setdiff(rownames(matRoe), pre_rns))
rns <- intersect(rns, rownames(matRoe))

tMR <- matRoe[rns, ]
tMC <- matCell[rns,]
tMR[tMR > threshold] = threshold
my.breaks.low <- seq(0, 1, by=0.01)
my.colors.low <- colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks.low))
my.breaks.mid <- seq(1.01, midRed, by=0.01)
my.colors.mid <- colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks.mid))
my.breaks.up <- seq(midRed + 0.01, threshold, by = (threshold - midRed)/99)
my.colors.up <- colorRampPalette(colors = c("#FD9AA0", "#550000"))(length(my.breaks.up))
my.breaks <- c(my.breaks.low, my.breaks.mid, my.breaks.up)
my.colors <- c(my.colors.low, my.colors.mid, my.colors.up)
widthi = ncol(tMR)/12 + 8
heighti = nrow(tMR)/7 + 8
pdf(file.path(figure_path, paste0("figure6_GSE169246_Roe.pdf")), width = widthi, height = heighti)
print(pheatmap(tMR,
               color = my.colors,
               breaks = my.breaks,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figure_path, paste0("figure6_GSE169246_Cell.pdf")), width = widthi, height = heighti)
print(pheatmap(tMC,
               show_colnames = T,
               show_rownames = T,
               display_numbers = T,
               number_format = "%.0f",
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()


cns <- c(
  "b_PDL1+Chemo_Pre-R",
  "b_PDL1+Chemo_Pre-NR",
  "b_PDL1+Chemo_Post-R",
  "b_PDL1+Chemo_Post-NR",
  "b_Chemo_Pre-R",
  "b_Chemo_Pre-NR",
  "b_Chemo_Post-R",
  "b_Chemo_Post-NR",
  "t_PDL1+Chemo_Pre-R",
  "t_PDL1+Chemo_Pre-NR",
  "t_PDL1+Chemo_Post-R",
  "t_PDL1+Chemo_Post-NR",
  "t_Chemo_Pre-R",
  "t_Chemo_Pre-NR",
  "t_Chemo_Post-R",
  "t_Chemo_Post-NR")

rns_lt <- allMD %>%
  group_by(predicted.celltype) %>%
  count %>%
  filter(n < 100) %>%
  pull(predicted.celltype)
tMR <- matRoe[setdiff(rns, rns_lt), cns]
tMC <- matCell[setdiff(rns, rns_lt), cns]
tMR[tMR > threshold] = threshold
my.breaks.low <- seq(0, 1, by=0.01)
my.colors.low <- colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks.low))
my.breaks.mid <- seq(1.01, midRed, by=0.01)
my.colors.mid <- colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks.mid))
my.breaks.up <- seq(midRed + 0.01, threshold, by = (threshold - midRed)/99)
my.colors.up <- colorRampPalette(colors = c("#FD9AA0", "#550000"))(length(my.breaks.up))
my.breaks <- c(my.breaks.low, my.breaks.mid, my.breaks.up)
my.colors <- c(my.colors.low, my.colors.mid, my.colors.up)
widthi = ncol(tMR)/12 + 6
heighti = nrow(tMR)/7 + 6
pdf(file.path(figure_path, paste0("figure6e_gt100_Roe.pdf")), width = widthi, height = heighti)
print(pheatmap(tMR,
               color = my.colors,
               breaks = my.breaks,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figure_path, paste0("figure6e_gt100_Cell.pdf")), width = widthi, height = heighti)
print(pheatmap(tMC,
               show_colnames = T,
               show_rownames = T,
               display_numbers = T,
               number_format = "%.0f",
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()

