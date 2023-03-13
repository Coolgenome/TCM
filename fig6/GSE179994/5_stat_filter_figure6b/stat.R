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

figure_path <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE179994/5_stat_filter_figure6b/")
if (!dir.exists(figure_path)) {
  dir.create(figure_path, recursive = T)
}
setwd(figure_path)

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


CD4_P_obj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE179994/2_extractTcell_proliferative/CD4ProlifSeuratObj_2022-10-31.rds")
CD4_P_meta <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE179994/4_mapping_filter_CD4_Proliferative/proliferative/meta_2022-10-31.tsv")
CD4_P_meta$predicted.celltype <- plyr::mapvalues(CD4_P_meta$predicted.celltype,
                                                 0:6,
                                                 c("P_CD4",
                                                   "P_CD4",
                                                   "P_CD4",
                                                   "P_CD4",
                                                   "P_CD4",
                                                   "P_CD4",
                                                   "P_Treg"))
CD4_P_meta <- CD4_P_meta %>%
  rename(Sample = sample, Patient = patient)

CD8_P_obj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE179994/2_extractTcell_proliferative/CD8ProlifSeuratObj_2022-10-31.rds")
CD8_P_meta <- as_tibble(CD8_P_obj@meta.data) %>%
  mutate(predicted.celltype = "P_CD8") %>%
  mutate(ID = Cells(CD8_P_obj)) %>% 
  rename(Sample = sample, Patient = patient)

CD8_meta <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE179994/MappingResult_filter/CD8/meta_2022-11-03.tsv") %>%
  rename(Sample = sample, Patient = patient)
CD8_meta$predicted.celltype <- plyr::mapvalues(CD8_meta$predicted.celltype,
                                             0:(length(CD8_CellType)-1),
                                             CD8_CellType)
CD8_obj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE179994/2_extractTcell_proliferative/CD8SeuratObj_2022-10-31.rds")
CD8_obj$predicted.celltype <- NA
CD8_obj$predicted.celltype <- CD8_meta$predicted.celltype[match(Cells(CD8_obj), CD8_meta$ID)]
Idents(CD8_obj) <- CD8_obj$predicted.celltype
pdf(file.path(getwd(), "CD8_predicted_celltype.pdf"))
DotPlot(CD8_obj, features = c("CD3D", "CD4", "CD8A", "HSPA1A", "HSPA1B", "MKI67")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        strip.text.y = element_text(angle = 0),
        strip.background = element_rect(colour=NA, fill=NA)) 
dev.off()


CD4_meta <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE179994/MappingResult_filter/CD4/meta_2022-11-03.tsv") %>%
rename(Sample = sample,
       Patient = patient)
CD4_meta$predicted.celltype <- plyr::mapvalues(CD4_meta$predicted.celltype,
                                               0:(length(CD4_CellType)-1),
                                               CD4_CellType)
CD4_obj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE179994/2_extractTcell_proliferative/CD4SeuratObj_2022-10-31.rds")
CD4_obj$predicted.celltype <- NA
CD4_obj$predicted.celltype <- CD4_meta$predicted.celltype[match(Cells(CD4_obj), CD4_meta$ID)]
Idents(CD4_obj) <- CD4_obj$predicted.celltype
pdf(file.path(getwd(), "CD4_predicted_celltype.pdf"))
DotPlot(CD4_obj, features = c("CD3D", "CD4", "CD8A", "HSPA1A", "HSPA1B", "MKI67")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        strip.text.y = element_text(angle = 0),
        strip.background = element_rect(colour=NA, fill=NA)) 
dev.off()

metaValues <- c("P1.pre", "P1.post.1", "P1.post.2", "P1.post.3", "P2.pre",
"P3.pre", "P4.pre", "P5.pre", "P6.pre", "P7.pre", "P8.pre", "P9.pre", "P10.pre",
"P10.post.1", "P11.pre", "P12.pre", "P13.post.1", "P13.post.2", "P13.pre",
"P14.pre", "P15.pre", "P16.pre", "P17.pre", "P18.pre", "P19.pre", "P19.post.1",
"P20.pre", "P21.pre", "P22.pre", "P23.pre", "P24.pre", "P25.pre", "P26.pre",
"P27.pre", "P28.pre", "P29.pre", "P29.post.1", "P30.pre", "P30.post.1",
"P33.pre", "P33.post.1", "P34.pre", "P35.post.1", "P35.pre", "P36.post.1",
"P37.post.1", "P38.post.1")

clinicTValues <- c("P1.pre.1", "P1.post.1", "P1.post.2", "P1.post.3", "P2",
"P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10.pre.1", "P10.post.1", "P11",
"P12", "P13.pre.1", "P13.post.1", "P13.post.3", "P14", "P15", "P16", "P17",
"P18", "P19.pre.1", "P19.post.1", "P20", "P21", "P22", "P23", "P24", "P25",
"P26", "P27", "P28", "P29.pre.1", "P29.post.1", "P30.pre.1", "P30.post.1",
"P33.pre.1", "P33.post.1", "P34", "P35.pre.1", "P35.post.1", "P36.post.1",
"P37.post.1", "P38.post.1" )

CD8_meta$Sample <- plyr::mapvalues(CD8_meta$Sample, metaValues, clinicTValues)
CD4_meta$Sample <- plyr::mapvalues(CD4_meta$Sample, metaValues, clinicTValues)

clinicT <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/data/GSE179994/ClinicData.txt") %>%
  rename(TumorType = `Tumor Type`,
         TreatmentHx = `Treatment Hx`,
         BiospySite = `Biospy Site`,
         SampleName = `Sample Name`)

CD8_meta$TumorType <- clinicT$TumorType[match(CD8_meta$Sample, clinicT$SampleName)]
CD8_meta$Treatment <- clinicT$Treatment[match(CD8_meta$Sample, clinicT$SampleName)]
CD8_meta$TreatmentHx <- clinicT$TreatmentHx[match(CD8_meta$Sample, clinicT$SampleName)]
CD8_meta$BiospySite <- clinicT$BiospySite[match(CD8_meta$Sample, clinicT$SampleName)]
CD8_meta$Response <- clinicT$Response[match(CD8_meta$Sample, clinicT$SampleName)]
CD8_meta$TreatmentHx[CD8_meta$TreatmentHx == "On treatment"] <- "On-treatment"
CD8_meta$Response <- plyr::mapvalues(CD8_meta$Response, c("Yes", "No"), c("R", "NR"))

CD4_meta$TumorType <- clinicT$TumorType[match(CD4_meta$Sample, clinicT$SampleName)]
CD4_meta$Treatment <- clinicT$Treatment[match(CD4_meta$Sample, clinicT$SampleName)]
CD4_meta$TreatmentHx <- clinicT$TreatmentHx[match(CD4_meta$Sample, clinicT$SampleName)]
CD4_meta$BiospySite <- clinicT$BiospySite[match(CD4_meta$Sample, clinicT$SampleName)]
CD4_meta$Response <- clinicT$Response[match(CD4_meta$Sample, clinicT$SampleName)]
CD4_meta$TreatmentHx[CD4_meta$TreatmentHx == "On treatment"] <- "On-treatment"
CD4_meta$Response <- plyr::mapvalues(CD4_meta$Response, c("Yes", "No"), c("R", "NR"))

CD8_P_meta$TumorType <- clinicT$TumorType[match(CD8_P_meta$Sample, clinicT$SampleName)]
CD8_P_meta$Treatment <- clinicT$Treatment[match(CD8_P_meta$Sample, clinicT$SampleName)]
CD8_P_meta$TreatmentHx <- clinicT$TreatmentHx[match(CD8_P_meta$Sample, clinicT$SampleName)]
CD8_P_meta$BiospySite <- clinicT$BiospySite[match(CD8_P_meta$Sample, clinicT$SampleName)]
CD8_P_meta$Response <- clinicT$Response[match(CD8_P_meta$Sample, clinicT$SampleName)]
CD8_P_meta$TreatmentHx[CD8_P_meta$TreatmentHx == "On treatment"] <- "On-treatment"
CD8_P_meta$Response <- plyr::mapvalues(CD8_P_meta$Response, c("Yes", "No"), c("R", "NR"))

CD4_P_meta$TumorType <- clinicT$TumorType[match(CD4_P_meta$Sample, clinicT$SampleName)]
CD4_P_meta$Treatment <- clinicT$Treatment[match(CD4_P_meta$Sample, clinicT$SampleName)]
CD4_P_meta$TreatmentHx <- clinicT$TreatmentHx[match(CD4_P_meta$Sample, clinicT$SampleName)]
CD4_P_meta$BiospySite <- clinicT$BiospySite[match(CD4_P_meta$Sample, clinicT$SampleName)]
CD4_P_meta$Response <- clinicT$Response[match(CD4_P_meta$Sample, clinicT$SampleName)]
CD4_P_meta$TreatmentHx[CD4_P_meta$TreatmentHx == "On treatment"] <- "On-treatment"
CD4_P_meta$Response <- plyr::mapvalues(CD4_P_meta$Response, c("Yes", "No"), c("R", "NR"))


cns_totalMeta <- c("ID", "Patient", "Sample", "predicted.celltype", "TumorType", "Treatment", "TreatmentHx", "BiospySite", "Response")

totalMeta <- bind_rows(CD8_meta[, cns_totalMeta], CD4_meta[, cns_totalMeta])
totalMeta <- bind_rows(totalMeta, CD4_P_meta[, cns_totalMeta])
totalMeta <- bind_rows(totalMeta, CD8_P_meta[, cns_totalMeta])

# Group A #####################################################################

totalMeta_A <- totalMeta %>% mutate(Site=BiospySite)
totalMeta_A$Site[totalMeta_A$Site %in% c("Left lung tumour", "Right lung tumour")] <- "Left/Right lung tumor"
totalMeta_A <- totalMeta_A %>%
    filter(Treatment == "Treatment-naive") %>%
    filter(TumorType %in% c("LUAD", "LUSC")) %>%
    filter(Site %in% c("Left/Right lung tumor", "LN metastasis"))
totalMeta_A$Site<- plyr::mapvalues(totalMeta_A$Site, c("Left/Right lung tumor", "LN metastasis"), c("Primary", "LN-Met"))
totalMeta_A <- totalMeta_A %>%
  mutate(group = paste0(Treatment, "-", Site))


CD8_obj$group <- totalMeta_A$group[match(Cells(CD8_obj), totalMeta_A$ID)]
Idents(CD8_obj) <- CD8_obj$group
temp_obj <- subset(CD8_obj, cells = Cells(CD8_obj)[!is.na(CD8_obj$group)])
pdf(file.path(getwd(), "CD8_A_group_dotplot.pdf"))
DotPlot(temp_obj, features = c("CD3D", "CD4", "CD8A", "CD8B", "HSPA1A", "HSPA1B")) +
 theme_classic() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        strip.text.y = element_text(angle = 0),
        strip.background = element_rect(colour=NA, fill=NA)) 
dev.off()

CD4_obj$group <- totalMeta_A$group[match(Cells(CD4_obj), totalMeta_A$ID)]
Idents(CD4_obj) <- CD4_obj$group
temp_obj <- subset(CD4_obj, cells = Cells(CD4_obj)[!is.na(CD4_obj$group)])
pdf(file.path(getwd(), "CD4_A_group_dotplot.pdf"))
DotPlot(temp_obj, features = c("CD3D", "CD4", "CD8A", "CD8B", "HSPA1A", "HSPA1B")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        strip.text.y = element_text(angle = 0),
        strip.background = element_rect(colour=NA, fill=NA)) 
dev.off()



###############################################################################
#                               get sample info                               #
###############################################################################

sampleMD <- totalMeta_A %>%
  group_by(Site, Sample) %>%
  dplyr::count() %>%
  ungroup() %>%
  group_by(Site) %>%
  dplyr::count()
write_tsv(sampleMD, file.path(figure_path,
                              paste0('NSCLC_Liu_treatmentNaive_sampleNum',
                                     "_", Sys.Date(), '.tsv')))

groupCellNum <- totalMeta_A %>%
  group_by(Site) %>%
  dplyr::count()
write_tsv(groupCellNum, file.path(figure_path,
                                  paste0('NSCLC_Liu_treatmentNaive_cellNum',
                                         "_", Sys.Date(), '.tsv')))

###############################################################################
#                            finish get sample info                           #
###############################################################################


# Group D #####################################################################

totalMeta_D <- totalMeta %>% mutate(Site=BiospySite)
totalMeta_D$Site[totalMeta_D$Site %in% c("Left lung tumour", "Right lung tumour")] <- "Left/Right lung tumor"
totalMeta_D <- totalMeta_D %>%
    filter(Site %in% c("Left/Right lung tumor")) %>%
    filter(TreatmentHx %in% c("Pre-treatment", "On-treatment")) %>%
    filter(Response %in% c("NR", "R"))
totalMeta_D$Site<- plyr::mapvalues(totalMeta_D$Site, c("Left/Right lung tumor"), c("Primary"))
totalMeta_D$TreatmentHx<- plyr::mapvalues(totalMeta_D$TreatmentHx, c("Pre-treatment", "On-treatment"), c("Pre", "On"))
totalMeta_D <- totalMeta_D %>%
    mutate(group = paste0("Primary-", TreatmentHx, "-", Response))

CD8_obj$group <- totalMeta_D$group[match(Cells(CD8_obj), totalMeta_D$ID)]
Idents(CD8_obj) <- CD8_obj$group
temp_obj <- subset(CD8_obj, cells = Cells(CD8_obj)[!is.na(CD8_obj$group)])
pdf(file.path(getwd(), "CD8_D_group_dotplot.pdf"))
DotPlot(temp_obj, features = c("CD3D", "CD4", "CD8A", "CD8B", "HSPA1A", "HSPA1B")) +
 theme_classic() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        strip.text.y = element_text(angle = 0),
        strip.background = element_rect(colour=NA, fill=NA)) 
dev.off()

CD4_obj$group <- totalMeta_D$group[match(Cells(CD4_obj), totalMeta_D$ID)]
Idents(CD4_obj) <- CD4_obj$group
temp_obj <- subset(CD4_obj, cells = Cells(CD4_obj)[!is.na(CD4_obj$group)])
pdf(file.path(getwd(), "CD4_D_group_dotplot.pdf"))
DotPlot(temp_obj, features = c("CD3D", "CD4", "CD8A", "CD8B", "HSPA1A", "HSPA1B")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        strip.text.y = element_text(angle = 0),
        strip.background = element_rect(colour=NA, fill=NA)) 
dev.off()

###############################################################################
#                               get sample info                               #
###############################################################################

sampleMD <- totalMeta_D %>%
  group_by(Response,
           Site,
           TreatmentHx,
           Sample) %>%
  dplyr::count() %>%
  ungroup() %>%
  group_by(Response,
           Site,
           TreatmentHx) %>%
  dplyr::count()
write_tsv(sampleMD, file.path(figure_path,
                              paste0('NSCLC_Liu_primary_sampleNum',
                                     "_", Sys.Date(), '.tsv')))

groupCellNum <- totalMeta_D %>%
  group_by(Response,
           Site,
           TreatmentHx) %>%
  dplyr::count()
write_tsv(groupCellNum, file.path(figure_path,
                                  paste0('NSCLC_Liu_primary_cellNum',
                                         "_", Sys.Date(), '.tsv')))

###############################################################################
#                            finish get sample info                           #
###############################################################################



# Group E #####################################################################

totalMeta_E <- totalMeta %>% mutate(Site=BiospySite)
totalMeta_E$Site[totalMeta_E$Site %in% c("Left lung tumour", "Right lung tumour")] <- "Left/Right lung tumor"
totalMeta_E <- totalMeta_E %>%
    filter(Site %in% c("LN metastasis")) %>%
    filter(TreatmentHx %in% c("Pre-treatment", "On-treatment")) %>%
    filter(Response %in% c("NR", "R"))
totalMeta_E$Site<- plyr::mapvalues(totalMeta_E$Site, c("LN metastasis"), c("LN-Met"))
totalMeta_E$TreatmentHx<- plyr::mapvalues(totalMeta_E$TreatmentHx, c("Pre-treatment", "On-treatment"), c("Pre", "On"))
totalMeta_E <- totalMeta_E %>%
    mutate(group = paste0("Metastatic-", TreatmentHx, "-", Response))


CD8_obj$group <- totalMeta_E$group[match(Cells(CD8_obj), totalMeta_E$ID)]
Idents(CD8_obj) <- CD8_obj$group
temp_obj <- subset(CD8_obj, cells = Cells(CD8_obj)[!is.na(CD8_obj$group)])
pdf(file.path(getwd(), "CD8_E_group_dotplot.pdf"))
DotPlot(temp_obj, features = c("CD3D", "CD4", "CD8A", "CD8B", "HSPA1A", "HSPA1B")) +
 theme_classic() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        strip.text.y = element_text(angle = 0),
        strip.background = element_rect(colour=NA, fill=NA)) 
dev.off()

CD4_obj$group <- totalMeta_E$group[match(Cells(CD4_obj), totalMeta_E$ID)]
Idents(CD4_obj) <- CD4_obj$group
temp_obj <- subset(CD4_obj, cells = Cells(CD4_obj)[!is.na(CD4_obj$group)])
pdf(file.path(getwd(), "CD4_E_group_dotplot.pdf"))
DotPlot(temp_obj, features = c("CD3D", "CD4", "CD8A", "CD8B", "HSPA1A", "HSPA1B")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        strip.text.y = element_text(angle = 0),
        strip.background = element_rect(colour=NA, fill=NA)) 
dev.off()


###############################################################################
#                               get sample info                               #
###############################################################################
sampleMD <- totalMeta_E %>%
  group_by(Response,
           Site,
           TreatmentHx,
           Sample) %>%
  dplyr::count() %>%
  ungroup() %>%
  group_by(Response,
           Site,
           TreatmentHx) %>%
  dplyr::count()
write_tsv(sampleMD, file.path(figure_path,
                              paste0('NSCLC_Liu_meta_sampleNum',
                                     "_", Sys.Date(), '.tsv')))

groupCellNum <- totalMeta_E %>%
  group_by(Response,
           Site,
           TreatmentHx) %>%
  dplyr::count()
write_tsv(groupCellNum, file.path(figure_path,
                                  paste0('NSCLC_Liu_meta_cellNum',
                                         "_", Sys.Date(), '.tsv')))

###############################################################################
#                            finish get sample info                           #
###############################################################################

totalMatRoe <- c()
totalMatCell <- c()

totalMeta_ADE <- bind_rows(totalMeta_A, totalMeta_D)
totalMeta_ADE <- bind_rows(totalMeta_ADE, totalMeta_E)


CD8_obj$group <- totalMeta_ADE$group[match(Cells(CD8_obj), totalMeta_ADE$ID)]
Idents(CD8_obj) <- CD8_obj$group
temp_obj <- subset(CD8_obj, cells = Cells(CD8_obj)[!is.na(CD8_obj$group)])
pdf(file.path(getwd(), "CD8_group_dotplot.pdf"))
DotPlot(temp_obj, features = c("CD3D", "CD4", "CD8A", "CD8B", "HSPA1A", "HSPA1B")) +
 theme_classic() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        strip.text.y = element_text(angle = 0),
        strip.background = element_rect(colour=NA, fill=NA)) 
dev.off()

CD4_obj$group <- totalMeta_ADE$group[match(Cells(CD4_obj), totalMeta_ADE$ID)]
Idents(CD4_obj) <- CD4_obj$group
temp_obj <- subset(CD4_obj, cells = Cells(CD4_obj)[!is.na(CD4_obj$group)])
pdf(file.path(getwd(), "CD4_group_dotplot.pdf"))
DotPlot(temp_obj, features = c("CD3D", "CD4", "CD8A", "CD8B", "HSPA1A", "HSPA1B")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        strip.text.y = element_text(angle = 0),
        strip.background = element_rect(colour=NA, fill=NA)) 
dev.off()

CD4_P_obj$group <- totalMeta_ADE$group[match(Cells(CD4_P_obj), totalMeta_ADE$ID)]
Idents(CD4_P_obj) <- CD4_P_obj$group
temp_obj <- subset(CD4_P_obj, cells = Cells(CD4_P_obj)[!is.na(CD4_P_obj$group)])
pdf(file.path(getwd(), "CD4_P_ADE_group_dotplot.pdf"))
DotPlot(temp_obj, features = c("CD3D", "CD4", "CD8A", "CD8B", "HSPA1A", "HSPA1B")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
    strip.text.y = element_text(angle = 0),
    strip.background = element_rect(colour = NA, fill = NA)
  )
dev.off()

CD8_P_obj$group <- totalMeta_ADE$group[match(Cells(CD8_P_obj), totalMeta_ADE$ID)]
Idents(CD8_P_obj) <- CD8_P_obj$group
temp_obj <- subset(CD8_P_obj, cells = Cells(CD8_P_obj)[!is.na(CD8_P_obj$group)])
pdf(file.path(getwd(), "CD8_P_ADE_group_dotplot.pdf"))
DotPlot(temp_obj, features = c("CD3D", "CD4", "CD8A", "CD8B", "HSPA1A", "HSPA1B")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
    strip.text.y = element_text(angle = 0),
    strip.background = element_rect(colour = NA, fill = NA)
  )
dev.off()




totalMeta_L_ADE <- list(A = totalMeta_A, D = totalMeta_D, E = totalMeta_E)

N <- dim(totalMeta)[1]
matRoe <- matrix(rep(0, length(unique(totalMeta_ADE$predicted.celltype)) *
                        length(unique(totalMeta_ADE$group))),
                 nrow = length(unique(totalMeta_ADE$predicted.celltype)),
                 ncol = length(unique(totalMeta_ADE$group)))

matCell <- matrix(rep(0, length(unique(totalMeta_ADE$predicted.celltype)) *
                        length(unique(totalMeta_ADE$group))),
                 nrow = length(unique(totalMeta_ADE$predicted.celltype)),
                 ncol = length(unique(totalMeta_ADE$group)))

rownames(matRoe) <- unique(totalMeta_ADE$predicted.celltype)
colnames(matRoe) <- unique(totalMeta_ADE$group)
rownames(matCell) <- unique(totalMeta_ADE$predicted.celltype)
colnames(matCell) <- unique(totalMeta_ADE$group)

ClusterT <- totalMeta %>%
  group_by(predicted.celltype) %>%
  dplyr::count()

for(tempGroup in unique(totalMeta_ADE$group)) {
  targetTcellMD <- totalMeta_ADE %>% filter(group == tempGroup)
  M <- dim(targetTcellMD)[1]
  for(gsci in unique(totalMeta_ADE$predicted.celltype)) {
    k <- ClusterT$n[ClusterT$predicted.celltype == gsci]
    if(gsci %in% targetTcellMD$predicted.celltype) {
      n <- dim(targetTcellMD %>% filter(predicted.celltype == gsci))[1]
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

cns <- c("Treatment-naive-Primary", "Treatment-naive-LN-Met", "Primary-Pre-R",
         "Metastatic-Pre-R", "Metastatic-On-R", "Metastatic-On-NR")
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
widthi = ncol(tMR)/12 + 4
heighti = nrow(tMR)/7 + 3
pdf(file.path(figure_path, paste0("figure6e_Roe.pdf")), width = widthi, height = heighti)
print(pheatmap(tMR,
               color = my.colors,
               breaks = my.breaks,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figure_path, paste0("figure6e_Cell.pdf")), width = widthi, height = heighti)
print(pheatmap(tMC,
               show_colnames = T,
               show_rownames = T,
               display_numbers = T,
               number_format = "%.0f",
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()


rns_gt <- totalMeta %>%
  group_by(predicted.celltype) %>%
  count %>%
  filter(n >= 100) %>%
  pull(predicted.celltype)
tMR <- matRoe[setdiff(rns_gt, rns),]
tMC <- matCell[setdiff(rns_gt, rns),]
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



