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

figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE179994/5_stat_filter_treatment-naive_violin/outs")
if(!dir.exists(figurePath)){
    dir.create(figurePath, recursive = T)
}
setwd(figurePath)

clinicT <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/data/GSE179994/ClinicData.txt") %>%
    rename(TumorType = `Tumor Type`,
           TreatmentHx = `Treatment Hx`,
           BiospySite = `Biospy Site`,
           SampleName = `Sample Name`)

CD8_meta <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE179994/MappingResult_filter/CD8/meta_2022-05-25.tsv") %>%
    rename(Sample = sample,
           Patient = patient)

CD4_meta <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE179994/MappingResult_filter/CD4/meta_2022-05-25.tsv") %>%
rename(Sample = sample,
       Patient = patient)

CD4Obj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE179994/TCells/ForMapping/CD4.rds")
CD8Obj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE179994/TCells/ForMapping/CD8.rds")

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

## # Group A #####################################################################
## CD8MetaA <- CD8_meta %>% mutate(Site=BiospySite)
## CD8MetaA$Site[CD8MetaA$Site %in% c("Left lung tumour", "Right lung tumour")] <- "Left/Right lung tumor"
## CD8MetaA <- CD8MetaA %>%
##     filter(Treatment == "Treatment-naive") %>%
##     filter(TumorType %in% c("LUAD", "LUSC")) %>%
##     filter(Site %in% c("Left/Right lung tumor", "LN metastasis", "Liver metastasis"))
## CD4MetaA <- CD4_meta %>% mutate(Site=BiospySite)
## CD4MetaA$Site[CD4MetaA$Site %in% c("Left lung tumour", "Right lung tumour")] <- "Left/Right lung tumor"
## CD4MetaA <- CD4MetaA %>%
##     filter(Treatment == "Treatment-naive") %>%
##     filter(TumorType %in% c("LUAD", "LUSC")) %>%
##     filter(Site %in% c("Left/Right lung tumor", "LN metastasis", "Liver metastasis"))
## metaL_A <- list(CD8 = CD8MetaA, CD4 = CD4MetaA)
## # Group B #####################################################################
## CD8MetaB <- CD8_meta %>% mutate(Site=BiospySite)
## CD8MetaB$Site[CD8MetaB$Site %in% c("Left lung tumour", "Right lung tumour")] <- "Left/Right lung tumor"
## CD8MetaB <- CD8MetaB %>%
##     filter(Treatment == "Treatment-naive") %>%
##     filter(TumorType %in% c("LUAD")) %>%
##     filter(Site %in% c("Left/Right lung tumor", "LN metastasis", "Liver metastasis"))
## CD4MetaB <- CD4_meta %>% mutate(Site=BiospySite)
## CD4MetaB$Site[CD4MetaB$Site %in% c("Left lung tumour", "Right lung tumour")] <- "Left/Right lung tumor"
## CD4MetaB <- CD4MetaB %>%
##     filter(Treatment == "Treatment-naive") %>%
##     filter(TumorType %in% c("LUAD")) %>%
##     filter(Site %in% c("Left/Right lung tumor", "LN metastasis", "Liver metastasis"))
## metaL_B <- list(CD8 = CD8MetaB, CD4 = CD4MetaB)
## # Group C #####################################################################
## CD8MetaC <- CD8_meta %>% mutate(Site=BiospySite)
## CD8MetaC$Site[CD8MetaC$Site %in% c("Left lung tumour", "Right lung tumour")] <- "Left/Right lung tumor"
## CD8MetaC <- CD8MetaC %>%
##     filter(Treatment == "Treatment-naive") %>%
##     filter(TumorType %in% c("LUSC")) %>%
##     filter(Site %in% c("Left/Right lung tumor", "LN metastasis", "Liver metastasis"))
## CD4MetaC <- CD4_meta %>% mutate(Site=BiospySite)
## CD4MetaC$Site[CD4MetaC$Site %in% c("Left lung tumour", "Right lung tumour")] <- "Left/Right lung tumor"
## CD4MetaC <- CD4MetaC %>%
##     filter(Treatment == "Treatment-naive") %>%
##     filter(TumorType %in% c("LUSC")) %>%
##     filter(Site %in% c("Left/Right lung tumor", "LN metastasis", "Liver metastasis"))
## metaL_C <- list(CD8 = CD8MetaC, CD4 = CD4MetaC)

# Group D #####################################################################

CD8MetaD <- CD8_meta %>% mutate(Site=BiospySite)
CD8MetaD$Site[CD8MetaD$Site %in% c("Left lung tumour", "Right lung tumour")] <- "Left/Right lung tumor"
CD8MetaD <- CD8MetaD %>%
    filter(Site %in% c("Left/Right lung tumor")) %>%
    filter(TreatmentHx %in% c("Pre-treatment", "On-treatment")) %>%
    filter(Response %in% c("NR", "R")) %>%
    mutate(group = paste0(TreatmentHx, "_", Response))

tempCD8Obj <- subset(CD8Obj, cells = CD8MetaD$cellid)
tempCD8Obj$group <- CD8MetaD$group[match(Cells(tempCD8Obj), CD8MetaD$cellid)]
Idents(tempCD8Obj) <- tempCD8Obj$group

CD4MetaD <- CD4_meta %>% mutate(Site=BiospySite)
CD4MetaD$Site[CD4MetaD$Site %in% c("Left lung tumour", "Right lung tumour")] <- "Left/Right lung tumor"
CD4MetaD <- CD4MetaD %>%
    filter(Site %in% c("Left/Right lung tumor")) %>%
    filter(TreatmentHx %in% c("Pre-treatment", "On-treatment")) %>%
    filter(Response %in% c("NR", "R")) %>%
    mutate(group = paste0(TreatmentHx, "_", Response))

tempCD4Obj <- subset(CD4Obj, cells = CD4MetaD$cellid)
tempCD4Obj$group <- CD4MetaD$group[match(Cells(tempCD4Obj), CD4MetaD$cellid)]
Idents(tempCD4Obj) <- tempCD4Obj$group

metaL_D <- list(CD8 = CD8MetaD, CD4 = CD4MetaD)
obj_D <- list(CD8 = tempCD8Obj, CD4 = tempCD4Obj)

# Group E #####################################################################
CD8MetaE <- CD8_meta %>% mutate(Site=BiospySite)
CD8MetaE$Site[CD8MetaE$Site %in% c("Left lung tumour", "Right lung tumour")] <- "Left/Right lung tumor"
CD8MetaE <- CD8MetaE %>%
    filter(Site %in% c("LN metastasis")) %>%
    filter(TreatmentHx %in% c("Pre-treatment", "On-treatment")) %>%
    filter(Response %in% c("NR", "R")) %>%
    mutate(group = paste0(TreatmentHx, "_", Response))

tempCD8Obj <- subset(CD8Obj, cells = CD8MetaE$cellid)
tempCD8Obj$group <- CD8MetaE$group[match(Cells(tempCD8Obj), CD8MetaE$cellid)]
Idents(tempCD8Obj) <- tempCD8Obj$group

CD4MetaE <- CD4_meta %>% mutate(Site=BiospySite)
CD4MetaE$Site[CD4MetaE$Site %in% c("Left lung tumour", "Right lung tumour")] <- "Left/Right lung tumor"
CD4MetaE <- CD4MetaE %>%
    filter(Site %in% c("LN metastasis")) %>%
    filter(TreatmentHx %in% c("Pre-treatment", "On-treatment")) %>%
    filter(Response %in% c("NR", "R")) %>%
    mutate(group = paste0(TreatmentHx, "_", Response))

tempCD4Obj <- subset(CD4Obj, cells = CD4MetaE$cellid)
tempCD4Obj$group <- CD4MetaE$group[match(Cells(tempCD4Obj), CD4MetaE$cellid)]
Idents(tempCD4Obj) <- tempCD4Obj$group

metaL_E <- list(CD8 = CD8MetaE, CD4 = CD4MetaE)
obj_E <- list(CD8 = tempCD8Obj, CD4 = tempCD4Obj)


hs_genes <- c("HSPA1A", "HSPA1B")
# plot ########################################################################


for(oi in 1:length(obj_D)){
    tempName <- names(obj_D[oi])
    tempObj <- obj_D[[oi]]
    Idents(tempObj) <- tempObj$group
    for(hg in hs_genes){
        pdf(file.path(getwd(), paste0(tempName,
                                      "_groupD_",
                                      "_gene_",
                                      hg,
                                      "_violinplot.pdf")), width = 10, height = 6)
        print(VlnPlot(tempObj, features = hg, pt.size = 0) +
              geom_boxplot(color = "black", fill = NA, outlier.shape = NA, width = 0.2) +
              ggpubr::stat_compare_means())
        dev.off()
    }

    targetMD <- tempObj@meta.data
    targetMD$HSPA1A <- 0
    targetMD$HSPA1A <- tempObj@assays$RNA@data["HSPA1A",]
    targetMD$HSPA1B <- 0
    targetMD$HSPA1B <- tempObj@assays$RNA@data["HSPA1B",]

    targetMD <- as_tibble(targetMD) %>%
        select(group, HSPA1A, HSPA1B) %>%
        gather(key = "Marker", value = "Expression", -group)

    g <- ggstatsplot::ggbetweenstats(
                          data = targetMD %>% filter(Marker == "HSPA1A"),
                          pairwise.display = "all",
                          p.adjust.method = "fdr",
                          x = group,
                          y = Expression)
    ggsave(file.path(figurePath, paste0(
                                     tempName,
                                     "_groupD_",
                                     "_gene_",
                                     "_HSPA1A_group_ggbetweenstats.pdf")),
           g, width = 210, height = 210, units = "mm")
    g <- ggstatsplot::ggbetweenstats(
                          data = targetMD %>% filter(Marker == "HSPA1B"),
                          pairwise.display = "all",
                          p.adjust.method = "fdr",
                          x = group,
                          y = Expression)
    ggsave(file.path(figurePath, paste0(
                                     tempName,
                                     "_groupD_",
                                     "_gene_",
                                     "_HSPA1B_group_ggbetweenstats.pdf")),
           g, width = 210, height = 210, units = "mm")
}

for(oi in 1:length(obj_E)){
    tempName <- names(obj_E[oi])
    tempObj <- obj_E[[oi]]
    Idents(tempObj) <- tempObj$group
    for(hg in hs_genes){
        pdf(file.path(getwd(), paste0(tempName,
                                      "_groupE_",
                                      "_gene_",
                                      hg,
                                      "_violinplot.pdf")), width = 10, height = 6)
        print(VlnPlot(tempObj, features = hg, pt.size = 0) +
              geom_boxplot(color = "black", fill = NA, outlier.shape = NA, width = 0.2) +
              ggpubr::stat_compare_means())
        dev.off()
    }


    targetMD <- tempObj@meta.data
    targetMD$HSPA1A <- 0
    targetMD$HSPA1A <- tempObj@assays$RNA@data["HSPA1A",]
    targetMD$HSPA1B <- 0
    targetMD$HSPA1B <- tempObj@assays$RNA@data["HSPA1B",]

    targetMD <- as_tibble(targetMD) %>%
        select(group, HSPA1A, HSPA1B) %>%
        gather(key = "Marker", value = "Expression", -group)

    g <- ggstatsplot::ggbetweenstats(
                          data = targetMD %>% filter(Marker == "HSPA1A"),
                          pairwise.display = "all",
                          p.adjust.method = "fdr",
                          x = group,
                          y = Expression)
    ggsave(file.path(figurePath, paste0(
                                     tempName,
                                     "_groupE_",
                                     "_gene_",
                                     "_HSPA1A_group_ggbetweenstats.pdf")),
           g, width = 210, height = 210, units = "mm")

    g <- ggstatsplot::ggbetweenstats(
                          data = targetMD %>% filter(Marker == "HSPA1B"),
                          pairwise.display = "all",
                          p.adjust.method = "fdr",
                          x = group,
                          y = Expression)
    ggsave(file.path(figurePath, paste0(
                                     tempName,
                                     "_groupE_",
                                     "_gene_",
                                     "_HSPA1B_group_ggbetweenstats.pdf")),
           g, width = 210, height = 210, units = "mm")
}
