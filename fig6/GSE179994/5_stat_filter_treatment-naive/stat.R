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

figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE179994/5_stat_filter_treatment-naive/outs")
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

# Group A #####################################################################

CD8MetaA <- CD8_meta %>% mutate(Site=BiospySite)
CD8MetaA$Site[CD8MetaA$Site %in% c("Left lung tumour", "Right lung tumour")] <- "Left/Right lung tumor"
CD8MetaA <- CD8MetaA %>%
    filter(Treatment == "Treatment-naive") %>%
    filter(TumorType %in% c("LUAD", "LUSC")) %>%
    filter(Site %in% c("Left/Right lung tumor", "LN metastasis", "Liver metastasis"))
CD4MetaA <- CD4_meta %>% mutate(Site=BiospySite)
CD4MetaA$Site[CD4MetaA$Site %in% c("Left lung tumour", "Right lung tumour")] <- "Left/Right lung tumor"
CD4MetaA <- CD4MetaA %>%
    filter(Treatment == "Treatment-naive") %>%
    filter(TumorType %in% c("LUAD", "LUSC")) %>%
    filter(Site %in% c("Left/Right lung tumor", "LN metastasis", "Liver metastasis"))
metaL_A <- list(CD8 = CD8MetaA, CD4 = CD4MetaA)

# Group B #####################################################################

CD8MetaB <- CD8_meta %>% mutate(Site=BiospySite)
CD8MetaB$Site[CD8MetaB$Site %in% c("Left lung tumour", "Right lung tumour")] <- "Left/Right lung tumor"
CD8MetaB <- CD8MetaB %>%
    filter(Treatment == "Treatment-naive") %>%
    filter(TumorType %in% c("LUAD")) %>%
    filter(Site %in% c("Left/Right lung tumor", "LN metastasis", "Liver metastasis"))
CD4MetaB <- CD4_meta %>% mutate(Site=BiospySite)
CD4MetaB$Site[CD4MetaB$Site %in% c("Left lung tumour", "Right lung tumour")] <- "Left/Right lung tumor"
CD4MetaB <- CD4MetaB %>%
    filter(Treatment == "Treatment-naive") %>%
    filter(TumorType %in% c("LUAD")) %>%
    filter(Site %in% c("Left/Right lung tumor", "LN metastasis", "Liver metastasis"))
metaL_B <- list(CD8 = CD8MetaB, CD4 = CD4MetaB)

# Group C #####################################################################

CD8MetaC <- CD8_meta %>% mutate(Site=BiospySite)
CD8MetaC$Site[CD8MetaC$Site %in% c("Left lung tumour", "Right lung tumour")] <- "Left/Right lung tumor"
CD8MetaC <- CD8MetaC %>%
    filter(Treatment == "Treatment-naive") %>%
    filter(TumorType %in% c("LUSC")) %>%
    filter(Site %in% c("Left/Right lung tumor", "LN metastasis", "Liver metastasis"))
CD4MetaC <- CD4_meta %>% mutate(Site=BiospySite)
CD4MetaC$Site[CD4MetaC$Site %in% c("Left lung tumour", "Right lung tumour")] <- "Left/Right lung tumor"
CD4MetaC <- CD4MetaC %>%
    filter(Treatment == "Treatment-naive") %>%
    filter(TumorType %in% c("LUSC")) %>%
    filter(Site %in% c("Left/Right lung tumor", "LN metastasis", "Liver metastasis"))
metaL_C <- list(CD8 = CD8MetaC, CD4 = CD4MetaC)


# Group D #####################################################################

CD8MetaD <- CD8_meta %>% mutate(Site=BiospySite)
CD8MetaD$Site[CD8MetaD$Site %in% c("Left lung tumour", "Right lung tumour")] <- "Left/Right lung tumor"
CD8MetaD <- CD8MetaD %>%
    filter(Site %in% c("Left/Right lung tumor")) %>%
    filter(TreatmentHx %in% c("Pre-treatment", "On-treatment")) %>%
    filter(Response %in% c("NR", "R")) %>%
    mutate(group = paste0(TreatmentHx, "_", Response))

CD4MetaD <- CD4_meta %>% mutate(Site=BiospySite)
CD4MetaD$Site[CD4MetaD$Site %in% c("Left lung tumour", "Right lung tumour")] <- "Left/Right lung tumor"
CD4MetaD <- CD4MetaD %>%
    filter(Site %in% c("Left/Right lung tumor")) %>%
    filter(TreatmentHx %in% c("Pre-treatment", "On-treatment")) %>%
    filter(Response %in% c("NR", "R")) %>%
    mutate(group = paste0(TreatmentHx, "_", Response))
metaL_D <- list(CD8 = CD8MetaD, CD4 = CD4MetaD)


# Group E #####################################################################
CD8MetaE <- CD8_meta %>% mutate(Site=BiospySite)
CD8MetaE$Site[CD8MetaE$Site %in% c("Left lung tumour", "Right lung tumour")] <- "Left/Right lung tumor"
CD8MetaE <- CD8MetaE %>%
    filter(Site %in% c("LN metastasis")) %>%
    filter(TreatmentHx %in% c("Pre-treatment", "On-treatment")) %>%
    filter(Response %in% c("NR", "R")) %>%
    mutate(group = paste0(TreatmentHx, "_", Response))

CD4MetaE <- CD4_meta %>% mutate(Site=BiospySite)
CD4MetaE$Site[CD4MetaE$Site %in% c("Left lung tumour", "Right lung tumour")] <- "Left/Right lung tumor"
CD4MetaE <- CD4MetaE %>%
    filter(Site %in% c("LN metastasis")) %>%
    filter(TreatmentHx %in% c("Pre-treatment", "On-treatment")) %>%
    filter(Response %in% c("NR", "R")) %>%
    mutate(group = paste0(TreatmentHx, "_", Response))
metaL_E <- list(CD8 = CD8MetaE, CD4 = CD4MetaE)


totalMD <- bind_rows(CD8_meta, CD4_meta)
TotalSampleCellNum <- totalMD %>%
    group_by(Patient, Sample) %>%
    count()


groupDEL <- list(groupD = metaL_D, groupE = metaL_E)
for(gi in 1:length(groupDEL)){
    groupName <- names(groupDEL[gi])
    tempGroupL <- groupDEL[[gi]]
    for(ggi in 1:length(tempGroupL)){
        tempName <- names(tempGroupL[ggi])
        tempMD <- tempGroupL[[ggi]]
        ClusterT <- tempMD %>%
            group_by(predicted.celltype) %>%
            dplyr::count()
        N <- dim(tempMD)[1]
        totalTargetTcellMD <- tempMD
        matRoe <- matrix(rep(0, length(unique(tempMD$predicted.celltype)) * length(unique(tempMD$group))),
                         nrow = length(unique(tempMD$predicted.celltype)),
                         ncol = length(unique(tempMD$group)))
        matCell <- matrix(rep(0, length(unique(tempMD$predicted.celltype)) * length(unique(tempMD$group))),
                          nrow = length(unique(tempMD$predicted.celltype)),
                          ncol = length(unique(tempMD$group)))
        rownames(matRoe) <- sort(unique(tempMD$predicted.celltype))
        colnames(matRoe) <- unique(tempMD$group)
        rownames(matCell) <- sort(unique(tempMD$predicted.celltype))
        colnames(matCell) <- unique(tempMD$group)
        ggplotdata <- c()
        for(tempGroup in unique(tempMD$group)){
            targetTcellMD <- totalTargetTcellMD %>% filter(group == tempGroup)
            M <- dim(targetTcellMD)[1]
            for(gsci in unique(tempMD$predicted.celltype)){
                k <- ClusterT$n[ClusterT$predicted.celltype == gsci]
                if(gsci %in% targetTcellMD$predicted.celltype){
                    n <- dim(targetTcellMD %>% filter(predicted.celltype == gsci))[1]
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
        matRoe <- matRoe[rowSums(matCell) > 100,]
        matRoe[matRoe>10] = 10
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
        pdf(file.path(figurePath, paste0(groupName, "-", tempName, "_response_roe_heatmap.pdf")), width = widthi, height = heighti)
        print(pheatmap(matRoe,
                       color = my.colors,
                       breaks = my.breaks,
                       show_colnames = T,
                       show_rownames = T,
                       cluster_cols = F,
                       cluster_rows = F,
                       border_color = F))
        dev.off()
        pdf(file.path(figurePath, paste0(groupName, "-", tempName, "_response_cell_heatmap.pdf")), width = widthi * 1.5, height = heighti)
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


groupABCL <- list(groupA = metaL_A, groupB = metaL_B, groupC = metaL_C)

for(gi in 1:length(groupABCL)){
    groupName <- names(groupABCL[gi])
    tempGroupL <- groupABCL[[gi]]
    for(ggi in 1:length(tempGroupL)){
        tempName <- names(tempGroupL[ggi])
        tempMD <- tempGroupL[[ggi]]
        ClusterT <- tempMD %>%
            group_by(predicted.celltype) %>%
            dplyr::count()
        N <- dim(tempMD)[1]
        totalTargetTcellMD <- tempMD
        matRoe <- matrix(rep(0, length(unique(tempMD$predicted.celltype)) * length(unique(tempMD$Site))),
                         nrow = length(unique(tempMD$predicted.celltype)),
                         ncol = length(unique(tempMD$Site)))
        matCell <- matrix(rep(0, length(unique(tempMD$predicted.celltype)) * length(unique(tempMD$Site))),
                          nrow = length(unique(tempMD$predicted.celltype)),
                          ncol = length(unique(tempMD$Site)))
        rownames(matRoe) <- sort(unique(tempMD$predicted.celltype))
        colnames(matRoe) <- unique(tempMD$Site)
        rownames(matCell) <- sort(unique(tempMD$predicted.celltype))
        colnames(matCell) <- unique(tempMD$Site)
        ggplotdata <- c()
        for(tempSite in unique(tempMD$Site)){
            targetTcellMD <- totalTargetTcellMD %>% filter(Site == tempSite)
            M <- dim(targetTcellMD)[1]
            for(gsci in unique(tempMD$predicted.celltype)){
                k <- ClusterT$n[ClusterT$predicted.celltype == gsci]
                if(gsci %in% targetTcellMD$predicted.celltype){
                    n <- dim(targetTcellMD %>% filter(predicted.celltype == gsci))[1]
                } else{
                    n <- 0
                }
                tempData <- tibble(Site = tempSite,
                                   Cluster = gsci,
                                   Roe = (n/M) / (k/N))
                ggplotdata <- bind_rows(ggplotdata, tempData)
                tempRowName <- as.character(gsci)
                tempColName <- tempSite
                matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
                matCell[tempRowName, tempColName] <- n
            }
        }
        matRoe <- matRoe[rowSums(matCell) > 100,]
        matRoe[matRoe>10] = 10
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
        pdf(file.path(figurePath, paste0(groupName, "-", tempName, "_response_roe_heatmap.pdf")), width = widthi, height = heighti)
        print(pheatmap(matRoe,
                       color = my.colors,
                       breaks = my.breaks,
                       show_colnames = T,
                       show_rownames = T,
                       cluster_cols = F,
                       cluster_rows = F,
                       border_color = F))
        dev.off()
        pdf(file.path(figurePath, paste0(groupName, "-", tempName, "_response_cell_heatmap.pdf")), width = widthi * 1.5, height = heighti)
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


