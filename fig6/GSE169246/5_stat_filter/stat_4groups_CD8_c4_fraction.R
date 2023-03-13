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

CD8_meta <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/MappingResult_filter/CD8/meta_2022-05-25.tsv")
CD8_meta$Sample <- stringr::str_extract(CD8_meta$ID, "(?<=^.{10,20}\\.).+")
CD8_meta$Patient <- stringr::str_extract(CD8_meta$Sample, "P\\d+")
CD8_meta$Tissue <- stringr::str_extract(CD8_meta$Sample, "\\w$")
CD8_meta$TumorTreatment <- stringr::str_extract(CD8_meta$Sample, "^[a-zA-Z]+")
CD8_meta$isResponse <- "NR"
CD8_meta$isResponse[CD8_meta$Patient %in% AllResponsePatients] <- "R"
CD8_meta$TreatmentType <- "PDL1+Chemo"
CD8_meta$TreatmentType[CD8_meta$Patient %in% TotalChemoPatients] <- "Chemo"

CD4_meta <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/MappingResult_filter/CD4/meta_2022-05-25.tsv")
CD4_meta$Sample <- stringr::str_extract(CD4_meta$ID, "(?<=^.{10,20}\\.).+")
CD4_meta$Patient <- stringr::str_extract(CD4_meta$Sample, "P\\d+")
CD4_meta$Tissue <- stringr::str_extract(CD4_meta$Sample, "\\w$")
CD4_meta$TumorTreatment <- stringr::str_extract(CD4_meta$Sample, "^[a-zA-Z]+")
CD4_meta$isResponse <- "NR"
CD4_meta$isResponse[CD4_meta$Patient %in% AllResponsePatients] <- "R"
CD4_meta$TreatmentType <- "PDL1+Chemo"
CD4_meta$TreatmentType[CD4_meta$Patient %in% TotalChemoPatients] <- "Chemo"

allMD <- list(CD8 = CD8_meta, CD4 = CD4_meta)

allT <- bind_rows(CD4_meta, CD8_meta)
allTSampleCN <- allT %>%
    group_by(Patient, Sample) %>%
    count

targetTreatments = c("Pre", "Post")
targetCellSet = c("AllT", "Split")
targetTreatmentType = c("PDL1+Chemo", "Chemo")


figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/5_stat_filter/outs4groups")
if(!dir.exists(figurePath)){
    dir.create(figurePath, recursive = T)
}
setwd(figurePath)


# CD8 fraction ################################################################

CD8C4Frac = CD8_meta %>%
    filter(Tissue == "t") %>% 
    filter(TumorTreatment != "Prog") %>% 
    group_by(Sample, TreatmentType, TumorTreatment, isResponse, predicted.celltype) %>%
    count %>%
    mutate(frac = n / dim(CD8_meta %>% filter(Tissue == "t"))[1]) %>%
    filter(predicted.celltype == 4) %>%
    mutate(group = paste0(TreatmentType, "-", isResponse))

groupOrder <- c("Chemo-NR", "Chemo-R", "PDL1+Chemo-NR", "PDL1+Chemo-R")
CD8C4Frac$group <- factor(CD8C4Frac$group, levels = groupOrder)

g <- ggstatsplot::ggbetweenstats(
                      data = CD8C4Frac %>% filter(TumorTreatment == "Pre"),
                      pairwise.display = "all",
                      p.adjust.method = "fdr",
                      x = group,
                      y = frac)
ggsave(file.path(figurePath, paste0("CD8C4_frac_Pre_group_ggbetweenstats.pdf")), g, width = 210, height = 210, units = "mm")

g <- ggstatsplot::ggbetweenstats(
                      data = CD8C4Frac %>% filter(TumorTreatment == "Post"),
                      pairwise.display = "all",
                      p.adjust.method = "fdr",
                      x = group,
                      y = frac)
ggsave(file.path(figurePath, paste0("CD8C4_frac_Post_group_ggbetweenstats.pdf")), g, width = 210, height = 210, units = "mm")


g <- ggplot(CD8C4Frac %>% filter(TumorTreatment != "Prog")) +
    geom_bar(aes(x = TreatmentType, y = frac, fill = TreatmentType), stat = "identity") +
    facet_grid(TumorTreatment ~ isResponse) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("CD8C4Frac.pdf", g)
 
