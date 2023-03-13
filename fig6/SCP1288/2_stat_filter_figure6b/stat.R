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

figure_path <- file.path("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/SCP1288/2_stat_filter_figure6b/")
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

CD4_meta <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/SCP1288/1_mapping_filter/CD4/meta_2022-11-05.tsv")
CD4_meta$predicted.celltype <- plyr::mapvalues(CD4_meta$predicted.celltype,
                                               0:(length(CD4_CellType) - 1),
                                               CD4_CellType)

CD4_meta %>%
  group_by(ICB_Exposed, ICB_Response, TKI_Exposed, InferCNV, organ__ontology_label) %>%
  count %>%
  as.data.frame

CD8_meta <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/SCP1288/1_mapping_filter/CD8/meta_2022-11-05.tsv")
CD8_meta$predicted.celltype <- plyr::mapvalues(CD8_meta$predicted.celltype,
                                               0:(length(CD8_CellType) - 1),
                                               CD8_CellType)

CD8_meta %>%
  group_by(ICB_Exposed, ICB_Response, TKI_Exposed, InferCNV, organ__ontology_label) %>%
  count %>%
  as.data.frame


total_meta <- bind_rows(CD4_meta, CD8_meta)

N <- dim(total_meta)[1]
matRoe <- matrix(rep(0, length(unique(total_meta$predicted.celltype)) *
                        length(unique(total_meta$ICB_Exposed))),
                 nrow = length(unique(total_meta$predicted.celltype)),
                 ncol = length(unique(total_meta$ICB_Exposed)))

matCell <- matrix(rep(0, length(unique(total_meta$predicted.celltype)) *
                        length(unique(total_meta$ICB_Exposed))),
                 nrow = length(unique(total_meta$predicted.celltype)),
                 ncol = length(unique(total_meta$ICB_Exposed)))

rownames(matRoe) <- unique(total_meta$predicted.celltype)
colnames(matRoe) <- unique(total_meta$ICB_Exposed)
rownames(matCell) <- unique(total_meta$predicted.celltype)
colnames(matCell) <- unique(total_meta$ICB_Exposed)

ClusterT <- total_meta %>%
  group_by(predicted.celltype) %>%
  dplyr::count()

for(tempGroup in unique(total_meta$ICB_Exposed)) {
  targetTcellMD <- total_meta %>% filter(ICB_Exposed == tempGroup)
  M <- dim(targetTcellMD)[1]
  for(gsci in unique(total_meta$predicted.celltype)) {
    k <- ClusterT$n[ClusterT$predicted.celltype == gsci]
    if(gsci %in% targetTcellMD$predicted.celltype) {
      n <- dim(targetTcellMD %>% filter(predicted.celltype == gsci))[1]
    } else {
      n <- 0
    }
    tempRowName <- gsci
    tempColName <- tempGroup
    matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
    matCell[tempRowName, tempColName] <- n
  }
}

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

tMR <- tMR[rowSums(tMC) > 100,]
tMC <- tMC[rowSums(tMC) > 100,]
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
pdf(file.path(figure_path, paste0("figure6h_Roe_1.pdf")), width = widthi, height = heighti)
print(pheatmap(tMR,
               color = my.colors,
               breaks = my.breaks,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figure_path, paste0("figure6h_gt100_Cell_1.pdf")), width = widthi, height = heighti)
print(pheatmap(tMC,
               show_colnames = T,
               show_rownames = T,
               display_numbers = T,
               number_format = "%.0f",
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()





total_meta <- bind_rows(CD4_meta, CD8_meta)
N <- dim(total_meta)[1]
ClusterT <- total_meta %>%
  group_by(predicted.celltype) %>%
  dplyr::count()
## total_meta %>%
##   group_by(ICB_Exposed, ICB_Response, TKI_Exposed, InferCNV, organ__ontology_label) %>%
##   count %>%
##   as.data.frame
total_meta <- total_meta %>%
  filter(ICB_Response %in% c("ICB_PD", "ICB_PR", "ICB_SD")) %>%
  filter(TKI_Exposed %in% "TKI")
total_meta$group <- plyr::mapvalues(total_meta$ICB_Response,
                                    c("ICB_PD", "ICB_PR", "ICB_SD"),
                                             c("NR", "R", "NR"))
matRoe <- matrix(rep(0, length(unique(total_meta$predicted.celltype)) *
                        length(unique(total_meta$group))),
                 nrow = length(unique(total_meta$predicted.celltype)),
                 ncol = length(unique(total_meta$group)))
matCell <- matrix(rep(0, length(unique(total_meta$predicted.celltype)) *
                        length(unique(total_meta$group))),
                 nrow = length(unique(total_meta$predicted.celltype)),
                 ncol = length(unique(total_meta$group)))
rownames(matRoe) <- unique(total_meta$predicted.celltype)
colnames(matRoe) <- unique(total_meta$group)
rownames(matCell) <- unique(total_meta$predicted.celltype)
colnames(matCell) <- unique(total_meta$group)
for(tempGroup in unique(total_meta$group)) {
  targetTcellMD <- total_meta %>% filter(group == tempGroup)
  M <- dim(targetTcellMD)[1]
  for(gsci in unique(total_meta$predicted.celltype)) {
    k <- ClusterT$n[ClusterT$predicted.celltype == gsci]
    if(gsci %in% targetTcellMD$predicted.celltype) {
      n <- dim(targetTcellMD %>% filter(predicted.celltype == gsci))[1]
    } else {
      n <- 0
    }
    tempRowName <- gsci
    tempColName <- tempGroup
    matRoe[tempRowName, tempColName] <- (n/M) / (k/N)
    matCell[tempRowName, tempColName] <- n
  }
}
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
tMR <- tMR[rowSums(tMC) > 100,]
tMC <- tMC[rowSums(tMC) > 100,]
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
pdf(file.path(figure_path, paste0("figure6h_Roe_2.pdf")), width = widthi, height = heighti)
print(pheatmap(tMR,
               color = my.colors,
               breaks = my.breaks,
               show_colnames = T,
               show_rownames = T,
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()
pdf(file.path(figure_path, paste0("figure6h_gt100_Cell_2.pdf")), width = widthi, height = heighti)
print(pheatmap(tMC,
               show_colnames = T,
               show_rownames = T,
               display_numbers = T,
               number_format = "%.0f",
               cluster_cols = F,
               cluster_rows = F,
               border_color = F))
dev.off()





CD4_obj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/SCP1288/0_merge/CD4.rds") %>%
  NormalizeData(verbose = T)
CD8_obj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/SCP1288/0_merge/CD8.rds") %>%
  NormalizeData(verbose = T)

CD4_obj$predicted.celltype <- CD4_meta$predicted.celltype[match(Cells(CD4_obj), CD4_meta$ID)]
targetCells <- Cells(CD4_obj)[
  CD4_obj$ICB_Response %in% c("ICB_PR", "ICB_PD", "ICB_SD") &
  CD4_obj$TKI_Exposed %in% c("TKI")]
temp_obj <- subset(CD4_obj, cells = targetCells)
temp_obj$group <- plyr::mapvalues(temp_obj$ICB_Response,
                                  c("ICB_PD", "ICB_PR", "ICB_SD"),
                                  c("NR", "R", "NR"))
Idents(temp_obj) <- temp_obj$group
NR_obj <- subset(temp_obj, idents = "NR")
Idents(NR_obj) <- NR_obj$predicted.celltype
p <- DotPlot(NR_obj, features = c("HSPA1A", "HSPA1B", "FOXP3"))
data_NR_CD4 <- p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]
data_NR_CD4$type <- "NR"
Idents(temp_obj) <- temp_obj$group
R_obj <- subset(temp_obj, idents = "R")
Idents(R_obj) <- R_obj$predicted.celltype
p <- DotPlot(R_obj, features = c("HSPA1A", "HSPA1B", "FOXP3"))
data_R_CD4 <- p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]
data_R_CD4$type <- "R"
data_CD4 <- bind_rows(data_NR_CD4, data_R_CD4)
data_CD4$data_type <- "CD4"


CD8_obj$predicted.celltype <- CD8_meta$predicted.celltype[match(Cells(CD8_obj), CD8_meta$ID)]
targetCells <- Cells(CD8_obj)[
  CD8_obj$ICB_Response %in% c("ICB_PR", "ICB_PD", "ICB_SD") &
  CD8_obj$TKI_Exposed %in% c("TKI")]
temp_obj <- subset(CD8_obj, cells = targetCells)
temp_obj$group <- plyr::mapvalues(temp_obj$ICB_Response,
                                             c("ICB_PD", "ICB_PR", "ICB_SD"),
                                             c("NR", "R", "NR"))
Idents(temp_obj) <- temp_obj$group
NR_obj <- subset(temp_obj, idents = "NR")
Idents(NR_obj) <- NR_obj$predicted.celltype
p <- DotPlot(NR_obj, features = c("HSPA1A", "HSPA1B", "FOXP3"))
data_NR_CD8 <- p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]
data_NR_CD8$type <- "NR"
Idents(temp_obj) <- temp_obj$group
R_obj <- subset(temp_obj, idents = "R")
Idents(R_obj) <- R_obj$predicted.celltype
p <- DotPlot(R_obj, features = c("HSPA1A", "HSPA1B", "FOXP3"))
data_R_CD8 <- p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]
data_R_CD8$type <- "R"
data_CD8 <- bind_rows(data_NR_CD8, data_R_CD8)
data_CD8$data_type <- "CD8"

data <- bind_rows(data_CD4, data_CD8) %>%
  filter(id %in% c("CD4_c4_Tstr",
                   "CD8_c4_Tstr",
                   "CD8_c1_Tex",
                   "CD4_c1_Treg",
                   "CD8_c2_Teff",
                   "CD8_c3_Tn",
                   "CD4_c2_Tn",
                   "CD8_c0_t-Teff",
                   "CD4_c0_Tcm",
                   "CD4_c5_Tctl"))
data$id <- factor(data$id, levels = c("CD4_c4_Tstr",
                                      "CD8_c4_Tstr",
                                      "CD8_c1_Tex",
                                      "CD4_c1_Treg",
                                      "CD8_c2_Teff",
                                      "CD8_c3_Tn",
                                      "CD4_c2_Tn",
                                      "CD8_c0_t-Teff",
                                      "CD4_c0_Tcm",
                                      "CD4_c5_Tctl"))
data$type <- factor(data$type, levels = c("R", "NR"))
plotx <- ggplot(data, aes(x = id, y = features.plot)) +
  geom_point(aes(fill = avg.exp.scaled, size = pct.exp),
             color = 'black',
             shape = 21,
             stroke = 0.01) +
  xlab("") + ylab("") +
  scale_fill_gradientn(
    colors = c("#5DBCFF", "#6DCCFF", "white", "#F494BE", "#F484AE")) +
  scale_size(range = c(0, 3.2), limits = c(0, 100),
             breaks = c(0, 20, 40, 60, 80, 100)) +
  scale_y_discrete(limits = rev) +
  facet_grid(. ~ type, space = "free") +
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
                             byrow = T,
                             override.aes = list(stroke = 0.4)),
         fill = guide_colourbar(title.position = "top", title.hjust = 0.5))
ggsave(file.path(figure_path, "CD8CD4_bubbleplot.pdf"),
       plotx,
       width = 210/2,
       height = 297/4,
       units = "mm")










CD4_meta$HSPA1A <- min(CD4_obj@assays$RNA@data["HSPA1A",])
CD4_meta$HSPA1A <- CD4_obj@assays$RNA@data["HSPA1A",][match(CD4_meta$ID, Cells(CD4_obj))]
CD4_meta$HSPA1B <- min(CD4_obj@assays$RNA@data["HSPA1B",])
CD4_meta$HSPA1B <- CD4_obj@assays$RNA@data["HSPA1B",][match(CD4_meta$ID, Cells(CD4_obj))]
CD4_meta$data <- "CD4"

CD8_meta$HSPA1A <- min(CD8_obj@assays$RNA@data["HSPA1A",])
CD8_meta$HSPA1A <- CD8_obj@assays$RNA@data["HSPA1A",][match(CD8_meta$ID, Cells(CD8_obj))]
CD8_meta$HSPA1B <- min(CD8_obj@assays$RNA@data["HSPA1B",])
CD8_meta$HSPA1B <- CD8_obj@assays$RNA@data["HSPA1B",][match(CD8_meta$ID, Cells(CD8_obj))]
CD8_meta$data <- "CD8"

totalMD <- bind_rows(CD4_meta, CD8_meta) 

# group 1 #####################################################################

targetMD_group1 <- totalMD %>%
  rename(sub_group = ICB_Exposed) %>% 
  select(data, sub_group, HSPA1A, HSPA1B) %>%
  gather(key = "Marker", value = "Expression", -data, -sub_group) %>%
  group_by(data, sub_group, Marker) %>%
  mutate(MD = median(Expression))

targetMD_group1$sub_group <- factor(targetMD_group1$sub_group, levels = c("NoICB", "ICB"))
g <- targetMD_group1 %>%
    ggplot() +
    geom_violin(aes(x = sub_group,
                    y = Expression,
                    fill = MD),
                lwd = 0.1,
                scale = "width") +
    stat_summary(aes(x = sub_group,
                     y = Expression),
                 fun = "mean",
                 geom = "point",
                 colour = "black",
                 size = 0.05) +
    scale_fill_gradient(low ="#6DCCFF",
                        high = "#F494BE") +
    facet_grid(data ~ Marker) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0))
ggsave(file.path(getwd(), "figure6h_violinplot_group1.pdf"), g, width = 4, height = 3)

# group 2 #####################################################################

totalMD <- bind_rows(CD4_meta, CD8_meta) 
targetMD_group2 <- totalMD %>%
  filter(ICB_Response %in% c("ICB_PD", "ICB_PR", "ICB_SD")) %>%
  filter(TKI_Exposed %in% "TKI")
targetMD_group2$sub_group <- plyr::mapvalues(targetMD_group2$ICB_Response,
                                             c("ICB_PD", "ICB_PR", "ICB_SD"),
                                             c("NR", "R", "NR"))
targetMD_group2 <- targetMD_group2 %>% 
  select(data, sub_group, HSPA1A, HSPA1B) %>%
  gather(key = "Marker", value = "Expression", -data, -sub_group) %>%
  group_by(data, sub_group, Marker) %>%
  mutate(MD = mean(Expression))

g <- targetMD_group2 %>%
  ggplot() +
  geom_violin(aes(x = sub_group,
                  y = Expression,
                  fill = MD),
              lwd = 0.1,
              scale = "width") +
  stat_summary(aes(x = sub_group,
                   y = Expression),
               fun = "mean",
               geom = "point",
               colour = "black",
               size = 0.05) +
  scale_fill_gradient(low ="#6DCCFF",
                      high = "#F494BE") +
  facet_grid(data ~ Marker) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0))
ggsave(file.path(getwd(), "figure6h_violinplot_group2.pdf"), g, width = 4, height = 3)
