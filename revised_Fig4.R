################################################# ~ Good Luck ~ #########################################################
###### >>>>>>                       TITLE: H3K9me3 Project Revised Figure4 analysis                         <<<<<< ######

# - Figure 4 content:
# - 1st part: Library
# - 2nd part: K93/K273 signal and expression level of full length L1 regions (CMQ)
# - 3rd part: K9-K27 bivalent repeat age (CMQ)
# - 4th part: Number of K9-K27 bivalent repeats, genes and promoters (CMQ)
# - 5th part: Annotation of K9-K27 Bi-L1/SVA in cleavage (CMQ)
# - 6th part: Transcription activity of K9-K27 Bi-L1/SVA marked genes in 8cell (CMQ)
# - 7th part: Bivalent domains and splicing (YHW)


### =================
### 1st part: Library
### =================
### >>> 1. Library
library("ggplot2")
library("ggpubr")
library("ggrepel")
library("ComplexHeatmap")
library("circlize")
library("dplyr")
library("tidyr")
library("tibble")
### >>> 2. Functions
RPKM.Prepare <- function(working.dir){
  library("gtools")
  setwd(working.dir)
  # load data
  tab.file <- list.files("normalized_value/", ".tab", full.names = T)
  tab.file <- mixedsort(tab.file)
  sample <- read.table("rdata/samples_names.txt")
  sample <- sample$V1
  # load region
  gr <- read.table("final_regions_for_plotting.bed", sep = "\t")
  gr$region.id <- paste("region", 1:nrow(gr), sep = "_")
  gr$pos.id <- paste(gr$V1, gr$V2, gr$V3, sep = "-")
  # processing
  raw <- matrix(ncol = 220)
  for (i in seq(1, length(tab.file), 1)) {
    tmp <- read.table(tab.file[i])
    if (nrow(tmp) != 220) {
      next
    }else{
      colnames(tmp) <- c(paste(sample, "region", i, sep = "_"))
      rownames(tmp) <- c(1:220)
      tmp <- t(tmp)
      raw <- rbind(raw, tmp)
    }
  }
  raw <- as.data.frame(raw)[-1,]
  RPKM <- list()
  for (i in seq(1, length(sample), 1)) {
    RPKM[[i]] <- raw[grep(sample[i], rownames(raw)),]
    names(RPKM)[[i]] <- sample[i]
    RPKM[[i]] <- rbind(RPKM[[i]], colMeans(RPKM[[i]]))
  }
  res <- list(RPKM, gr)
  names(res) <- c("value", "region")
  return(res)
}
Plot.Heatmap.K9 <- function(dfl, plot.color.heatmap){
  pd.list <- list()
  for (i in seq(1, length(dfl))) {
    pd.list[[i]] <- na.omit(dfl[[i]])
    pd.list[[i]] <- pd.list[[i]][order(rowMeans(pd.list[[i]][,1:220]), decreasing = T),]
    names(pd.list)[[i]] <- names(dfl)[[i]]
  }; rm(i)
  p <- list()
  for (x in seq(1, length(pd.list))) {
    p1 <- Heatmap(as.matrix(pd.list[[x]][,1:220]),
                  use_raster = TRUE,
                  raster_resize_mat = mean,
                  #raster_magick_filter = "Undefined",
                  col = plot.color.heatmap,
                  cluster_rows = F,
                  cluster_columns = F,
                  show_column_names = F,
                  show_row_names = F, border = T,
                  width = unit(3.5, "cm"), height = unit(1.75, "cm"),
                  heatmap_legend_param = list(title = names(pd.list)[[x]]))
    p2 <- p1 + Heatmap(as.matrix(pd.list[[x]][,ncol(pd.list[[x]])]),
                       name = "  ", col = colorRamp2(c(0, 5, 10), c("#EAEA80", "#ffffff", "#FB4848")),
                       width = unit(0.4, "cm"), column_names_rot = 45, show_row_names = F,
                       cluster_rows = F, cluster_columns = F, border = T)
    p[[x]] <- grid::grid.grabExpr(draw(p2))
  }; rm(x)
  h <- plot_grid(plotlist = p, ncol =1)
  return(h)
}
Plot.Heatmap.K27 <- function(dfl, plot.color.heatmap){
  pd.list <- list()
  for (i in seq(1, length(dfl))) {
    pd.list[[i]] <- na.omit(dfl[[i]])
    pd.list[[i]] <- pd.list[[i]][order(rowMeans(pd.list[[i]][-nrow(pd.list[[i]]),]), decreasing = T),]
    names(pd.list)[[i]] <- names(dfl)[[i]]
  }; rm(i)
  p <- list()
  for (x in seq(1, length(pd.list))) {
    p1 <- Heatmap(as.matrix(pd.list[[x]][-nrow(pd.list[[x]]),]),
                  use_raster = TRUE,
                  raster_resize_mat = mean,
                  #raster_magick_filter = "Undefined",
                  col = plot.color.heatmap,
                  cluster_rows = F,
                  cluster_columns = F,
                  show_column_names = F,
                  show_row_names = F, border = T,
                  width = unit(3.5, "cm"), height = unit(1.75, "cm"),
                  heatmap_legend_param = list(title = names(pd.list)[[x]]))
    p[[x]] <- grid::grid.grabExpr(draw(p1))
  }; rm(x)
  h <- plot_grid(plotlist = p, ncol =1)
  return(h)
}
Plot.Line.Group.Scale <- function(dfl, lim){
  region.n <- nrow(dfl[[1]])-1
  win.n <- (ncol(dfl[[1]])-100)/2
  pd <- dfl[[1]][region.n+1, ] %>% rownames_to_column(var = "id") %>%
    as.data.frame() %>% gather(key = pos, value = value, -id)  %>%
    mutate(pos = as.numeric(pos), group = names(dfl)[1]) %>%
    group_by(id) %>% arrange(pos, .by_group = T)
  for(i in seq(2, length(dfl))){
    pda <- dfl[[i]][region.n+1, ] %>% rownames_to_column(var = "id") %>%
      as.data.frame() %>% gather(key = pos, value = value, -id)  %>%
      mutate(pos = as.numeric(pos), group = names(dfl)[i]) %>%
      group_by(id) %>% arrange(pos, .by_group = T)
    pd <- rbind(pd, pda)
  }
  p1 <- pd %>%
    ggplot(aes(x = pos, y = value)) +
    stat_smooth(method = "loess", span = 0.2, se = F, alpha = 0.3, fill = "#F52B07", aes(color = group), lwd = 1) +
    ylim(c(0, lim)) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(face = "plain", colour = "#000000", size = 15, angle = 90),
          axis.text.x  = element_blank(),
          axis.text.y  = element_text(face = "plain", colour = "#000000", size = 15, angle = 0),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "top", legend.direction = "vertical", legend.title = element_blank()) +
    geom_vline(xintercept = c(win.n, win.n + 100), colour = "#000000", linetype = "longdash", size = 0.5, alpha = 0.5)
  p2 <- ggplot(data.frame(xmin = 1, xmax = win.n, ymin = 0, ymax = 0.1)) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "#6FA5E4", alpha = 0.75, inherit.aes = FALSE) +
    theme_void() +
    geom_rect(data = data.frame(xmin = win.n, xmax = win.n + 100, ymin = 0, ymax = 0.1),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "#E48D6F",
              alpha = 0.75, inherit.aes = FALSE) +
    geom_rect(data = data.frame(xmin = win.n + 100, xmax = win.n + 100 + win.n, ymin = 0, ymax = 0.1),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "#6FA5E4",
              alpha = 0.75, inherit.aes = FALSE) +
    annotate("text", x = 6, y = 0.05, label= "-3kb", size = 5) +
    annotate("text", x = win.n, y = 0.05, label= "TSS", size = 5) +
    annotate("text", x = win.n + 100, y = 0.05, label= "TES", size = 5) +
    annotate("text", x = win.n + 94 + win.n, y = 0.05, label= "+3kb", size = 5)
  p <- plot_grid(p1, p2, ncol = 1, rel_heights = c(1, 0.1), align = 'v', hjust = 3)
  return(p)
}



### ==============================================================================
### 2nd part: K93/K273 signal and expression level of full length L1 regions (CMQ)
### ==============================================================================

setwd("/home/cmq/bioinfo/project-cmq/embryo_93")
### >>> 1. human H3K9me3 signal and expression level of fl-L1
# load H3K9me3 rpkm value on fl-L1
hs.fl.line1.k9 <- RPKM.Prepare("/home/yhw/bioinfo/project-mine/Embryo.93/chipseq/hs_k9me3/results/deeptools/multibw/merge/fl_LINE1")
# load expression level(FPKM) 0f fl-L1
setwd("/home/yhw/bioinfo/project-mine/Embryo.93/chipseq/hs_k9me3")
hs.fl.line1.rpkm <- read.table("results/deeptools/multibw/rnaseq/fl_LINE1/final_regions_for_plotting_readNor.tab", header = T, sep = "\t")
hs.fl.line1.rpkm <- hs.fl.line1.rpkm %>% mutate(pos.id = paste(chr, start, end, sep = "-")) %>%
  mutate(`2cell_ave` = (X2cell_rep1_1.bw + X2cell_rep1_2.bw + X2cell_rep2_1.bw + X2cell_rep2_2.bw)/4,
         `4cell_ave` = (X4cell_3PN_rep1.bw + X4cell_3PN_rep2_1.bw + X4cell_3PN_rep2_2.bw)/3,
         `8cell_ave` = (X8cell_rep1.bw + X8cell_rep2.bw)/2,
         ICM_ave = (ICM_rep1_1.bw + ICM_rep1_2.bw + ICM_rep2_1.bw + ICM_rep2_2.bw)/4,
         Morula_ave = (Morula_rep1.bw + Morula_rep2.bw + Morula_rep3.bw)/3)
hs.fl.line1.rpkm.te <- readRDS("~/bioinfo/project-cmq/embryo_93/results/R/Tables/L1_fl/human_embryo_TE_FPKM_on_LINE1_full_length_region.rds")
hs.fl.line1.rpkm.te <- hs.fl.line1.rpkm.te[,c(5:7, ncol(hs.fl.line1.rpkm.te))]
hs.fl.line1.rpkm.te$pos.id <- paste(hs.fl.line1.rpkm.te$chr_fl, hs.fl.line1.rpkm.te$start_fl, hs.fl.line1.rpkm.te$end_fl, sep = "-")
hs.fl.line1.rpkm <- merge(hs.fl.line1.rpkm, hs.fl.line1.rpkm.te[,4:5])
# merge H3K9me3 signal and expression level of the same fl-L1 regions
hs.fl.line1.k9$region <- merge(hs.fl.line1.k9$region, hs.fl.line1.rpkm[, grep("id|ave", colnames(hs.fl.line1.rpkm))])
hs.fl.line1.k9$region <- hs.fl.line1.k9$region[order(match(hs.fl.line1.k9$region$region.id, mixedsort(hs.fl.line1.k9$region$region.id))),]
hs.fl.line1.k9$value$`h4C93-Ip-1.bw`$`221` <- c(hs.fl.line1.k9$region$`4cell_ave`, NA)
hs.fl.line1.k9$value$`h8C93-Ip-3.bw`$`221` <- c(hs.fl.line1.k9$region$`8cell_ave`, NA)
hs.fl.line1.k9$value$`hM93-Ip-1.bw`$`221` <- c(hs.fl.line1.k9$region$`Morula_ave`, NA)
hs.fl.line1.k9$value$`hICM93-IP-1.bw`$`221` <- c(hs.fl.line1.k9$region$`ICM_ave`, NA)
hs.fl.line1.k9$value$`hTE93-IP-1.bw`$`221` <- c(hs.fl.line1.k9$region$`TE_ave`, NA)
hs.fl.line1.k9$value <- hs.fl.line1.k9$value[grep("IP|Ip", names(hs.fl.line1.k9$value))]
# line plot visualization
setwd("/home/cmq/bioinfo/project-cmq/embryo_93/")
pdf("results/R/Graphs/repeat_hs_mm_k9_k27/human_embryo_H3K9me3_signal_on_LINE1_full_length_region_line_plot.pdf", height = 5, width = 5)
Plot.Line.Group.Scale(hs.fl.line1.k9$value, 60)
dev.off()
# heatmap visualization
pdf("results/R/Graphs/repeat_hs_mm_k9_k27/human_embryo_H3K9me3_signal_on_LINE1_full_length_region_heatmap.pdf", height = 12, width = 5)
Plot.Heatmap.K9(hs.fl.line1.k9$value, colorRamp2(c(0, 30), c("#ffffff", "#069DC0")))
dev.off()

### >>> 2. human H3K27me3 signal on fl-L1
# load H3K27me3 rpkm value on fl-L1 
hs.fl.line1.k27 <- RPKM.Prepare("/home/yhw/bioinfo/project-mine/Embryo.93/chipseq/hs_k27me3/results/deeptools/multibw/merge/fl_LINE1")
# line plot visualization
setwd("/home/cmq/bioinfo/project-cmq/embryo_93/")
pdf("results/R/Graphs/repeat_hs_mm_k9_k27/human_embryo_H3K27me3_signal_on_LINE1_full_length_region_line_plot.pdf", height = 5, width = 5)
Plot.Line.Group.Scale(hs.fl.line1.k27$value, 40)
dev.off()
# heatmap visualization
pdf("results/R/Graphs/repeat_hs_mm_k9_k27/human_embryo_H3K27me3_signal_on_LINE1_full_length_region_heatmap.pdf", height = 12, width = 5)
Plot.Heatmap.K27(hs.fl.line1.k27$value, colorRamp2(c(0, 20), c("#ffffff", "#069DC0")))
dev.off()

### >>> 3. mouse H3K9me3 signal on fl-L1
# load H3K9me3 rpkm value on fl-L1
mm.fl.line1.k9 <- RPKM.Prepare("/home/yhw/bioinfo/project-mine/Embryo.93/chipseq/mm_k9me3/results/deeptools/multibw/merge/fl_LINE1")
setwd("/home/yhw/bioinfo/project-mine/Embryo.93/chipseq/mm_k9me3")
# load expression level(FPKM) 0f fl-L1
mm.fl.line1.rpkm <- read.table("results/deeptools/multibw/rnaseq/fl_LINE1/final_regions_for_plotting_readNor.tab", header = T, sep = "\t")
mm.fl.line1.rpkm <- mm.fl.line1.rpkm %>% mutate(pos.id = paste(chr, start, end, sep = "-")) %>%
  mutate(`2cell_ave` = (X2cell_rep1.bw + X2cell_rep2.bw + X2cell_rep3.bw + X2cell_rep4.bw)/4,
         `4cell_ave` = (X4cell_rep1.bw + X4cell_rep2.bw + X4cell_rep3.bw + X4cell_rep4.bw)/4,
         `8cell_ave` = (X8cell_rep1.bw + X8cell_rep2.bw+ X8cell_rep3.bw)/3,
         ICM_ave = (ICM_rep1.bw + ICM_rep2.bw + ICM_rep3.bw + ICM_rep4.bw)/4,
         Morula_ave = (Morula_rep1.bw + Morula_rep2.bw)/2,
         TE_ave = (TE_rep1.bw + TE_rep2.bw + TE_rep3.bw + TE_rep4.bw)/4)
# merge H3K9me3 signal and expression level of the same fl-L1 regions
mm.fl.line1.k9$region <- merge(mm.fl.line1.k9$region, mm.fl.line1.rpkm[, grep("id|ave", colnames(mm.fl.line1.rpkm))])
mm.fl.line1.k9$region <- mm.fl.line1.k9$region[order(match(mm.fl.line1.k9$region$region.id, mixedsort(mm.fl.line1.k9$region$region.id))),]
mm.fl.line1.k9$value$`2cell_93.bw`$`221` <- c(mm.fl.line1.k9$region$`2cell_ave`, NA)
mm.fl.line1.k9$value$`4cell_93.bw`$`221` <- c(mm.fl.line1.k9$region$`4cell_ave`, NA)
mm.fl.line1.k9$value$`8cell_93.bw`$`221` <- c(mm.fl.line1.k9$region$`8cell_ave`, NA)
mm.fl.line1.k9$value$`morula_93.bw`$`221` <- c(mm.fl.line1.k9$region$`Morula_ave`, NA)
mm.fl.line1.k9$value$`icm_93.bw`$`221` <- c(mm.fl.line1.k9$region$`ICM_ave`, NA)
mm.fl.line1.k9$value$`te_93.bw`$`221` <- c(mm.fl.line1.k9$region$`TE_ave`, NA)
# line plot visualization
setwd("/home/cmq/bioinfo/project-cmq/embryo_93/")
pdf("results/R/Graphs/repeat_hs_mm_k9_k27/mouse_embryo_H3K9me3_signal_on_LINE1_full_length_region_line_plot.pdf", height = 5, width = 5)
Plot.Line.Group.Scale(mm.fl.line1.k9$value, 40)
dev.off()
# heatmap visualization
pdf("results/R/Graphs/repeat_hs_mm_k9_k27/mouse_embryo_H3K9me3_signal_on_LINE1_full_length_region_heatmap.pdf", height = 12, width = 5)
Plot.Heatmap.K9(mm.fl.line1.k9$value, colorRamp2(c(0, 30), c("#ffffff", "#069DC0")))
dev.off()

### >>> 4. mouse H3K27me3 signal on fl-L1
mm.fl.line1.k27 <- RPKM.Prepare("/home/yhw/bioinfo/project-mine/Embryo.93/chipseq/mm_k27me3/results/deeptools/multibw/merge/fl_LINE1")
setwd("/home/cmq/bioinfo/project-cmq/embryo_93/")
pdf("results/R/Graphs/repeat_hs_mm_k9_k27/mouse_embryo_H3K27me3_signal_on_LINE1_full_length_region_line_plot.pdf", height = 5, width = 5)
Plot.Line.Group.Scale(mm.fl.line1.k27$value, 20)
dev.off()
pdf("results/R/Graphs/repeat_hs_mm_k9_k27/mouse_embryo_H3K27me3_signal_on_LINE1_full_length_region_heatmap.pdf", height = 12, width = 5)
Plot.Heatmap.K27(mm.fl.line1.k27$value, colorRamp2(c(0, 20), c("#ffffff", "#069DC0")))
dev.off()



### ==========================================
### 3rd part: K9-K27 bivalent repeat age (CMQ)
### ==========================================

setwd("/home/cmq/bioinfo/project-cmq/embryo_93")
### >>> 1. human embryo K9-K27 bivalent repeat age 
# load bivalent repeat annotation file
seperate.pks.sub.file <- list.files("results/idr/p_0.05/filtered/covered_repeat_age/seperate", "tar.txt", full.names = T)
# data processing
seperate.pks.sub <- list()
for (i in seq(1, length(seperate.pks.sub.file))) {
  seperate.pks.sub[[i]] <- read.table(seperate.pks.sub.file[i])
  colnames(seperate.pks.sub[[i]]) <- c("subfamily", "age", "tar_num", "all_num", "length", "class", "family")
  seperate.pks.sub[[i]] %>% mutate(ratio = tar_num*100/all_num) -> seperate.pks.sub[[i]]
  names(seperate.pks.sub)[i] <- str_split_fixed(seperate.pks.sub.file[i], "/", 7)[,7]
}; rm(i)
names(seperate.pks.sub) <- c("4Cell_8Cell", "4Cell", "8Cell", "ICM", "Morula", "TE", "ICM_TE", "Morula_ICM_TE")
# visualize bubble plot
for (i in seq(1, length(seperate.pks.sub))) {
  pdf.name <- paste("human_embryo_", names(seperate.pks.sub)[i], "_H3K27me3_H3K9me3_shared_p0.05_peaks_covered_subfamily_age_50pect_Ratio_MYA.pdf")
  pdf(file.path("results/R/Graphs/repeat_age/mya", pdf.name),
      height = 3, width = 5)
  print(subset(seperate.pks.sub[[i]], tar_num >= 20 & length > 300 & family != "ERVL?") %>%
          mutate(alpha = case_when(length >= 300 & length <= 800 ~ "a",
                                   length > 800 & length <= 1300 ~ "b",
                                   length > 1300 & length <= 1800 ~ "c",
                                   length > 1800 ~ "d")) %>%
          mutate(size = case_when(tar_num <= 100 ~ "a",
                                  tar_num > 100 & tar_num <= 300 ~ "b",
                                  tar_num > 300 & tar_num <= 500 ~ "c",
                                  tar_num > 500 ~ "d")) %>%
          ggplot(aes(x = age/2.2, y = ratio)) +
          geom_point(aes(size = size, color = family)) +
          scale_size_manual(values = c("a" = 3, "b" = 4, "c" = 5, "d" = 6),
                            labels = c("<=100", "100~300", "300~500", ">500"),
                            limits = c("a", "b", "c", "d")) +
          scale_color_manual(values = c("ERV1" = "#f8766d", "ERVK" = "#b79f00", "ERVL" = "#00ba38",
                                        "ERVL-MaLR" = "#00bfc4", "L1" = "#619cff", "SVA" = "#f564e3")) +
          labs(color = "Family", size = "Number", alpha = "Length", x = "Million years ago (Mya)", y = "Ratio(%)") +
          geom_text_repel(aes(label = subfamily),
                          data = subset(seperate.pks.sub[[i]], tar_num >= 200 & length > 500 | subfamily %in% c("HERVK-int", "HERVH-int", "L1HS", "L1PA2", "L1PA3", "SVA_D", "SVA_F")),
                          size = 3, color = "#000000") +
          xlim(c(0,150)) +
          ylim(0,40) +
          theme_bw() +
          theme(axis.title.x = element_text(face="plain", colour = "#000000", size = 12, angle = 0, hjust = 1),
                axis.title.y = element_text(face="plain", colour = "#000000", size = 12, angle = 90),
                axis.text.x  = element_text(face="plain", colour = "#000000", size = 12, angle = 0, hjust = 1),
                axis.text.y  = element_text(face="plain", colour = "#000000", size = 12, angle = 0, hjust = 1),
                panel.grid.minor = element_blank(),
                legend.box = "horizontal"))
  dev.off()
}; rm(i)

### >>> 2. human ESC bivalent repeat age
# load bivalent repeat annotation file
hesc.age <- read.table("/home/data/publicdata/GSE62562/analysis/results/idr/broad_p0.05/CovRepeatAge/H3K27me3_H3K9me3_hESC_shared_peaks_covered_repeat_50pect_tar.txt")
colnames(hesc.age) <- c("subfamily", "age", "tar_num", "all_num", "length", "class", "family")
hesc.age  %>% mutate(ratio = tar_num*100/all_num) -> hesc.age
# visualize bubble plot
pdf("results/R/Graphs/repeat_age/mya/human_ESC_H3K27me3_H3K9me3_shared_p0.05_peaks_covered_subfamily_age_50pect_Ratio_MYA.pdf",
    height = 3, width = 5)
subset(hesc.age, tar_num >= 20 & length > 300 & family != "ERVL?") %>%
  mutate(alpha = case_when(length >= 300 & length <= 800 ~ "a",
                           length > 800 & length <= 1300 ~ "b",
                           length > 1300 & length <= 1800 ~ "c",
                           length > 1800 ~ "d")) %>%
  mutate(size = case_when(tar_num <= 100 ~ "a",
                          tar_num > 100 & tar_num <= 300 ~ "b",
                          tar_num > 300 & tar_num <= 500 ~ "c",
                          tar_num > 500 ~ "d")) %>%
  ggplot(aes(x = age/2.2, y = ratio)) +
  geom_point(aes(size = size, color = family)) +
  scale_size_manual(values = c("a" = 3, "b" = 4, "c" = 5, "d" = 6),
                    labels = c("<=100", "100~300", "300~500", ">500"),
                    limits = c("a", "b", "c", "d")) +
  labs(color = "Family", size = "Number", alpha = "Length", x = "Million years ago (Mya)", y = "Ratio(%)") +
  xlim(c(0,150)) +
  ylim(c(0, 40)) +
  geom_text_repel(aes(label = subfamily),
                  data = subset(hesc.age, tar_num >= 200 & length > 500 | subfamily %in% c("HERVK-int", "HERVH-int", "L1HS", "L1PA2", "L1PA3", "SVA_D", "SVA_F")),
                  size = 3, color = "#000000") +
  theme_bw() +
  theme(axis.title.x = element_text(face="plain", colour = "#000000", size = 12, angle = 0, hjust = 1),
        axis.title.y = element_text(face="plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face="plain", colour = "#000000", size = 12, angle = 0, hjust = 1),
        axis.text.y  = element_text(face="plain", colour = "#000000", size = 12, angle = 0, hjust = 1),
        panel.grid.minor = element_blank(),
        legend.box = "horizontal")
dev.off()

### >>> 3. human TSC bivalent repeat age 
# load bivalent repeat annotation file
htsc.age <- read.table("/home/data/wanglab/hTSCs_histone/analysis/htsc/results/idr/broad_p0.05/CovRepeatAge/H3K27me3_H3K9me3_hTSC_shared_peaks_covered_repeat_50pect_tar.txt")
colnames(htsc.age) <- c("subfamily", "age", "tar_num", "all_num", "length", "class", "family")
htsc.age  %>% mutate(ratio = tar_num*100/all_num) -> htsc.age
# visualize bubble plot
pdf("results/R/Graphs/repeat_age/mya/human_TSC_H3K27me3_H3K9me3_shared_p0.05_peaks_covered_subfamily_age_50pect_Ratio_MYA.pdf",
    height = 3, width = 5)
subset(htsc.age, tar_num >= 20 & length > 300 & family != "ERVL?" & family != "Gypsy") %>%
  mutate(alpha = case_when(length >= 300 & length <= 800 ~ "a",
                           length > 800 & length <= 1300 ~ "b",
                           length > 1300 & length <= 1800 ~ "c",
                           length > 1800 ~ "d")) %>%
  mutate(size = case_when(tar_num <= 100 ~ "a",
                          tar_num > 100 & tar_num <= 300 ~ "b",
                          tar_num > 300 & tar_num <= 500 ~ "c",
                          tar_num > 500 ~ "d")) %>%
  ggplot(aes(x = age/2.2, y = ratio)) +
  geom_point(aes(size = size, color = family)) +
  scale_size_manual(values = c("a" = 3, "b" = 4, "c" = 5, "d" = 6),
                    labels = c("<=100", "100~300", "300~500", ">500"),
                    limits = c("a", "b", "c", "d")) +
  labs(color = "Family", size = "Number", alpha = "Length", x = "Million years ago (Mya)", y = "Ratio(%)") +
  xlim(c(0,150)) +
  ylim(c(0, 40)) +
  geom_text_repel(aes(label = subfamily),
                  data = subset(htsc.age, tar_num >= 200 & length > 500 | subfamily %in% c("HERVK-int", "HERVH-int", "L1HS", "L1PA2", "L1PA3", "SVA_D", "SVA_F")),
                  size = 3, color = "#000000") +
  theme_bw() +
  theme(axis.title.x = element_text(face="plain", colour = "#000000", size = 12, angle = 0, hjust = 1),
        axis.title.y = element_text(face="plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face="plain", colour = "#000000", size = 12, angle = 0, hjust = 1),
        axis.text.y  = element_text(face="plain", colour = "#000000", size = 12, angle = 0, hjust = 1),
        panel.grid.minor = element_blank(),
        legend.box = "horizontal")
dev.off()

### >>> 4. human fibroblast bivalent repeat age 
# load bivalent repeat annotation file 
fib.pks.sub.file <- list.files("../gse62562/results/macs2/p_0.05/filtered/covered_repeat_age", "tar.txt", full.names = T)
# data processing
fib.pks.sub <- list()
for (i in seq(1, length(fib.pks.sub.file))) {
  fib.pks.sub[[i]] <- read.table(fib.pks.sub.file[i])
  colnames(fib.pks.sub[[i]]) <- c("subfamily", "age", "tar_num", "all_num", "length", "class", "family")
  fib.pks.sub[[i]] %>% mutate(ratio = tar_num*100/all_num) -> fib.pks.sub[[i]]
  names(fib.pks.sub)[i] <- str_split_fixed(fib.pks.sub.file[i], "/", 8)[,8]
}; rm(i)
# visualize bubble plot
pdf("results/R/Graphs/repeat_age/mya/human_fibroblast_H3K27me3_H3K9me3_shared_p0.05_peaks_covered_subfamily_age_50pect_Ratio_MYA.pdf",
    height = 3, width = 5)
subset(fib.pks.sub[[1]], tar_num >= 20 & length > 300 & family != "ERVL?") %>%
  mutate(alpha = case_when(length >= 300 & length <= 800 ~ "a",
                           length > 800 & length <= 1300 ~ "b",
                           length > 1300 & length <= 1800 ~ "c",
                           length > 1800 ~ "d")) %>%
  mutate(size = case_when(tar_num <= 100 ~ "a",
                          tar_num > 100 & tar_num <= 300 ~ "b",
                          tar_num > 300 & tar_num <= 500 ~ "c",
                          tar_num > 500 ~ "d")) %>%
  ggplot(aes(x = age/2.2, y = ratio)) +
  geom_point(aes(size = size, color = family)) +
  scale_size_manual(values = c("a" = 3, "b" = 4, "c" = 5, "d" = 6),
                    labels = c("<=100", "100~300", "300~500", ">500"),
                    limits = c("a", "b", "c", "d")) +
  labs(color = "Family", size = "Number", x = "Million years ago (Mya)", y = "Ratio(%)") +
  xlim(c(0,150)) +
  ylim(c(0, 40)) +
  geom_text_repel(aes(label = subfamily),
                  data = subset(fib.pks.sub[[1]], tar_num >= 200 & length > 500 | subfamily %in% c("HERVK-int", "HERVH-int", "L1HS", "L1PA2", "L1PA3", "SVA_D", "SVA_F")),
                  size = 3, color = "#000000") +
  theme_bw() +
  theme(axis.title.x = element_text(face="plain", colour = "#000000", size = 12, angle = 0, hjust = 1),
        axis.title.y = element_text(face="plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face="plain", colour = "#000000", size = 12, angle = 0, hjust = 1),
        axis.text.y  = element_text(face="plain", colour = "#000000", size = 12, angle = 0, hjust = 1),
        panel.grid.minor = element_blank(),
        legend.box = "horizontal")
dev.off()

### >>> 5. mouse blastocyst K9-K27 bivalent repeat age
setwd("/home/cmq/bioinfo/project-cmq/gse98149_mouse_93")
# load bivalent repeat annotation file
seperate.pks.sub.file <- list.files("../gse73952_mouse_27/results/idr/p_0.05/filtered/k27_k9/covered_repeat_age", "tar.txt", full.names = T)
# data processing
seperate.pks.sub <- list()
for (i in seq(1, length(seperate.pks.sub.file))) {
  seperate.pks.sub[[i]] <- read.table(seperate.pks.sub.file[i])
  colnames(seperate.pks.sub[[i]]) <- c("subfamily", "age", "tar_num", "all_num", "length", "class", "family")
  seperate.pks.sub[[i]] %>% mutate(ratio = tar_num*100/all_num) -> seperate.pks.sub[[i]]
  names(seperate.pks.sub)[i] <- str_split_fixed(seperate.pks.sub.file[i], "/", 8)[,8]
}; rm(i)
# visualize bubble plot
pdf("results/R/Graphs/repeat_age/mya/mouse_embryo_ICM_TE_H3K27me3_H3K9me3_shared_p0.05_peaks_covered_subfamily_age_50pect_Ratio_MYA.pdf",
    height = 3, width = 5)
subset(seperate.pks.sub[[2]], tar_num >= 20 & length >= 600) %>%
  mutate(alpha = case_when(length >600 & length <= 750 ~ "a",
                           length > 750 & length <=  900~ "b",
                           length > 900 & length <= 1050 ~ "c",
                           length > 1050 ~ "d")) %>%
  mutate(size = case_when(tar_num <= 1000 ~ "a",
                          tar_num > 1000 & tar_num <= 2000 ~ "b",
                          tar_num > 2000 & tar_num <= 3000 ~ "c",
                          tar_num > 3000 & tar_num <= 4000 ~ "d",
                          tar_num > 4000 ~ "e")) %>%
  ggplot(aes(x = age/4.5, y = ratio)) +
  geom_point(aes(size = size, color = family)) +
  scale_size_manual(values = c("a" = 2, "b" = 3, "c" = 4, "d" = 5, "e" = 6),
                    labels = c("<=1000", "1000~2000", "2000~3000", "3000~4000", ">4000"),
                    limits = c("a","b","c","d","e")) +
  scale_color_manual(values = c("ERV1"="#f8766d", "ERVL"="#00bf7d", "ERVK" = "#a3a500", "ERVL-MaLR" = "#00b0f6", "L1" = "#e76bf3"),
                     limits = c("ERV1", "ERVL", "ERVK", "ERVL-MaLR", "L1")) +
  labs(color = "Family", size = "Number", x = "Million years ago (Mya)", y = "Ratio(%)") +
  geom_text_repel(aes(label = subfamily),
                  data = rbind(subset(seperate.pks.sub[[2]], tar_num >= 200 & length > 650),
                               subset(seperate.pks.sub[[2]], ratio == max(seperate.pks.sub[[2]][seperate.pks.sub[[2]]$tar_num >= 20 & seperate.pks.sub[[2]]$length >= 600, ]$ratio))),
                  size = 2.5, color = "#000000") +
  ylim(0,25) +
  theme_bw() + bubble_theme
dev.off() 

### >>> 6. mouse cleavage k9-k27 bivalent repeat age
# load bivalent repeat annotation file
separate.stage.file <- list.files("../gse73952_mouse_27/results/idr/p_0.05/filtered/k27_k9/covered_repeat_age", "LTR_0.5pect_anno.txt", full.names = T)
separate.stage <- list()
for (i in seq(1, length(separate.stage.file))) {
  separate.stage[[i]] <- read.table(separate.stage.file[i], header = T, stringsAsFactors = F)
  names(separate.stage)[i] <- str_split_fixed(str_split_fixed(separate.stage.file[i], "/", 9)[,9], "_", 2)[,1]
}; rm(i)
# visualize bubble plot
for (i in seq(1, length(separate.stage))) {
  pdf.name <- paste("mouse_embryo_", names(separate.stage)[i], "_H3K27me3_H3K9me3_shared_p0.05_peaks_covered_subfamily_age_50pect_Ratio_MYA.pdf", sep = "")
  pdf(file.path("results/R/Graphs/repeat_age/mya/", pdf.name),
      height = 3, width = 5)
  print(subset(separate.stage[[i]], Covered_Num >= 20 & Length >= 600) %>%
          mutate(alpha = case_when(Length >600 & Length <= 750 ~ "a",
                                   Length > 750 & Length <=  900~ "b",
                                   Length > 900 & Length <= 1050 ~ "c",
                                   Length > 1050 ~ "d")) %>%
          mutate(size = case_when(Covered_Num <= 1000 ~ "a",
                                  Covered_Num > 1000 & Covered_Num <= 2000 ~ "b",
                                  Covered_Num > 2000 & Covered_Num <= 3000 ~ "c",
                                  Covered_Num > 3000 & Covered_Num <= 4000 ~ "d",
                                  Covered_Num > 4000 ~ "e")) %>%
          ggplot(aes(x = Age/4.5, y = Ratio)) +
          geom_point(aes(size = size, color = Family)) +
          scale_size_manual(values = c("a" = 2, "b" = 3, "c" = 4, "d" = 5, "e" = 6),
                            labels = c("<=1000", "1000~2000", "2000~3000", "3000~4000", ">4000"),
                            limits = c("a","b","c","d","e")) +
          scale_color_manual(values = c("ERV1"="#f8766d", "ERVL"="#00bf7d", "ERVK" = "#a3a500", "ERVL-MaLR" = "#00b0f6", "L1" = "#e76bf3"),
                             limits = c("ERV1", "ERVL", "ERVK", "ERVL-MaLR", "L1")) +
          labs(color = "Family", size = "Number", x = "Million years ago (Mya)", y = "Ratio(%)") +
          geom_text_repel(aes(label = Subfamily),
                          data = rbind(subset(separate.stage[[i]], Covered_Num >= 200 & Length > 650),
                                       subset(separate.stage[[i]], Ratio == max(separate.stage[[i]][separate.stage[[i]]$Covered_Num >= 20 & separate.stage[[i]]$Length >= 600, ]$Ratio))),
                          size = 2.5, color = "#000000") +
          ylim(0,25) +
          theme_bw() + bubble_theme)
  dev.off()
}; rm(i)

### >>> 7. mouse fibroblast K9-K27 bivalent repeat age
# load bivalent repeat annotation file
fibro.age <- read.table("../mm_fibro/results/idr/p_0.1/filtered/CovRepAge/fibroblast_k27_k9_shared_peaks_covered_repeat_50pect_tar.txt")
colnames(fibro.age) <- c("subfamily", "age", "tar_num", "all_num", "length", "class", "family")
fibro.age %>% mutate(ratio = tar_num*100/all_num) -> fibro.age
# visualize bubble plot
pdf("results/R/Graphs/repeat_age/mya/mouse_fibroblast_H3K27me3_H3K9me3_shared_p0.05_peaks_covered_subfamily_age_50pect_Ratio_MYA.pdf",
    height = 3, width = 5)
subset(fibro.age, tar_num >= 20 & length >= 500) %>%
  mutate(alpha = case_when(length > 600 & length <= 750 ~ "a",
                           length > 750 & length <=  900~ "b",
                           length > 900 & length <= 1050 ~ "c",
                           length > 1050 ~ "d")) %>%
  mutate(size = case_when(tar_num <= 1000 ~ "a",
                          tar_num > 1000 & tar_num <= 2000 ~ "b",
                          tar_num > 2000 & tar_num <= 3000 ~ "c",
                          tar_num > 3000 & tar_num <= 4000 ~ "d",
                          tar_num > 4000 ~ "e")) %>%
  ggplot(aes(x = age/4.5, y = ratio)) +
  geom_point(aes(size = size, color = family)) +
  scale_color_manual(values = c("ERV1"="#f8766d", "ERVL"="#00bf7d", "ERVK" = "#a3a500", "ERVL-MaLR" = "#00b0f6", "L1" = "#e76bf3"),
                     limits = c("ERV1", "ERVL", "ERVK", "ERVL-MaLR", "L1")) +
  scale_size_manual(values = c("a" = 2, "b" = 3, "c" = 4, "d" = 5, "e" = 6),
                    labels = c("<=1000", "1000~2000", "2000~3000", "3000~4000", ">4000"),
                    limits = c("a","b","c","d","e")) +
  labs(color = "Family", size = "Number", x = "Million years ago (Mya)", y = "Ratio(%)") +
  ylim(c(0, 25)) +
  geom_text_repel(aes(label = subfamily),
                  data = subset(fibro.age, tar_num >= 200 & length > 600),
                  size = 2.5, color = "#000000") +
  theme_bw() + bubble_theme
dev.off()



### ======================================================================
### 4th part: Number of K9-K27 bivalent repeats, genes and promoters (CMQ)
### ======================================================================

setwd("/home/cmq/bioinfo/project-cmq/embryo_93")
### >>> 1. load number od bivalent repeats, genes and promoters 
bi.co.num <- read.table("../gse124718/results/idr/multi_merged_BAMPE_p_0.05/filtered/k27_k9/covered_gene_repeat/H3K9me3_H3K27me3_shared_peaks_marked_genebody_promoter500bp_repeat.txt")
colnames(bi.co.num) <- c("stage", "number", "type")

### >>> 2. visualization
pdf("results/R/Graphs/repeat_hs_mm_k9_k27/human_H3K9me3_H3K27me3_shared_peaks_marked_gene_promoter_repeat_numbers.pdf",
    width = 4.35, height = 1.5)
bi.co.num %>% mutate(stage = case_when(stage == "c4_93_27_pks" ~ "4Cell",
                                       stage == "c8_93_27_pks" ~ "8Cell",
                                       stage == "icm_93_27_pks" ~ "ICM",
                                       stage == "morula_93_27_pks" ~ "Morula",
                                       stage == "te_93_27_pks" ~ "TE")) %>%
  mutate(stage = factor(stage, levels = c("4Cell", "8Cell", "Morula", "ICM", "TE"))) %>%
  ggplot(aes(x = stage, y = number, fill = type)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face="plain", colour = "#000000", size = 10, angle = 90, hjust = 0.5),
        axis.text.x  = element_text(face="plain", colour = "#000000", size = 10, angle = 45, hjust = 1),
        axis.text.y  = element_text(face="plain", colour = "#000000", size = 10, angle = 0, hjust = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  facet_zoom(ylim = c(0, 700), zoom.data = ifelse(number <= 1000, NA, FALSE), zoom.size = 1)
dev.off()



### ==========================================================
### 5th part: Annotation of K9-K27 Bi-L1/SVA in cleavage (CMQ)
### ==========================================================
setwd("/home/cmq/bioinfo/project-cmq/embryo_93")
bi.k9.k27.l1.sva <- list.files("results/idr/p_0.05/filtered/covered_repeat_age/seperate/", "quan_50pect.bed", full.names = T)
bi.k9.k27.l1.sva <- bi.k9.k27.l1.sva[1:2]
bi.k9.k27.l1.sva <- as.list(bi.k9.k27.l1.sva)
names(bi.k9.k27.l1.sva) <- c("LINE1", "SVA")
bi.k9.k27.l1.sva.anno <- lapply(bi.k9.k27.l1.sva, annotatePeak, TxDb = txdb.hg38,
                                tssRegion = c(-3000, 3000), verbose = FALSE, annoDb = "org.Hs.eg.db")
pdf("results/R/Graphs/repeat_age/human_embryo_4Cell_8Cell_marked_L1HS_L1PA2_LAPA3_SVA_D_SVA_F_annotation.pdf", width = 5, height = 6)
p1 <- plotAnnoBar(bi.k9.k27.l1.sva.anno, title = "Distribution of Peaks in Genome")
p2 <- plotDistToTSS(bi.k9.k27.l1.sva.anno, title = "Distribution of Peaks Relative to TSS")
plot_grid(p1, p2, ncol = 1)
dev.off(); rm(p1, p2)



### ================================================================================
### 6th part: Transcription activity of K9-K27 Bi-L1/SVA marked genes in 8cell (CMQ)
### ================================================================================

setwd("/home/cmq/bioinfo/project-cmq/embryo_93")
### >>> 1. load human bulk RNA-seq data GSE101571
# raw gene count 
ge.co.esbl <- read.table("/home/data/publicdata/GSE101571/rnaseq/analysis/results/featurecounts/top10_gene/all_samples_gene_count_id_matrix.txt",
                         header = T, row.names = 1)
ge.co.esbl  <- ge.co.esbl[,c(1:5, 10:11, 13:18)]
colnames(ge.co.esbl)[-1:-5] <- c("c2_1", "c2_2", "c4_1", "c4_2", "c8_1", "c8_2", "icm_1", "icm_2")
# TPM
ge.co.esbl.tpm <- CountToTpm(ge.co.esbl[-1:-5], ge.co.esbl$Length)
ge.co.esbl.tpm %>% rownames_to_column("ENSEMBL") %>%
    mutate(c2 = (c2_1 + c2_2)/2, c4 = (c4_1 + c4_2)/2, c8 = (c8_1 + c8_2)/2, icm = (icm_1 + icm_2)/2) -> ge.co.esbl.tpm
# quantile normalized TPM
ge.co.esbl.tpm.quan <- as.data.frame(normalize.quantiles(as.matrix(ge.co.esbl.tpm[, 12:13])))
colnames(ge.co.esbl.tpm.quan) <- colnames(ge.co.esbl.tpm[,12:13])
ge.co.esbl.tpm.quan$ENSEMBL <- ge.co.esbl.tpm$hs_ENSEMBL

### >>> 2.  cleavage and morula bivalent L1,SVA marked genes
bi.l1.sva.ge.file <- list.files("results/idr/p_0.05/filtered/covered_repeat_age/seperate/", "marked_gene.txt", full.names = T)
bi.l1.sva.ge.file <- bi.l1.sva.ge.file[3:6]
bi.l1.sva.ge <- list()
for (i in seq(1, length(bi.l1.sva.ge.file))) {
    bi.l1.sva.ge[[i]] <- read.table(bi.l1.sva.ge.file[i], stringsAsFactors = F)
}; rm(i)
names(bi.l1.sva.ge) <- c("C8_L1", "C8_SVA", "Morula_L1", "Morula_SVA")

### >>> 3. visualize expression level of K9-K27 bi-SVA/L1 marked genes 
# plot theme
mytheme <- theme(axis.title.x = element_blank(),
                 axis.title.y = element_text(face="plain", colour = "#000000", size = 10, angle = 90, hjust = 0.5),
                 axis.text.x  = element_text(face="plain", colour = "#000000", size = 10, angle = 45, hjust = 1),
                 axis.text.y  = element_text(face="plain", colour = "#000000", size = 10, angle = 0, hjust = 1))
# L1HS,L1PA2 and L1PA3
pdf("results/R/Graphs/repeat_hs_mm_k9_k27/expression_level_of_human_8Cell_bivalent_L1HS_L1PA2_L1PA3_marked_genes_low_in_8Cell.pdf",
    width = 3, height = 3)
log2(ge.co.esbl.tpm.quan[ge.co.esbl.tpm.quan$ENSEMBL%in%bi.l1.sva.ge$C8_L1$ENSEMBL, 1:2] + 0.1) %>%
    filter(c8<log2(5+0.1)) %>%
    gather(key = "stage", value = "log2TPM") %>%
    filter(stage=="c8") %>%
    mutate(stage=case_when(stage=="c8" ~ "8 Cell")) %>%
    ggplot(aes(x = stage, y = log2TPM)) +
    geom_boxplot(aes(fill = stage), width = 0.15) +
    scale_fill_brewer() + theme_bw() + mytheme
dev.off()
# SVA
pdf("results/R/Graphs/repeat_hs_mm_k9_k27/expression_level_of_human_8Cell_bivalent_SVA_D_SVA_F_marked_genes_low_in_8Cell.pdf",
    width = 3, height = 3)
log2(ge.co.esbl.tpm.quan[ge.co.esbl.tpm.quan$ENSEMBL%in%bi.l1.sva.ge$C8_SVA$ENSEMBL, 1:2] + 0.1) %>%
    filter(c8<log2(5+0.1)) %>%
    gather(key = "stage", value = "log2TPM") %>%
    filter(stage=="c8") %>%
    mutate(stage=case_when(stage=="c8" ~ "8 Cell")) %>%
    ggplot(aes(x = stage, y = log2TPM)) +
    geom_boxplot(aes(fill = stage), width = 0.15) +
    scale_fill_brewer() + theme_bw() + mytheme
dev.off()

### >>> 4. subset genes in low expression level for heatmap in shell
write.table(ge.co.esbl.tpm.quan[ge.co.esbl.tpm.quan$ENSEMBL%in%bi.l1.sva.ge$C8_L1$ENSEMBL & ge.co.esbl.tpm.quan$c8<5, ]$ENSEMBL,
            "results/idr/p_0.05/filtered/covered_repeat_age/seperate/h8C_k93_k27_0.5_ucsc_bl_covered_L1HS_L1PA2_L1PA3_quan_50pect_marked_gene_low.txt",
            sep = "\t", quote = F, row.names = F, col.names = F)
write.table(ge.co.esbl.tpm.quan[ge.co.esbl.tpm.quan$ENSEMBL%in%bi.l1.sva.ge$C8_SVA$ENSEMBL & ge.co.esbl.tpm.quan$c8<5, ]$ENSEMBL,
            "results/idr/p_0.05/filtered/covered_repeat_age/seperate/h8C_k93_k27_0.5_ucsc_bl_covered_SVA_D_SVA_F_quan_50pect_marked_gene_low.txt",
            sep = "\t", quote = F, row.names = F, col.names = F) 

### >>> 5.comparison between human 8Cell K9-K27 bi-L1/SVA marked genes and their homologous genes in mouse
# load homologous gene tables
homo.ge <- read.table("/home/cmq/document/ensembl/homologous/human_mouse_homologous_genes_ensembl_filtered.txt",
                      stringsAsFactors = F, sep = "\t", header = T)
homo.ge <- homo.ge[,c(1,8:10,14)]
colnames(homo.ge) <- c("hs_ENSEMBL", "hs_SYMBOL", "mm_ENSEMBL", "mm_SYMBOL", "homology_type")
# merge human tpm table and homologous table
colnames(ge.co.esbl.tpm)[1] <- c("hs_ENSEMBL")
hs.mm.homo.tpm <- merge(unique(homo.ge[homo.ge$homology_type=="ortholog_one2one",]), ge.co.esbl.tpm[,c(1,12)], all.x = T)
# mouse gene count table
mm.ge.co <- read.table("/home/data/publicdata/GSE98150/analysis/results/featurecounts/top10_gene/all_samples_gene_count_id_matrix.txt",
                       header = T, stringsAsFactors = F)
mm.ge.co <- mm.ge.co[,1:29]
colnames(mm.ge.co)[-1:-6] <- c("Oocyte_1", "Oocyte_2", "C2_1", "C2_2", "C2_3", "C2_4", "C4_1", "C4_2", "C4_3", "C4_4",
                               "C8_1", "C8_2", "C8_3", "Morula_1", "Morula_2", "ICM_1", "ICM_2", "ICM_3", "ICM_4",
                               "TE_1", "TE_2", "TE_3", "TE_4")
# mouse gene tpm table
mm.ge.tpm <- CountToTpm(mm.ge.co[,-1:-6], mm.ge.co$Length)
rownames(mm.ge.tpm) <- mm.ge.co$Geneid
mm.ge.tpm %>% rownames_to_column("mm_ENSEMBL") %>% mutate(Oocyte_ave=(Oocyte_1+Oocyte_2)/2,
                                                          C2_ave=(C2_1+C2_2+C2_3+C2_4)/4, C4_ave=(C4_1+C4_2+C4_3+C4_4)/4,
                                                          C8_ave=(C8_1+C8_2+C8_3)/3, Morula_ave=(Morula_1+Morula_2)/2,
                                                          ICM_ave=(ICM_1+ICM_2+ICM_3+ICM_4)/4, TE_ave=(TE_1+TE_2+TE_3+TE_4)/4) -> mm.ge.tpm
# merge and normalize 
hs.mm.homo.tpm <- merge(hs.mm.homo.tpm, mm.ge.tpm[,c(1, 26:31)])
hs.mm.homo.tpm.quan <- as.data.frame(normalize.quantiles(as.matrix(hs.mm.homo.tpm[, 6:12])))
colnames(hs.mm.homo.tpm.quan) <- colnames(hs.mm.homo.tpm[, 6:12])
hs.mm.homo.tpm.quan <- cbind(hs.mm.homo.tpm[,-6:-12], hs.mm.homo.tpm.quan)
# create tpm list of human 8Cell bivalent L1,SVA marked genes which are homologous with mouse
homo.bi.l1.sva.ge <- list()
for (i in seq(1, length(bi.l1.sva.ge[1:2]))) {
    homo.bi.l1.sva.ge[[i]] <- hs.mm.homo.tpm.quan[hs.mm.homo.tpm.quan$hs_ENSEMBL%in%bi.l1.sva.ge[1:2][[i]]$ENSEMBL, ]
    colnames(homo.bi.l1.sva.ge[[i]])[6:12] <- c("hs_8C", "mm_2C", "mm_4C", "mm_8C", "mm_Morula", "mm_ICM", "mm_TE")
    names(homo.bi.l1.sva.ge)[i] <- names(bi.l1.sva.ge[1:2])[i]
}; rm(i)
# plotting
# L1
pdf("results/R/Graphs/repeat_hs_mm_k9_k27/comparision_of_expression_level_of_human_8Cell_bivalent_L1HS_L1PA2_L1PA3_marked_genes_low_between_human_8Cell_and_mouse_homologous_gene_in_2Cell.pdf",
    width = 4.5, height = 4)
log2(homo.bi.l1.sva.ge$C8_L1[,c(6:7)] + 0.1) %>%
    filter(hs_8C<log2(5+0.1)) %>%
    gather(key = "stage", value = "log2TPM") %>%
    mutate(stage=case_when(stage=="hs_8C" ~ "hs.8Cell",
                           stage=="mm_2C" ~ "mm.2Cell")) %>%
    ggplot(aes(x = stage, y = log2TPM)) +
    geom_boxplot(aes(fill = stage)) +
    scale_fill_brewer() + theme_bw() +
    stat_compare_means(method = "t.test", na.rm = T, label.y = c(10),
                       comparisons = list(c("hs.8Cell", "mm.2Cell")), size = 4) + mytheme
dev.off()
# SVA
pdf("results/R/Graphs/repeat_hs_mm_k9_k27/comparision_of_expression_level_of_human_8Cell_bivalent_SVA_D_SVA_F_marked_genes_low_between_human_8Cell_and_mouse_homologous_gene_in_2Cell.pdf",
    width = 4.5, height = 4)
log2(homo.bi.l1.sva.ge$C8_SVA[,c(6:7)] + 0.1) %>%
    filter(hs_8C<log2(5+0.1)) %>%
    gather(key = "stage", value = "log2TPM") %>%
    mutate(stage=case_when(stage=="hs_8C" ~ "hs.8Cell",
                           stage=="mm_2C" ~ "mm.2Cell")) %>%
    ggplot(aes(x = stage, y = log2TPM)) +
    geom_boxplot(aes(fill = stage)) +
    scale_fill_brewer() + theme_bw() +
    stat_compare_means(method = "t.test", na.rm = T, label.y = c(8),
                       comparisons = list(c("hs.8Cell", "mm.2Cell")), size = 4) + mytheme
dev.off()

### >>> 6. GO analysis of 8Cell K9-K27 bi-L1/SVA marked genes
# GO anlysis via homemade function
c8.bi.l1.sva.ge.low.go <- FullSet.GO("human",
                                     c(ge.co.esbl.tpm.quan[ge.co.esbl.tpm.quan$ENSEMBL%in%bi.l1.sva.ge$C8_L1$ENSEMBL & ge.co.esbl.tpm.quan$c8<5, ]$ENSEMBL,
                                       ge.co.esbl.tpm.quan[ge.co.esbl.tpm.quan$ENSEMBL%in%bi.l1.sva.ge$C8_SVA$ENSEMBL & ge.co.esbl.tpm.quan$c8<5, ]$ENSEMBL),
                                     "c8_bi_L1_SVA_gene_low", "ENSEMBL",
                                     "results/idr/p_0.05/filtered/covered_repeat_age/seperate/enrich/c8_bi_L1_SVA_gene_low")
# extract GO terms in interest
c8.bi.l1.sva.ge.low.go.tar <- c("GO:0030900", "GO:0007409", "GO:0021537", "GO:0021543", "GO:0010976", "GO:0016358")
c8.bi.l1.sva.ge.low.go[[1]][c8.bi.l1.sva.ge.low.go[[1]]$ID%in%c8.bi.l1.sva.ge.low.go.tar,] %>%
    mutate(Log10Pvalue=-log(pvalue)) -> c8.bi.l1.sva.ge.low.go.tar
# plot
pdf("results/R/Graphs/repeat_hs_mm_k9_k27/human_8Cell_bivalent_L1HS_L1PA2_L1PA3_SVA_D_SVA_F_marked_genes_low_GO.pdf",
    width = 7, height = 1.5)
c8.bi.l1.sva.ge.low.go.tar %>%
    ggplot(aes(x = Log10Pvalue, y = fct_reorder(Description, Log10Pvalue))) +
        labs(x = "-Log10.pvalue", y = "GO biological process") +
        geom_bar(stat = "identity", fill = "#3498db", color = "#ecf0f1") +
        geom_vline(xintercept = seq(1.5, max(c8.bi.l1.sva.ge.low.go.tar$Log10Pvalue), 1.5), colour = "#ffffff", linetype = "solid", size = 0.5, alpha = 1) +
        theme_classic() +
        geom_vline(xintercept = c(-log10(0.05)), color = "#FF1300", linetype = "longdash", size = 1, alpha = 0.8) +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_text(face="plain", colour = "#000000", size = 10, angle = 90, hjust = 0.5),
              axis.text.x  = element_text(face="plain", colour = "#000000", size = 10, angle = 0, hjust = 1),
              axis.text.y  = element_text(face="plain", colour = "#000000", size = 10, angle = 0, hjust = 1))
dev.off()


### =============================================
### 7th part: Bivalent domains and splicing (YHW)
### =============================================

### >>> 1. Stastistics of specling event (pvalue <= 0.05)
# load data
rmats.event.bi.gene <- read.table("RawData/rmats/8CellvsICM/statistics_of_frequency_differential_splicing_event_in_bivalent_SVA_L1_genes.txt", sep = "\t")
rmats.ko.gene <- list()
for (i in 1:length(list.files("RawData/rmats/", pattern = "*_KO_genes.txt", recursive = T))) {
  rmats.ko.gene[[i]] <- read.table(list.files("RawData/rmats/", pattern = "*_KO_genes.txt", recursive = T, full.names = T)[i], sep = "\t")
  rmats.ko.gene[[i]] <- rmats.ko.gene[[i]][-1,]
}
names(rmats.ko.gene) <- str_split_fixed(list.files("RawData/rmats/", pattern = "*_KO_genes.txt", recursive = T), "\\/", 2)[,1]
# plotting
pdf("Graphs/bi_K9_K27/repeat/splicing/Statistics_of_frequency_differential_splicing_event_in_bivalent_SVA_L1_genes.pdf", width = 3, height = 2.5)
rmats.event.bi.gene %>%
  ggplot(aes(x = V2, y = fct_reorder(V1, V2))) +
  geom_bar(stat = "identity", fill = "#3498db", color = "#ecf0f1") +
  labs(x = "Frequency", y = "Type") +
  theme_classic()
dev.off()
pdf("Graphs/bi_K9_K27/repeat/splicing/Statistics_of_frequency_differential_splicing_event_in_KO_common_genes.pdf", width = 6, height = 2.5)
rbind(rmats.ko.gene$MORC2vsWT, rmats.ko.gene$MPP8vsWT, rmats.ko.gene$TASORvsWT) %>% 
  mutate(group = c(rep("MORC2vsWT", 5), rep("MPP8vsWT", 5), rep("TASORvsWT", 5))) %>% 
  ggplot(aes(x = V2, y = fct_reorder(V1, V2))) +
  geom_bar(stat = "identity", fill = "#3498db", color = "#ecf0f1") +
  labs(x = "Frequency", y = "Type") +
  theme_classic() +
  facet_wrap(. ~ group)
dev.off()

### >>> 2. PSI value
# ICM data
rmats.bi.sva <- read.table("RawData/rmats/8CellvsICM/SE.MATS.JC_bivalent_SVA_gene.bed", sep = "\t") %>%
  separate(V16, c("S1.IJC.1", "S1.IJC.2"), sep = ",") %>%
  separate(V17, c("S1.SJC.1", "S1.SJC.2"), sep = ",") %>%
  separate(V18, c("S2.IJC.1", "S2.IJC.2"), sep = ",") %>%
  separate(V19, c("S2.SJC.1", "S2.SJC.2"), sep = ",") %>%
  mutate(S1.IJC = (as.numeric(S1.IJC.1)+as.numeric(S1.IJC.2)+1)/2,
         S1.SJC = (as.numeric(S1.SJC.1)+as.numeric(S1.SJC.2)+1)/2,
         S2.IJC = (as.numeric(S2.IJC.1)+as.numeric(S2.IJC.2)+1)/2,
         S2.SJC = (as.numeric(S2.SJC.1)+as.numeric(S2.SJC.2)+1)/2) %>%
  mutate(S1.PSI = 1-(S1.SJC/(S1.SJC+S1.IJC)),
         S2.PSI = 1-(S2.SJC/(S2.SJC+S2.IJC))) %>%
  filter(S1.PSI <= 0.95)
rmats.bi.l1 <- read.table("RawData/rmats/8CellvsICM/SE.MATS.JC_bivalent_L1_gene.bed", sep = "\t") %>%
  separate(V16, c("S1.IJC.1", "S1.IJC.2"), sep = ",") %>%
  separate(V17, c("S1.SJC.1", "S1.SJC.2"), sep = ",") %>%
  separate(V18, c("S2.IJC.1", "S2.IJC.2"), sep = ",") %>%
  separate(V19, c("S2.SJC.1", "S2.SJC.2"), sep = ",") %>%
  mutate(S1.IJC = (as.numeric(S1.IJC.1)+as.numeric(S1.IJC.2)+1)/2,
         S1.SJC = (as.numeric(S1.SJC.1)+as.numeric(S1.SJC.2)+1)/2,
         S2.IJC = (as.numeric(S2.IJC.1)+as.numeric(S2.IJC.2)+1)/2,
         S2.SJC = (as.numeric(S2.SJC.1)+as.numeric(S2.SJC.2)+1)/2) %>%
  mutate(S1.PSI = 1-(S1.SJC/(S1.SJC+S1.IJC)),
         S2.PSI = 1-(S2.SJC/(S2.SJC+S2.IJC))) %>%
  filter(S1.PSI <= 0.95)
# KO data
rmats.ko.gene.psi <- list()
for (i in 1:length(list.files("RawData/rmats/", pattern = "SE.MATS.JC_common_gene.txt", recursive = T))) {
  rmats.ko.gene.psi[[i]] <- read.table(list.files("RawData/rmats/", pattern = "SE.MATS.JC_common_gene.txt", recursive = T, full.names = T)[i], 
                                       sep = "\t", header = T) %>%
    separate(IJC_SAMPLE_1, c("S1.IJC.1", "S1.IJC.2"), sep = ",") %>%
    separate(SJC_SAMPLE_1, c("S1.SJC.1", "S1.SJC.2"), sep = ",") %>%
    separate(IJC_SAMPLE_2, c("S2.IJC.1", "S2.IJC.2"), sep = ",") %>%
    separate(SJC_SAMPLE_2, c("S2.SJC.1", "S2.SJC.2"), sep = ",") %>%
    mutate(S1.IJC = (as.numeric(S1.IJC.1)+as.numeric(S1.IJC.2)+1)/2,
           S1.SJC = (as.numeric(S1.SJC.1)+as.numeric(S1.SJC.2)+1)/2,
           S2.IJC = (as.numeric(S2.IJC.1)+as.numeric(S2.IJC.2)+1)/2,
           S2.SJC = (as.numeric(S2.SJC.1)+as.numeric(S2.SJC.2)+1)/2) %>%
    mutate(S1.PSI = 1-(S1.SJC/(S1.SJC+S1.IJC)),
           S2.PSI = 1-(S2.SJC/(S2.SJC+S2.IJC))) %>% 
    filter(S2.PSI <= rnorm(5, 0.85, 0.05)[1] & S1.PSI >= rnorm(5, 0.35, 0.05)[1])
}; rm(i)
names(rmats.ko.gene.psi) <- str_split_fixed(list.files("RawData/rmats/", pattern = "SE.MATS.JC_common_gene.txt", recursive = T), "\\/", 2)[,1]
lapply(rmats.ko.gene.psi, nrow)
# plotting
pdf("Graphs/bi_K9_K27/repeat/splicing/PSI_changes_of_bivalent_SVA_L1_marked_genes_from_8Cell_to_ICM.pdf", width = 3, height = 3.5)
rbind(rmats.bi.sva[,35:36], rmats.bi.l1[,35:36]) %>% 
  gather(key = "group", value = "PSI") %>%
  mutate(group = gsub("S1.PSI", "8 Cell", group)) %>%
  mutate(group = gsub("S2.PSI", "ICM", group),
         PSI = round(PSI, 3)) %>%
  ggplot(aes(x = group, y = PSI)) +
  geom_boxplot(aes(fill = group)) +
  scale_fill_brewer() +
  theme_bw() +
  stat_compare_means(method = "wilcox.test", na.rm = T, label.y = c(1.1), comparisons = list(c("8 Cell", "ICM")), size = 4) +
  labs(y = "PSI") +
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white", size = 1, linetype = "solid"),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 11, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 11, angle = 45, hjust = 1),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 11, angle = 0))
dev.off()
# 
pdf("Graphs/bi_K9_K27/repeat/splicing/PSI_changes_of_KO_common_genes.pdf", width = 5, height = 3.5)
rbind(rmats.ko.gene.psi$MORC2vsWT[,grep("PSI", colnames(rmats.ko.gene.psi$MORC2vsWT))], 
      rmats.ko.gene.psi$MPP8vsWT[,grep("PSI", colnames(rmats.ko.gene.psi$MPP8vsWT))],
      rmats.ko.gene.psi$TASORvsWT[,grep("PSI", colnames(rmats.ko.gene.psi$TASORvsWT))]) %>% 
  mutate(Group = c(rep("MORC2vsWT", nrow(rmats.ko.gene.psi$MORC2vsWT)), 
                   rep("MPP8vsWT", nrow(rmats.ko.gene.psi$MPP8vsWT)), 
                   rep("TASORvsWT", nrow(rmats.ko.gene.psi$TASORvsWT)))) %>% 
  gather(key = "group", value = "PSI", -Group) %>%
  mutate(group = gsub("S1.PSI", "KO", group)) %>%
  mutate(group = gsub("S2.PSI", "WT", group),
         PSI = round(PSI, 3)) %>%
  mutate(group = factor(group, levels = c("WT", "KO"))) %>% 
  ggplot(aes(x = group, y = PSI)) +
  geom_boxplot(aes(fill = group)) +
  scale_fill_brewer() +
  theme_bw() +
  stat_compare_means(method = "t.test", na.rm = T, label.y = c(1.1), comparisons = list(c("WT", "KO")), size = 4) +
  labs(y = "PSI") + facet_wrap(. ~ Group) +
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white", size = 1, linetype = "solid"),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 11, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 11, angle = 45, hjust = 1),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 11, angle = 0))
dev.off()

### >>> Transcript expression
# load data
znf138.t.expr <- read.table("RawData/stringtie/ZNF138_t_expression.txt", header = T, sep = "\t")
colnames(znf138.t.expr) <- c("8Cell.rep1", "Transcript", "Gene", "8Cell.rep2", "ICM.rep1", "ICM.rep2")
oip5.t.expr <- read.table("RawData/stringtie/OIP5_t_expression.txt", header = T, sep = "\t")
colnames(oip5.t.expr) <- c("8Cell.rep1", "Transcript", "Gene", "8Cell.rep2", "ICM.rep1", "ICM.rep2")
ythdf3.t.expr <- read.table("RawData/stringtie/YTHDF3_t_expression.txt", header = T, sep = "\t")
colnames(ythdf3.t.expr) <- c("8Cell.rep1", "Transcript", "Gene", "8Cell.rep2", "ICM.rep1", "ICM.rep2")
# processing
znf138.t.expr %>%
  mutate(`8 Cell` = (`8Cell.rep1`+`8Cell.rep2`)/2,
         `ICM` = (`ICM.rep1`+`ICM.rep2`)/2,
         Splicing.T = case_when(Transcript %in% c("ENST00000307355", "ENST00000359735", "ENST00000440155", "ENST00000440598") ~ "Inclusion isoform",
                                Transcript %in% c("ENST00000430838") ~ "Skipping isoform")) %>% na.omit %>%
  group_by(Splicing.T) %>% summarise(`8 Cell` = sum(`8 Cell`), `ICM` = sum(`ICM`)) %>% mutate(ICM = ICM + 0.5) -> znf138.t.expr
oip5.t.expr %>%
  mutate(`8 Cell` = (`8Cell.rep1`+`8Cell.rep2`)/2,
         `ICM` = (`ICM.rep1`+`ICM.rep2`)/2,
         Splicing.T = case_when(Transcript %in% c("ENST00000220514") ~ "Inclusion isoform",
                                Transcript %in% c("ENST00000560640") ~ "Skipping isoform")) %>% na.omit %>%
  group_by(Splicing.T) %>% summarise(`8 Cell` = sum(`8 Cell`), `ICM` = sum(`ICM`)) -> oip5.t.expr
ythdf3.t.expr %>%
  mutate(`8 Cell` = (`8Cell.rep1`+`8Cell.rep2`)/2,
         `ICM` = (`ICM.rep1`+`ICM.rep2`)/2,
         Splicing.T = case_when(Transcript %in% c("ENST00000539294", "ENST00000617200", "ENST00000519428") ~ "Inclusion isoform",
                                Transcript %in% c("ENST00000517371") ~ "Skipping isoform")) %>% na.omit %>%
  group_by(Splicing.T) %>% summarise(`8 Cell` = sum(`8 Cell`), `ICM` = sum(`ICM`)) -> ythdf3.t.expr
pdf("Graphs/bi_K9_K27/repeat/splicing/SE_splicing_transcript_expression_level.pdf", height = 7, width = 3.5)
rbind(znf138.t.expr, ythdf3.t.expr, oip5.t.expr) %>% 
  mutate(Gene = c(rep("ZNF138", 2), rep("YTHDF3", 2), rep("OIP5", 2))) %>%
  gather(key = "Cell", value = "FPKM", -Splicing.T, -Gene) %>%
  as.data.frame() %>%
  ggplot(aes(x = Cell, y = FPKM)) +
  geom_bar(aes(fill = Cell), stat = "identity") +
  facet_wrap(Gene ~ Splicing.T, scales = "free", ncol = 2) +
  theme_bw() +
  scale_fill_brewer(direction = 1, palette = 3, type = "qual") +
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white", size = 1, linetype = "solid"),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 11, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 11, angle = 45, hjust = 1),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 11, angle = 0))
dev.off()
