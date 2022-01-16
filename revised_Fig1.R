################################################# ~ Good Luck ~ #########################################################
###### >>>>>>                       TITLE: H3K9me3 Project Revised Figure1 analysis                         <<<<<< ######

# - Figure 1 content:
# - 1st part: Global setting
# - 2nd part: QC: Hierarchical clustering (CMQ)
# - 3th part: QC: correlation analysis between technical replicates (YHW)
# - 4rd part: QC: correlation analysis between large and low cell number (YHW)
# - 5th part: Peaks enrichment analysis (CMQ)
# - 6th part: K93 and K273 Bivalent or univalent promoters (CMQ)
# - 7th part: Peaks genomic distribution (YHW)
# - 8th part: Expression level of repeats in embryo
# - 9th part: Define stage-specific peaks (YHW)



### ========================
### 1st part: Global setting
### ========================
setwd("/home/cmq/bioinfo/project-cmq/embryo_93")
# plot
library("ggplot2")
library("ggpubr")
library("ggrepel")
library("ComplexHeatmap")
library("circlize")
library("dendextend")
# data processing
library("dplyr")
library("tidyr")
library("tibble")
# math
library("amap")
# normalization
library("limma")
#save.image("scripts/revised_Fig2.RData")
#load("scripts/revised_Fig2.RData")



### ===========================================
### 2nd part: QC: Hierarchical clustering (CMQ)
### ==========================================

### >>> 1. load human embryo H3K9me3 global signal table
h93.emb.glo <- read.table("results/R/RawData/deeptools/coverage/GobalCoverage/gobal_readNor.tab", header = T, stringsAsFactors = F)
colnames(h93.emb.glo)[4:19] <- c("h4C93_Input_1", "h4C93_Input_2", "h4C93_Ip_1", "h4C93_Ip_2",
                                 "h8C93_Input_2", "h8C93_Ip_3", "h8C93_Ip_4",
                                 "hICM93_Input", "hICM93_Ip_1", "hICM93_Ip_2",
                                 "hM93_Input", "hM93_Ip_1", "hM93_Ip_2",
                                 "hTE93_Input", "hTE93_Ip_1", "hTE93_Ip_2")
h93.emb.glo$name <- paste(h93.emb.glo$chr, ":", h93.emb.glo$start, "-", h93.emb.glo$end, sep = "")

### >>> 2. load hTSCs and hESCs global signal table
h93.cell.glo <- read.table("/home/cmq/bioinfo/project-cmq/cell_line_93/results/R/Rawdata/deeptools/coverage/GobalCoverage/gobal_readNor.tab",
                           header = T, stringsAsFactors = F)
colnames(h93.cell.glo)[4:9] <- c("TSC_Ip_1", "TSC_Ip_2", "ESC_Input_1", "ESC_Ip_1", "TSC_Input_1", "TSC_Input_2")
h93.cell.glo$name <- paste(h93.cell.glo$chr, ":", h93.cell.glo$start, "-", h93.cell.glo$end, sep = "")

### >>> 3. merge embryo and cell line H3K9me3 signal
h93.glo <- merge(h93.emb.glo, h93.cell.glo, by = "name")

### >>> 4. clustering and visualization
# calculate LFC of RPKM
(h93.glo[, -c(1:4, 21:23)] + 0.1) %>%
  mutate(h4C1 = log2(h4C93_Ip_1/h4C93_Input_1), h4C2 = log2(h4C93_Ip_2/h4C93_Input_2),
         h8C3 = log2(h8C93_Ip_3/h8C93_Input_2), h8C4 = log2(h8C93_Ip_4/h8C93_Input_2),
         hICM1 = log2(hICM93_Ip_1/hICM93_Input), hICM2 = log2(hICM93_Ip_2/hICM93_Input),
         Morula1 = log2(hM93_Ip_1/hM93_Input), Morula2 = log2(hM93_Ip_2/hM93_Input),
         TE1 = log2(hTE93_Ip_1/hTE93_Input), TE2 = log2(hTE93_Ip_2/hTE93_Input),
         ESC = log2(ESC_Ip_1/ESC_Input_1), TSC = log2(TSC_m/TSC_Input_m)) -> h93.glo.dend
h93.glo.dend <- h93.glo.dend[, 27:38]
# remove batch effect between embryo and cell line
h93.glo.dend.rb <- removeBatchEffect(h93.glo.dend, batch = c(rep("A", 10), "B", "C"))
# clustering
dend <- h93.glo.dend.rb[,-ncol(h93.glo.dend.rb)] %>% t %>% scale %>% hcluster(method = "pearson", link = "complete") %>% as.dendrogram
labels(dend) <- c("4Cell-1", "4Cell-2", "8Cell-1", "8Cell-2", "ESC", "ICM-1", "ICM-2", "TE-1", "TE-2", "Morula-1", "Morula-2")
# plot
pdf("results/R/Graphs/correlation/human_embryo_H3K9me3_all_stage_global_dendrogram.pdf", width = 4, height = 6)
dend %>%
  # control the number, color and width of branches
  set("branches_k_color", value = c("#E3BF03", "#437ECE", "#E87A2A", "#B430BE", "#089498"), k = 5) %>% set("branches_lwd", 1.8) %>%
  # control the size and color of branch labels
  set("labels_cex", 1) %>% set("labels_col", value = c("#E3BF03", "#437ECE", "#E87A2A", "#B430BE", "#089498"), k = 5) %>%
  set("hang_leaves", 0.3) %>%
  # control size, shape and color of leaves
  set("leaves_cex", 1) %>% set("leaves_pch", 19) %>%
  set("leaves_col", c("#E3BF03")) %>%
  plot(horiz = TRUE, axes = T)
dend %>% rect.dendrogram(k = 5, horiz = TRUE, border = 8, lty = 5, lwd = 1)
dev.off()



### =====================================================================
### 3th part: QC: correlation analysis between technical replicates (YHW)
### =====================================================================

### >>> 1. correlation of H3K9me3 signal on all stage merged peaks between biological replicates
# load H3K9me3 signal (RPKM)
h93.pks.file <- list.files("results/deeptools/coverage/PeaksCoverage", ".tab", full.names = T)
h93.pks.file <- h93.pks.file[1]
h93.pks <- read.table(h93.pks.file, header = T, stringsAsFactors = F)
colnames(h93.pks)[-1:-3] <- c("h4C93_Input_1", "h4C93_Input_2", "h4C93_Ip_1", "h4C93_Ip_2",
                              "h8C93_Input_2", "h8C93_Ip_3", "h8C93_Ip_4",
                              "hICM93_Input", "hICM93_Ip_1", "hICM93_Ip_2",
                              "hM93_Input", "hM93_Ip_1", "hM93_Ip_2",
                              "hTE93_Input", "hTE93_Ip_1", "hTE93_Ip_2")
# visualization via density plot
pdf("results/R/Graphs/correlation/peaks/human_embryo_H3K9me3_all_stage_replicates_on_all_peaks_log2rpkm_cor_density.pdf", height = 3.5, width = 3)
plot.Lab.palette <- colorRampPalette(c("#ffffff", "#0B397D", "#30999A"))
sample <- c("h4C93_Ip", "h8C93_Ip", "hM93_Ip", "hICM93_Ip", "hTE93_Ip")
stage <- c("4 Cell", "8 Cell", "Morula", "ICM", "TE")
for (x in seq(1, length(sample))){
  print(grep(sample[x], colnames(h93.pks), value = T))
  smoothScatter(log2(h93.pks[, grep(sample[x], colnames(h93.pks), value = T)[1]] + 0.1),
                log2(h93.pks[, grep(sample[x], colnames(h93.pks), value = T)[2]] + 0.1),
                colramp = plot.Lab.palette, cex = 0.5, pch = 16,
                nrpoints = 60, col = "#494949", main = stage[x],
                xlab = "Replicates 1", ylab = "Replicates 2", xlim = c(0,15), ylim = c(0,15))
}; rm(x)
dev.off()
# print correlation coefficient 
for (x in seq(1, length(sample))){
  print(grep(sample[x], colnames(h93.pks), value = T))
  print(cor.test(h93.pks[, grep(sample[x], colnames(h93.pks), value = T)[1]],
                 h93.pks[, grep(sample[x], colnames(h93.pks), value = T)[2]],
                 method = "pearson"))
}; rm(x)

### >>> 2. correlation between H3K9me3 and ATAC-seq signal on H3K9me3 and ATAC-seq merged peaks
# load RPKM tables
rpkm.atac.file <- list.files("results/deeptools/coverage/diff_his/h3k9me3_atac", ".tab", full.names = T)
rpkm.atac <- list()
for (i in seq(1, length(rpkm.atac.file))) {
  rpkm.atac[[i]] <- read.table(rpkm.atac.file[i], header = T, stringsAsFactors = F)
  colnames(rpkm.atac[[i]])[-1:-3] <- c("atac_4Cell_1", "atac_4Cell_2", "atac_8Cell_1", "atac_8Cell_2", "atac_ESC_1", "atac_ESC_2",
                                       "k93_4Cell_input", "k93_4Cell", "k93_8Cell_input", "k93_8Cell", "k93_ICM_input", "k93_ICM",
                                       "k93_Morula_input", "k93_Morula", "k93_TE_input", "k93_TE", "atac_ICM_1", "atac_ICM_2",
                                       "k93_ESC500_input", "k93_ESC50_input", "k93_ESC500", "k93_ESC50", "atac_TE_1", "atac_TE_2")
  rpkm.atac[[i]] %>% mutate(atac_4Cell_ave = (atac_4Cell_1 + atac_4Cell_2)/2, atac_8Cell_ave = (atac_8Cell_1 + atac_8Cell_2)/2,
                            atac_ICM_ave = (atac_ICM_1 + atac_ICM_2)/2, atac_TE_ave = (atac_TE_1 + atac_TE_2)/2,
                            atac_ESC50_ave = (atac_ESC_1 + atac_ESC_2)/2, atac_ESC500_ave = (atac_ESC_1 + atac_ESC_2)/2) %>%
    mutate(k93_4Cell_LFC = log2((k93_4Cell+0.1)/(k93_4Cell_input+0.1)), k93_8Cell_LFC = log2((k93_8Cell+0.1)/(k93_8Cell_input+0.1)),
           k93_ICM_LFC = log2((k93_ICM+0.1)/(k93_ICM_input+0.1)), k93_TE_LFC = log2((k93_TE+0.1)/(k93_TE_input+0.1)),
           k93_ESC50_LFC = log2((k93_ESC50+0.1)/(k93_ESC50_input+0.1)), k93_ESC500_LFC = log2((k93_ESC500+0.1)/(k93_ESC500_input+0.1))) %>%
    mutate(k93_4Cell_FC = (k93_4Cell+0.1)/(k93_4Cell_input+0.1), k93_8Cell_FC = (k93_8Cell+0.1)/(k93_8Cell_input+0.1),
           k93_ICM_FC = (k93_ICM+0.1)/(k93_ICM_input+0.1), k93_TE_FC = (k93_TE+0.1)/(k93_TE_input+0.1),
           k93_ESC50_FC = (k93_ESC50+0.1)/(k93_ESC50_input+0.1), k93_ESC500_FC = (k93_ESC500+0.1)/(k93_ESC500_input+0.1)) %>%
    mutate(Type = case_when(end-start >= 2000 ~ "long", end-start < 2000 ~ "short")) -> rpkm.atac[[i]]
  
  names(rpkm.atac)[i] <- str_split_fixed(str_split_fixed(rpkm.atac.file[i], "/", 6)[,6], "_", 4)[,3]
}; rm(i)
# correlation analysis
rpkm.atac.cor <- matrix(ncol = 3, nrow = 6)
colnames(rpkm.atac.cor) <- c("Stage", "Pearson", "Spearman")
for (i in seq(1, length(rpkm.atac))) {
  tmp <- rpkm.atac[[i]][rpkm.atac[[i]]$Type=="long",grep(paste(names(rpkm.atac)[i], "_", sep = ""), colnames(rpkm.atac[[i]]))]
  rpkm.atac.cor[i,1] <- paste(grep("_ave", colnames(tmp), value = T), grep("_LFC", colnames(tmp), value = T), sep = "_vs_")
  rpkm.atac.cor[i,2] <- as.numeric(cor.test(tmp[,grep("_ave", colnames(tmp))], tmp[,grep("_LFC", colnames(tmp))], method = "pearson")$estimate)
  rpkm.atac.cor[i,3] <- as.numeric(cor.test(tmp[,grep("_ave", colnames(tmp))], tmp[,grep("_LFC", colnames(tmp))], method = "spearman")$estimate)
}; rm(i, tmp)
# plot
rpkm.atac.cor <- as.data.frame(rpkm.atac.cor)
rpkm.atac.cor$Stage <- str_split_fixed(rpkm.atac.cor$Stage, "_", 3)[,2]
rpkm.atac.cor %>% gather(key = "Method", value = "Coefficient", -Stage) %>%
  mutate(Coefficient = as.numeric(Coefficient)) %>%
  mutate(Stage = factor(Stage, levels = c("4Cell", "8Cell", "ICM", "TE", "ESC50", "ESC500")),
         Method = factor(Method, levels = c("Spearman", "Pearson"))) -> rpkm.atac.cor
pdf("results/R/Graphs/correlation/diff_his/human_embryo_correlation_of_signal_of_H3K9me3_vs_ATAC_seq_noTE_on_peaks.pdf",
    height = 2.5, width = 3)
rpkm.atac.cor %>% filter(!Stage == "TE") %>%
  ggplot(aes(x = Stage, y = Coefficient)) +
  geom_bar(stat = "identity", fill = "#0B488A") + theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(face="plain", colour = "#000000", size = 10, angle = 90),
        axis.text.x  = element_text(face="plain", colour = "#000000", size = 10, angle = 45),
        axis.text.y  = element_text(face="plain", colour = "#000000", size = 10, angle = 0)) +
  facet_grid(. ~ Method)
dev.off()

### >>> 3. correlation between H3K9me3 and LiCAT-seq signal on H3K9me3 and LiCAT-seq merged peaks
# load rpkm table
rpkm.licat.file <- list.files("results/deeptools/coverage/diff_his/h3k9me3_licat_2", ".tab", full.names = T)
rpkm.licat <- list()
for (i in seq(1, length(rpkm.licat.file))) {
  rpkm.licat[[i]] <- read.table(rpkm.licat.file[i], header = T, stringsAsFactors = F)
  rpkm.licat[[i]] <- rpkm.licat[[i]][,-c(17,19)]
  colnames(rpkm.licat[[i]])[-1:-3] <- c("k93_4Cell_input", "k93_4Cell", "k93_8Cell_input", "k93_8Cell", "k93_ICM_input", "k93_ICM",
                                        "k93_Morula_input", "k93_Morula", "k93_TE_input", "k93_TE", "k93_ESC500n_input", "k93_ESC50n_input",
                                        "k93_ESC500n", "k93_ESC50n", "licat_ESC", "licat_Morula", "licat_ICM", "licat_4Cell", "licat_8Cell", "licat_TE")
  rpkm.licat[[i]] %>%
    mutate(licat_ESC500n = licat_ESC, licat_ESC50n = licat_ESC) %>%
    mutate(k93_4Cell_LFC = log2((k93_4Cell+0.1)/(k93_4Cell_input+0.1)), k93_8Cell_LFC = log2((k93_8Cell+0.1)/(k93_8Cell_input+0.1)),
           k93_ICM_LFC = log2((k93_ICM+0.1)/(k93_ICM_input+0.1)), k93_TE_LFC = log2((k93_TE+0.1)/(k93_TE_input+0.1)),
           k93_Morula_LFC = log2((k93_Morula+0.1)/(k93_Morula_input+0.1)),
           k93_ESC50n_LFC = log2((k93_ESC50n+0.1)/(k93_ESC50n_input+0.1)), k93_ESC500n_LFC = log2((k93_ESC500n+0.1)/(k93_ESC500n_input+0.1))) %>%
    mutate(k93_4Cell_FC = (k93_4Cell+0.1)/(k93_4Cell_input+0.1), k93_8Cell_FC = (k93_8Cell+0.1)/(k93_8Cell_input+0.1),
           k93_Morula_FC = (k93_Morula+0.1)/(k93_Morula_input+0.1),
           k93_ICM_FC = (k93_ICM+0.1)/(k93_ICM_input+0.1), k93_TE_FC = (k93_TE+0.1)/(k93_TE_input+0.1),
           k93_ESC50n_FC = (k93_ESC50n+0.1)/(k93_ESC50n_input+0.1), k93_ESC500n_FC = (k93_ESC500n+0.1)/(k93_ESC500n_input+0.1)) %>%
    mutate(Type = case_when(end-start >= 2000 ~ "long", end-start < 2000 ~ "short")) -> rpkm.licat[[i]]
}; rm(i)
names(rpkm.licat) <- c("4Cell", "8Cell", "ESC50n", "ESC500n", "ICM", "Morula", "TE")
# correlation analysis
rpkm.licat.cor <- matrix(ncol = 3, nrow = 7)
colnames(rpkm.licat.cor) <- c("Stage", "Pearson", "Spearman")
for (i in seq(1, length(rpkm.licat))) {
  tmp <- rpkm.licat[[i]][rpkm.licat[[i]]$Type=="long",grep(names(rpkm.licat)[i], colnames(rpkm.licat[[i]]))]
  rpkm.licat.cor[i,1] <- paste(grep("licat", colnames(tmp), value = T), grep("_LFC", colnames(tmp), value = T), sep = "_vs_")
  rpkm.licat.cor[i,2] <- as.numeric(cor.test(tmp[,grep("licat", colnames(tmp))], tmp[,grep("_LFC", colnames(tmp))], method = "pearson")$estimate)
  rpkm.licat.cor[i,3] <- as.numeric(cor.test(tmp[,grep("licat", colnames(tmp))], tmp[,grep("_LFC", colnames(tmp))], method = "spearman")$estimate)
}; rm(i, tmp)
# plot
rpkm.licat.cor <- as.data.frame(rpkm.licat.cor)
rpkm.licat.cor$Stage <- str_split_fixed(rpkm.licat.cor$Stage, "_", 3)[,2]
rpkm.licat.cor %>% gather(key = "Method", value = "Coefficient", -Stage) %>%
  mutate(Coefficient = as.numeric(Coefficient)) %>%
  mutate(Stage = factor(Stage, levels = c("4Cell", "8Cell", "Morula", "ICM", "TE", "ESC50n", "ESC500n")),
         Method = factor(Method, levels = c("Spearman", "Pearson"))) -> rpkm.licat.cor
pdf("results/R/Graphs/correlation/diff_his/human_embryo_correlation_of_signal_of_H3K9me3_vs_LiCAT_seq_on_peaks.pdf",
    height = 2.5, width = 3)
rpkm.licat.cor %>% ggplot(aes(x = Stage, y = Coefficient)) +
  geom_bar(stat = "identity", fill = "#0B488A") + theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(face="plain", colour = "#000000", size = 10, angle = 90),
        axis.text.x  = element_text(face="plain", colour = "#000000", size = 10, angle = 45),
        axis.text.y  = element_text(face="plain", colour = "#000000", size = 10, angle = 0)) +
  facet_grid(. ~ Method)
dev.off()



### ==========================================================================
### 4rd part: QC: correlation analysis between large and low cell number (YHW)
### ==========================================================================

### >>> 1. Load data
hesc.ln.global <- read.table("RawData/deeptools/multibw/hESC_low_num_global_raw/raw_bam_readNor.tab", sep = "\t", header = T)
colnames(hesc.ln.global)[-1:-3] <- c("Input.large", "IP.large", "Input.public1", "Input.public2", "IP.public1", "IP.public2",
                                     str_split_fixed(colnames(hesc.ln.global)[-1:-9], "_", 2)[, 1])

### >>> 2. Correlation analysis
cor.test(hesc.ln.global[, -1:-3]$IP50, 
         hesc.ln.global[, -1:-3]$IP.large, 
         method = "pearson")
cor.test(log2(hesc.ln.global[, -1:-3]$IP500+0.1),
         log2(hesc.ln.global[, -1:-3]$IP.large+0.1), 
         method = "pearson")
pdf("Graphs/dynamic_change/correlation/Correlation_H3K9me3_signal_between_large_and_small_cell_num_50_vs_bulk.pdf", height = 4.5, width = 4)
smoothScatter(log2(hesc.ln.global[, -1:-3]$IP50+0.1), log2(hesc.ln.global[, -1:-3]$IP.large+0.1), 
              colramp = colorRampPalette(c("#ffffff", "#0B397D", "#30999A")), 
              cex = 0.75, pch = 16, nrpoints = 100, col = "#494949", main = "50 vs Bulk", 
              xlab = "Group 1", ylab = "Group 2", xlim = c(0, 12), ylim = c(0, 12))
dev.off()
pdf("Graphs/dynamic_change/correlation/Correlation_H3K9me3_signal_between_large_and_small_cell_num_500_vs_bulk.pdf", height = 4.5, width = 4)
smoothScatter(log2(hesc.ln.global[, -1:-3]$IP500+0.1), log2(hesc.ln.global[, -1:-3]$IP.large+0.1), 
              colramp = colorRampPalette(c("#ffffff", "#0B397D", "#30999A")), 
              cex = 0.75, pch = 16, nrpoints = 100, col = "#494949", main = "500 vs Bulk", 
              xlab = "Group 1", ylab = "Group 2", xlim = c(0, 12), ylim = c(0, 12))
dev.off()
png("Graphs/Correlation_H3K9me3_signal_between_large_and_small_cell_num.png", res = 600, height = 3, width = 3, units = "in")
log2(hesc.ln.global[, -1:-3] + 0.1) %>% 
  ggplot(aes(x = IP500, y = IP.large)) + 
  geom_point(size = 1.5, color = "#636e72", alpha = 0.3, shape = 16) +
  geom_smooth(method = "lm", se = TRUE, colour = "#c44569") +
  labs(x = "Replicates 1", y = "Replicates 2") +
  stat_cor(method = "pearson", p.accuracy = 0.0001, size = 5) + 
  theme_classic() + ggtitle(label = "") + 
  labs(x = "hESCs H3K9me3 (500 cell)", y = "hESCs H3K9me3 (10e7 cell)") +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0))
dev.off()

### >>> 3. Signal of captured and missed peaks
hs.k9.cap.mis <- list()
for (i in seq(1, length(list.files("RawData/bedtools/quality_control", ".bed")))) {
  tmp <- read.table(list.files("RawData/bedtools/quality_control", ".bed", full.names = T)[i], sep = "\t", header = F)
  tmp$id <- paste(tmp$V1, ":", tmp$V2, "-", tmp$V3, sep = "")
  hs.k9.cap.mis[[i]] <- tmp
}; rm(i, tmp)
names(hs.k9.cap.mis) <- list.files("RawData/bedtools/quality_control", ".bed")
hesc.large.pks <- read.table("RawData/deeptools/multibw/hESCs_p0.01pks_raw/hESCs_large_cell_number_p0.01_readNor.tab", 
                             sep = "\t", header = T)
colnames(hesc.large.pks)[-1:-3] <- c("Input.large", "IP.large", "Input.public1", "Input.public2", "IP.public1", "IP.public2",
                                     str_split_fixed(colnames(hesc.ln.global)[-1:-9], "_", 2)[, 1])
hesc.large.pks$id <- paste(hesc.large.pks$chr, ":", hesc.large.pks$start, "-", hesc.large.pks$end, sep = "")
hesc.large.pks %>% mutate(type.500cell = case_when(id %in% hs.k9.cap.mis$hESCs_large_cell_cap_peak_by_500cell.bed$id ~ "500cell cap",
                                                   id %in% hs.k9.cap.mis$hESCs_large_cell_mis_peak_by_500cell.bed$id ~ "500cell mis"),
                          type.50cell = case_when(id %in% hs.k9.cap.mis$hESCs_large_cell_cap_peak_by_50cell.bed$id ~ "50cell cap",
                                                  id %in% hs.k9.cap.mis$hESCs_large_cell_mis_peak_by_50cell.bed$id ~ "50cell mis")) -> hesc.large.pks
p1 <- hesc.large.pks[, c(8, 15, 16)] %>%
  ggplot(aes(x = type.500cell, y = log2(IP.public1+0.1))) +
  geom_boxplot(aes(fill = type.500cell), width = 0.6, outlier.size = 0.3) +
  theme_bw() +
  scale_fill_manual(values = c(`500cell cap`= "#00AFBB", `500cell mis` = "#E7B800")) +
  theme(legend.position = "top",
        strip.background = element_rect(colour = "black", fill = "white", size = 0.75, linetype = "solid"),
        strip.text.x = element_text(face = "plain", size = 11, colour = "black", angle = 0),
        strip.text.y = element_text(face = "plain", size = 11, colour = "black", angle = 0),
        panel.grid.minor = element_blank(),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 11, angle = 45),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 11),
        axis.title  = element_text(face = "plain", colour = "#000000", size = 13)) +
  labs(y = "log2(FPKM+1)", x = "Groups") +
  stat_compare_means(method = "t.test", na.rm = T, label.y = 9, comparisons = list(c("500cell cap", "500cell mis")), size = 4)
p2 <- hesc.large.pks[, c(8, 15, 16)] %>%
  ggplot(aes(x = type.50cell, y = log2(IP.public1+0.1))) +
  geom_boxplot(aes(fill = type.50cell), width = 0.6, outlier.size = 0.3) +
  theme_bw() +
  scale_fill_manual(values = c(`50cell cap`= "#00AFBB", `50cell mis` = "#E7B800")) +
  theme(legend.position = "top",
        strip.background = element_rect(colour = "black", fill = "white", size = 0.75, linetype = "solid"),
        strip.text.x = element_text(face = "plain", size = 11, colour = "black", angle = 0),
        strip.text.y = element_text(face = "plain", size = 11, colour = "black", angle = 0),
        panel.grid.minor = element_blank(),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 11, angle = 45),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 11),
        axis.title  = element_text(face = "plain", colour = "#000000", size = 13)) +
  labs(y = "log2(FPKM+1)", x = "Groups") +
  stat_compare_means(method = "t.test", na.rm = T, label.y = 9, comparisons = list(c("50cell cap", "50cell mis")), size = 4)
pdf("Graphs/dynamic_change/statistics/Captured_and_missed_peaks_FPKM_relative_known_and_large_cell_data.pdf", height = 4, width = 4)
plot_grid(p1, p2, cols = 2)
dev.off()

### >>> 4. Ratio of captured peaks
pdf("Graphs/dynamic_change/statistics/Captured_peaks_ratio_relative_known_and_large_cell_data.pdf", height = 3.5, width = 4.5)
data.frame(cell = c(rep("hESCs 500 cell", 3), rep("hESCs 50 cell", 3), 
                    rep("Embryo 8 cell H3K9me3", 3), rep("Embryo 8 cell H3K27ac", 3)),
           peak.type = rep(c("total peaks", "gene peaks", "repeat peaks"), 4),
           histone = c(rep("H3K9me3", 6), rep("H3K9me3", 3), rep("H3K27ac", 3)), 
           cap.ratio = c(12436/14630, 886/1098, 11543/13554, 
                         10843/14630, 792/1098, 10564/13554,
                         8151/13587, 896/1451, 6448/10549,
                         48971/76729, 5526/8377, 17168/31042)*100) %>%
  mutate(peak.type = factor(peak.type, levels = c("total peaks", "gene peaks", "repeat peaks"))) %>%
  filter(cell %in% c("hESCs 50 cell", "hESCs 500 cell")) %>%
  ggplot(aes(x = cell, y = cap.ratio, fill = peak.type)) +
  geom_col(position = position_dodge(0.75), width = 0.6) +
  scale_fill_manual(values = c(`total peaks` = "#00AFBB", `gene peaks` = "#E7B800", `repeat peaks` = "#FC4E08")) +
  ylim(c(0, 100)) + theme_bw() + 
  #facet_grid(. ~ histone, scales = "free") +
  theme(strip.background = element_rect(colour = "black", fill = "white", size = 0.75, linetype = "solid"),
        strip.text.x = element_text(face = "plain", size = 11, colour = "black", angle = 0),
        strip.text.y = element_text(face = "plain", size = 11, colour = "black", angle = 0),
        panel.grid.minor = element_blank(),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 11, angle = 45),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 11),
        axis.title  = element_text(face = "plain", colour = "#000000", size = 13)) +
  labs(y = "Ratio (percentage %)", x = "Groups")
dev.off()



### =========================================
### 5th part: Peaks enrichment analysis (CMQ)
### =========================================

### >>> 1. H3K9me3 peaks
# load regionR results from shell
# repeat
re.cl.en50 <- readRDS("results/R/RawData/RegionSetEnrich/p_0.05/class_50_enrichment_analysis_data.rds")
names(re.cl.en50) <- c("4 Cell", "8 Cell", "ICM", "Morula", "TE")

for(i in seq(1, length(re.cl.en50))){
  re.cl.en50[[i]]$FoldChange <- re.cl.en50[[i]]$Observed/re.cl.en50[[i]]$Expected
  used <- rownames(re.cl.en50[[i]])%in%c("LINE", "SINE", "LTR", "DNA", "Satellite")
  re.cl.en50[[i]] <- re.cl.en50[[i]][used, -5]
}; rm(i);rm(used)
re.cl.en50 <- do.call(rbind, re.cl.en50)
re.cl.en50 <- cbind(re.cl.en50, str_split_fixed(rownames(re.cl.en50), "\\.", 2))
colnames(re.cl.en50)[6:7] <- c("Stage", "Type")
# promoter
pro.en50 <- readRDS("results/R/RawData/RegionSetEnrich/p_0.05/Homo_sapiens.GRCh38.97.ens2ucsc.gene.promoter.1k_50_enrichment_analysis_data.rds")
names(pro.en50) <- c("4 Cell", "8 Cell", "ICM", "Morula", "TE")
for(i in seq(1, length(pro.en50))){
  pro.en50[[i]] <- pro.en50[[i]][-5, 1]
}
pro.en50 <- as.data.frame(do.call(rbind, pro.en50))
colnames(pro.en50) <- c("Observed", "Expected", "Zscore", "Pvalue")
pro.en50 %>% rownames_to_column("Stage") %>% mutate(FoldChange = Observed/Expected, Type = "Promoter") -> pro.en50
rownames(pro.en50) <- paste(pro.en50$Stage, pro.en50$Type, sep = ".")
# genebody
ge.en50 <- readRDS("results/R/RawData/RegionSetEnrich/p_0.05/Homo_sapiens.GRCh38.97.ens2ucsc.gene.body_50_enrichment_analysis_data.rds")
names(ge.en50) <- c("4 Cell", "8 Cell", "ICM", "Morula", "TE")

for(i in seq(1, length(ge.en50))){
  ge.en50[[i]] <- ge.en50[[i]][-5, 1]
}
ge.en50 <- as.data.frame(do.call(rbind, ge.en50))
colnames(ge.en50) <- c("Observed", "Expected", "Zscore", "Pvalue")
ge.en50 %>% rownames_to_column("Stage") %>% mutate(FoldChange = Observed/Expected, Type = "Gene") -> ge.en50
rownames(ge.en50) <- paste(ge.en50$Stage, ge.en50$Type, sep = ".")
# combine all stages and types
all.en <- rbind(pro.en50, re.cl.en50, ge.en50)
# visualization
all.en$Type <- factor(all.en$Type, levels = c("LTR", "LINE", "SINE", "DNA", "Satellite", "Promoter", "Gene"))
empty_bar <- 3
to_add <- data.frame(matrix(NA, empty_bar*nlevels(all.en$Type), ncol(all.en)))
colnames(to_add) <- colnames(all.en)
to_add$Type <- rep(levels(all.en$Type), each = empty_bar)
all.en <- rbind(all.en, to_add)
all.en <- all.en[order(match(as.character(all.en$Stage), c("4 Cell", "8 Cell", "Morula", "ICM", "TE", "NA"))), ]
all.en <- all.en[order(match(as.character(all.en$Type), c("LTR", "LINE", "SINE", "DNA", "Satellite", "Promoter", "Gene"))), ]
all.en$id <- seq(1, nrow(all.en))
label_data <- all.en[all.en$Type%in%c("Promoter", "LINE"),]
number_of_bar <- nrow(all.en)
angle <- 90 - 360 * (label_data$id-0.5)/number_of_bar
label_data$hjust <- ifelse(angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)
base_data <- all.en %>%
  group_by(Type) %>%
  summarize(start=min(id), end=max(id) - empty_bar) %>%
  rowwise() %>%
  mutate(title=mean(c(start, end)))
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]
pdf("results/R/Graphs/region_set_enrich/Human_Embryo_H3K9me3_RegionSetEnrichment_Analysis_Circular_BarChart_p_0.05.pdf",
    width = 6, height = 6)
ggplot(all.en, aes(x = as.factor(id), y = FoldChange)) +
  geom_bar(aes(fill = Type), stat = "identity") +
  labs(title = "H3K9me3") +
  geom_segment(data=grid_data, aes(x = end, y = 2, xend = start, yend = 2), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
  geom_segment(data=grid_data, aes(x = end, y = 1.5, xend = start, yend = 1.5), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
  geom_segment(data=grid_data, aes(x = end, y = 1, xend = start, yend = 1), colour = "#212F3D", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
  geom_segment(data=grid_data, aes(x = end, y = 0.5, xend = start, yend = 0.5), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
  annotate("text", x = rep(max(all.en$id), 4), y = c(0.5, 1, 1.5, 2), label = c("0.5", "1.0", "1.5", "2.0") ,
           color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  geom_hline(yintercept = 1, colour = "#000000", linetype = "longdash", size = 0.25, alpha = 0.8) +
  ylim(-3, 7) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(0, 10), "cm"),
    plot.title = element_text(hjust = 0.5, size = 15)) +
  coord_polar(start = 0) +
  geom_text(data = label_data, aes(x = id, y = FoldChange + 0.5, label = Stage, hjust = hjust),
            color = "black", fontface = "bold", alpha = 1, size = 2, angle = label_data$angle, inherit.aes = FALSE) +
  geom_segment(data = base_data, aes(x = start, y = -0.5, xend = end, yend = -0.5),
               colour = "black", alpha = 1, size = 1, inherit.aes = FALSE) +
  geom_text(data = base_data, aes(x = title, y = -1.5, label = Type), colour = "black",
            alpha = 0.8, size = 2.5, fontface = "bold", inherit.aes = FALSE)
dev.off()

### >>> 2. H3K27me3 peaks
# load regionR results from shell
# repeat
re.cl.27.en50 <- readRDS("/home/cmq/bioinfo/project-cmq/gse124718/results/idr/multi_merged_BAMPE_p_0.05/filtered/class_50_enrichment_analysis_data.rds")
names(re.cl.27.en50) <- c("2 Cell", "4 Cell", "8 Cell", "ICM", "TE")
for(i in seq(1, length(re.cl.27.en50))){
  re.cl.27.en50[[i]]$FoldChange <- re.cl.27.en50[[i]]$Observed/re.cl.27.en50[[i]]$Expected
  used <- rownames(re.cl.27.en50[[i]])%in%c("LINE", "SINE", "LTR", "DNA", "Satellite")
  re.cl.27.en50[[i]] <- re.cl.27.en50[[i]][used, -5]
}; rm(i);rm(used)
re.cl.27.en50 <- do.call(rbind, re.cl.27.en50)
re.cl.27.en50 <- cbind(re.cl.27.en50, str_split_fixed(rownames(re.cl.27.en50), "\\.", 2))
colnames(re.cl.27.en50)[6:7] <- c("Stage", "Type")
# promoter
pro.27.en50 <- readRDS("/home/cmq/bioinfo/project-cmq/gse124718/results/idr/multi_merged_BAMPE_p_0.05/filtered/Homo_sapiens.GRCh38.97.ens2ucsc.gene.promoter.1k_50_enrichment_analysis_data.rds")
names(pro.27.en50) <- c("2 Cell", "4 Cell", "8 Cell", "ICM", "TE")
as.data.frame(do.call(rbind, pro.27.en50))[,-5] %>% rownames_to_column("Stage") %>%
  mutate(FoldChange = as.numeric(Observed)/as.numeric(Expected), Type = "Promoter") -> pro.27.en50
rownames(pro.27.en50) <- paste(pro.27.en50$Stage, pro.27.en50$Type, sep = ".")
# gene body
ge.27.en50 <- readRDS("/home/cmq/bioinfo/project-cmq/gse124718/results/idr/multi_merged_BAMPE_p_0.05/filtered/Homo_sapiens.GRCh38.97.ens2ucsc.gene.body_50_enrichment_analysis_data.rds")
names(ge.27.en50) <- c("2 Cell", "4 Cell", "8 Cell", "ICM", "TE")
as.data.frame(do.call(rbind, ge.27.en50))[,-5] %>% rownames_to_column("Stage") %>%
  mutate(FoldChange = as.numeric(Observed)/as.numeric(Expected), Type = "Gene") -> ge.27.en50
rownames(ge.27.en50) <- paste(ge.27.en50$Stage, ge.27.en50$Type, sep = ".")
# repeat (morula)
re.cl.27.mo.en50 <- readRDS("/home/yhw/bioinfo/project-mine/Embryo.93/chipseq/gse123023/results/idr/broad_p0.05/filtered/class_50_enrichment_analysis_data.rds")
names(re.cl.27.mo.en50) <- c("Morula")
re.cl.27.mo.en50[[1]]$FoldChange <- re.cl.27.mo.en50[[1]]$Observed/re.cl.27.mo.en50[[1]]$Expected
re.cl.27.mo.en50[[1]] <- re.cl.27.mo.en50[[1]][rownames(re.cl.27.mo.en50[[1]])%in%c("LINE", "SINE", "LTR", "DNA", "Satellite"), -5]
re.cl.27.mo.en50 <- do.call(rbind, re.cl.27.mo.en50)
re.cl.27.mo.en50 <- cbind(re.cl.27.mo.en50, str_split_fixed(rownames(re.cl.27.mo.en50), "\\.", 2))
colnames(re.cl.27.mo.en50)[6:7] <- c("Stage", "Type")
# promoter (morula)
pro.27.mo.en50 <- readRDS("/home/yhw/bioinfo/project-mine/Embryo.93/chipseq/gse123023/results/idr/broad_p0.05/filtered/Homo_sapiens.GRCh38.97.ens2ucsc.gene.promoter.1k_50_enrichment_analysis_data.rds")
names(pro.27.mo.en50) <- c("Morula")
as.data.frame(do.call(rbind, pro.27.mo.en50))[,-5] %>% rownames_to_column("Stage") %>%
  mutate(FoldChange = as.numeric(Observed)/as.numeric(Expected), Type = "Promoter") -> pro.27.mo.en50
rownames(pro.27.mo.en50) <- paste(pro.27.mo.en50$Stage, pro.27.mo.en50$Type, sep = ".")
# gene body (morula)
ge.27.mo.en50 <- readRDS("/home/yhw/bioinfo/project-mine/Embryo.93/chipseq/gse123023/results/idr/broad_p0.05/filtered/Homo_sapiens.GRCh38.97.ens2ucsc.gene.body_50_enrichment_analysis_data.rds")
names(ge.27.mo.en50) <- c("Morula")
as.data.frame(do.call(rbind, ge.27.mo.en50))[,-5] %>% rownames_to_column("Stage") %>%
  mutate(FoldChange = as.numeric(Observed)/as.numeric(Expected), Type = "Gene") -> ge.27.mo.en50
rownames(ge.27.mo.en50) <- paste(ge.27.mo.en50$Stage, ge.27.mo.en50$Type, sep = ".")
# combine all stages and types
all.27.en50 <- rbind(re.cl.27.en50, pro.27.en50, ge.27.en50, re.cl.27.mo.en50, pro.27.mo.en50, ge.27.mo.en50)
# visualization
all.27.en50$Type <- factor(all.27.en50$Type, levels = c("LTR", "LINE", "SINE", "DNA", "Satellite", "Promoter", "Gene"))
empty_bar <- 3
to_add <- data.frame(matrix(NA, empty_bar*nlevels(all.27.en50$Type), ncol(all.27.en50)))
colnames(to_add) <- colnames(all.27.en50)
to_add$Type <- rep(levels(all.27.en50$Type), each = empty_bar)
all.27.en50 <- rbind(all.27.en50, to_add)
all.27.en50 <- all.27.en50[order(match(as.character(all.27.en50$Stage), c("2 Cell", "4 Cell", "8 Cell", "Morula", "ICM", "TE", "NA"))), ]
all.27.en50 <- all.27.en50[order(match(as.character(all.27.en50$Type), c("LTR", "LINE", "SINE", "DNA", "Satellite", "Promoter", "Gene"))), ]
all.27.en50$id <- seq(1, nrow(all.27.en50))
label_data <- all.27.en50[all.27.en50$Type%in%c("Promoter", "LINE"),]
number_of_bar <- nrow(all.27.en50)
angle <- 90 - 360 * (label_data$id-0.5)/number_of_bar
label_data$hjust <- ifelse(angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)
base_data <- all.27.en50 %>%
  group_by(Type) %>%
  summarize(start=min(id), end=max(id) - empty_bar) %>%
  rowwise() %>%
  mutate(title=mean(c(start, end)))
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]
pdf("results/R/Graphs/region_set_enrich/Human_Embryo_H3K27me3_RegionSetEnrichment_Analysis_Circular_BarChart_p_0.05.pdf",
    width = 8, height = 6)
ggplot(all.27.en50, aes(x = as.factor(id), y = FoldChange)) +
  geom_bar(aes(fill = Type), stat = "identity") +
  labs(title = "H3K27me3") +
  geom_segment(data=grid_data, aes(x = end, y = 2, xend = start, yend = 2), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
  geom_segment(data=grid_data, aes(x = end, y = 1.5, xend = start, yend = 1.5), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
  geom_segment(data=grid_data, aes(x = end, y = 1, xend = start, yend = 1), colour = "#212F3D", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
  geom_segment(data=grid_data, aes(x = end, y = 0.5, xend = start, yend = 0.5), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
  annotate("text", x = rep(max(all.27.en50$id), 4), y = c(0.5, 1, 1.5, 2), label = c("0.5", "1.0", "1.5", "2.0") ,
           color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  geom_hline(yintercept = 1, colour = "#000000", linetype = "longdash", size = 0.25, alpha = 0.8) +
  ylim(-3, 5) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(0, 30), "cm"),
    plot.title = element_text(hjust = 0.5, size = 15)) +
  coord_polar(start = 0) +
  geom_text(data = label_data, aes(x = id, y = FoldChange + 0.5, label = Stage, hjust = hjust),
            color = "black", fontface = "bold", alpha = 1, size = 2.5, angle = label_data$angle, inherit.aes = FALSE) +
  geom_segment(data = base_data, aes(x = start, y = -0.5, xend = end, yend = -0.5),
               colour = "black", alpha = 1, size = 1, inherit.aes = FALSE) +
  geom_text(data = base_data, aes(x = title, y = -1.5, label = Type), colour = "black",
            alpha = 0.8, size = 2, fontface = "bold", inherit.aes = FALSE)
dev.off()



### ============================================================
### 6th part: K93 and K273 Bivalent or univalent promoters (CMQ)
### ============================================================

### >>> 1. load H3K9me3 and H3K27me3 peaks marked genes(promoters)
# H3K27me3
k27.pks.pro.file <- list.files("../gse124718/results/idr/multi_merged_BAMPE_p_0.05/filtered/peaks_covered_promoter",
                               "covered_promoter500bp_25pect.bed", full.names = T)
k27.pks.pro <- list()
for (i in seq(1, length(k27.pks.pro.file))) {
  k27.pks.pro[[i]] <- read.table(k27.pks.pro.file[i])
  k27.pks.pro[[i]] <- str_split_fixed(k27.pks.pro[[i]][,4], "/", 2)[,1]
  names(k27.pks.pro)[i] <- str_split_fixed(str_split_fixed(k27.pks.pro.file[i], "/", 8)[,8], "_", 2)[,1]
}; rm(i)
names(k27.pks.pro) <- c("Morula", "C2", "C4", "C8", "ICM", "TE")
# H3K9me3
k9.pks.pro.file <- list.files("results/idr/p_0.05/filtered/peaks_covered_promoter", "covered_promoter500bp_25pect.bed", full.names = T)
k9.pks.pro <- list()
for (i in seq(1, length(k9.pks.pro.file))) {
  k9.pks.pro[[i]] <- read.table(k9.pks.pro.file[i])
  k9.pks.pro[[i]] <- str_split_fixed(k9.pks.pro[[i]][,4], "/", 2)[,1]
  names(k9.pks.pro)[i] <- str_split_fixed(str_split_fixed(k9.pks.pro.file[i], "/", 6)[,6], "_", 2)[,1]
}; rm(i)
names(k9.pks.pro) <- c("C4", "C8", "ICM", "Morula", "TE")

### >>> 2. define bivalent or univalent modified genes

no.k9.k27 <- list(C4 = setdiff(setdiff(str_split_fixed(ge.id$V4, "/", 2)[,1], unique(k9.pks.pro$C4)), unique(k27.pks.pro$C4)),
                  C8 = setdiff(setdiff(str_split_fixed(ge.id$V4, "/", 2)[,1], unique(k9.pks.pro$C8)), unique(k27.pks.pro$C8)),
                  Morula = setdiff(setdiff(str_split_fixed(ge.id$V4, "/", 2)[,1], unique(k9.pks.pro$Morula)),unique(k27.pks.pro$Morula)),
                  ICM = setdiff(setdiff(str_split_fixed(ge.id$V4, "/", 2)[,1], unique(k9.pks.pro$ICM)), unique(k27.pks.pro$ICM)),
                  TE = setdiff(setdiff(str_split_fixed(ge.id$V4, "/", 2)[,1], unique(k9.pks.pro$TE)), unique(k27.pks.pro$TE)))

k9.k27.bivalent <- list(C4 = intersect(unique(k27.pks.pro$C4), unique(k9.pks.pro$C4)),
                        C8 = intersect(unique(k27.pks.pro$C8), unique(k9.pks.pro$C8)),
                        Morula = intersect(unique(k27.pks.pro$Morula), unique(k9.pks.pro$Morula)),
                        ICM = intersect(unique(k27.pks.pro$ICM), unique(k9.pks.pro$ICM)),
                        TE = intersect(unique(k27.pks.pro$TE), unique(k9.pks.pro$TE)))

k9.univalent <- list(C4 = setdiff(unique(k9.pks.pro$C4), unique(k27.pks.pro$C4)),
                     C8 = setdiff(unique(k9.pks.pro$C8), unique(k27.pks.pro$C8)),
                     Morula = setdiff(unique(k9.pks.pro$Morula), unique(k27.pks.pro$Morula)),
                     ICM = setdiff(unique(k9.pks.pro$ICM), unique(k27.pks.pro$ICM)),
                     TE = setdiff(unique(k9.pks.pro$TE), unique(k27.pks.pro$TE)))

k27.univalent <- list(C4 = setdiff(unique(k27.pks.pro$C4), unique(k9.pks.pro$C4)),
                      C8 = setdiff(unique(k27.pks.pro$C8), unique(k9.pks.pro$C8)),
                      Morula = setdiff(unique(k27.pks.pro$Morula), unique(k9.pks.pro$Morula)),
                      ICM = setdiff(unique(k27.pks.pro$ICM), unique(k9.pks.pro$ICM)),
                      TE = setdiff(unique(k27.pks.pro$TE), unique(k9.pks.pro$TE)))

### >>> 3. create table for plot
# 4Cell to 8Cell
c4.to.c8 <- data.frame(c8_no = c(length(intersect(no.k9.k27$C4, no.k9.k27$C8)),
                                 length(intersect(k9.univalent$C4, no.k9.k27$C8)),
                                 length(intersect(k27.univalent$C4, no.k9.k27$C8)),
                                 length(intersect(k9.k27.bivalent$C4, no.k9.k27$C8))),
                       c8_k9 = c(length(intersect(no.k9.k27$C4, k9.univalent$C8)),
                                 length(intersect(k9.univalent$C4, k9.univalent$C8)),
                                 length(intersect(k27.univalent$C4, k9.univalent$C8)),
                                 length(intersect(k9.k27.bivalent$C4, k9.univalent$C8))),
                       c8_k27 = c(length(intersect(no.k9.k27$C4, k27.univalent$C8)),
                                  length(intersect(k9.univalent$C4, k27.univalent$C8)),
                                  length(intersect(k27.univalent$C4, k27.univalent$C8)),
                                  length(intersect(k9.k27.bivalent$C4, k27.univalent$C8))),
                       c8_k9_k27 = c(length(intersect(no.k9.k27$C4, k9.k27.bivalent$C8)),
                                     length(intersect(k9.univalent$C4, k9.k27.bivalent$C8)),
                                     length(intersect(k27.univalent$C4, k9.k27.bivalent$C8)),
                                     length(intersect(k9.k27.bivalent$C4, k9.k27.bivalent$C8))))
rownames(c4.to.c8) <- c("c4_no", "c4_k9", "c4_k27", "c4_k9_k27")
c4.to.c8 %>%
  rownames_to_column %>%
  gather(key = 'key', value = 'value', -rowname) -> c4.to.c8.tr 
# 8Cell to Morula
c8.to.mo <- data.frame(mo_no_icm = c(length(intersect(no.k9.k27$C8, no.k9.k27$Morula))/2,
                                     length(intersect(k9.univalent$C8, no.k9.k27$Morula))/2,
                                     length(intersect(k27.univalent$C8, no.k9.k27$Morula))/2,
                                     length(intersect(k9.k27.bivalent$C8, no.k9.k27$Morula))/2),
                       mo_k9_icm = c(length(intersect(no.k9.k27$C8, k9.univalent$Morula))/2,
                                     length(intersect(k9.univalent$C8, k9.univalent$Morula))/2,
                                     length(intersect(k27.univalent$C8, k9.univalent$Morula))/2,
                                     length(intersect(k9.k27.bivalent$C8, k9.univalent$Morula))/2),
                       mo_k27_icm = c(length(intersect(no.k9.k27$C8, k27.univalent$Morula))/2,
                                      length(intersect(k9.univalent$C8, k27.univalent$Morula))/2,
                                      length(intersect(k27.univalent$C8, k27.univalent$Morula))/2,
                                      length(intersect(k9.k27.bivalent$C8, k27.univalent$Morula))/2),
                       mo_k9_k27_icm = c(length(intersect(no.k9.k27$C8, k9.k27.bivalent$Morula))/2,
                                         length(intersect(k9.univalent$C8, k9.k27.bivalent$Morula))/2,
                                         length(intersect(k27.univalent$C8, k9.k27.bivalent$Morula))/2,
                                         length(intersect(k9.k27.bivalent$C8, k9.k27.bivalent$Morula))/2),
                       mo_no_te = c(length(intersect(no.k9.k27$C8, no.k9.k27$Morula))/2,
                                    length(intersect(k9.univalent$C8, no.k9.k27$Morula))/2,
                                    length(intersect(k27.univalent$C8, no.k9.k27$Morula))/2,
                                    length(intersect(k9.k27.bivalent$C8, no.k9.k27$Morula))/2),
                       mo_k9_te = c(length(intersect(no.k9.k27$C8, k9.univalent$Morula))/2,
                                    length(intersect(k9.univalent$C8, k9.univalent$Morula))/2,
                                    length(intersect(k27.univalent$C8, k9.univalent$Morula))/2,
                                    length(intersect(k9.k27.bivalent$C8, k9.univalent$Morula))/2),
                       mo_k27_te = c(length(intersect(no.k9.k27$C8, k27.univalent$Morula))/2,
                                     length(intersect(k9.univalent$C8, k27.univalent$Morula))/2,
                                     length(intersect(k27.univalent$C8, k27.univalent$Morula))/2,
                                     length(intersect(k9.k27.bivalent$C8, k27.univalent$Morula))/2),
                       mo_k9_k27_te = c(length(intersect(no.k9.k27$C8, k9.k27.bivalent$Morula))/2,
                                        length(intersect(k9.univalent$C8, k9.k27.bivalent$Morula))/2,
                                        length(intersect(k27.univalent$C8, k9.k27.bivalent$Morula))/2,
                                        length(intersect(k9.k27.bivalent$C8, k9.k27.bivalent$Morula))/2))
rownames(c8.to.mo) <- c("c8_no", "c8_k9", "c8_k27", "c8_k9_k27")
c8.to.mo %>%
  rownames_to_column %>%
  gather(key = 'key', value = 'value', -rowname) -> c8.to.mo.tr
# Morula to ICM
mo.to.icm <- data.frame(icm_no = c(length(intersect(no.k9.k27$Morula, no.k9.k27$ICM)),
                                   length(intersect(k9.univalent$Morula, no.k9.k27$ICM)),
                                   length(intersect(k27.univalent$Morula, no.k9.k27$ICM)),
                                   length(intersect(k9.k27.bivalent$Morula, no.k9.k27$ICM))),
                        icm_k9 = c(length(intersect(no.k9.k27$Morula, k9.univalent$ICM)),
                                   length(intersect(k9.univalent$Morula, k9.univalent$ICM)),
                                   length(intersect(k27.univalent$Morula, k9.univalent$ICM)),
                                   length(intersect(k9.k27.bivalent$Morula, k9.univalent$ICM))),
                        icm_k27 = c(length(intersect(no.k9.k27$Morula, k27.univalent$ICM)),
                                   length(intersect(k9.univalent$Morula, k27.univalent$ICM)),
                                   length(intersect(k27.univalent$Morula, k27.univalent$ICM)),
                                   length(intersect(k9.k27.bivalent$Morula, k27.univalent$ICM))),
                        icm_k9_k27 = c(length(intersect(no.k9.k27$Morula, k9.k27.bivalent$ICM)),
                                       length(intersect(k9.univalent$Morula, k9.k27.bivalent$ICM)),
                                       length(intersect(k27.univalent$Morula, k9.k27.bivalent$ICM)),
                                       length(intersect(k9.k27.bivalent$Morula, k9.k27.bivalent$ICM))))

rownames(mo.to.icm) <- c("mo_no_icm", "mo_k9_icm", "mo_k27_icm", "mo_k9_k27_icm")
mo.to.icm %>%
  rownames_to_column %>%
  gather(key = 'key', value = 'value', -rowname) -> mo.to.icm.tr
# Morula to TE
mo.to.te <- data.frame(te_no = c(length(intersect(no.k9.k27$Morula, no.k9.k27$TE)),
                                 length(intersect(k9.univalent$Morula, no.k9.k27$TE)),
                                 length(intersect(k27.univalent$Morula, no.k9.k27$TE)),
                                 length(intersect(k9.k27.bivalent$Morula, no.k9.k27$TE))),
                       te_k9 = c(length(intersect(no.k9.k27$Morula, k9.univalent$TE)),
                                 length(intersect(k9.univalent$Morula, k9.univalent$TE)),
                                 length(intersect(k27.univalent$Morula, k9.univalent$TE)),
                                 length(intersect(k9.k27.bivalent$Morula, k9.univalent$TE))),
                       te_k27 = c(length(intersect(no.k9.k27$Morula, k27.univalent$TE)),
                                  length(intersect(k9.univalent$Morula, k27.univalent$TE)),
                                  length(intersect(k27.univalent$Morula, k27.univalent$TE)),
                                  length(intersect(k9.k27.bivalent$Morula, k27.univalent$TE))),
                       te_k9_k27 = c(length(intersect(no.k9.k27$Morula, k9.k27.bivalent$TE)),
                                     length(intersect(k9.univalent$Morula, k9.k27.bivalent$TE)),
                                     length(intersect(k27.univalent$Morula, k9.k27.bivalent$TE)),
                                     length(intersect(k9.k27.bivalent$Morula, k9.k27.bivalent$TE))))

rownames(mo.to.te) <- c("mo_no_te", "mo_k9_te", "mo_k27_te", "mo_k9_k27_te")
mo.to.te %>%
  rownames_to_column %>%
  gather(key = 'key', value = 'value', -rowname) -> mo.to.te.tr
# combine all stage and create data in specific format
sankey.data <- rbind(c4.to.c8.tr, c8.to.mo.tr, mo.to.icm.tr, mo.to.te.tr)
colnames(sankey.data) <- c("source", "target", "value")
sankey.nodes <- data.frame(name=c(as.character(sankey.data$source), as.character(sankey.data$target)) %>% unique())
sankey.data$IDsource=match(sankey.data$source, sankey.nodes$name)-1
sankey.data$IDtarget=match(sankey.data$target, sankey.nodes$name)-1

sankey.data$link_group <- NA
sankey.data[grep("no", sankey.data$source),]$link_group <- "no"
sankey.data[grep("k9", sankey.data$source),]$link_group <- "k9"
sankey.data[grep("k27", sankey.data$source),]$link_group <- "k27"
sankey.data[grep("k9_k27", sankey.data$source),]$link_group <- "k9_k27"

sankey.nodes$node_group <- NA
sankey.nodes[grep("no", sankey.nodes$name),]$node_group <- "a"
sankey.nodes[grep("k9", sankey.nodes$name),]$node_group <- "b"
sankey.nodes[grep("k27", sankey.nodes$name),]$node_group <- "c"
sankey.nodes[grep("k9_k27", sankey.nodes$name),]$node_group <- "d"

### >>> 4. plot
# nodes color
sankey.color ='d3.scaleOrdinal() .domain(["c4_no", "c4_k9","c4_k27", "c4_k9_k27", "c8_no", "c8_k9", "c8_k27", "c8_k9_k27", 
"mo_no_icm", "mo_k9_icm", "mo_k27_icm", "mo_k9_k27_icm", "mo_no_te", "mo_k9_te", "mo_k27_te", "mo_k9_k27_te", 
"icm_no", "icm_k9", "icm_k27", "icm_k9_k27", "te_no", "te_k9", "te_k27", "te_k9_k27"]) 
.range(["#b5ccea", "#ffa2a0", "#0880BA", "#48be86", "#b5ccea", "#ffa2a0", "#0880BA", "#48be86", "#b5ccea", "#ffa2a0", "#0880BA", "#48be86", 
"#b5ccea", "#ffa2a0", "#0880BA", "#48be86", "#b5ccea", "#ffa2a0", "#0880BA", "#48be86", "#b5ccea", "#ffa2a0", "#0880BA", "#48be86"])'
# connection color
connection.color = 'd3.scaleOrdinal() .domain(["no", "k9", "k27", "k9_k27", "a", "b", "c", "d"]) 
.range(["#C8D5F6", "#ffc185", "#cab7d9", "#9ddd90", 
"#0880BA", "#ff8b25", "#9652AF", "#40a940"])'
# plot
sankey.p <- sankeyNetwork(Links = sankey.data, Nodes = sankey.nodes,
                          Source = "IDsource", Target = "IDtarget",
                          Value = "value", NodeID = "name", iterations = 0,
                          sinksRight=FALSE, colourScale=connection.color, nodeWidth=15, fontSize=12, nodePadding=10,
                          width = 350, height = 250, fontFamily ="Arial", LinkGroup="link_group", NodeGroup = "node_group")
sankey.p
saveWidget(sankey.p, file = "results/R/Graphs/sankey_plots/human_embryo_H3K9me3_H3K27me3_p0.05_peaks_covered_promoter500bp_sankey_plot_revised.html")

### >>> 5. GO analysis of genes marked by established and disappeared H3K9me3 and H3K27me3 peaks from morula to icm
# morula to icm gained k9 (univalent)
mo.to.icm.gain.k9.uni <- intersect(k9.univalent$ICM, no.k9.k27$Morula)
mo.to.icm.gain.k9.uni.go <- FullSet.GO("human", mo.to.icm.gain.k9.uni, "mo_to_icm_gain_k9_uni",
                                       "ENSEMBL", "results/R/Tables/GO/bivalent/mo_to_icm_gain_k9_uni")
# morula to icm gained k27 (univalent)
mo.to.icm.gain.k27.uni <- intersect(k27.univalent$ICM, no.k9.k27$Morula)
mo.to.icm.gain.k27.uni.go <- FullSet.GO("human", mo.to.icm.gain.k27.uni, "mo_to_icm_gain_k27_uni",
                                       "ENSEMBL", "results/R/Tables/GO/bivalent/mo_to_icm_gain_k27_uni")
# morula to icm gained k9,k27 (bivalent)
mo.to.icm.gain.k9.k27.bi <- intersect(k9.k27.bivalent$ICM, no.k9.k27$Morula)
mo.to.icm.gain.k9.k27.bi.go <- FullSet.GO("human", mo.to.icm.gain.k9.k27.bi, "mo_to_icm_gain_k9_k27_bi",
                                          "ENSEMBL", "results/R/Tables/GO/bivalent/mo_to_icm_gain_k9_k27_bi")

# morula to icm losed k9 (univalent)
mo.to.icm.lose.k9.uni <- intersect(k9.univalent$Morula, no.k9.k27$ICM)
mo.to.icm.lose.k9.uni.go <- FullSet.GO("human", mo.to.icm.lose.k9.uni, "mo_to_icm_lose_k9_uni",
                                       "ENSEMBL", "results/R/Tables/GO/bivalent/mo_to_icm_lose_k9_uni")
# morula to icm losed k27 (univalent)
mo.to.icm.lose.k27.uni <- intersect(k27.univalent$Morula, no.k9.k27$ICM)
mo.to.icm.lose.k27.uni.go <- FullSet.GO("human", mo.to.icm.lose.k27.uni, "mo_to_icm_lose_k27_uni",
                                        "ENSEMBL", "results/R/Tables/GO/bivalent/mo_to_icm_lose_k27_uni")
# morula to icm losed k9,k27 (bivalent)
mo.to.icm.lose.k9.k27.bi <- intersect(k9.k27.bivalent$Morula, no.k9.k27$ICM)
mo.to.icm.lose.k9.k27.bi.go <- FullSet.GO("human", mo.to.icm.lose.k9.k27.bi, "mo_to_icm_lose_k9_k27_bi",
                                          "ENSEMBL", "results/R/Tables/GO/bivalent/mo_to_icm_lose_k9_k27_bi")

### >>> 6. visualize GO analysis results
# GO term for plot
tar.term.sc.pla <- c("GO:0001892", "GO:0060706", "GO:0001890", "GO:0019827", "GO:0072089", "GO:0048864")
# extract GO term in interest and add annotation
mo.to.icm.gain.k9.uni.go[[1]][tar.term.sc.pla, c(1, 2, 3, 4, 5)] %>% mutate(modi = "H3K9me3", type = "Established") -> mo.to.icm.gain.k9.uni.go.tar
mo.to.icm.gain.k27.uni.go[[1]][tar.term.sc.pla, c(1, 2, 3, 4, 5)] %>% mutate(modi = "H3K27me3", type = "Established") -> mo.to.icm.gain.k27.uni.go.tar
mo.to.icm.gain.k9.k27.bi.go[[1]][tar.term.sc.pla, c(1, 2, 3, 4, 5)] %>% mutate(modi = "Bivalent", type = "Established") -> mo.to.icm.gain.k9.k27.bi.go.tar
mo.to.icm.lose.k9.uni.go[[1]][tar.term.sc.pla, c(1, 2, 3, 4, 5)] %>% mutate(modi = "H3K9me3", type = "Disappeared") -> mo.to.icm.lose.k9.uni.go.tar
mo.to.icm.lose.k27.uni.go[[1]][tar.term.sc.pla, c(1, 2, 3, 4, 5)] %>% mutate(modi = "H3K27me3", type = "Disappeared") -> mo.to.icm.lose.k27.uni.go.tar
rbind(na.omit(mo.to.icm.lose.k27.uni.go.tar),
      data.frame(ID="GO:0060706", Description="cell differentiation involved in embryonic placenta development",
                 GeneRatio=0, BgRatio=0, pvalue=1, modi="H3K27me3", type="Disappeared")) -> mo.to.icm.lose.k27.uni.go.tar
mo.to.icm.lose.k9.k27.bi.go[[1]][tar.term.sc.pla, c(1, 2, 3, 4, 5)] %>% mutate(modi = "Bivalent", type = "Disappeared") -> mo.to.icm.lose.k9.k27.bi.go.tar
rbind(na.omit(mo.to.icm.lose.k9.k27.bi.go.tar),
      data.frame(ID=c("GO:0001892", "GO:0060706", "GO:0001890", "GO:0048864"),
                 Description=c("embryonic placenta development", "cell differentiation involved in embryonic placenta development",
                               "placenta development", "stem cell development"),
                 GeneRatio=rep(0, 4), BgRatio=rep(0, 4), pvalue=rep(1, 4),
                 modi = rep("Bivalent", 4), type = rep("Disappeared", 4))) -> mo.to.icm.lose.k9.k27.bi.go.tar
# data processing
mo.to.icm.go.tar <- rbind(mo.to.icm.gain.k9.uni.go.tar, mo.to.icm.gain.k27.uni.go.tar, mo.to.icm.gain.k9.k27.bi.go.tar,
                          mo.to.icm.lose.k9.uni.go.tar, mo.to.icm.lose.k27.uni.go.tar, mo.to.icm.lose.k9.k27.bi.go.tar)
mo.to.icm.go.tar %>% mutate(GeneRation.n = sapply(GeneRatio, function(x) eval(parse(text=as.character(x)))),
                            BgRatio.n = sapply(BgRatio, function(x) eval(parse(text=as.character(x))))) %>%
  mutate(Ratio = GeneRation.n/BgRatio.n)  %>%
  mutate(Sig = case_when(pvalue >= 0.05 ~ "Non Significant", pvalue < 0.05 ~ "Significant")) -> mo.to.icm.go.tar
mo.to.icm.go.tar[mo.to.icm.go.tar$Ratio == "NaN",10] <- 0
mo.to.icm.go.tar[mo.to.icm.go.tar$ID == "GO:0060706", 2] <- "cell differentiation involved in \n embryonic placenta development"
# plot
pdf("results/R/Graphs/GO/morula_to_icm/flows_of_modification_on_human_embryo_from_morula_to_icm.pdf",
    width = 9.5, height = 5)
mo.to.icm.go.tar %>% mutate(modi = factor(modi, levels = c("H3K9me3", "H3K27me3", "Bivalent")),
                            type = factor(type, levels = c("Established", "Disappeared"))) %>%
  ggplot(aes(x = Ratio, y = Description)) +
  geom_point(aes(color = Sig), size = 4) +
  facet_grid(type ~ modi) +
  scale_color_manual(name = "Significance",
                     values = c("Significant" = "#FF6800",
                                "Non Significant" = "#8F8F8F")) +
  theme_bw() +
  xlab("Observed/Expected") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(colour = "black", fill = "white", size = 1, linetype = "solid"),
        strip.text = element_text(face="plain", colour = "#000000", size = 12, angle = 0, hjust = 0.5),
        axis.title.x = element_text(face="plain", colour = "#000000", size = 12, angle = 0, hjust = 0.5),
        axis.title.y = element_blank(),
        axis.text.x  = element_text(face="plain", colour = "#000000", size = 12, angle = 0, hjust = 1),
        axis.text.y  = element_text(face="plain", colour = "#000000", size = 12, angle = 0, hjust = 1),
        legend.direction = "horizontal",
        legend.position = "bottom")
dev.off()



### ==========================================
### 7th part: Peaks genomic distribution (YHW)
### ==========================================

### >>> 1. Load peaks
# human
hs.k9me3.pks.file <- list.files("RawData/peaks/hs_k9", pattern = ".bed", full.names = T)
hs.k9me3.pks.file <- as.list(hs.k9me3.pks.file)
names(hs.k9me3.pks.file) <- c("4 Cell", "8 Cell", "ICM", "Morula", "TE", "ICM-gained specific", "TE-gained specific")
hs.k9me3.pks <- list()
for (i in seq(1, length(hs.k9me3.pks.file))) {
  hs.k9me3.pks[[i]] <- read.table(hs.k9me3.pks.file[[i]])[, 1:3]
  hs.k9me3.pks[[i]] <- hs.k9me3.pks[[i]][sample(seq(1, nrow(hs.k9me3.pks[[i]])), 5000), ]
  hs.k9me3.pks[[i]] <- toGRanges(hs.k9me3.pks[[i]])
  hs.k9me3.pks[[i]] <- filterChromosomes(hs.k9me3.pks[[i]], organism = "hg38")
}
names(hs.k9me3.pks) <- names(hs.k9me3.pks.file)
# mouse
mm.k9me3.pks.file <- list.files("RawData/peaks/mm_k9", pattern = ".bed", full.names = T)
mm.k9me3.pks.file <- as.list(mm.k9me3.pks.file)
names(mm.k9me3.pks.file) <- c("ICM-gained specific", "TE-gained specific")
# human satellite
hs.sate <- read.table("RawData/peaks/GRCh38_RepeatMasker_Satellite_region.bed", sep = "\t") %>% 
  separate(col = "V4", into = c("class", "family", "subfamily", "locus"), sep = ":")

### >>> 2. Load annotation data
txdb.mm10 <- TxDb.Mmusculus.UCSC.mm10.knownGene
txdb.hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene 

### >>> 3. Annotate peaks
hs.k9me3.pks.anno <- lapply(hs.k9me3.pks.file, annotatePeak, 
                            TxDb = txdb.hg38, tssRegion = c(-3000, 3000), verbose = FALSE, annoDb = "org.Hs.eg.db")
mm.k9me3.pks.anno <- lapply(mm.k9me3.pks.file, annotatePeak, 
                            TxDb = txdb.mm10, tssRegion = c(-3000, 3000), verbose = FALSE, annoDb = "org.Mm.eg.db")
install.packages("BSgenome.Hsapiens.UCSC.hg38")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library("BSgenome.Hsapiens.UCSC.hg38")
all.regions <- hs.k9me3.pks[1:5]
names(all.regions) <- NULL
joined.regions <- do.call(c, all.regions)
pdf("Graphs/quality_control/peaks_distribution/human_H3K9me3_all_stage_peaks_distribution_on_chromosome.pdf", height = 5, width = 5)
kp <- plotKaryotype(plot.type = 2)
kpPlotRegions(kp, data = getGenomeAndMask("hg38")$mask, r0 = 0, r1 = 1, col = "lightgray", border = "lightgray")
invisible(lapply(seq_len(5), function(i) {
  kpPlotRegions(kp, all.regions[[i]],
                r0=(1/5)*(i-1), r1=(1/5)*i)
}))
kpPlotCoverage(kp, data = joined.regions, data.panel = 2, col = "#e12293")
dev.off(); rm(all.regions, joined.regions, kp)

### >>> 4. Visualization
pdf("Graphs/quality_control/peaks_distribution/human_H3K9me3_all_stage_peaks_annotation.pdf", width = 6, height = 10)
p1 <- plotAnnoBar(hs.k9me3.pks.anno[1:5], title = "Distribution of Peaks in Genome")
p2 <- plotDistToTSS(hs.k9me3.pks.anno[1:5], title = "Distribution of Peaks Relative to TSS")
plot_grid(p1, p2, ncol = 1)
dev.off(); rm(p1, p2)
pdf("Graphs/diff_between_ICM_and_TE/human_H3K9me3_peaks_annotation.pdf", width = 6, height = 10)
p1 <- plotAnnoBar(hs.k9me3.pks.anno[-1:-5], title = "Distribution of Peaks in Genome")
p2 <- plotDistToTSS(hs.k9me3.pks.anno[-1:-5], title = "Distribution of Peaks Relative to TSS")
plot_grid(p1, p2, ncol = 1)
dev.off(); rm(p1, p2)
pdf("Graphs/diff_between_ICM_and_TE/mouse_H3K9me3_peaks_annotation.pdf", width = 6, height = 10)
p1 <- plotAnnoBar(mm.k9me3.pks.anno, title = "Distribution of Peaks in Genome")
p2 <- plotDistToTSS(mm.k9me3.pks.anno, title = "Distribution of Peaks Relative to TSS")
plot_grid(p1, p2, ncol = 1)
dev.off(); rm(p1, p2)



### ===============================================
### 8th part: Expression level of repeats in embryo
### ===============================================

### >>> 1. load E-MTAB-3929 gene count
ge.co <- readRDS("/home/cmq/bioinfo/project-ymz/scrna/cell_2016/results/R/Rawdata/count/2016_cell_sce.rds")
# gene count matrix
ge.counts <- as.data.frame(assay(ge.co, "counts"))
# metadata of gene
ge.fe.meta <- rowData(ge.co)
#subset a set of mitochondrial genes
ge.mt <- ge.fe.meta[grep("^MT-", ge.fe.meta$Symbol),]

### >>> 2. load cell metadata
# annotation of cells 
cell.anno <- readRDS("results/R/Rawdata/count/cell_anno.rds") 
# metadata of cell
ge.cell.meta <- as.data.frame(colData(ge.co))

### >>> 3. create seurat object
ge.sar <- CreateSeuratObject(counts = ge.counts, meta.data = ge.cell.meta)
ge.sar[["pect.mt"]] <- PercentageFeatureSet(ge.sar, features = match(ge.mt$Geneid, rownames(ge.sar@assays$RNA@counts)))
ge.sar <- subset(ge.sar, subset = nFeature_RNA > 3000 & nFeature_RNA < 20000 & nCount_RNA < 10000000 & pect.mt < 30)

### >>> 4. subset E5 cells and define cell identity
ge.sar.e5 <- subset(ge.sar, stage == "E5") %>% NormalizeData(ge.sar.e5) %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% ScaleData(features = rownames(ge.sar.e5))
ge.sar.e5 <- RunPCA(ge.sar.e5, features = VariableFeatures(object = ge.sar.e5))
ElbowPlot(ge.sar.e5)
ge.sar.e5 <- FindNeighbors(ge.sar.e5, dims = 1:10)
ge.sar.e5 <- FindClusters(ge.sar.e5, resolution = 0.3)
ge.sar.e5 <- RunUMAP(ge.sar.e5, dims = 1:10) %>% RunTSNE(dims = 1:10)

### >>> 5. visualize lineage marker genes for cell identity definition
VlnPlot(ge.sar.e5, features = c("ENSG00000111704", "ENSG00000204531", "ENSG00000181449",
                                "ENSG00000107485", "ENSG00000179348", "ENSG00000165556"), flip = T, stack = T)
VlnPlot(ge.sar.e5, features = c("ENSG00000066135","ENSG00000127663","ENSG00000107077","ENSG00000186280",
                                "ENSG00000235268","ENSG00000255855","ENSG00000143379",
                                "ENSG00000136169","ENSG00000101945","ENSG00000152455","ENSG00000276043"),
        ncol = 3, flip = T, stack = T)

### >>> 6. define TE cells according to CDX2 expression level
cell.anno.e5 <- as.data.frame(ge.sar.e5@meta.data)
cdx2.tpm.E5c0c2 <- as.data.frame(assay(ge.co, "TPM"))["ENSG00000165556", cell.anno.e5[cell.anno.e5$seurat_clusters %in% c("0","2"), ]$cell_name]

### >>> 7. expression level of repeat elements 
# load raw count of repeat elements
repeat.co.2016 <- readRDS("cell_2016/results/count/all_TEcounts.rds")
repeat.co.2016.e345 <- repeat.co.2016[,c(grep("E3.", colnames(repeat.co.2016)),
                                         grep("E4.", colnames(repeat.co.2016)), grep("E5.", colnames(repeat.co.2016)))]
# load repeat annotation
repeat.anno <- read.table("/home/cx/reference_genome/ucsc/hg38/dna_annotation/GRCh38_RepeatMasker_annotation_repeat_element.saf",
                          stringsAsFactors = F)
colnames(repeat.anno) <- c("id", "chr", "start", "end", "strand")
repeat.co.2016.e345 <- merge(repeat.co.2016.e345, repeat.anno)
rownames(repeat.co.2016.e345) <- repeat.co.2016.e345$id
repeat.co.2016.e345$length <- repeat.co.2016.e345$end - repeat.co.2016.e345$start
# calculate TPM of repeat elements
repeat.tpm.2016.e345 <- CountToTpm(repeat.co.2016.e345[,2:494], repeat.co.2016.e345$length)
repeat.tpm.2016.e345 <- repeat.tpm.2016.e345 %>% rownames_to_column("id") %>% merge(repeat.tpm.2016.e345, repeat.anno)
repeat.tpm.2016.e345$pos.id <- paste(repeat.tpm.2016.e345$chr, repeat.tpm.2016.e345$start, repeat.tpm.2016.e345$end, sep = ":")
repeat.tpm.2016.e345 <- separate(repeat.tpm.2016.e345, "id", c("class", "family", "name", "locus"), sep = "\\.", remove = FALSE)
repeat.tpm.2016.e345 <- repeat.tpm.2016.e345[rowSums(repeat.tpm.2016.e345[,6:498]) > 0, ]

### >>> 8. differential expression analysis of repeat elements between ICM and TE
# load function of edgeR
edgeR <- function(count, meta, g1, g2, lfc, sig, dir){
  library("edgeR")
  dir.create(dir, recursive = T)
  # obtain samples index
  g1.cols <- grep(paste("^", g1, "$", sep = ""), meta$CellType)
  g2.cols <- grep(paste("^", g2, "$", sep = ""), meta$CellType)
  # re-organize count data and meta data
  tar.meta <- meta[c(g1.cols, g2.cols), ]
  print(table(tar.meta$CellType))
  tar.count <- count[, c(g1.cols, g2.cols)]
  # differential expression analysis
  dge <- DGEList(counts = tar.count, samples = tar.meta)
  dge <- calcNormFactors(dge)
  cdr <- scale(colMeans(dge$counts > 0))
  design <- model.matrix(~ cdr + dge$samples$CellType)
  dge <- estimateDisp(dge, design = design)
  fit <- glmQLFit(dge, design = design)
  qlf <- glmQLFTest(fit)
  tt <- topTags(qlf, n = Inf)
  degs.up <- subset(tt$table, logFC >= lfc & PValue <= sig)
  degs.down <- subset(tt$table, logFC <= -lfc & PValue <= sig)
  degs.all <- subset(tt$table, abs(logFC) >= lfc & PValue <= sig)
  degs.all.res <- list(tt$table, degs.all, degs.up, degs.down)
  write.table(tt$table, paste(dir, "/edgeR_t1_", g1, "_Vs_", g2, "_all_gene_results.txt", sep = ""), sep = "\t", quote = F)
  write.table(degs.all, paste(dir, "/edgeR_t2_", g1, "_Vs_", g2, "_all_sig_results.txt", sep = ""), sep = "\t", quote = F)
  write.table(degs.up, paste(dir, "/edgeR_t3_", g1, "_Vs_", g2, "_up_sig_gene_fc", lfc, "_sig", sig, ".txt", sep = ""), sep = "\t", quote = F)
  write.table(degs.down, paste(dir, "/edgeR_t3_", g1, "_Vs_", g2, "_down_sig_gene_fc", lfc, "_sig", sig, ".txt", sep = ""), sep = "\t", quote = F)
  return(degs.all.res)
}
# data preprocessing
repeat.co.2016.e345 <- repeat.co.2016[,c(grep("E3.", colnames(repeat.co.2016), value = T),
                                         grep("E4.", colnames(repeat.co.2016), value = T),
                                         colnames(cdx2.tpm.E5c0c2[, cdx2.tpm.E5c0c2 > 20]),
                                         cell.anno.e5[cell.anno.e5$seurat_clusters == 3, ]$cell_name)]
repeat.co.2016.e345 <- repeat.co.2016.e345[rowSums(repeat.co.2016.e345) > 30,]
# cell metadata preparing
meta.e345 <- data.frame(cell_name = colnames(repeat.co.2016.e345))
meta.e345 %>% mutate(CellType = case_when(cell_name %in% grep("E3.", colnames(repeat.co.2016.e345), value = T) ~ "8Cell",
                                          cell_name %in% grep("E4.", colnames(repeat.co.2016.e345), value = T) ~ "Morula",
                                          cell_name %in% cell.anno.e5[cell.anno.e5$seurat_clusters == 3, ]$cell_name ~ "ICM",
                                          cell_name %in% colnames(cdx2.tpm.E5c0c2[, cdx2.tpm.E5c0c2 > 20]) ~ "TE")) -> meta.e345
meta.e345$CellId <- ave(meta.e345$CellType, meta.e345$CellType, FUN = function(i) paste0(i, '_', seq_along(i)))
rownames(meta.e345) <- meta.e345$CellId
if(all(colnames(repeat.co.2016.e345)==meta.e345$cell_name)){
  colnames(repeat.co.2016.e345) <- meta.e345$CellId
}
# run edger
icm.vs.te.repeat <- edgeR(repeat.co.2016.e345, meta.e345, "ICM", "TE", 1, 0.05,
                          "/home/cmq/bioinfo/project-cmq/embryo_93/results/R/Tables/sc_2016/edgeR/repeat/hs_icm_vs_te_repeat")

### >>> 9. expression level of ICM specific ERVL-MaLR and ERV1
# extract target repeat
re <- c(grep("ERV1\\.", rownames(subset(icm.vs.te.repeat[[2]], logFC < -1.5 & PValue < 0.05)), value = T),
        grep("ERVL-MaLR", rownames(subset(icm.vs.te.repeat[[2]], logFC < -1.5 & PValue < 0.05)), value = T))
# extract target repeat TPM in 8Cell, Morula, ICM and TE
icm.spe.l1l2ervl.tpm <- repeat.tpm.2016.e345[repeat.tpm.2016.e345$id %in% re, c("id", "class", "family", "name", "locus", "pos.id", cell.anno.e5[cell.anno.e5$seurat_clusters == 3, ]$cell_name,
                                                                                                                                    colnames(cdx2.tpm.E5c0c2[, cdx2.tpm.E5c0c2 > 20]),
                                                                                                                                    grep("E[3-4].", colnames(repeat.tpm.2016.e345), value = T))]

colnames(icm.spe.l1l2ervl.tpm)[colnames(icm.spe.l1l2ervl.tpm) %in% cell.anno.e5[cell.anno.e5$seurat_clusters == 3, ]$cell_name] <- "ICM"
colnames(icm.spe.l1l2ervl.tpm)[colnames(icm.spe.l1l2ervl.tpm) %in% colnames(cdx2.tpm.E5c0c2[, cdx2.tpm.E5c0c2 > 20])] <- "TE"
colnames(icm.spe.l1l2ervl.tpm)[grep("E3.", colnames(icm.spe.l1l2ervl.tpm))] <- "8Cell"
colnames(icm.spe.l1l2ervl.tpm)[grep("E4.", colnames(icm.spe.l1l2ervl.tpm))] <- "Morula"
rownames(icm.spe.l1l2ervl.tpm) <- icm.spe.l1l2ervl.tpm$id
# calculate sum of TPM in different repeat family
icm.spe.l1l2ervl.tpm <- rbind(icm.spe.l1l2ervl.tpm,
                              colSums(icm.spe.l1l2ervl.tpm[grep("ERVL-MaLR", rownames(icm.spe.l1l2ervl.tpm)), ]),
                              colSums(icm.spe.l1l2ervl.tpm[grep("ERV1\\.", rownames(icm.spe.l1l2ervl.tpm)), ]))
rownames(icm.spe.l1l2ervl.tpm)[186] <- "ERVL-MaLR"
rownames(icm.spe.l1l2ervl.tpm)[187] <- "ERV1"
# transform data form and visualization
icm.spe.l1l2ervl.tpm[186:187,] %>% rownames_to_column("family") %>% gather(key = "cell", value = "tpm", - family) %>%
  mutate(stage = str_split_fixed(cell, "\\.", 2)[,1]) -> icm.spe.l1l2ervl.tpm

pdf("/home/cmq/bioinfo/project-cmq/embryo_93/results/R/Graphs/icm_vs_te/hs_ICM_specific_expressed_ERVL-MaLR_ERV1_violin_plot.pdf",
    width = 8, height = 4)
icm.spe.l1l2ervl.tpm %>% mutate(stage=factor(stage, levels = c("8Cell", "Morula", "ICM", "TE"))) %>%
  ggplot(aes(x = stage, y = log2(tpm + 0.1), fill = stage)) +
  geom_violin() +
  facet_wrap(. ~ family,scales = "free") +
  ylab("log2(Sum of TPM)") +
  theme_bw() +
  theme(axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.title.x = element_blank(),
        axis.text.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.text.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        legend.text = element_text(face = "plain", colour = "#000000", size = 8, angle = 0),
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white")) +
  stat_compare_means(method = "t.test", na.rm = T, label.y = c(12, 13, 14),
                     comparisons = list(c("8Cell", "TE"), c("Morula", "TE"), c("ICM", "TE")), size = 4)
dev.off()

### >>> 10. expression level of TE specific L1 and L2
# extract target repeat  
re <- c(grep("\\.L1\\.", rownames(subset(icm.vs.te.repeat[[2]], logFC > 1.5 & PValue < 0.05)), value = T),
        grep("\\.L2\\.", rownames(subset(icm.vs.te.repeat[[2]], logFC > 1.5 & PValue < 0.05)), value = T))
# extract target repeat TPM in 8Cell, Morula, ICM and TE
te.spe.l1l2ervl.tpm <- repeat.tpm.2016.e345[repeat.tpm.2016.e345$id %in% re, c("id", "class", "family", "name", "locus", "pos.id",
                                                                             cell.anno.e5[cell.anno.e5$seurat_clusters == 3, ]$cell_name,
                                                                             colnames(cdx2.tpm.E5c0c2[, cdx2.tpm.E5c0c2 > 20]),
                                                                             grep("E[3-4].", colnames(repeat.tpm.2016.e345), value = T))]
colnames(te.spe.l1l2ervl.tpm)[colnames(te.spe.l1l2ervl.tpm) %in% cell.anno.e5[cell.anno.e5$seurat_clusters == 3, ]$cell_name] <- "ICM"
colnames(te.spe.l1l2ervl.tpm)[colnames(te.spe.l1l2ervl.tpm) %in% colnames(cdx2.tpm.E5c0c2[, cdx2.tpm.E5c0c2 > 20])] <- "TE"
colnames(te.spe.l1l2ervl.tpm)[grep("E3.", colnames(te.spe.l1l2ervl.tpm))] <- "8Cell"
colnames(te.spe.l1l2ervl.tpm)[grep("E4.", colnames(te.spe.l1l2ervl.tpm))] <- "Morula"
rownames(te.spe.l1l2ervl.tpm) <- te.spe.l1l2ervl.tpm$id
# calculate mean and variance of TPM in order to select specific repeats
te.spe.l1l2ervl.tpm$TE_ave <- rowMeans(te.spe.l1l2ervl.tpm[,55:142])
te.spe.l1l2ervl.tpm$TE_var <- rowVars(as.matrix(te.spe.l1l2ervl.tpm[,55:142]))
te.spe.l1l2ervl.tpm <- te.spe.l1l2ervl.tpm[te.spe.l1l2ervl.tpm$TE_ave > 1.5 & te.spe.l1l2ervl.tpm$TE_var < 60, c(-1:-6, -360:-361)]
# calculate sum of TPM in different repeat family
te.spe.l1l2ervl.tpm <- rbind(te.spe.l1l2ervl.tpm,
                             colSums(te.spe.l1l2ervl.tpm[grep("\\.L1\\.", rownames(te.spe.l1l2ervl.tpm)), ]),
                             colSums(te.spe.l1l2ervl.tpm[grep("\\.L2\\.", rownames(te.spe.l1l2ervl.tpm)), ]))
rownames(te.spe.l1l2ervl.tpm)[330] <- "L1"
rownames(te.spe.l1l2ervl.tpm)[331] <- "L2"
# transform data form and visualization
te.spe.l1l2ervl.tpm[330:331,] %>% rownames_to_column("family") %>% gather(key = "cell", value = "tpm", - family) %>%
  mutate(stage = str_split_fixed(cell, "\\.", 2)[,1]) -> te.spe.l1l2ervl.tpm

pdf("/home/cmq/bioinfo/project-cmq/embryo_93/results/R/Graphs/icm_vs_te/hs_TE_specific_expressed_L1_L2_violin_plot.pdf",
    width = 8, height = 4)
te.spe.l1l2ervl.tpm %>% mutate(stage=factor(stage, levels = c("8Cell", "Morula", "ICM", "TE"))) %>%
  ggplot(aes(x = stage, y = log2(tpm + 0.1), fill = stage)) +
  geom_violin() +
  facet_wrap(. ~ family,scales = "free") +
  ylab("log2(Sum of TPM)") +
  theme_bw() +
  theme(axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.title.x = element_blank(),
        axis.text.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.text.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        legend.text = element_text(face = "plain", colour = "#000000", size = 8, angle = 0),
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white")) +
  stat_compare_means(method = "t.test", na.rm = T, label.y = c(12, 13, 14),
                     comparisons = list(c("8Cell", "ICM"), c("Morula", "ICM"), c("ICM", "TE")), size = 4)
dev.off()



### ===========================================
### 9th part: Define stage-specific peaks (YHW)
### ===========================================

### >>> 1. Load data (ChIP-seq signal on large H3K9me3 domains)
k9.domain <- read.table("RawData/deeptools/multibw/all_stage_large_domain_3kb/H3K9me3_all_stage_large_domain_3kb_readNor.tab", header = T)
colnames(k9.domain)[-1:-3] <- c("4 Cell input", "4 Cell", "8 Cell input", "8 Cell",
                                "ICM input", "ICM", "Morula input", "Morula", "TE input", "TE")
k9.domain %>% mutate(`4 Cell Ratio` = (`4 Cell`+0.1)/(`4 Cell input`+0.1),
                     `8 Cell Ratio` = (`8 Cell`+0.1)/(`8 Cell input`+0.1),
                     `Morula Ratio` = (Morula+0.1)/(`Morula input`+0.1),
                     `ICM Ratio` = (ICM+0.1)/(`ICM input`+0.1),
                     `TE Ratio` = (TE+0.1)/(`TE input`+0.1),
                     pos.id = paste(chr, ":", start, "-", end, sep = "")) -> k9.domain
k9.domain <- k9.domain[!duplicated(k9.domain$pos.id),]
rownames(k9.domain) <- k9.domain$pos.id

### >>> 2. Quantile-normalization
k9.domain.quan <- as.data.frame(normalize.quantiles(as.matrix(k9.domain[, grep("Ratio", colnames(k9.domain))])), row.names = rownames(k9.domain))
colnames(k9.domain.quan) <- colnames(k9.domain)[grep("Ratio", colnames(k9.domain))]
boxplot(log2(k9.domain.quan+0.1))

### >>> 3. Identify stage-specific domains
library(foreach)
library(doParallel)
Entropy <- function(df, ncore){
  registerDoParallel(ncore)
  entropy <- function(vec){
    pg <- vec/sum(vec)
    pg <- pg[pg>0]
    hg <- sum(-(pg)*log2(pg))
    q <- hg - log2(pg)
    res <- list(hg, q)
    names(res) <- c("Hg", "Qg")
    return(res)
  }
  r <- foreach(row = seq(1, nrow(df)), .combine = rbind) %dopar% { entropy(df[row, ]) }
  rownames(r) <- rownames(df)
  return(r)
}
k9.domain.quan.entro <- Entropy(k9.domain.quan, 20)
k9.domain.quan.entro <- as.data.frame(k9.domain.quan.entro)
k9.domain.quan.entro$Hg <- unlist(k9.domain.quan.entro$Hg)
summary(k9.domain.quan.entro$Hg)
table(k9.domain.quan.entro$Hg>2.26)

### >>> 4. Visualization
pdf("Graphs/stage_specfic_3kb_domains/H3K9me3_stage_specfic_3kb_domains_fold_change_heatmap.pdf", height = 10, width = 10)
Heatmap(as.matrix(log2(k9.domain.quan[k9.domain.quan.entro$Hg<1.5,]+0.1)), name = "Log2(Enrichment)",
        col = colorRamp2(c(-1, 2, 5), c("#015a9e", "#ffffff", "#c60d0d")),
        cluster_rows = T, show_row_dend = T, row_dend_side = "left", row_dend_width = unit(1, "cm"), 
        clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
        cluster_columns = F, show_column_dend = T, column_dend_height = unit(1, "cm"), column_dend_side = "top",
        clustering_distance_columns = "euclidean", clustering_method_columns = "complete",
        show_column_names = T, column_names_rot = 45, column_names_side = "bottom", column_title = "H3K9me3",
        show_row_names = F, row_names_rot = 0, row_names_side = "right", row_split = 5, border = T,
        width = unit(6, "cm"), height = unit(8, "cm"), use_raster = T,
        top_annotation = HeatmapAnnotation(DataType = anno_block(gp = gpar(fill = c("#FB9A99")),
                                                                 labels = c("H3K9me3"),
                                                                 labels_gp = gpar(col = "white", fontsize = 10),
                                                                 height = unit(4, "mm"))))
dev.off()

### >>> 5. Output specific locus id
write.table(as.data.frame(rownames(k9.domain.quan[k9.domain.quan.entro$Hg<1.5,])[row_order(p)[[1]]]), 
            "Tables/specific_domains/8cell_specific_3kb_domains.txt", 
            row.names = F, col.names = F, quote = F)
write.table(as.data.frame(rownames(k9.domain.quan[k9.domain.quan.entro$Hg<1.5,])[row_order(p)[[2]]]), 
            "Tables/specific_domains/te_specific_3kb_domains.txt", 
            row.names = F, col.names = F, quote = F)
write.table(as.data.frame(rownames(k9.domain.quan[k9.domain.quan.entro$Hg<1.5,])[row_order(p)[[3]]]), 
            "Tables/specific_domains/morula_specific_3kb_domains.txt", 
            row.names = F, col.names = F, quote = F)
write.table(as.data.frame(rownames(k9.domain.quan[k9.domain.quan.entro$Hg<1.5,])[row_order(p)[[4]]]), 
            "Tables/specific_domains/4cell_specific_3kb_domains.txt", 
            row.names = F, col.names = F, quote = F)
write.table(as.data.frame(rownames(k9.domain.quan[k9.domain.quan.entro$Hg<1.5,])[row_order(p)[[5]]]), 
            "Tables/specific_domains/icm_specific_3kb_domains.txt", 
            row.names = F, col.names = F, quote = F)
write.table(as.data.frame(rownames(k9.domain.quan[k9.domain.quan.entro$Hg>2.26,])), 
            "Tables/specific_domains/no_specific_3kb_domains.txt", 
            row.names = F, col.names = F, quote = F)

### >>> 6. Calculate the number and type of marked element
subfamily.anno <- read.table("Tables/specific_domains/subfamily_anno.txt", sep = "\t")
colnames(subfamily.anno) <- c("total", "subfamily")
stage.spe.domain.marked <- list()
for (i in seq(1, length(list.files("Tables/specific_domains", pattern = "repeat.txt")))) {
  tmp <- read.table(list.files("Tables/specific_domains", pattern = "repeat.txt", full.names = T)[i], sep = "\t")
  tmp$pie <- NA
  tmp$pie[grep("ENSG", tmp$V4)] <- "gene"
  tmp$pie[-grep("ENSG", tmp$V4)] <- "repeat"
  tmp <- cbind(tmp, str_split_fixed(tmp$V4, ":", 4))
  colnames(tmp) <- c("chr", "start", "end", "id", "type", "class", "family", "subfamily", "locus")
  tmp[grep("ENSG", tmp$id), 6:9] <- rep("gene", 4)
  stage.spe.domain.marked[[i]] <- tmp
}; rm(i, tmp)
names(stage.spe.domain.marked) <- str_split_fixed(list.files("Tables/specific_domains", pattern = "repeat.txt"), "_", 2)[, 1]
library("ggplot2")
library("moonBook")
library("webr")
for (i in seq(1, length(stage.spe.domain.marked))) {
  pd.family <- as.data.frame(table(stage.spe.domain.marked[[i]]$subfamily)) %>% 
    full_join(subfamily.anno, by = c("Var1" = "subfamily")) %>%
    filter(total > 200, Freq > 50) %>%
    mutate(ratio = Freq/total) %>% filter(Var1 != "gene") %>% top_n(30, wt = ratio)
  pdf(paste("Graphs/", names(stage.spe.domain.marked)[i], "_specific_3kb_domain_marked_gene_repeat.pdf", sep = ""), height = 5, width = 5)
  PieDonut(subset(stage.spe.domain.marked[[i]], subfamily %in% c(as.character(pd.family$Var1), "gene")),
           aes(pies = family, donuts = subfamily), addDonutLabel = TRUE, showRatioDonut = F,
           showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.01), color = "#ffffff",
           ratioByGroup = F, labelposition = 1, pieLabelSize = 2, donutLabelSize = 2)
  dev.off()
}; rm(i)
# GO
stage.spe.domain.marked.no.go <- FullSet.GO("human", grep("ENSG", str_split_fixed(stage.spe.domain.marked$no$id, ":", 2)[,1], value = T), 
                                            "all_stage_shared_H3K9me3_domains_covered_genes", "ENSEMBL", "Tables/GO/all_stage_shared_H3K9me3_domains")
pdf("Graphs/stage_specfic_3kb_domains/no_specific_3kb_domain_marked_gene_GO.pdf", height = 2.5, width = 6)
subset(stage.spe.domain.marked.no.go[[1]][, 1:7], ID %in% c("GO:0048562", "GO:0048333", "GO:0001707", "GO:0048568", "GO:0048704", "GO:0048562")) %>%
  mutate(Log10Pvalue = -log10(pvalue)) %>%
  ggplot(aes(x = Log10Pvalue, y = fct_reorder(Description, Log10Pvalue))) + 
  labs(x = "-Log10.pvalue", y = "GO biological process") +
  geom_bar(stat = "identity", fill = "#3498db", color = "#ecf0f1", width = 0.75) + 
  theme_classic() + 
  theme(axis.title.x = element_text(face="plain", colour = "#1D1D1E", size = 10, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#1D1D1E", size = 10, angle = 90), 
        axis.text.x  = element_text(face="plain", colour = "#1D1D1E", size = 10, angle = 0),
        axis.text.y  = element_text(face="plain", colour = "#1D1D1E", size = 10, angle = 0))
dev.off()

### >>> 7. Calculate the number and type of marked element (specific peaks)
stage.spe.pks.marked <- list()
for (i in seq(1, length(list.files("Tables/specific_peaks", pattern = "repeat.txt")))) {
  tmp <- read.table(list.files("Tables/specific_peaks", pattern = "repeat.txt", full.names = T)[i], sep = "\t")
  tmp$pie <- NA
  tmp$pie[grep("ENSG", tmp$V4)] <- "gene"
  tmp$pie[-grep("ENSG", tmp$V4)] <- "repeat"
  tmp <- cbind(tmp, str_split_fixed(tmp$V4, ":", 4))
  colnames(tmp) <- c("chr", "start", "end", "id", "type", "class", "family", "subfamily", "locus")
  tmp[grep("ENSG", tmp$id), 6:9] <- rep("gene", 4)
  stage.spe.pks.marked[[i]] <- tmp
}; rm(i, tmp)
names(stage.spe.pks.marked) <- str_split_fixed(list.files("Tables/specific_peaks", pattern = "repeat.txt"), "_", 2)[, 1]
library("ggplot2")
library("moonBook")
library("webr")
for (i in seq(1, length(stage.spe.pks.marked))) {
  pd.family <- as.data.frame(table(stage.spe.pks.marked[[i]]$subfamily)) %>% 
    full_join(subfamily.anno, by = c("Var1" = "subfamily")) %>%
    filter(total > 200, Freq > 50) %>%
    mutate(ratio = Freq/total) %>% filter(Var1 != "gene") %>% top_n(30, wt = ratio)
  pdf(paste("Graphs/", names(stage.spe.pks.marked)[i], "_specific_3kb_peaks_marked_gene_repeat.pdf", sep = ""), height = 5, width = 5)
  PieDonut(subset(stage.spe.pks.marked[[i]], subfamily %in% c(as.character(pd.family$Var1), "gene")),
           aes(pies = family, donuts = subfamily), addDonutLabel = TRUE, showRatioDonut = F,
           showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.01), color = "#ffffff",
           ratioByGroup = F, labelposition = 1, pieLabelSize = 2, donutLabelSize = 2)
  dev.off()
}; rm(i)

### >>> 8. GO analysis of H3K9me3 stage specific peaks marked genes
pks.co.ge.file <- list.files("/home/yhw/bioinfo/project-mine/Embryo.93/chipseq/hs_k9me3/results/bedtools/specific_large_k9_3kb_domain", 
                             pattern = "pcg10.txt", full.names = T)
FullSet.GO <- edit(FullSet.GO)
clusterProfiler.GO <- edit(clusterProfiler.GO)
pks.co.ge <- list()
pks.co.ge.go <- list()
for (i in seq(1, length(pks.co.ge.file))) {
  pks.co.ge[[i]] <- read.table(pks.co.ge.file[i], header = F)
  pks.co.ge[[i]] <- str_split_fixed(pks.co.ge[[i]]$V4, "/", 2)[,1]
  bname <- str_split_fixed(str_split_fixed(pks.co.ge.file[i], "/", 12)[,12], "_", 2)[,1]
  od <- paste("results/R/Tables/GO/stage_specific_pks", bname, sep = "/")
  pks.co.ge.go[[i]] <- FullSet.GO("human", pks.co.ge[[i]], bname, "ENSEMBL", od)
  names(pks.co.ge.go)[i] <- bname
}; rm(bname, od, i)
