# =============
# Configuration
# =============
setwd("/home/cmq/bioinfo/project-cmq/embryo_93")
# plot
library("ggplot2")
library("ggpubr")
library("ggrepel")
library("ComplexHeatmap")
library("circlize")
# data processing
library("dplyr")
library("tidyr")
library("tibble")


# ================================================================
# ICM H3K9me3 and H3K4me3 bivalent peaks marked repeats annotation 
# ================================================================

### >>> 1. load K9-K4 Bi-SVA marked repeats annotation files
c8.icm.93.43.file <- list.files("results/tables/naive_icm", "Retroposon_0.5pect_anno.txt", full.names = T) 
c8.icm.93.43.pks <- list()
for (i in seq(1, length(c8.icm.93.43.file))) {
        c8.icm.93.43.pks[[i]] <- read.table(c8.icm.93.43.file[i])
        colnames(c8.icm.93.43.pks[[i]]) <- c8.icm.93.43.pks[[i]][nrow(c8.icm.93.43.pks[[i]]),]
        c8.icm.93.43.pks[[i]] <- c8.icm.93.43.pks[[i]][-nrow(c8.icm.93.43.pks[[i]]),]
        c8.icm.93.43.pks[[i]] %>% mutate(Age = as.numeric(Age), Covered_Num = as.numeric(Covered_Num), 
                                         All_Num = as.numeric(All_Num), Length = as.numeric(Length), 
                                         Ratio = as.numeric(Ratio)) -> c8.icm.93.43.pks[[i]]
        names(c8.icm.93.43.pks)[i] <- str_split_fixed(c8.icm.93.43.file[i], "/", 4)[,4]
}; rm(i)

### >>> 2. visualization via bubble plot
pdf("results/R/Graphs/repeat_age/mya/human_embryo_ICM_H3K9me3_H3K4me3_shared_p0.05_p0.01_peaks_covered_subfamily_age_50pect_Ratio_MYA.pdf", 
   height = 3, width = 5)
subset(c8.icm.93.43.pks[[3]], Covered_Num >= 5 & Length > 50 & Family != "ERVL?") %>% 
        mutate(alpha = case_when(Length >= 50 & Length <= 400 ~ "a", 
                                 Length > 400 & Length <= 600 ~ "b", 
                                 Length > 600 & Length <= 800 ~ "c", 
                                 Length > 800 ~ "d")) %>% 
        mutate(size = case_when(Covered_Num <= 50 ~ "a", 
                                Covered_Num > 50 & Covered_Num <= 100 ~ "b",
                                Covered_Num > 100 & Covered_Num <= 200 ~ "c")) %>%  
        ggplot(aes(x = Age/2.2, y = Ratio)) +
        geom_point(aes(size = size, color = Family)) + 
        scale_size_manual(values = c("a" = 3, "b" = 5, "c" = 6), 
                          labels = c("<=50", "50~100", "100~200"), 
                          limits = c("a", "b", "c")) + 
        labs(color = "Family", size = "Number", x = "Million years ago (Mya)", y = "Ratio(%)") + 
        #xlim(c(-2, 320)) +
        geom_text_repel(aes(label = Subfamily), 
                        data = subset(c8.icm.93.43.pks[[3]], Covered_Num >= 10 & Length > 50 & Family != "ERVL?"), 
                        size = 3, color = "#000000") + 
        #ylim(0,40) + 
        theme_bw() + 
        theme(axis.title.x = element_text(face="plain", colour = "#000000", size = 12, angle = 0, hjust = 1),
              axis.title.y = element_text(face="plain", colour = "#000000", size = 12, angle = 90),  
              axis.text.x  = element_text(face="plain", colour = "#000000", size = 12, angle = 0, hjust = 1),
              axis.text.y  = element_text(face="plain", colour = "#000000", size = 12, angle = 0, hjust = 1),
              panel.grid.minor = element_blank(), 
              legend.box = "horizontal")
dev.off()


# ======================================
# GO analysis of K9-K4 Bi-SVA near genes 
# ======================================

### >>> 1. load GREAT results
sva.near.ge.go.file <- list.files("results/tables/naive_icm/enrich/", "tar.txt", full.names = T)
sva.near.ge.go <- list()
for (i in seq(1,length(sva.near.ge.go.file))) {
    sva.near.ge.go[[i]] <- read.table(sva.near.ge.go.file[i], sep = "\t", quote = "")
    colnames(sva.near.ge.go[[i]]) <- c("ID", "Description", "PValue")
    sva.near.ge.go[[i]]$Log10PValue <- -log10(sva.near.ge.go[[i]]$PValue)
}; rm(i)
names(sva.near.ge.go) <- c("HERVH", "Bivalent", "Non-Bivalent")

### >>> 2. extract GO term in interest
sva.near.ge.go.tar <- left_join(sva.near.ge.go$Bivalent, sva.near.ge.go$`Non-Bivalent`, by = "ID")[, c(2,4,7)]
sva.near.ge.go.tar$Log10PValue.HERVH <- c(rep(NA, 8), sva.near.ge.go$HERVH$Log10PValue)
sva.near.ge.go.tar$Log10PValue.LTR5_Hs <- rep(NA, 9)
sva.near.ge.go.tar$Log10PValue.LTR7B <- rep(NA, 9)

### >>> 3. visualization
pdf("results/R/Graphs/icm_naive/human_embryo_H3K9me3_H3K4me3_bivalent_SVA_nearby_genes_in_ICM_GO_related_to_DNA_repair.pdf", 
    height = 3.5, width = 8)
sva.near.ge.go.tar %>%
    gather(key = "Desc", value = "pvalue", -Description.x) %>%
    mutate(Type = case_when(Desc == "Log10PValue.x" ~ "Bi.SVA",
                            Desc == "Log10PValue.y" ~ "Non.Bi.SVA", 
                            Desc == "Log10PValue.HERVH" ~ "Bi.HERVH-int",
                            Desc == "Log10PValue.LTR5_Hs" ~ "Bi.LTR5_Hs",
                            Desc == "Log10PValue.LTR7B" ~ "Bi.LTR7B"),
           Log10PValue = case_when(pvalue %in% NA ~ 0.1,
                                   ! pvalue %in% NA ~ pvalue)) %>%
    filter(Description.x != "DNA repair") %>%
    mutate(Type = factor(Type, levels = c("Bi.SVA", "Bi.HERVH-int", "Non.Bi.SVA", "Bi.LTR5_Hs", "Bi.LTR7B"))) %>%
    ggplot(aes(x = Type, y = fct_reorder(Description.x, Log10PValue, .desc = F), color = Type)) + 
    geom_point(aes(size = Log10PValue)) + 
    scale_size_continuous(range = c(1.3,8)) + 
    theme_classic() + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(), 
          axis.text.x  = element_text(face="plain", colour = "#000000", size = 10, angle = 45, hjust = 1),
          axis.text.y  = element_text(face="plain", colour = "#000000", size = 10, angle = 0, hjust = 1))
dev.off()



