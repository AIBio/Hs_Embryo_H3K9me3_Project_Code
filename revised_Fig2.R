################################################# ~ Good Luck ~ #########################################################
###### >>>>>>                       TITLE: H3K9me3 Project Revised Figure2 analysis                         <<<<<< ######

# - Figure 2 content:
# - 1st part: Global setting
# - 2nd part: Save Data
# - 3rd part: 8Cell Enhancer Marked Repeat Annotation
# - 4th part: Plot the results of TOBIAS
# - 5th part: 8Cell Enhancer Motif Analysis Visualization
# - 6th part: REs and PEs targeted genes
# - 7th Part: REs and PEs Enrichment analysis on repeats
# - 8th Part: Detect selection on REs and PEs
# - 9th Part: Hic data analysis



### ========================
### 1st part: Global setting
### ========================
setwd("/home/cmq/bioinfo/project-cmq/embryo_93")
#save.image("scripts/revised_Fig2.RData")
#load("scripts/revised_Fig2.RData")



### ===================
### 2nd part: Save Data
### ===================
.libPaths("/home/cmq/software/anaconda3/envs/rs.p37.r4111/lib/R/library")
library(devtools)
library(BiocManager)
### >>> 1. CRAN packages
cran.pks <- c("forcats", "dplyr", "tidyr", "tidyverse", "stringr", "ggplot2", "ggrepel", "circlize", "cowplot", "ggpubr",
              "gkmSVM")
for (pks in cran.pks) {
  if (!require(pks, character.only = TRUE)) {
    install.packages(pks)
    library(pks, character.only = T)
  } else {
    library(pks, character.only = T)
  }
}
### >>> 2. Bioconductor packages
bioc.pks <- c("ComplexHeatmap", "BSgenome.Hsapiens.UCSC.hg38.masked", "GenomicAlignments", "IRanges")
for (pks in bioc.pks) {
  if (!require(pks, character.only = TRUE)) {
    BiocManager::install(pks)
    library(pks, character.only = T)
  } else {
    library(pks, character.only = T)
  }
}
### >>> 3. GitHub packages
if (!require("webr", character.only = TRUE)) {
  devtools::install_github("cardiomoon/webr")
  library("webr", character.only = T)
} else {
  library("webr", character.only = T)
}
if (!require("moonBook", character.only = TRUE)) {
  devtools::install_github("cardiomoon/moonBook")
  library("moonBook", character.only = T)
} else {
  library("moonBook", character.only = T)
}
### >>> 4. Functions
CountToTpm <- function(count, length){
  for (i in seq(1, ncol(count))){
    numer <- log(count[, i]) - log(length)
    denom <- log(sum(exp(numer)))
    count[, i] <- exp(numer - denom + log(1e6))
  }
  return(count)
}
DbiIDtrans <- function(genelist, intype, outtype, species){
  # Load the required packages
  # load the packages
  library("org.Hs.eg.db")
  library("org.Mm.eg.db")
  library("org.Rn.eg.db")
  library("AnnotationDbi")
  # Judge the species
  if (species == "human"){
    database <- "org.Hs.eg.db"
    database_short <- "hsa"
  } else if (species == "mouse"){
    database <- "org.Mm.eg.db"
    database_short <- "mmu"
  } else if (species == "rat"){
    database <-  "org.Rn.eg.db"
    database_short <- "rno"
  } else {
    print("You may need to define a new species dataset in the function")
  }
  # Transforming(keytypes list: keytypes(org.Hs.eg.db))
  if (species == "human"){
    newgenelist <- mapIds(x = org.Hs.eg.db, keys = unique(as.character(genelist)),
                          column = outtype, keytype = intype, multiVals = "first")
  } else if (species == "mouse"){
    newgenelist <- mapIds(x = org.Mm.eg.db, keys = unique(as.character(genelist)),
                          column = outtype, keytype = intype, multiVals = "first")
  } else if (species == "rat"){
    newgenelist <- mapIds(x = org.Rn.eg.db, keys = unique(as.character(genelist)),
                          column = outtype, keytype = intype, multiVals = "first")
  } else {
    print("You may need to define a new species dataset in the function")
  }
  # Retuen the results
  return(newgenelist)
}
ParseSNPdensity <- function(working.dir){
  library("dplyr")
  library("gtools")
  old <- getwd()
  setwd(working.dir)
  # load data
  count.file <- list.files("counts/", ".bed", full.names = T)
  count.file <- mixedsort(count.file)
  print(count.file)
  # load region
  gr <- read.table("final_regions_for_count_SNP.bed", sep = "\t")
  gr$region.id <- paste("region", 1:nrow(gr), sep = "_")
  gr$pos.id <- paste(gr$V1, ":", gr$V2, "-", gr$V3, sep = "")
  # processing
  raw <- matrix(ncol = 41)
  for (i in seq(1, length(count.file), 1)) {
    tmp <- read.table(count.file[i])
    if (nrow(tmp) != 41) {
      next  
    }else{
      colnames(tmp) <- c("chr", "start", "end", "count1", "count2", "total", "ratio")
      rownames(tmp) <- paste("bins", 1:41, sep = "_")
      tmp$ratio <- as.numeric(tmp$ratio)
      tmp <- t(tmp)
      raw <- rbind(raw, tmp[7,])  
    }
  }
  raw <- raw %>% na.omit() %>% as.data.frame(row.names = paste("region", 1:length(count.file), sep = "_"))
  raw <- sapply(raw, as.numeric) %>% as.data.frame(row.names = rownames(raw))
  res <- list(raw, gr)
  names(res) <- c("value", "region")
  setwd(old)
  return(res)
}
ParseRate <- function(working.dir){
  library("dplyr")
  library("gtools")
  old <- getwd()
  setwd(working.dir)
  # load data
  count.file <- list.files("rate/", ".bed", full.names = T)
  count.file <- mixedsort(count.file)
  print(count.file)
  # load region
  gr <- read.table("final_regions_for_count_SNP.bed", sep = "\t")
  gr$region.id <- paste("region", 1:nrow(gr), sep = "_")
  gr$pos.id <- paste(gr$V1, ":", gr$V2, "-", gr$V3, sep = "")
  # processing
  raw <- matrix(ncol = 41)
  for (i in seq(1, length(count.file), 1)) {
    tmp <- read.table(count.file[i])
    if (nrow(tmp) != 41) {
      next  
    }else{
      colnames(tmp) <- c("chr", "start", "end", "rate")
      rownames(tmp) <- paste("bins", 1:41, sep = "_")
      tmp$rate <- as.numeric(tmp$rate)
      tmp <- t(tmp)
      raw <- rbind(raw, tmp[4,])  
    }
  }
  raw <- raw %>% na.omit() %>% as.data.frame(row.names = paste("region", 1:length(count.file), sep = "_"))
  raw <- sapply(raw, as.numeric) %>% as.data.frame(row.names = rownames(raw))
  res <- list(raw, gr)
  names(res) <- c("value", "region")
  setwd(old)
  return(res)
}



### =================================================
### 3rd part: 8Cell Enhancer Marked Repeat Annotation
### =================================================

### >>> 1. Set output directory
outdir <- file.path(getwd(), "results/R/Graphs/revised_Fig2")
if (! dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}

### >>> 2. load enhancer annotation files
ehc.anno.file <- list.files("results/bedtools/revised_8c_enhancer/AnnoRepeat", "Retroposon_0.5pect_anno.txt", full.names = T)
ehc.anno <- list()
for (i in seq(1, length(ehc.anno.file))) {
  ehc.anno[[i]] <- read.table(ehc.anno.file[i], stringsAsFactors = F, header = T)
  names(ehc.anno)[i] <- str_split_fixed(str_split_fixed(ehc.anno.file[i], "/", 5)[,5], "_covered", 2)[,1]
}; rm(i)

### >>> 3. visualization via donut plot
for (i in seq(1, length(ehc.anno))) {
  pd <- subset(ehc.anno[[i]], ! Class %in% c("Satellite", "Simple_repeat", "tRNA", "Low_complexity") & 
               Covered_Num >= 10 & All_Num >= 30) %>% top_n(30, wt = Ratio)
  pd <- data.frame(Subfamily = rep(pd$Subfamily, pd$Covered_Num), 
                   Family = rep(pd$Family, pd$Covered_Num))
  pdf(paste(outdir, "/PieDonut_plot_to_show_", names(ehc.anno)[i], "_marked_repeat_anno_top30_subfamily.pdf", sep = ""), 
      height = 5, width = 5)
  PieDonut(pd, aes(pies = Family, donuts = Subfamily), addDonutLabel = TRUE, showRatioDonut = T,
           showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.01), color = "#ffffff",
           ratioByGroup = F, labelposition = 1, pieLabelSize = 2, donutLabelSize = 2)
  dev.off()
}; rm(i, pd)

### >>> 4. Contribution by repeats or non-repeats
# - length of overlapped regions between enhancers and repeats >= 200
# - PE non-repeat 8731 repeat 14809
# - RE93 non-repeat 2897 repeat 10982
# - RE273 non-repeat 563 repeat 1444
pdf(file.path(outdir, "PieDonut_plot_to_show_ratio_of_enhancers_derived_by_repeats_or_not.pdf"), height = 4.5, width = 5)
data.frame(Type = rep(c("Repeats", "Non.repeats"), 3),
           Group = c(rep("PEs", 2), rep("REs.K93", 2), rep("REs.K273", 2)),
           Num = c(14809, 8731, 10982, 2897, 1444, 563),
           Ratio = c(14809/(14809+8731), 8731/(14809+8731), 10982/(10982+2897), 2897/(10982+2897), 1444/(1444+563), 563/(1444+563))) %>%
  mutate(Label = as.character(paste(round(Ratio, 2)*100, "%", sep = ""))) %>% 
  ggplot(aes(fill = Type, y = Ratio, x = Group)) + 
  geom_bar(aes(color = Type), position = "stack", stat = "identity", width = 0.75) +
  geom_text(aes(label = Label), position = position_stack(vjust= 0.5), colour = "#000000", size = 4) +
  coord_polar(theta = "y") +
  theme_bw() +
  labs(x = "Groups", y = "Ratio") +
  theme(axis.title.x = element_text(face="plain", colour = "#1D1D1E", size = 15, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#1D1D1E", size = 15, angle = 90), 
        axis.text.x  = element_text(face="plain", colour = "#1D1D1E", size = 13, angle = 0),
        axis.text.y  = element_text(face="plain", colour = "#1D1D1E", size = 13, angle = 0),
        panel.grid.minor = element_blank())
dev.off()

# testing
pd %>% arrange(desc(Covered_Num)) %>% 
  ggplot(aes(x = fct_reorder(Subfamily, Covered_Num, .desc = T), y = Covered_Num)) + 
  geom_bar(stat = "identity")



### ====================================
### 4th part: Plot the results of TOBIAS
### ====================================

### >>> 1. Set output directory
outdir <- file.path(getwd(), "results/R/Graphs/revised_Fig2")
if (! dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}

### >>> 2. Prepare data
# - gene expression
pd.expr <- gse101571.ge.tpm
pd.expr %>% rownames_to_column(var = "id") %>% 
  mutate(name = DbiIDtrans(genelist = id, intype = "ENSEMBL", outtype = "SYMBOL", species = "human"),
         GV.Oocyte = (GV.oocyte.1+GV.oocyte.2)/2, 
         MII.Oocyte = (MII.oocyte.1+MII.oocyte.2)/2,
         `2Cell` = (`2Cell.1`+`2Cell.1`)/2, 
         `4Cell` = (`4Cell.1`+`4Cell.1`)/2,
         `8Cell` = (`8Cell.1`+`8Cell.1`)/2, 
         ICM = (ICM.1+ICM.2)/2) %>% na.omit() -> pd.expr
pd.expr <- pd.expr[,15:21]
# - co-factors of K9-factors: SETDB1, SETDB2, SUV39H1, SUV39H2, UHRF1, CBX1, CBX5, TRIM28
hs.k9.binds <- read.table("/home/cmq/bioinfo/project-cmq/embryo_93/metadata/K9_factor_list_binding_proteins.txt", sep = "\t")
colnames(hs.k9.binds) <- c("Interactor.A", "Interactor.B", "Experimental.System", "Experimental.Type", "Publication")
hs.k9.eraser.binds <- read.table("/home/cmq/bioinfo/project-cmq/embryo_93/metadata/K9_erasers_list_binding_proteins.txt", sep = "\t")
colnames(hs.k9.eraser.binds) <- c("Interactor.A", "Interactor.B", "Experimental.System", "Experimental.Type", "Publication")
# - List of human genes encoding for KRAB-containing proteins
hs.krab.gene <- read.table("/home/cmq/bioinfo/project-cmq/embryo_93/metadata/pone.0056721.s007.txt", sep = "\t", header = T)


### >>> 2. Reprogrammed enhancers: H3K9me3
# all enhancers
res.k9 <- read.table("/home/cmq/bioinfo/project-cmq/embryo_93/results/tobias/8cell_reprogrammed_enhancer_H3K9me3/BIND_output/bindetect_results.txt",
                     sep = "\t", header = T)
res.k9 %>% 
  mutate(group = case_when((X4cell_8cell_change <= -0.2 & -log10(X4cell_8cell_pvalue) >= 50) | 
                             (X4cell_8cell_change < 0 & -log10(X4cell_8cell_pvalue) >= 100) ~ "8C most enriched",
                           (X4cell_8cell_change <= 0 ) & (-log10(X4cell_8cell_pvalue) >= 50 & -log10(X4cell_8cell_pvalue) < 100 ) ~ "8C less enriched",
                           (X4cell_8cell_change >= 0.2 & -log10(X4cell_8cell_pvalue) >= 50) | 
                             (X4cell_8cell_change > 0 & -log10(X4cell_8cell_pvalue) >= 100) ~ "4C most enriched",
                           (X4cell_8cell_change > 0 ) & (-log10(X4cell_8cell_pvalue) >= 50 & -log10(X4cell_8cell_pvalue) < 100 ) ~ "4C less enriched",
                           -log10(X4cell_8cell_pvalue) < 50 ~ "Low enriched")) -> pd
pd <- merge(pd, pd.expr, by = "name") %>% filter(`8Cell` >= 5 | `4Cell` >= 5)
Reduce(intersect, list(a = pd$name,
                       b = unique(c(hs.k9.binds$Interactor.A, hs.k9.binds$Interactor.B)),
                       c = hs.krab.gene$Gene.symbol))
Reduce(intersect, list(a = pd$name,
                       b = unique(c(hs.k9.eraser.binds$Interactor.A, hs.k9.eraser.binds$Interactor.B))))
pdf(file.path(outdir, "Volcano_plot_to_show_differential_binding_score_between_8cell_4cell_in_REs_H3K9me3.pdf"), height = 4.5, width = 6)
pd %>% 
  ggplot(aes(x = -X4cell_8cell_change, y = -log10(X4cell_8cell_pvalue))) +
  geom_point(aes(fill = group, size = group), shape = 21, color = "#000000") +
  scale_size_manual(values = c(`8C most enriched` = 2.5, `8C less enriched` = 2, 
                               `4C most enriched` = 2.5, `4C less enriched` = 2, 
                               `Low enriched` = 1.5)) +
  scale_fill_manual(values = c(`8C most enriched` = "#F53C3C", `8C less enriched` = "#F58282", 
                               `4C most enriched` = "#0ACEDF", `4C less enriched` = "#78D4DC", 
                               `Low enriched` = "#a5a5a5")) +
  labs(x = "Differential binding score", y = "-log10(pvalue)", fill = "Group", size = "Group") +
  geom_label_repel(aes(label = name),
                   data = subset(pd, name %in% c("Dux", "DUXA", "NFYA", "OTX2", "YY1", "ZSCAN4", "ZKSCAN5", "ZNF263", "ZNF682", "ZNF684")),
                   box.padding = 0.2, point.padding = 0.2, segment.color = "#130f40",
                   label.size = 0.2, fill = "White", size = 3) +
  theme_bw() + 
  theme(axis.title.x = element_text(face="plain", colour = "#1D1D1E", size = 15, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#1D1D1E", size = 15, angle = 90), 
        axis.text.x  = element_text(face="plain", colour = "#1D1D1E", size = 13, angle = 0),
        axis.text.y  = element_text(face="plain", colour = "#1D1D1E", size = 13, angle = 0),
        panel.grid.minor = element_blank())
dev.off()
# SVA enhancers
res.k9 <- read.table("/home/cmq/bioinfo/project-cmq/embryo_93/results/tobias/8cell_reprogrammed_enhancer_H3K9me3_derived_by_SVA/BIND_output/bindetect_results.txt",
                     sep = "\t", header = T)
res.k9 %>% 
  mutate(group = case_when((X4cell_8cell_change <= -0.2 & -log10(X4cell_8cell_pvalue) >= 50) | 
                             (X4cell_8cell_change < 0 & -log10(X4cell_8cell_pvalue) >= 100) ~ "8C most enriched",
                           (X4cell_8cell_change <= 0 ) & (-log10(X4cell_8cell_pvalue) >= 50 & -log10(X4cell_8cell_pvalue) < 100 ) ~ "8C less enriched",
                           (X4cell_8cell_change >= 0.2 & -log10(X4cell_8cell_pvalue) >= 50) | 
                             (X4cell_8cell_change > 0 & -log10(X4cell_8cell_pvalue) >= 100) ~ "4C most enriched",
                           (X4cell_8cell_change > 0 ) & (-log10(X4cell_8cell_pvalue) >= 50 & -log10(X4cell_8cell_pvalue) < 100 ) ~ "4C less enriched",
                           -log10(X4cell_8cell_pvalue) < 50 ~ "Low enriched")) -> pd
pd <- merge(pd, pd.expr, by = "name") %>% filter(`8Cell` >= 5 | `4Cell` >= 5) %>% filter(abs(X4cell_8cell_change) <= 2)
Reduce(intersect, list(a = pd$name,
                       b = unique(c(hs.k9.binds$Interactor.A, hs.k9.binds$Interactor.B)),
                       c = hs.krab.gene$Gene.symbol))
pdf(file.path(outdir, "Volcano_plot_to_show_differential_binding_score_between_8cell_4cell_in_SVA_REs_H3K9me3.pdf"), height = 4.5, width = 6)
pd %>% 
  ggplot(aes(x = -X4cell_8cell_change, y = -log10(X4cell_8cell_pvalue))) +
  geom_point(aes(fill = group, size = group), shape = 21, color = "#000000") +
  scale_size_manual(values = c(`8C most enriched` = 2.5, `8C less enriched` = 2, 
                               `4C most enriched` = 2.5, `4C less enriched` = 2, 
                               `Low enriched` = 1.5)) +
  scale_fill_manual(values = c(`8C most enriched` = "#F53C3C", `8C less enriched` = "#F58282", 
                               `4C most enriched` = "#0ACEDF", `4C less enriched` = "#78D4DC", 
                               `Low enriched` = "#a5a5a5")) +
  labs(x = "Differential binding score", y = "-log10(pvalue)", fill = "Group", size = "Group") +
  geom_label_repel(aes(label = name),
                   data = subset(pd, name %in% c("Dux", "DUXA", "NFYA", "OTX2", "YY1", "ZSCAN4", "ZKSCAN5", "ZNF263", "ZNF682", "ZNF684")),
                   box.padding = 0.2, point.padding = 0.2, segment.color = "#130f40",
                   label.size = 0.2, fill = "White", size = 3) +
  theme_bw() + 
  theme(axis.title.x = element_text(face="plain", colour = "#1D1D1E", size = 15, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#1D1D1E", size = 15, angle = 90), 
        axis.text.x  = element_text(face="plain", colour = "#1D1D1E", size = 13, angle = 0),
        axis.text.y  = element_text(face="plain", colour = "#1D1D1E", size = 13, angle = 0),
        panel.grid.minor = element_blank())
dev.off()


### >>> 3. Reprogrammed enhancers: H3K27me3
# - all enhancers
res.k27 <- read.table("/home/cmq/bioinfo/project-cmq/embryo_93/results/tobias/8cell_reprogrammed_enhancer_H3K27me3/BIND_output/bindetect_results.txt",
                      sep = "\t", header = T)
res.k27 %>% 
  mutate(group = case_when((X4cell_8cell_change <= -0.2 & -log10(X4cell_8cell_pvalue) >= 50) | 
                             (X4cell_8cell_change < 0 & -log10(X4cell_8cell_pvalue) >= 100) ~ "8C most enriched",
                           (X4cell_8cell_change <= 0 ) & (-log10(X4cell_8cell_pvalue) >= 50 & -log10(X4cell_8cell_pvalue) < 100 ) ~ "8C less enriched",
                           (X4cell_8cell_change >= 0.2 & -log10(X4cell_8cell_pvalue) >= 50) | 
                             (X4cell_8cell_change > 0 & -log10(X4cell_8cell_pvalue) >= 100) ~ "4C most enriched",
                           (X4cell_8cell_change > 0 ) & (-log10(X4cell_8cell_pvalue) >= 50 & -log10(X4cell_8cell_pvalue) < 100 ) ~ "4C less enriched",
                           -log10(X4cell_8cell_pvalue) < 50 ~ "Low enriched")) -> pd
pd <- merge(pd, pd.expr, by = "name") %>% filter(`8Cell` >= 5 | `4Cell` >= 5)
Reduce(intersect, list(a = pd$name,
                       b = unique(c(hs.k9.binds$Interactor.A, hs.k9.binds$Interactor.B)),
                       c = hs.krab.gene$Gene.symbol))
pdf(file.path(outdir, "Volcano_plot_to_show_differential_binding_score_between_8cell_4cell_in_REs_H3K27me3.pdf"), height = 4.5, width = 6)
pd %>% 
  ggplot(aes(x = -X4cell_8cell_change, y = -log10(X4cell_8cell_pvalue))) +
  geom_point(aes(fill = group, size = group), shape = 21, color = "#000000") +
  scale_size_manual(values = c(`8C most enriched` = 2.5, `8C less enriched` = 2, 
                               `4C most enriched` = 2.5, `4C less enriched` = 2, 
                               `Low enriched` = 1.5)) +
  scale_fill_manual(values = c(`8C most enriched` = "#F53C3C", `8C less enriched` = "#F58282", 
                               `4C most enriched` = "#0ACEDF", `4C less enriched` = "#78D4DC", 
                               `Low enriched` = "#a5a5a5")) +
  labs(x = "Differential binding score", y = "-log10(pvalue)", fill = "Group", size = "Group") +
  geom_label_repel(aes(label = name),
                   data = subset(pd, name %in% c("Dux", "DUXA", "NFYA", "OTX2", "YY1", "ZSCAN4", "ZKSCAN5", "ZNF263", "ZNF682", "ZNF684")),
                   box.padding = 0.2, point.padding = 0.2, segment.color = "#130f40",
                   label.size = 0.2, fill = "White", size = 3) +
  theme_bw() + 
  theme(axis.title.x = element_text(face="plain", colour = "#1D1D1E", size = 15, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#1D1D1E", size = 15, angle = 90), 
        axis.text.x  = element_text(face="plain", colour = "#1D1D1E", size = 13, angle = 0),
        axis.text.y  = element_text(face="plain", colour = "#1D1D1E", size = 13, angle = 0),
        panel.grid.minor = element_blank())
dev.off()
# - SVA enhancers
res.k27 <- read.table("/home/cmq/bioinfo/project-cmq/embryo_93/results/tobias/8cell_reprogrammed_enhancer_H3K27me3_derived_by_SVA/BIND_output/bindetect_results.txt",
                      sep = "\t", header = T)
res.k27 %>% 
  mutate(group = case_when((X4cell_8cell_change <= -0.2 & -log10(X4cell_8cell_pvalue) >= 50) | 
                             (X4cell_8cell_change < 0 & -log10(X4cell_8cell_pvalue) >= 100) ~ "8C most enriched",
                           (X4cell_8cell_change <= 0 ) & (-log10(X4cell_8cell_pvalue) >= 50 & -log10(X4cell_8cell_pvalue) < 100 ) ~ "8C less enriched",
                           (X4cell_8cell_change >= 0.2 & -log10(X4cell_8cell_pvalue) >= 50) | 
                             (X4cell_8cell_change > 0 & -log10(X4cell_8cell_pvalue) >= 100) ~ "4C most enriched",
                           (X4cell_8cell_change > 0 ) & (-log10(X4cell_8cell_pvalue) >= 50 & -log10(X4cell_8cell_pvalue) < 100 ) ~ "4C less enriched",
                           -log10(X4cell_8cell_pvalue) < 50 ~ "Low enriched")) -> pd
pd <- merge(pd, pd.expr, by = "name") %>% filter(`8Cell` >= 5 | `4Cell` >= 5) %>% filter(abs(X4cell_8cell_change) <= 2)
Reduce(intersect, list(a = pd$name,
                       b = unique(c(hs.k9.binds$Interactor.A, hs.k9.binds$Interactor.B)),
                       c = hs.krab.gene$Gene.symbol))
pdf(file.path(outdir, "Volcano_plot_to_show_differential_binding_score_between_8cell_4cell_in_SVA_REs_H3K27me3.pdf"), height = 4.5, width = 6)
pd %>% 
  ggplot(aes(x = -X4cell_8cell_change, y = -log10(X4cell_8cell_pvalue))) +
  geom_point(aes(fill = group, size = group), shape = 21, color = "#000000") +
  scale_size_manual(values = c(`8C most enriched` = 2.5, `8C less enriched` = 2, 
                               `4C most enriched` = 2.5, `4C less enriched` = 2, 
                               `Low enriched` = 1.5)) +
  scale_fill_manual(values = c(`8C most enriched` = "#F53C3C", `8C less enriched` = "#F58282", 
                               `4C most enriched` = "#0ACEDF", `4C less enriched` = "#78D4DC", 
                               `Low enriched` = "#a5a5a5")) +
  labs(x = "Differential binding score", y = "-log10(pvalue)", fill = "Group", size = "Group") +
  geom_label_repel(aes(label = name),
                   data = subset(pd, name %in% c("Dux", "DUXA", "NFYA", "OTX2", "YY1", "ZSCAN4", "ZKSCAN5", "ZNF263", "ZNF682", "ZNF684")),
                   box.padding = 0.2, point.padding = 0.2, segment.color = "#130f40",
                   label.size = 0.2, fill = "White", size = 3) +
  theme_bw() + 
  theme(axis.title.x = element_text(face="plain", colour = "#1D1D1E", size = 15, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#1D1D1E", size = 15, angle = 90), 
        axis.text.x  = element_text(face="plain", colour = "#1D1D1E", size = 13, angle = 0),
        axis.text.y  = element_text(face="plain", colour = "#1D1D1E", size = 13, angle = 0),
        panel.grid.minor = element_blank())
dev.off()

### >>> 4. Primed enhancers
# - all enhancers
pes <- read.table("/home/cmq/bioinfo/project-cmq/embryo_93/results/tobias/8cell_primed_enhancer/BIND_output/bindetect_results.txt",
                  sep = "\t", header = T)
pes %>% 
  mutate(group = case_when((X4cell_8cell_change <= -0.2 & -log10(X4cell_8cell_pvalue) >= 50) | 
                             (X4cell_8cell_change < 0 & -log10(X4cell_8cell_pvalue) >= 100) ~ "8C most enriched",
                           (X4cell_8cell_change <= 0 ) & (-log10(X4cell_8cell_pvalue) >= 50 & -log10(X4cell_8cell_pvalue) < 100 ) ~ "8C less enriched",
                           (X4cell_8cell_change >= 0.2 & -log10(X4cell_8cell_pvalue) >= 50) | 
                             (X4cell_8cell_change > 0 & -log10(X4cell_8cell_pvalue) >= 100) ~ "4C most enriched",
                           (X4cell_8cell_change > 0 ) & (-log10(X4cell_8cell_pvalue) >= 50 & -log10(X4cell_8cell_pvalue) < 100 ) ~ "4C less enriched",
                           -log10(X4cell_8cell_pvalue) < 50 ~ "Low enriched")) -> pd
pd <- merge(pd, pd.expr, by = "name") %>% filter(`8Cell` >= 5 | `4Cell` >= 5)
Reduce(intersect, list(a = pd$name,
                       b = unique(c(hs.k9.binds$Interactor.A, hs.k9.binds$Interactor.B)),
                       c = hs.krab.gene$Gene.symbol))
pdf(file.path(outdir, "Volcano_plot_to_show_differential_binding_score_between_8cell_4cell_in_PEs.pdf"), height = 4.5, width = 6)
pd %>% 
  ggplot(aes(x = -X4cell_8cell_change, y = -log10(X4cell_8cell_pvalue))) +
  geom_point(aes(fill = group, size = group), shape = 21, color = "#000000") +
  scale_size_manual(values = c(`8C most enriched` = 2.5, `8C less enriched` = 2, 
                               `4C most enriched` = 2.5, `4C less enriched` = 2, 
                               `Low enriched` = 1.5)) +
  scale_fill_manual(values = c(`8C most enriched` = "#F53C3C", `8C less enriched` = "#F58282", 
                               `4C most enriched` = "#0ACEDF", `4C less enriched` = "#78D4DC", 
                               `Low enriched` = "#a5a5a5")) +
  labs(x = "Differential binding score", y = "-log10(pvalue)", fill = "Group", size = "Group") +
  geom_label_repel(aes(label = name),
                   data = subset(pd, name %in% c("Dux", "DUXA", "NFYA", "OTX2", "YY1", "ZSCAN4", "ZKSCAN5", "ZNF263", "ZNF682", "ZNF684")),
                   box.padding = 0.2, point.padding = 0.2, segment.color = "#130f40",
                   label.size = 0.2, fill = "White", size = 3) +
  theme_bw() + 
  theme(axis.title.x = element_text(face="plain", colour = "#1D1D1E", size = 15, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#1D1D1E", size = 15, angle = 90), 
        axis.text.x  = element_text(face="plain", colour = "#1D1D1E", size = 13, angle = 0),
        axis.text.y  = element_text(face="plain", colour = "#1D1D1E", size = 13, angle = 0),
        panel.grid.minor = element_blank())
dev.off()
# - SVA enhancers
pes <- read.table("/home/cmq/bioinfo/project-cmq/embryo_93/results/tobias/8cell_primed_enhancer_derived_by_SVA/BIND_output/bindetect_results.txt",
                  sep = "\t", header = T)
pes %>% 
  mutate(group = case_when((X4cell_8cell_change <= -0.2 & -log10(X4cell_8cell_pvalue) >= 50) | 
                             (X4cell_8cell_change < 0 & -log10(X4cell_8cell_pvalue) >= 100) ~ "8C most enriched",
                           (X4cell_8cell_change <= 0 ) & (-log10(X4cell_8cell_pvalue) >= 50 & -log10(X4cell_8cell_pvalue) < 100 ) ~ "8C less enriched",
                           (X4cell_8cell_change >= 0.2 & -log10(X4cell_8cell_pvalue) >= 50) | 
                             (X4cell_8cell_change > 0 & -log10(X4cell_8cell_pvalue) >= 100) ~ "4C most enriched",
                           (X4cell_8cell_change > 0 ) & (-log10(X4cell_8cell_pvalue) >= 50 & -log10(X4cell_8cell_pvalue) < 100 ) ~ "4C less enriched",
                           -log10(X4cell_8cell_pvalue) < 50 ~ "Low enriched")) -> pd
pd <- merge(pd, pd.expr, by = "name") %>% filter(`8Cell` >= 5 | `4Cell` >= 5) %>% filter(abs(X4cell_8cell_change) <= 2)
Reduce(intersect, list(a = pd$name,
                       b = unique(c(hs.k9.binds$Interactor.A, hs.k9.binds$Interactor.B)),
                       c = hs.krab.gene$Gene.symbol))
pdf(file.path(outdir, "Volcano_plot_to_show_differential_binding_score_between_8cell_4cell_in_SVA_PEs.pdf"), height = 4.5, width = 6)
pd %>% 
  ggplot(aes(x = -X4cell_8cell_change, y = -log10(X4cell_8cell_pvalue))) +
  geom_point(aes(fill = group, size = group), shape = 21, color = "#000000") +
  scale_size_manual(values = c(`8C most enriched` = 2.5, `8C less enriched` = 2, 
                               `4C most enriched` = 2.5, `4C less enriched` = 2, 
                               `Low enriched` = 1.5)) +
  scale_fill_manual(values = c(`8C most enriched` = "#F53C3C", `8C less enriched` = "#F58282", 
                               `4C most enriched` = "#0ACEDF", `4C less enriched` = "#78D4DC", 
                               `Low enriched` = "#a5a5a5")) +
  labs(x = "Differential binding score", y = "-log10(pvalue)", fill = "Group", size = "Group") +
  geom_label_repel(aes(label = name),
                   data = subset(pd, name %in% c("Dux", "DUXA", "NFYA", "OTX2", "YY1", "ZSCAN4", "ZKSCAN5", "ZNF263", "ZNF682", "ZNF684")),
                   box.padding = 0.2, point.padding = 0.2, segment.color = "#130f40",
                   label.size = 0.2, fill = "White", size = 3) +
  theme_bw() + 
  theme(axis.title.x = element_text(face="plain", colour = "#1D1D1E", size = 15, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#1D1D1E", size = 15, angle = 90), 
        axis.text.x  = element_text(face="plain", colour = "#1D1D1E", size = 13, angle = 0),
        axis.text.y  = element_text(face="plain", colour = "#1D1D1E", size = 13, angle = 0),
        panel.grid.minor = element_blank())
dev.off()



### =====================================================
### 5th part: 8Cell Enhancer Motif Analysis Visualization
### =====================================================

### >>> 1. Set output directory
outdir <- file.path(getwd(), "results/R/Graphs/revised_Fig2")
if (! dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}

### >>> 2. load table for plot
# motif pvalue
ehc.motif <- read.table("results/homer/findMotifs/revised_8c_enhancer/motif_log10pvalue.txt", header = T)
rownames(ehc.motif) <-ehc.motif$SYMBOL

# gene expression table
# single cell GSE36552
ge.tpm.gse36552 <- read.table("results/R/Tables/reprogrammed/enrich/hs_embryo_sc_rna_seq_gene_tpm.txt", sep = "\t", header = T)
colnames(ge.tpm.gse36552)
# single cell GSE44183
ge.co.gse44183 <- read.table("/home/data/publicdata/GSE44183/analysis/human/results/featurecounts/top10_gene/all_samples_gene_count_name_matrix.txt", 
                             header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
colnames(ge.co.gse44183)[-1:-5] <- c("blood", "oocyte_1", "oocyte_2", "oocyte_3", "pronucleus_1", "pronucleus_2", "pronucleus_3", 
                                     "zygote_1", "zygote_2", "2cell_1", "2cell_2", "2cell_3", "4cell_1", "4cell_2", "4cell_3", "4cell_4", 
                                     "8cell_1", "8cell_2", "8cell_3", "8cell_4", "8cell_5", "8cell_6", "8cell_7", "8cell_8", "8cell_9", 
                                     "8cell_10", "8cell_11", "morula_1", "morula_2", "morula_3") 
ge.tpm.gse44183 <- CountToTpm(ge.co.gse44183[,-1:-5], ge.co.gse44183$Length)
ge.tpm.gse44183 %>% rownames_to_column("SYMBOL") %>% 
  mutate(Oocyte = rowMeans(ge.tpm.gse44183[,grep("oocyte_", colnames(ge.tpm.gse44183))]),
         Zygote = rowMeans(ge.tpm.gse44183[,grep("zygote_", colnames(ge.tpm.gse44183))]),
         `2Cell` = rowMeans(ge.tpm.gse44183[,grep("2cell_", colnames(ge.tpm.gse44183))]), 
         `4Cell` = rowMeans(ge.tpm.gse44183[,grep("4cell_", colnames(ge.tpm.gse44183))]),
         `8Cell` = rowMeans(ge.tpm.gse44183[,grep("8cell_", colnames(ge.tpm.gse44183))]),
         `Morula` = rowMeans(ge.tpm.gse44183[,grep("morula_", colnames(ge.tpm.gse44183))])) -> ge.tpm.gse44183
ge.tpm.gse44183 <- ge.tpm.gse44183[,c(1,32:37)]
colnames(ge.tpm.gse44183)

# bulk GSE101571
ge.co.gse101571 <- read.table("/home/data/publicdata/GSE101571/rnaseq/analysis/results/featurecounts/top10_gene_mul/all_samples_gene_count_name_matrix.txt", 
                              row.names = 1, header = T, stringsAsFactors = F, sep = "\t")
ge.co.gse101571 <- ge.co.gse101571[,c(1:5, 10:11, 13:18)]
colnames(ge.co.gse101571)[-1:-5] <- c("2Cell_1", "2Cell_2", "4Cell_1", "4Cell_2", "8Cell_1", "8Cell_2", "ICM_1", "ICM_2")
ge.tpm.gse101571 <- CountToTpm(ge.co.gse101571[,-1:-5], ge.co.gse101571$Length)
ge.tpm.gse101571 %>% mutate(`2Cell_ave`=rowMeans(ge.tpm.gse101571[,1:2]), `4Cell_ave`=rowMeans(ge.tpm.gse101571[,3:4]), 
                            `8Cell_ave`=rowMeans(ge.tpm.gse101571[,5:6]), ICM_ave=rowMeans(ge.tpm.gse101571[,7:8])) -> ge.tpm.gse101571
ge.tpm.gse101571$SYMBOL <- rownames(ge.tpm.gse101571)
colnames(ge.tpm.gse101571)

ge.tpm.gse44183[ge.tpm.gse44183$SYMBOL=="DUX4",]
ge.tpm.gse36552[ge.tpm.gse36552$SYMBOL%in%"DUX4",]
ge.tpm.gse101571["DUX4",]

### >>> 3. merge gene expression level and motif pvalue table
ehc.motif <- merge(ehc.motif, ge.tpm.gse36552[,c(3:6, 8)]) %>% merge(ge.tpm.gse44183[,c(1,4:7)]) %>% merge(ge.tpm.gse101571[,9:13])
colnames(ehc.motif)
write.table(ehc.motif, file.path(outdir, "8cell_enhancer_enriched_motif.txt"), sep = "\t", quote = F, row.names = F, col.names = T)

### >>> 4. visualization via bubble plot
pd <- ehc.motif[,c(1:4, 15)]
pdf(file.path(outdir, "Bubble_plot_to_show_8cell_enhancer_enriched_motif.pdf"), height = 5.25, width = 4)
pd %>% gather(key=group, value = log10P, -SYMBOL, -`8Cell_ave`) %>% 
  mutate(group=str_split_fixed(group, "_", 2)[,1]) %>% 
  mutate(group=case_when(group=="pe" ~ "PEs", 
                         group=="re93" ~ "H3K9me3 REs",
                         group=="re273" ~ "H3K27me3 REs")) %>% 
  mutate(group=factor(group, levels = c("PEs", "H3K9me3 REs", "H3K27me3 REs"))) %>% 
  mutate(size=case_when(log10P<=2 ~ "a", (log10P>2&log10P<=50) ~ "b", (log10P>50&log10P<=200) ~ "c", (log10P>200&log10P<=400) ~ "d", log10P>400 ~ "e")) %>% 
  mutate(color=case_when(`8Cell_ave`<1 ~ "a", (`8Cell_ave`>=1&`8Cell_ave`<10) ~ "b", (`8Cell_ave`>=10&`8Cell_ave`<20) ~ "c", 
                         (`8Cell_ave`>=20&`8Cell_ave`<50) ~ "d", (`8Cell_ave`>=50) ~ "e")) %>% 
  ggplot(aes(x=group, y=SYMBOL)) + 
  geom_point(aes(size = size, color = color), alpha = 1, shape = 16) + 
  scale_size_manual(values = c(a = 1, b = 4, c = 5.5, d = 7, e = 8.5), labels = c("<=2", "2~50", "50~200", "200~400", ">400"), name = "-log10(Pvalue)") + 
  scale_color_manual(values = c(a = "#99a3a4", b = "#a9cce3", c = "#5499c7", d = "#2471a3", e = "#21618c"), 
                     labels = c("<1", "1~10", "10~20", "20~50", ">=50"), name = "8Cell GSE101571\nTPM") + 
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x  = element_text(face="plain", colour = "#000000", size = 12, angle = 45, hjust = 1),
        axis.text.y  = element_text(face="plain", colour = "#000000", size = 12, angle = 0, hjust = 1),
        panel.grid = element_blank())
dev.off()
library(ComplexHeatmap)
### >>> 5. visualization via heatmap
pd <- ehc.motif
pd[pd$pe < -log10(0.05),]$pe <- NA
pd[pd$re273 < -log10(0.05),]$re273 <- NA
rownames(pd) <- pd$SYMBOL
pdf(file.path(outdir, "Heatmap_to_show_8cell_enhancer_enriched_motif.pdf"), height = 6, width = 8)
h1 <- Heatmap(as.matrix(pd[,6:8]), 
              #col = colorRamp2(c(0, 10, 20, 30, 40), c("#481b6d", "#1e9b8a", "#63A97D", "#72cf56", "#f7e620")), 
              #col = colorRampPalette(c("#481b6d", "#1e9b8a", "#63A97D", "#72cf56", "#f7e620"))(20),
              col = colorRamp2(c(0, 20, 40), c("#1078d4", "#ffffff", "#f03039")), 
              name = "GSE36552 TPM",
              # clustering setting
              cluster_rows = T, show_row_dend = F, row_dend_side = "left", row_dend_width = unit(1, "cm"), 
              clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
              cluster_columns = F, show_column_dend = T, column_dend_height = unit(1, "cm"), column_dend_side = "top",
              clustering_distance_columns = "euclidean", clustering_method_columns = "complete",
              # names, rotation and position of columns and row names
              show_column_names = T, column_names_rot = 45, column_names_side = "bottom",
              show_row_names = T, row_names_rot = 0, row_names_side = "left",
              # size of heatmap, also including heatmap_height and heatmap_width
              width = unit(2, "cm"), height = unit(8, "cm"),
              # control row names size
              row_names_max_width = max_text_width(rownames(pd), gp = gpar(fontsize = 10)), 
              top_annotation = HeatmapAnnotation(Sample = c("4Cell", "8Cell", "Morula"),  
                                                 simple_anno_size = unit(0.2, "cm"), border = T, show_annotation_name = F, 
                                                 col = list(Sample = c("4Cell" = "#1B9E77", 
                                                                       "8Cell" = "#D95F02", "Morula" = "#7570B3"))))
h2 <- Heatmap(as.matrix(pd[,10:12]), 
              #col = colorRamp2(c(0, 10, 20, 30, 40), c("#481b6d", "#1e9b8a", "#63A97D", "#72cf56", "#f7e620")),
              #col = colorRampPalette(c("#481b6d", "#1e9b8a", "#63A97D", "#72cf56", "#f7e620"))(20),
              col = colorRamp2(c(0, 20, 40), c("#1078d4", "#ffffff", "#f03039")), 
              name = "GSE44183 TPM",
              # clustering setting
              cluster_rows = T, show_row_dend = F, row_dend_side = "left", row_dend_width = unit(1, "cm"), 
              clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
              cluster_columns = F, show_column_dend = T, column_dend_height = unit(1, "cm"), column_dend_side = "top",
              clustering_distance_columns = "euclidean", clustering_method_columns = "complete",
              # names, rotation and position of columns and row names
              show_column_names = T, column_names_rot = 45, column_names_side = "bottom",
              show_row_names = F, row_names_rot = 0, row_names_side = "left",
              # size of heatmap, also including heatmap_height and heatmap_width
              width = unit(2, "cm"), height = unit(8, "cm"),
              # control row names size
              top_annotation = HeatmapAnnotation(Sample = c("4Cell", "8Cell", "Morula"),  
                                                 simple_anno_size = unit(0.2, "cm"), border = T, show_annotation_name = F, 
                                                 col = list(Sample = c("4Cell" = "#1B9E77", 
                                                                       "8Cell" = "#D95F02", "Morula" = "#7570B3"))))

h3 <- Heatmap(as.matrix(pd[,14:16]), 
              #col = colorRamp2(c(0, 10, 20, 30, 40), c("#481b6d", "#1e9b8a", "#63A97D", "#72cf56", "#f7e620")),
              #col = colorRampPalette(c("#481b6d", "#1e9b8a", "#63A97D", "#72cf56", "#f7e620"))(20),
              col = colorRamp2(c(0, 20, 40), c("#1078d4", "#ffffff", "#f03039")), 
              name = "GSE101571 TPM",
              # clustering setting
              cluster_rows = T, show_row_dend = F, row_dend_side = "left", row_dend_width = unit(1, "cm"), 
              clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
              cluster_columns = F, show_column_dend = T, column_dend_height = unit(1, "cm"), column_dend_side = "top",
              clustering_distance_columns = "euclidean", clustering_method_columns = "complete",
              # names, rotation and position of columns and row names
              show_column_names = T, column_names_rot = 45, column_names_side = "bottom",
              show_row_names = F, row_names_rot = 0, row_names_side = "left",
              # size of heatmap, also including heatmap_height and heatmap_width
              width = unit(2, "cm"), height = unit(8, "cm"),
              # control row names size
              top_annotation = HeatmapAnnotation(Sample = c("4Cell", "8Cell", "ICM"),  
                                                 simple_anno_size = unit(0.2, "cm"), border = T, show_annotation_name = F, 
                                                 col = list(Sample = c("4Cell" = "#1B9E77", 
                                                                       "8Cell" = "#D95F02", "ICM" = "#e3e306"))))

h4 <- Heatmap(pd[,2], 
              name = "PEs", 
              col = colorRamp2(c(0, 20, 40, 80, 200), c("#99a3a4", "#a9cce3", "#5499c7", "#2471a3", "#21618c")),
              width = unit(0.5, "cm"), column_names_rot = 45, na_col = "#B2B5B5", show_row_names = F, show_column_names = T, 
              cluster_rows = F, cluster_columns = F, border = T)
h5 <- Heatmap(pd[,3], 
              name = "93 REs", 
              col = colorRamp2(c(0, 20, 40, 80, 200), c("#99a3a4", "#a9cce3", "#5499c7", "#2471a3", "#21618c")),
              width = unit(0.5, "cm"), column_names_rot = 45, na_col = "#B2B5B5", show_row_names = F, show_column_names = T, 
              cluster_rows = F, cluster_columns = F, border = T)

h1 + h2 + h3 + h4 + h5 + Heatmap(pd[,4], 
                                 name = "273 REs", 
                                 col = colorRamp2(c(0, 20, 40, 80, 200), c("#99a3a4", "#a9cce3", "#5499c7", "#2471a3", "#21618c")), 
                                 width = unit(0.5, "cm"), column_names_rot = 45, na_col = "#B2B5B5", show_row_names = F, show_column_names = T, 
                                 cluster_rows = F, cluster_columns = F, border = T)
dev.off(); rm(h1, h2, h3, h4, h5, pd)



### ====================================
### 6th part: REs and PEs targeted genes
### ====================================

### >>> 1. Set output directory
outdir <- file.path(getwd(), "results/R/Graphs/revised_Fig2")
if (! dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}

### >>> 2. Load major ZGA genes
hs.major.zga <- read.table("metadata/hs_zga_gene_pc_lncrna_ensembl.bed") %>% separate(V4, c("EnsemblID", "Symbol"), sep = ":")


### >>> 3. Load genes in loops
hs.gene.in.loops <- list()
for (i in seq(1, length(list.files("results/bedtools/revised_8c_enhancer", "*full.txt")))) {
  hs.gene.in.loops[[i]] <- read.table(list.files("results/bedtools/revised_8c_enhancer", "*full.txt", full.names = T)[[i]], sep = "\t") %>% 
    separate(V16, c("EnsemblID", "Symbol"), sep = ":") %>% unique()
  hs.gene.in.loops[[i]] <- subset(hs.gene.in.loops[[i]], EnsemblID %in% hs.major.zga$EnsemblID)
  hs.gene.in.loops[[i]] <- hs.gene.in.loops[[i]][!duplicated(hs.gene.in.loops[[i]]$EnsemblID),]
}; rm(i)
names(hs.gene.in.loops) <- gsub("\\_full.txt", "", gsub("genes_interacted_with_", "", list.files("results/bedtools/revised_8c_enhancer", "*full.txt")))


### >>> 4. Classification
# - 2cell
c2.k93 <- setdiff(hs.gene.in.loops$RE_k93_enhancers_in_2cell$EnsemblID, 
                  c(hs.gene.in.loops$PE_enhancers_in_2cell$EnsemblID, hs.gene.in.loops$RE_k273_enhancers_in_2cell$EnsemblID)) %>% unique()
c2.k273 <- setdiff(hs.gene.in.loops$RE_k273_enhancers_in_2cell$EnsemblID, 
                   c(hs.gene.in.loops$PE_enhancers_in_2cell$EnsemblID, hs.gene.in.loops$RE_k93_enhancers_in_2cell$EnsemblID)) %>% unique()
c2.pe <- setdiff(hs.gene.in.loops$PE_enhancers_in_2cell$EnsemblID, 
                 c(hs.gene.in.loops$RE_k93_enhancers_in_2cell$EnsemblID, hs.gene.in.loops$RE_k273_enhancers_in_2cell$EnsemblID)) %>% unique()
c2.k93.k273.pe <- intersect(hs.gene.in.loops$RE_k93_enhancers_in_2cell$EnsemblID,
                            intersect(hs.gene.in.loops$PE_enhancers_in_2cell$EnsemblID, hs.gene.in.loops$RE_k273_enhancers_in_2cell$EnsemblID)) %>% unique()
c2.k93.k273 <- setdiff(intersect(hs.gene.in.loops$RE_k93_enhancers_in_2cell$EnsemblID, hs.gene.in.loops$RE_k273_enhancers_in_2cell$EnsemblID),
                       c2.k93.k273.pe) %>% unique()
c2.k93.pe <- setdiff(intersect(hs.gene.in.loops$RE_k93_enhancers_in_2cell$EnsemblID, hs.gene.in.loops$PE_enhancers_in_2cell$EnsemblID),
                     c2.k93.k273.pe) %>% unique()
c2.k273.pe <- setdiff(intersect(hs.gene.in.loops$RE_k273_enhancers_in_2cell$EnsemblID, hs.gene.in.loops$PE_enhancers_in_2cell$EnsemblID),
                      c2.k93.k273.pe) %>% unique()
# - 8cell
c8.k93 <- setdiff(hs.gene.in.loops$RE_k93_enhancers_in_8cell$EnsemblID, 
                  c(hs.gene.in.loops$PE_enhancers_in_8cell$EnsemblID, hs.gene.in.loops$RE_k273_enhancers_in_8cell$EnsemblID)) %>% unique()
c8.k273 <- setdiff(hs.gene.in.loops$RE_k273_enhancers_in_8cell$EnsemblID, 
                   c(hs.gene.in.loops$PE_enhancers_in_8cell$EnsemblID, hs.gene.in.loops$RE_k93_enhancers_in_8cell$EnsemblID)) %>% unique()
c8.pe <- setdiff(hs.gene.in.loops$PE_enhancers_in_8cell$EnsemblID, 
                 c(hs.gene.in.loops$RE_k93_enhancers_in_8cell$EnsemblID, hs.gene.in.loops$RE_k273_enhancers_in_8cell$EnsemblID)) %>% unique()
c8.k93.k273.pe <- intersect(hs.gene.in.loops$RE_k93_enhancers_in_8cell$EnsemblID,
                            intersect(hs.gene.in.loops$PE_enhancers_in_8cell$EnsemblID, hs.gene.in.loops$RE_k273_enhancers_in_8cell$EnsemblID)) %>% unique()
c8.k93.k273 <- setdiff(intersect(hs.gene.in.loops$RE_k93_enhancers_in_8cell$EnsemblID, hs.gene.in.loops$RE_k273_enhancers_in_8cell$EnsemblID),
                       c8.k93.k273.pe) %>% unique()
c8.k93.pe <- setdiff(intersect(hs.gene.in.loops$RE_k93_enhancers_in_8cell$EnsemblID, hs.gene.in.loops$PE_enhancers_in_8cell$EnsemblID),
                     c8.k93.k273.pe) %>% unique()
c8.k273.pe <- setdiff(intersect(hs.gene.in.loops$RE_k273_enhancers_in_8cell$EnsemblID, hs.gene.in.loops$PE_enhancers_in_8cell$EnsemblID),
                      c8.k93.k273.pe) %>% unique()
# - classify ZGA genes
hs.major.zga <- hs.major.zga %>%
  mutate(c2 = case_when(EnsemblID %in% c2.pe ~ "PEs",
                        EnsemblID %in% c2.k93 ~ "REs.K9",
                        EnsemblID %in% c2.k273 ~ "REs.K273"),
         c8 = case_when(EnsemblID %in% c8.pe ~ "PEs",
                        EnsemblID %in% c8.k93 ~ "REs.K9",
                        EnsemblID %in% c8.k273 ~ "REs.K273"))
hs.major.zga$c2[is.na(hs.major.zga$c2)] <- "None"
hs.major.zga$c8[is.na(hs.major.zga$c8)] <- "None"
# - Plotting the classification
pdf("Graphs/reprogrammed/hic/2cell_hic_loops_with_major_ZGA_genes_classification_by_enhancer.pdf", height = 5, width = 5)
PieDonut(hs.major.zga,
         r0 = getOption("PieDonut.r0", 0.3), r1 = getOption("PieDonut.r1", 0.7), r2 = getOption("PieDonut.r2", 0.9),
         aes(pies = c2), addDonutLabel = TRUE, showRatioDonut = F,
         showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.01), color = "#ffffff",
         ratioByGroup = F, labelposition = 1, pieLabelSize = 3, donutLabelSize = 3)
dev.off()
pdf("Graphs/reprogrammed/hic/8cell_hic_loops_with_major_ZGA_genes_classification_by_enhancer.pdf", height = 5, width = 5)
PieDonut(hs.major.zga,
         r0 = getOption("PieDonut.r0", 0.3), r1 = getOption("PieDonut.r1", 0.7), r2 = getOption("PieDonut.r2", 0.9),
         aes(pies = c8), addDonutLabel = TRUE, showRatioDonut = F,
         showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.01), color = "#ffffff",
         ratioByGroup = F, labelposition = 1, pieLabelSize = 3, donutLabelSize = 3)
dev.off()

### >>> 5. Make a statistic 
# - all enhancers
enhancer.stat <- data.frame(Stage = c(rep("2C", 3), rep("8C", 3)),
                            Group = rep(c("REs.K93", "REs.K273", "PEs"), 2),
                            Num = c(length(unique(hs.gene.in.loops$c8_RE_k93_enhancers_in_2cell$EnsemblID)),
                                    length(unique(hs.gene.in.loops$c8_RE_k273_enhancers_in_2cell$EnsemblID)),
                                    length(unique(hs.gene.in.loops$c8_PE_enhancers_in_2cell$EnsemblID)),
                                    length(unique(hs.gene.in.loops$c8_RE_k93_enhancers_in_8cell$EnsemblID)),
                                    length(unique(hs.gene.in.loops$c8_RE_k273_enhancers_in_8cell$EnsemblID)),
                                    length(unique(hs.gene.in.loops$c8_PE_enhancers_in_8cell$EnsemblID))))
pdf(file.path(outdir, "Bar_plot_to_show_changes_of_ZGA_genes_intersected_with_enhancers.pdf"), height = 3.5, width = 4)
enhancer.stat %>% 
  ggplot(aes(x = Stage, y = Num/nrow(hs.major.zga))) +
  geom_bar(aes(fill = Stage), width = 0.75, stat = "identity") +
  facet_wrap(. ~ Group) +
  labs(x = "Stage", y = "Ratio in ZGA genes") +
  theme_bw() +
  scale_fill_manual(values = c(`8C` = "#F44C4C", `2C` = "#0B88CF")) +
  theme(axis.title.x = element_text(face="plain", colour = "#1D1D1E", size = 15, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#1D1D1E", size = 15, angle = 90), 
        axis.text.x  = element_text(face="plain", colour = "#1D1D1E", size = 13, angle = 0),
        axis.text.y  = element_text(face="plain", colour = "#1D1D1E", size = 13, angle = 0),
        panel.grid.minor = element_blank(), legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white", size = 0.5, linetype = "solid"),
        strip.text.x = element_text(face = "plain", size = 13, colour = "black", angle = 0),
        strip.text.y = element_text(face = "plain", size = 13, colour = "black", angle = 0))
dev.off()
# - SVA enhancers
enhancer.stat <- data.frame(Stage = c(rep("2C", 3), rep("8C", 3)),
                            Group = rep(c("REs.K93", "REs.K273", "PEs"), 2),
                            Num = c(length(unique(hs.gene.in.loops$sva_c8_RE_k93_enhancers_in_2cell$EnsemblID)),
                                    length(unique(hs.gene.in.loops$sva_c8_RE_k273_enhancers_in_2cell$EnsemblID)),
                                    length(unique(hs.gene.in.loops$sva_c8_PE_enhancers_in_2cell$EnsemblID)),
                                    length(unique(hs.gene.in.loops$sva_c8_RE_k93_enhancers_in_8cell$EnsemblID)),
                                    length(unique(hs.gene.in.loops$sva_c8_RE_k273_enhancers_in_8cell$EnsemblID)),
                                    length(unique(hs.gene.in.loops$sva_c8_PE_enhancers_in_8cell$EnsemblID))))
pdf(file.path(outdir, "Bar_plot_to_show_changes_of_ZGA_genes_intersected_with_SVA_enhancers.pdf"), height = 3.5, width = 4)
enhancer.stat %>% 
  ggplot(aes(x = Stage, y = Num/nrow(hs.major.zga))) +
  geom_bar(aes(fill = Stage), width = 0.75, stat = "identity") +
  facet_wrap(. ~ Group) +
  labs(x = "Stage", y = "Ratio in ZGA genes") +
  theme_bw() +
  scale_fill_manual(values = c(`8C` = "#F44C4C", `2C` = "#0B88CF")) +
  theme(axis.title.x = element_text(face="plain", colour = "#1D1D1E", size = 15, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#1D1D1E", size = 15, angle = 90), 
        axis.text.x  = element_text(face="plain", colour = "#1D1D1E", size = 13, angle = 0),
        axis.text.y  = element_text(face="plain", colour = "#1D1D1E", size = 13, angle = 0),
        panel.grid.minor = element_blank(), legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white", size = 0.5, linetype = "solid"),
        strip.text.x = element_text(face = "plain", size = 13, colour = "black", angle = 0),
        strip.text.y = element_text(face = "plain", size = 13, colour = "black", angle = 0))
dev.off()


### >>> 6. PEs/REs-nearby genes
# - load expr data
gse101571.ge.count <- read.table("results/featurecount/GSE101571_all_samples_gene_count_id_matrix.txt", header = T, row.names = 1)
colnames(gse101571.ge.count)[-1:-5] <- c("GV.oocyte.1", "GV.oocyte.2", "MII.oocyte.1", "MII.oocyte.2", "2Cell.1", "2Cell.2", "2Cell.3",
                                         "4Cell.1", "4Cell.2", "8Cell.1", "8Cell.2", "ICM.1", "ICM.2")
keep <- rowSums(gse101571.ge.count[,-1:-5]) >= 10
gse101571.ge.count <- gse101571.ge.count[keep, ]; rm(keep)
gse101571.ge.tpm <- CountToTpm(gse101571.ge.count[,-1:-5], gse101571.ge.count$Length)
# - make a plot (enhancer-targeted genes)
pdf(file.path(outdir, "Box_plot_to_show_changes_of_genes_expr_intersected_with_enhancers.pdf"), height = 3.5, width = 10)
p1 <- gse101571.ge.tpm[hs.gene.in.loops$c8_RE_k93_enhancers_in_8cell$EnsemblID, ] %>%
  mutate(GV.Oocyte = (GV.oocyte.1+GV.oocyte.2)/2, MII.Oocyte = (MII.oocyte.1+MII.oocyte.2)/2,
         `2Cell` = (`2Cell.1`+`2Cell.1`)/2, `4Cell` = (`4Cell.1`+`4Cell.1`)/2,
         `8Cell` = (`8Cell.1`+`8Cell.1`)/2, ICM = (ICM.1+ICM.2)/2) %>%
  gather(key = "Celltype", value = "TPM") %>%
  filter(Celltype %in% c("MII.Oocyte", "2Cell", "4Cell", "8Cell", "ICM")) %>%
  mutate(Celltype = factor(Celltype, levels = c("MII.Oocyte", "2Cell", "4Cell", "8Cell", "ICM"))) %>%
  ggplot(aes(x = Celltype, y = log2(TPM+0.1))) +
  geom_boxplot(aes(fill = Celltype)) +
  scale_fill_brewer() + theme_bw() + ggtitle("REs.K93 interacted genes") +
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank(),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 11, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 11, angle = 45, hjust = 1),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 11, angle = 0)) +
  stat_compare_means(method = "t.test", na.rm = T, label.y = c(14), comparisons = list(c("4Cell", "8Cell")), size = 4)
p2 <- gse101571.ge.tpm[hs.gene.in.loops$c8_RE_k273_enhancers_in_8cell$EnsemblID, ] %>%
  mutate(GV.Oocyte = (GV.oocyte.1+GV.oocyte.2)/2, MII.Oocyte = (MII.oocyte.1+MII.oocyte.2)/2,
         `2Cell` = (`2Cell.1`+`2Cell.1`)/2, `4Cell` = (`4Cell.1`+`4Cell.1`)/2,
         `8Cell` = (`8Cell.1`+`8Cell.1`)/2, ICM = (ICM.1+ICM.2)/2) %>%
  gather(key = "Celltype", value = "TPM") %>%
  filter(Celltype %in% c("MII.Oocyte", "2Cell", "4Cell", "8Cell", "ICM")) %>%
  mutate(Celltype = factor(Celltype, levels = c("MII.Oocyte", "2Cell", "4Cell", "8Cell", "ICM"))) %>%
  ggplot(aes(x = Celltype, y = log2(TPM+0.1))) +
  geom_boxplot(aes(fill = Celltype)) +
  scale_fill_brewer() + theme_bw() + ggtitle("REs.K273 interacted genes") +
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank(),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 11, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 11, angle = 45, hjust = 1),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 11, angle = 0)) +
  stat_compare_means(method = "t.test", na.rm = T, label.y = c(14), comparisons = list(c("4Cell", "8Cell")), size = 4)
p3 <- gse101571.ge.tpm[hs.gene.in.loops$c8_PE_enhancers_in_8cell$EnsemblID, ] %>%
  mutate(GV.Oocyte = (GV.oocyte.1+GV.oocyte.2)/2, MII.Oocyte = (MII.oocyte.1+MII.oocyte.2)/2,
         `2Cell` = (`2Cell.1`+`2Cell.1`)/2, `4Cell` = (`4Cell.1`+`4Cell.1`)/2,
         `8Cell` = (`8Cell.1`+`8Cell.1`)/2, ICM = (ICM.1+ICM.2)/2) %>%
  gather(key = "Celltype", value = "TPM") %>%
  filter(Celltype %in% c("MII.Oocyte", "2Cell", "4Cell", "8Cell", "ICM")) %>%
  mutate(Celltype = factor(Celltype, levels = c("MII.Oocyte", "2Cell", "4Cell", "8Cell", "ICM"))) %>%
  ggplot(aes(x = Celltype, y = log2(TPM+0.1))) +
  geom_boxplot(aes(fill = Celltype)) +
  scale_fill_brewer() + theme_bw() + ggtitle("PEs interacted genes") +
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank(),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 11, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 11, angle = 45, hjust = 1),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 11, angle = 0)) +
  stat_compare_means(method = "t.test", na.rm = T, label.y = c(14), comparisons = list(c("4Cell", "8Cell")), size = 4)
plot_grid(p1, p2, p3, ncol = 3)
dev.off(); rm(p1, p2, p3)
# - make a plot (SVA enhancer-targeted genes)
pdf(file.path(outdir, "Box_plot_to_show_changes_of_genes_expr_intersected_with_SVA_enhancers.pdf"), height = 3.5, width = 10)
p1 <- gse101571.ge.tpm[hs.gene.in.loops$sva_c8_RE_k93_enhancers_in_8cell$EnsemblID, ] %>%
  mutate(GV.Oocyte = (GV.oocyte.1+GV.oocyte.2)/2, MII.Oocyte = (MII.oocyte.1+MII.oocyte.2)/2,
         `2Cell` = (`2Cell.1`+`2Cell.1`)/2, `4Cell` = (`4Cell.1`+`4Cell.1`)/2,
         `8Cell` = (`8Cell.1`+`8Cell.1`)/2, ICM = (ICM.1+ICM.2)/2) %>%
  gather(key = "Celltype", value = "TPM") %>%
  filter(Celltype %in% c("MII.Oocyte", "2Cell", "4Cell", "8Cell", "ICM")) %>%
  mutate(Celltype = factor(Celltype, levels = c("MII.Oocyte", "2Cell", "4Cell", "8Cell", "ICM"))) %>%
  ggplot(aes(x = Celltype, y = log2(TPM+0.1))) +
  geom_boxplot(aes(fill = Celltype)) +
  scale_fill_brewer() + theme_bw() + ggtitle("REs.K93 interacted genes") +
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank(),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 11, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 11, angle = 45, hjust = 1),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 11, angle = 0)) +
  stat_compare_means(method = "t.test", na.rm = T, label.y = c(14), comparisons = list(c("4Cell", "8Cell")), size = 4)
p2 <- gse101571.ge.tpm[hs.gene.in.loops$sva_c8_RE_k273_enhancers_in_8cell$EnsemblID, ] %>%
  mutate(GV.Oocyte = (GV.oocyte.1+GV.oocyte.2)/2, MII.Oocyte = (MII.oocyte.1+MII.oocyte.2)/2,
         `2Cell` = (`2Cell.1`+`2Cell.1`)/2, `4Cell` = (`4Cell.1`+`4Cell.1`)/2,
         `8Cell` = (`8Cell.1`+`8Cell.1`)/2, ICM = (ICM.1+ICM.2)/2) %>%
  gather(key = "Celltype", value = "TPM") %>%
  filter(Celltype %in% c("MII.Oocyte", "2Cell", "4Cell", "8Cell", "ICM")) %>%
  mutate(Celltype = factor(Celltype, levels = c("MII.Oocyte", "2Cell", "4Cell", "8Cell", "ICM"))) %>%
  ggplot(aes(x = Celltype, y = log2(TPM+0.1))) +
  geom_boxplot(aes(fill = Celltype)) +
  scale_fill_brewer() + theme_bw() + ggtitle("REs.K273 interacted genes") +
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank(),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 11, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 11, angle = 45, hjust = 1),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 11, angle = 0)) +
  stat_compare_means(method = "t.test", na.rm = T, label.y = c(14), comparisons = list(c("4Cell", "8Cell")), size = 4)
p3 <- gse101571.ge.tpm[hs.gene.in.loops$sva_c8_PE_enhancers_in_8cell$EnsemblID, ] %>%
  mutate(GV.Oocyte = (GV.oocyte.1+GV.oocyte.2)/2, MII.Oocyte = (MII.oocyte.1+MII.oocyte.2)/2,
         `2Cell` = (`2Cell.1`+`2Cell.1`)/2, `4Cell` = (`4Cell.1`+`4Cell.1`)/2,
         `8Cell` = (`8Cell.1`+`8Cell.1`)/2, ICM = (ICM.1+ICM.2)/2) %>%
  gather(key = "Celltype", value = "TPM") %>%
  filter(Celltype %in% c("MII.Oocyte", "2Cell", "4Cell", "8Cell", "ICM")) %>%
  mutate(Celltype = factor(Celltype, levels = c("MII.Oocyte", "2Cell", "4Cell", "8Cell", "ICM"))) %>%
  ggplot(aes(x = Celltype, y = log2(TPM+0.1))) +
  geom_boxplot(aes(fill = Celltype)) +
  scale_fill_brewer() + theme_bw() + ggtitle("PEs interacted genes") +
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank(),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 11, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 11, angle = 45, hjust = 1),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 11, angle = 0)) +
  stat_compare_means(method = "t.test", na.rm = T, label.y = c(14), comparisons = list(c("4Cell", "8Cell")), size = 4)
plot_grid(p1, p2, p3, ncol = 3)
dev.off(); rm(p1, p2, p3)

### >>> 7. GO enrichment analysis of targeted genes
# - REs H3K9me3
write.table(na.omit(unique(hs.gene.in.loops$c8_RE_k93_enhancers_in_8cell$Symbol)), 
            "/home/cmq/bioinfo/project-cmq/embryo_93/results/R/Tables/GO/targeted_genes/REs_K9_targeted_genes.txt", quote = F, row.names = F, col.names = F)
FullSet.GO(species = "human", genelist = unique(hs.gene.in.loops$c8_RE_k93_enhancers_in_8cell$EnsemblID), basename = "REs_K9", 
           genetype = "ENSEMBL", outdir = "/home/cmq/bioinfo/project-cmq/embryo_93/results/R/Tables/GO/targeted_genes")
# - REs H3K27me3
write.table(na.omit(unique(hs.gene.in.loops$c8_RE_k273_enhancers_in_8cell$Symbol)), 
            "/home/cmq/bioinfo/project-cmq/embryo_93/results/R/Tables/GO/targeted_genes/REs_K27_targeted_genes.txt", quote = F, row.names = F, col.names = F)
FullSet.GO(species = "human", genelist = unique(hs.gene.in.loops$c8_RE_k273_enhancers_in_8cell$EnsemblID), basename = "REs_K27", 
           genetype = "ENSEMBL", outdir = "/home/cmq/bioinfo/project-cmq/embryo_93/results/R/Tables/GO/targeted_genes")
# - PEs
write.table(na.omit(unique(hs.gene.in.loops$c8_PE_enhancers_in_8cell$Symbol)), 
            "/home/cmq/bioinfo/project-cmq/embryo_93/results/R/Tables/GO/targeted_genes/PEs_targeted_genes.txt", quote = F, row.names = F, col.names = F)
FullSet.GO(species = "human", genelist = unique(hs.gene.in.loops$c8_PE_enhancers_in_8cell$EnsemblID), basename = "PEs", 
           genetype = "ENSEMBL", outdir = "/home/cmq/bioinfo/project-cmq/embryo_93/results/R/Tables/GO/targeted_genes")
# - REs H3K9me3 SVA
write.table(na.omit(unique(hs.gene.in.loops$sva_c8_RE_k93_enhancers_in_8cell$Symbol)), 
            "/home/cmq/bioinfo/project-cmq/embryo_93/results/R/Tables/GO/targeted_genes/REs_K9_SVA_targeted_genes.txt", quote = F, row.names = F, col.names = F)
FullSet.GO(species = "human", genelist = unique(hs.gene.in.loops$sva_c8_RE_k93_enhancers_in_8cell$EnsemblID), basename = "REs_K9_SVA", 
           genetype = "ENSEMBL", outdir = "/home/cmq/bioinfo/project-cmq/embryo_93/results/R/Tables/GO/targeted_genes")
# - REs H3K27me3 SVA
write.table(na.omit(unique(hs.gene.in.loops$sva_c8_RE_k273_enhancers_in_8cell$Symbol)), 
            "/home/cmq/bioinfo/project-cmq/embryo_93/results/R/Tables/GO/targeted_genes/REs_K27_SVA_targeted_genes.txt", quote = F, row.names = F, col.names = F)
FullSet.GO(species = "human", genelist = unique(hs.gene.in.loops$sva_c8_RE_k273_enhancers_in_8cell$EnsemblID), basename = "REs_K27_SVA", 
           genetype = "ENSEMBL", outdir = "/home/cmq/bioinfo/project-cmq/embryo_93/results/R/Tables/GO/targeted_genes")
# - PEs SVA
write.table(na.omit(unique(hs.gene.in.loops$sva_c8_PE_enhancers_in_8cell$Symbol)), 
            "/home/cmq/bioinfo/project-cmq/embryo_93/results/R/Tables/GO/targeted_genes/PEs_SVA_targeted_genes.txt", quote = F, row.names = F, col.names = F)
FullSet.GO(species = "human", genelist = unique(hs.gene.in.loops$sva_c8_PE_enhancers_in_8cell$EnsemblID), basename = "PEs_SVA", 
           genetype = "ENSEMBL", outdir = "/home/cmq/bioinfo/project-cmq/embryo_93/results/R/Tables/GO/targeted_genes")
# - Visualize the GO enrichment results
targeted.go <- list()
for (file in list.files("results/R/Tables/GO/targeted_genes", "*BP.txt", full.names = T)) {
  index <- grep(file, list.files("results/R/Tables/GO/targeted_genes", "*BP.txt", full.names = T))
  targeted.go[[index]] <- read.table(file, sep = "\t", header = T)
  rownames(targeted.go[[index]]) <- targeted.go[[index]]$Description
}
names(targeted.go) <- gsub("_GO_BP.txt", "", list.files("results/R/Tables/GO/targeted_genes", "*BP.txt"))
for (i in 1:length(targeted.go)) {
  pdf(file.path(outdir, paste("Bar_plot_to_show_GO_term_of_", names(targeted.go)[i], "_regulated_genes.pdf", sep = "")), height = 4, width = 8)
  print(targeted.go[[i]][1:5,] %>% 
          mutate(Log10Pvalue = -log10(pvalue)) %>%
          ggplot(aes(x = Log10Pvalue, y = fct_reorder(Description, Log10Pvalue))) + 
          labs(x = "-Log10.pvalue", y = "GO biological process") +
          geom_bar(stat = "identity", fill = "#3498db", color = "#ecf0f1", width = 0.75) + 
          theme_classic() +
          geom_vline(xintercept = c(-log10(0.05)), color = "#FF1300", linetype = "longdash", size = 1, alpha = 0.8) + 
          theme(axis.title.x = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 0),
                axis.title.y = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 90, hjust = 1), 
                axis.text.x  = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 0, hjust = 1),
                axis.text.y  = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 0, hjust = 1)))
  dev.off()
}

# - Visualize the KO-phenotype enrichment results
targeted.phenotype <- list()
for (file in list.files("results/R/Tables/GO/targeted_genes", "*table.txt", full.names = T)) {
  index <- grep(file, list.files("results/R/Tables/GO/targeted_genes", "*table.txt", full.names = T))
  targeted.phenotype[[index]] <- read.table(file, sep = "\t", header = T)
  rownames(targeted.phenotype[[index]]) <- targeted.phenotype[[index]]$Term
}
names(targeted.phenotype) <- gsub("_Level_4_2021", "",
                                  gsub("_table.txt", "",
                                       gsub("targeted_genes_", "", list.files("results/R/Tables/GO/targeted_genes", "*table.txt"))))
# - REs H3K9me3
pdf(file.path(outdir, "Bar_plot_to_show_KO_phenotype_enrichment_of_REs_K9_regulated_genes.pdf"), height = 4, width = 8)
targeted.phenotype$REs_K9_MGI_Mammalian_Phenotype[c("embryonic growth arrest MP:0001730",
                                                    "failure to gastrulate MP:0001696",
                                                    "absent inner cell mass MP:0004964",
                                                    "abnormal blastocyst morphology MP:0004957",
                                                    "failure of zygotic cell division MP:0003406"),] %>% 
  mutate(Log10Pvalue = -log10(P.value)) %>%
  ggplot(aes(x = Log10Pvalue, y = fct_reorder(Term, Log10Pvalue))) + 
  labs(x = "-Log10.pvalue", y = "MGI KO phenotype") +
  geom_bar(stat = "identity", fill = "#3498db", color = "#ecf0f1", width = 0.75) + 
  theme_classic() +
  geom_vline(xintercept = c(-log10(0.05)), color = "#FF1300", linetype = "longdash", size = 1, alpha = 0.8) + 
  theme(axis.title.x = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 90), 
        axis.text.x  = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 0, hjust = 1),
        axis.text.y  = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 0, hjust = 1))
dev.off()
pdf(file.path(outdir, "Bar_plot_to_show_KO_phenotype_enrichment_of_REs_K9_SVA_regulated_genes.pdf"), height = 4, width = 8)
targeted.phenotype$REs_K9_SVA_MGI_Mammalian_Phenotype[c("failure to gastrulate MP:0001696",
                                                        "embryonic lethality, complete penetrance MP:0011092",
                                                        "abnormal blastocyst morphology MP:0004957",
                                                        "embryonic lethality prior to organogenesis MP:0013292",
                                                        "embryonic lethality before implantation, incomplete penetrance MP:0011104"),] %>% 
  mutate(Log10Pvalue = -log10(P.value)) %>%
  ggplot(aes(x = Log10Pvalue, y = fct_reorder(Term, Log10Pvalue))) + 
  labs(x = "-Log10.pvalue", y = "MGI KO phenotype") +
  geom_bar(stat = "identity", fill = "#3498db", color = "#ecf0f1", width = 0.75) + 
  theme_classic() +
  geom_vline(xintercept = c(-log10(0.05)), color = "#FF1300", linetype = "longdash", size = 1, alpha = 0.8) + 
  theme(axis.title.x = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 90), 
        axis.text.x  = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 0, hjust = 1),
        axis.text.y  = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 0, hjust = 1))
dev.off()
# - PEs
pdf(file.path(outdir, "Bar_plot_to_show_KO_phenotype_enrichment_of_PEs_regulated_genes.pdf"), height = 4, width = 8)
targeted.phenotype$PEs_MGI_Mammalian_Phenotype[c("embryonic lethality prior to organogenesis MP:0013292",
                                                 "embryonic lethality, complete penetrance MP:0011092",
                                                 "embryonic growth arrest MP:0001730",
                                                 "failure to gastrulate MP:0001696",
                                                 "abnormal blastocyst morphology MP:0004957"),] %>% 
  mutate(Log10Pvalue = -log10(P.value)) %>%
  ggplot(aes(x = Log10Pvalue, y = fct_reorder(Term, Log10Pvalue))) + 
  labs(x = "-Log10.pvalue", y = "MGI KO phenotype") +
  geom_bar(stat = "identity", fill = "#3498db", color = "#ecf0f1", width = 0.75) + 
  theme_classic() +
  geom_vline(xintercept = c(-log10(0.05)), color = "#FF1300", linetype = "longdash", size = 1, alpha = 0.8) + 
  theme(axis.title.x = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 90), 
        axis.text.x  = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 0, hjust = 1),
        axis.text.y  = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 0, hjust = 1))
dev.off()
pdf(file.path(outdir, "Bar_plot_to_show_KO_phenotype_enrichment_of_PEs_SVA_regulated_genes.pdf"), height = 4, width = 8)
targeted.phenotype$PEs_SVA_MGI_Mammalian_Phenotype[c("decreased embryo size MP:0001698",
                                                     "failure of primitive streak formation MP:0001693",
                                                     "embryonic lethality prior to organogenesis MP:0013292",
                                                     "failure to gastrulate MP:0001696",
                                                     "embryonic lethality, complete penetrance MP:0011092"),] %>% 
  mutate(Log10Pvalue = -log10(P.value)) %>%
  ggplot(aes(x = Log10Pvalue, y = fct_reorder(Term, Log10Pvalue))) + 
  labs(x = "-Log10.pvalue", y = "MGI KO phenotype") +
  geom_bar(stat = "identity", fill = "#3498db", color = "#ecf0f1", width = 0.75) + 
  theme_classic() +
  geom_vline(xintercept = c(-log10(0.05)), color = "#FF1300", linetype = "longdash", size = 1, alpha = 0.8) + 
  theme(axis.title.x = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 90), 
        axis.text.x  = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 0, hjust = 1),
        axis.text.y  = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 0, hjust = 1))
dev.off()



### ====================================================
### 7th Part: REs and PEs Enrichment analysis on repeats
### ====================================================

### >>> 1. Set output directory
outdir <- file.path(getwd(), "results/R/Graphs/revised_Fig2")
if (! dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}

### >>> 2. Load major ZGA genes
ehc.rr <- readRDS("results/bedtools/revised_8c_enhancer/RegionR/name_25_enrichment_analysis_data.rds")
names(ehc.rr) <- c("PE", "RE273", "RE93")
for(i in seq(1, length(ehc.rr))){
  ehc.rr[[i]]$FoldChange <- ehc.rr[[i]]$Observed/ehc.rr[[i]]$Expected
  ehc.rr[[i]] <- ehc.rr[[i]][, -5]
}; rm(i)
ehc.rr <- do.call(rbind, ehc.rr)
ehc.rr <- cbind(ehc.rr, str_split_fixed(rownames(ehc.rr), "\\.", 2))
colnames(ehc.rr)[6:7] <- c("Region", "Name")


### >>> 3. Make a plot
# PE
pd <- subset(ehc.rr, Pvalue<=0.05 & Region=="PE") %>% 
  arrange(desc(FoldChange)) %>% 
  filter(Zscore > 0) %>% 
  mutate(Group = case_when(FoldChange >= 10 | Zscore >= 30 ~ "Most enriched",
                           FoldChange < 10 | Zscore < 30 ~ "Less enriched"))
pdf(file.path(outdir, "Scatter_plot_to_show_enrichment_analysis_results_of_repeats_in_PEs.pdf"), height = 3.5, width = 5)
pd %>% 
  ggplot(aes(x = FoldChange, y = Zscore)) +
  geom_point(aes(size = Observed, fill = Group), shape = 21) +
  scale_fill_manual(values = c(`Most enriched` = "#F53C3C", `Less enriched` = "#30D2DF")) +
  theme_bw() +
  labs(x = "Observed / Expected", y = "Enrichment score", fill = "Group", size = "Expected value") +
  geom_label_repel(aes(label = Name),
                   data = subset(pd, FoldChange >= 10 | Zscore >= 30),
                   box.padding = 0.2, point.padding = 0.2, segment.color = "#130f40",
                   label.size = 0.5, fill = NA, size = 3) +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 13, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 13, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 13, angle = 0),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 13, angle = 0),
        panel.grid.minor = element_blank())
dev.off()
# REs H3K9me3
pd <- subset(ehc.rr, Pvalue<=0.05 & Region=="RE93") %>% 
  arrange(desc(FoldChange)) %>% 
  filter(Zscore > 0) %>% 
  mutate(Group = case_when(FoldChange >= 10 | Zscore >= 30 ~ "Most enriched",
                           FoldChange < 10 | Zscore < 30 ~ "Less enriched"))
pdf(file.path(outdir, "Scatter_plot_to_show_enrichment_analysis_results_of_repeats_in_REs_H3K9me3.pdf"), height = 3.5, width = 5)
pd %>% 
  ggplot(aes(x = FoldChange, y = Zscore)) +
  geom_point(aes(size = Observed, fill = Group), shape = 21) +
  scale_fill_manual(values = c(`Most enriched` = "#F53C3C", `Less enriched` = "#30D2DF")) +
  labs(x = "Observed / Expected", y = "Enrichment score", fill = "Group", size = "Expected value") +
  theme_bw() +
  geom_label_repel(aes(label = Name),
                   data = subset(pd, FoldChange >= 10 | Zscore >= 30),
                   box.padding = 0.2, point.padding = 0.2, segment.color = "#130f40",
                   label.size = 0.5, fill = NA, size = 3) +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 13, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 13, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 13, angle = 0),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 13, angle = 0),
        panel.grid.minor = element_blank())
dev.off()
# REs H3K27me3
pd <- subset(ehc.rr, Pvalue<=0.05 & Region=="RE273") %>% 
  arrange(desc(FoldChange)) %>% 
  filter(Zscore > 0) %>% 
  mutate(Group = case_when(FoldChange >= 10 | Zscore >= 30 ~ "Most enriched",
                           FoldChange < 10 | Zscore < 30 ~ "Less enriched"))
pdf(file.path(outdir, "Scatter_plot_to_show_enrichment_analysis_results_of_repeats_in_REs_H3K27me3.pdf"), height = 3.5, width = 5)
pd %>% 
  ggplot(aes(x = FoldChange, y = Zscore)) +
  geom_point(aes(size = Observed, fill = Group), shape = 21) +
  scale_fill_manual(values = c(`Most enriched` = "#F53C3C", `Less enriched` = "#30D2DF")) +
  theme_bw() +
  labs(x = "Observed / Expected", y = "Enrichment score", fill = "Group", size = "Expected value") +
  geom_label_repel(aes(label = Name),
                   data = subset(pd, FoldChange >= 10 | Zscore >= 30),
                   box.padding = 0.2, point.padding = 0.2, segment.color = "#130f40",
                   label.size = 0.5, fill = NA, size = 3) +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 13, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 13, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 13, angle = 0),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 13, angle = 0),
        panel.grid.minor = element_blank())
dev.off()



### =========================================
### 8th Part: Detect selection on REs and PEs
### =========================================

### >>> 1. Generate null sequences (negative set)
outdir <- file.path(getwd(), "results/R/Graphs/revised_Fig2")
if (! dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}


### >>> 2. Generate null sequences (negative set)
# - load genome
hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
# - 
c8.enhancer <- "/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer/8cell_enhancer.bed"
c8.enhancer.neg <- "/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer/8cell_enhancer_negSet.bed"
posSet.fa <- "/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer/8cell_enhancer_posSet.fa"
negSet.fa <- "/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer/8cell_enhancer_negSet.fa"
genNullSeqs(c8.enhancer, 
            outputBedFN = c8.enhancer.neg, 
            outputPosFastaFN = posSet.fa,
            outputNegFastaFN = negSet.fa, 
            xfold = 1, 
            repeat_match_tol = 0.02, 
            GC_match_tol = 0.02, 
            length_match_tol = 0.02, 
            batchsize = 5000, 
            nMaxTrials = 20, 
            genome = hg38)


### >>> 3. Distribution of deltaSVMs
# load deltaSVMs data
pe.re.svms <- list()
for (files in list.files("/home/yhw/bioinfo/project-mine/Embryo.93/evo/revised_8c_enhancer", "*output.txt", full.names = T)) {
  index <- grep(files, list.files("/home/yhw/bioinfo/project-mine/Embryo.93/evo/revised_8c_enhancer", "*output.txt", full.names = T))
  pe.re.svms[[index]] <- read.table(files)
  names(pe.re.svms)[index] <- gsub("_output_in_enhancers.txt", "", gsub("detalSVM_", "", basename(files)))
}; rm(files, index)
# load control
svms.pvalue <- list()
for (files in list.files("/home/yhw/bioinfo/project-mine/Embryo.93/evo/revised_8c_enhancer", "*(0|1).output$", full.names = T)) {
  index <- grep(files, list.files("/home/yhw/bioinfo/project-mine/Embryo.93/evo/revised_8c_enhancer", "*(0|1).output$", full.names = T))
  svms.pvalue[[index]] <- read.table(files)
  names(svms.pvalue)[index] <- gsub("_random_contrl", "", gsub(".processed.output", "", basename(files)))
}; rm(files, index)
# - plotting deltaSVMs (REs H3K9me3)
pdf(file.path(outdir, "Evo_selection_deltaSVMs_distribution_in_REs_H3K9me3.pdf"), width = 4, height = 3.5)
sample.n <- sample(1:length(pe.re.svms$REs_H3K9me3_output.txt$V2), length(svms.pvalue$REs_H3K9me3_10000.output$V2))
data.frame(deltasvm = c(-pe.re.svms$REs_H3K9me3_output.txt$V2[sample.n], 
                        -svms.pvalue$REs_H3K9me3_10000.output$V2),
           Group = c(rep("REs.K9", length(pe.re.svms$REs_H3K9me3_output.txt$V2[sample.n])),
                     rep("Ctrl", length(svms.pvalue$REs_H3K9me3_10000.output$V2)))) %>% 
  ggplot(aes(deltasvm)) + 
  geom_histogram(aes(fill = Group), binwidth = 0.5, color = "#000000", alpha = 0.8, lwd = 0.35) +
  scale_fill_manual(values = c("REs.K9" = "#df233d", "Ctrl" = "#3093dc")) +
  geom_vline(xintercept = c(0), colour = "#d3df23", linetype = "longdash", size = 0.5, alpha = 1) +
  annotate("text", x = 10, y = 9000, label = paste("N = ", sum(-pe.re.svms$REs_H3K9me3_output.txt$V2[sample.n]>=0), sep = ""), color = "#df233d", size = 5) +
  annotate("text", x = -10, y = 9000, label = paste("N = ", sum(-pe.re.svms$REs_H3K9me3_output.txt$V2[sample.n]<0), sep = ""), color = "#df233d", size = 5) +
  annotate("text", x = 10, y = 6000, label = paste("N = ", sum(-svms.pvalue$REs_H3K9me3_10000.output$V2>=0), sep = ""), color = "#3093dc", size = 5) +
  annotate("text", x = -10, y = 6000, label = paste("N = ", sum(-svms.pvalue$REs_H3K9me3_10000.output$V2<0), sep = ""), color = "#3093dc", size = 5) +
  labs(x = "deltaSVMs", y = "Frequency") + theme_bw() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        panel.grid.minor = element_blank(), legend.position = "top")
dev.off()
# - plotting deltaSVMs (REs H3K27me3)
pdf(file.path(outdir, "Evo_selection_deltaSVMs_distribution_in_REs_H3K27me3.pdf"), width = 4, height = 3.5)
sample.n <- sample(1:length(pe.re.svms$REs_H3K27me3_output.txt$V2), length(svms.pvalue$REs_H3K27me3_10000.output$V2))
data.frame(deltasvm = c(-pe.re.svms$REs_H3K27me3_output.txt$V2[sample.n], 
                        -svms.pvalue$REs_H3K27me3_10000.output$V2),
           Group = c(rep("REs.K27", length(pe.re.svms$REs_H3K27me3_output.txt$V2[sample.n])),
                     rep("Ctrl", length(svms.pvalue$REs_H3K27me3_10000.output$V2)))) %>% 
  ggplot(aes(deltasvm)) + 
  geom_histogram(aes(fill = Group), binwidth = 0.5, color = "#000000", alpha = 0.8, lwd = 0.35) +
  scale_fill_manual(values = c("REs.K27" = "#df233d", "Ctrl" = "#3093dc")) +
  geom_vline(xintercept = c(0), colour = "#d3df23", linetype = "longdash", size = 0.5, alpha = 1) +
  annotate("text", x = 10, y = 3000, label = paste("N = ", sum(-pe.re.svms$REs_H3K27me3_output.txt$V2[sample.n]>=0), sep = ""), color = "#df233d", size = 5) +
  annotate("text", x = -10, y = 3000, label = paste("N = ", sum(-pe.re.svms$REs_H3K27me3_output.txt$V2[sample.n]<0), sep = ""), color = "#df233d", size = 5) +
  annotate("text", x = 10, y = 2000, label = paste("N = ", sum(-svms.pvalue$REs_H3K27me3_10000.output$V2>=0), sep = ""), color = "#3093dc", size = 5) +
  annotate("text", x = -10, y = 2000, label = paste("N = ", sum(-svms.pvalue$REs_H3K27me3_10000.output$V2<0), sep = ""), color = "#3093dc", size = 5) +
  labs(x = "deltaSVMs", y = "Frequency") + theme_bw() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        panel.grid.minor = element_blank(), legend.position = "top")
dev.off()
# - plotting deltaSVMs (PEs)
pdf(file.path(outdir, "Evo_selection_deltaSVMs_distribution_in_PEs.pdf"), width = 4, height = 3.5)
sample.n <- sample(1:length(pe.re.svms$PEs_output.txt$V2), length(svms.pvalue$PEs_10000.output$V2))
data.frame(deltasvm = c(-pe.re.svms$PEs_output.txt$V2[sample.n], 
                        -svms.pvalue$PEs_10000.output$V2),
           Group = c(rep("PEs", length(pe.re.svms$PEs_output.txt$V2[sample.n])),
                     rep("Ctrl", length(svms.pvalue$PEs_10000.output$V2)))) %>% 
  ggplot(aes(deltasvm)) + 
  geom_histogram(aes(fill = Group), binwidth = 0.5, color = "#000000", alpha = 0.8, lwd = 0.35) +
  scale_fill_manual(values = c("PEs" = "#df233d", "Ctrl" = "#3093dc")) +
  geom_vline(xintercept = c(0), colour = "#d3df23", linetype = "longdash", size = 0.5, alpha = 1) +
  annotate("text", x = 10, y = 1500, label = paste("N = ", sum(-pe.re.svms$PEs_output.txt$V2[sample.n]>=0), sep = ""), color = "#df233d", size = 5) +
  annotate("text", x = -10, y = 1500, label = paste("N = ", sum(-pe.re.svms$PEs_output.txt$V2[sample.n]<0), sep = ""), color = "#df233d", size = 5) +
  annotate("text", x = 10, y = 1000, label = paste("N = ", sum(-svms.pvalue$PEs_10000.output$V2>=0), sep = ""), color = "#3093dc", size = 5) +
  annotate("text", x = -10, y = 1000, label = paste("N = ", sum(-svms.pvalue$PEs_10000.output$V2<0), sep = ""), color = "#3093dc", size = 5) +
  labs(x = "deltaSVMs", y = "Frequency") + theme_bw() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        panel.grid.minor = element_blank(), legend.position = "top")
dev.off()
# - plotting pvalue (REs H3K9me3)
svms.pvalue <- list()
for (files in list.files("/home/yhw/bioinfo/project-mine/Embryo.93/evo/revised_8c_enhancer", "*processed.output$", full.names = T)) {
  index <- grep(files, list.files("/home/yhw/bioinfo/project-mine/Embryo.93/evo/revised_8c_enhancer", "*processed.output$", full.names = T))
  svms.pvalue[[index]] <- read.table(files) %>% group_by(V7) %>% 
      summarise(pvalue = 1-V16)
  names(svms.pvalue)[index] <- gsub("_random_contrl", "", gsub(".processed.output", "", basename(files)))
}; rm(files, index)
pdf(file.path(outdir, "Evo_selection_pvalue_for_deltaSVMs_distribution_in_REs_H3K9me3.pdf"), width = 5, height = 3.5)
pos.ratio <- paste("Negative sites:\nRatio (pvalue <= 0.05) = ", 
                   round(sum(svms.pvalue$REs_H3K9me3_10000$pvalue<=0.05)/
                           length(svms.pvalue$REs_H3K9me3_10000$pvalue), 2), sep = "")
svms.pvalue$REs_H3K9me3_10000 %>% 
  ggplot(aes(x = pvalue)) + 
  geom_histogram(binwidth = 0.01, color = "#000000", fill = "#3093dc", lwd = 0.35) +
  geom_vline(xintercept = 0.05, colour = "#f32d3c", linetype = "longdash", size = 0.7, alpha = 1) +
  annotate("text", x = 0.5, y = 4000, label = pos.ratio, color = "#000000", size = 5) +
  labs(x = "Pvalue", y = "Frequency") + theme_bw() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        panel.grid.minor = element_blank())
dev.off(); rm(pos.ratio)
# - plotting pvalue (REs H3K27me3)
pdf(file.path(outdir, "Evo_selection_pvalue_for_deltaSVMs_distribution_in_REs_H3K27me3.pdf"), width = 5, height = 3.5)
pos.ratio <- paste("Negative sites:\nRatio (pvalue <= 0.05) = ", 
                   round(sum(svms.pvalue$REs_H3K27me3_10000$pvalue<=0.05)/
                           length(svms.pvalue$REs_H3K27me3_10000$pvalue), 2), sep = "")
svms.pvalue$REs_H3K27me3_10000 %>% 
  ggplot(aes(x = pvalue)) + 
  geom_histogram(binwidth = 0.01, color = "#000000", fill = "#3093dc", lwd = 0.35) +
  geom_vline(xintercept = 0.05, colour = "#f32d3c", linetype = "longdash", size = 0.7, alpha = 1) +
  annotate("text", x = 0.5, y = 1500, label = pos.ratio, color = "#000000", size = 5) +
  labs(x = "Pvalue", y = "Frequency") + theme_bw() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        panel.grid.minor = element_blank())
dev.off(); rm(pos.ratio)
# - plotting pvalue (PEs)
pdf(file.path(outdir, "Evo_selection_pvalue_for_deltaSVMs_distribution_in_PEs.pdf"), width = 5, height = 3.5)
pos.ratio <- paste("Negative sites:\nRatio (pvalue <= 0.05) = ", 
                   round(sum(svms.pvalue$PEs_10000$pvalue<=0.05)/
                           length(svms.pvalue$PEs_10000$pvalue), 2), sep = "")
svms.pvalue$PEs_10000 %>% 
  ggplot(aes(x = pvalue)) + 
  geom_histogram(binwidth = 0.01, color = "#000000", fill = "#3093dc", lwd = 0.35) +
  geom_vline(xintercept = 0.05, colour = "#f32d3c", linetype = "longdash", size = 0.7, alpha = 1) +
  annotate("text", x = 0.5, y = 500, label = pos.ratio, color = "#000000", size = 5) +
  labs(x = "Pvalue", y = "Frequency") + theme_bw() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        panel.grid.minor = element_blank())
dev.off(); rm(pos.ratio)


### >>> 4. Allele frequency
# - load deltaSVMs
pe.re.svms <- list()
for (files in list.files("/home/yhw/bioinfo/project-mine/Embryo.93/evo/revised_8c_enhancer", "*output_in_enhancers.txt", full.names = T)) {
  index <- grep(files, list.files("/home/yhw/bioinfo/project-mine/Embryo.93/evo/revised_8c_enhancer", "*output_in_enhancers.txt", full.names = T))
  pe.re.svms[[index]] <- read.table(files)
  names(pe.re.svms)[index] <- gsub("_output_in_enhancers.txt", "", gsub("detalSVM_", "", basename(files)))
}; rm(files, index)
# - load allele frequency
pe.re.alle <- list()
for (files in list.files("/home/yhw/bioinfo/project-mine/Embryo.93/evo/revised_8c_enhancer", "*allele_freq*", full.names = T)) {
  index <- grep(files, list.files("/home/yhw/bioinfo/project-mine/Embryo.93/evo/revised_8c_enhancer", "*allele_freq*", full.names = T))
  pe.re.alle[[index]] <- read.table(files)
  pe.re.alle[[index]]$pos.id <- paste(pe.re.alle[[index]]$V1, ":", pe.re.alle[[index]]$V2, "-", pe.re.alle[[index]]$V3, sep = "")
  names(pe.re.alle)[index] <- gsub("allele_freq_", "", str_split_fixed(basename(files), "\\.", 2)[,1])
}; rm(files, index)
# plotting for all population
pd <- pe.re.svms$REs_H3K9me3 %>% group_by(V7) %>% 
  summarise(snvs.count = length(V5),
            pos.count = sum(V14 > 0),
            pos.ratio = sum(V14 > 0)/length(V14)) %>% as.data.frame() %>% 
  mutate(group = case_when(pos.ratio >= 0.50 ~ "Positive sites",
                           pos.ratio <= 0.25 ~ "Negative sites")) %>% na.omit()
table(pd$group)
pd <- pe.re.alle$REs_H3K9me3_in_ALL %>% 
  filter(V16 > 1000) %>% group_by(V7) %>% 
  summarise(value = mean(V16/V14)) %>% as.data.frame() %>% 
  mutate(group = case_when(V7 %in% subset(pd, group == "Positive sites")$V7 ~ "Positive sites",
                           V7 %in% subset(pd, group == "Negative sites")$V7 ~ "Negative sites")) %>% na.omit() %>% 
  mutate(group = factor(group, levels = c("Positive sites", "Negative sites")))
pdf(file.path(outdir, "Evo_selection_allele_frequency_in_ALL_population_in_REs_H3K9me3_positive_and_negative_sites.pdf"), width = 2.5, height = 3.5)
do.call("rbind", replicate(1, pd, simplify = FALSE)) %>% 
  ggplot(aes(x = group, y = value)) +
  geom_boxplot(aes(fill = group), width  = 0.6) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  labs(x = "", y = "Ratio (allele frequeny / population size)") + theme_bw() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        legend.position = "none", panel.grid.minor = element_blank()) +
  stat_compare_means(method = "t.test", na.rm = T, label.y = 0.15, 
                     comparisons = list(c("Positive sites", "Negative sites")), size = 4)
dev.off()
pd <- pe.re.svms$REs_H3K27me3 %>% group_by(V7) %>% 
  summarise(snvs.count = length(V5),
            pos.count = sum(V14 > 0),
            pos.ratio = sum(V14 > 0)/length(V14)) %>% as.data.frame() %>% 
  mutate(group = case_when(pos.ratio >= 0.50 ~ "Positive sites",
                           pos.ratio <= 0.25 ~ "Negative sites")) %>% na.omit()
table(pd$group)
pd <- pe.re.alle$REs_H3K27me3_in_ALL %>% 
  filter(V16 > 1000) %>% group_by(V7) %>% 
  summarise(value = mean(V16/V14)) %>% as.data.frame() %>% 
  mutate(group = case_when(V7 %in% subset(pd, group == "Positive sites")$V7 ~ "Positive sites",
                           V7 %in% subset(pd, group == "Negative sites")$V7 ~ "Negative sites")) %>% na.omit() %>% 
  mutate(group = factor(group, levels = c("Positive sites", "Negative sites")))
pdf(file.path(outdir, "Evo_selection_allele_frequency_in_ALL_population_in_REs_H3K27me3_positive_and_negative_sites.pdf"), width = 2.5, height = 3.5)
do.call("rbind", replicate(1, pd, simplify = FALSE)) %>% 
  ggplot(aes(x = group, y = value)) +
  geom_boxplot(aes(fill = group), width  = 0.6) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  labs(x = "", y = "Ratio (allele frequeny / population size)") + theme_bw() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        legend.position = "none", panel.grid.minor = element_blank()) +
  stat_compare_means(method = "t.test", na.rm = T, label.y = 0.15, 
                     comparisons = list(c("Positive sites", "Negative sites")), size = 4)
dev.off()
pd <- pe.re.svms$PEs %>% group_by(V7) %>% 
  summarise(snvs.count = length(V5),
            pos.count = sum(V14 > 0),
            pos.ratio = sum(V14 > 0)/length(V14)) %>% as.data.frame() %>% 
  mutate(group = case_when(pos.ratio >= 0.50 ~ "Positive sites",
                           pos.ratio <= 0.25 ~ "Negative sites")) %>% na.omit()
table(pd$group)
pd <- pe.re.alle$PEs_in_ALL %>% 
  filter(V16 > 1000) %>% group_by(V7) %>% 
  summarise(value = mean(V16/V14)) %>% as.data.frame() %>% 
  mutate(group = case_when(V7 %in% subset(pd, group == "Positive sites")$V7 ~ "Positive sites",
                           V7 %in% subset(pd, group == "Negative sites")$V7 ~ "Negative sites")) %>% na.omit() %>% 
  mutate(group = factor(group, levels = c("Positive sites", "Negative sites")))
pdf(file.path(outdir, "Evo_selection_allele_frequency_in_ALL_population_in_PEs_positive_and_negative_sites.pdf"), width = 2.5, height = 3.5)
do.call("rbind", replicate(1, pd, simplify = FALSE)) %>% 
  ggplot(aes(x = group, y = value)) +
  geom_boxplot(aes(fill = group), width  = 0.6) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  labs(x = "", y = "Ratio (allele frequeny / population size)") + theme_bw() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        legend.position = "none", panel.grid.minor = element_blank()) +
  stat_compare_means(method = "t.test", na.rm = T, label.y = 0.15, 
                     comparisons = list(c("Positive sites", "Negative sites")), size = 4)
dev.off()
# plotting for CHB
pd <- pe.re.svms$REs_H3K9me3 %>% group_by(V7) %>% 
  summarise(snvs.count = length(V5),
            pos.count = sum(V14 > 0),
            pos.ratio = sum(V14 > 0)/length(V14)) %>% as.data.frame() %>% 
  mutate(group = case_when(pos.ratio >= 0.50 ~ "Positive sites",
                           pos.ratio <= 0.25 ~ "Negative sites")) %>% na.omit()
table(pd$group)
pd <- pe.re.alle$REs_H3K9me3_in_CHB %>% 
  filter(V16 > 50) %>% group_by(V7) %>% 
  summarise(value = mean(V16/V14)) %>% as.data.frame() %>% 
  mutate(group = case_when(V7 %in% subset(pd, group == "Positive sites")$V7 ~ "Positive sites",
                           V7 %in% subset(pd, group == "Negative sites")$V7 ~ "Negative sites")) %>% na.omit() %>% 
  mutate(group = factor(group, levels = c("Positive sites", "Negative sites")))
pdf(file.path(outdir, "Evo_selection_allele_frequency_in_Chinese_CHB_population_in_REs_H3K9me3_positive_and_negative_sites.pdf"), width = 2.5, height = 3.5)
do.call("rbind", replicate(1, pd, simplify = FALSE)) %>% 
  ggplot(aes(x = group, y = value)) +
  geom_boxplot(aes(fill = group), width  = 0.6) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  labs(x = "", y = "Ratio (allele frequeny / population size)") + theme_bw() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        legend.position = "none", panel.grid.minor = element_blank()) +
  stat_compare_means(method = "t.test", na.rm = T, label.y = 0.15, 
                     comparisons = list(c("Positive sites", "Negative sites")), size = 4)
dev.off()
pd <- pe.re.svms$REs_H3K27me3 %>% group_by(V7) %>% 
  summarise(snvs.count = length(V5),
            pos.count = sum(V14 > 0),
            pos.ratio = sum(V14 > 0)/length(V14)) %>% as.data.frame() %>% 
  mutate(group = case_when(pos.ratio >= 0.50 ~ "Positive sites",
                           pos.ratio <= 0.25 ~ "Negative sites")) %>% na.omit()
table(pd$group)
pd <- pe.re.alle$REs_H3K27me3_in_CHB %>% 
  filter(V16 > 50) %>% group_by(V7) %>% 
  summarise(value = mean(V16/V14)) %>% as.data.frame() %>% 
  mutate(group = case_when(V7 %in% subset(pd, group == "Positive sites")$V7 ~ "Positive sites",
                           V7 %in% subset(pd, group == "Negative sites")$V7 ~ "Negative sites")) %>% na.omit() %>% 
  mutate(group = factor(group, levels = c("Positive sites", "Negative sites")))
pdf(file.path(outdir, "Evo_selection_allele_frequency_in_Chinese_CHB_population_in_REs_H3K27me3_positive_and_negative_sites.pdf"), width = 2.5, height = 3.5)
do.call("rbind", replicate(1, pd, simplify = FALSE)) %>% 
  ggplot(aes(x = group, y = value)) +
  geom_boxplot(aes(fill = group), width  = 0.6) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  labs(x = "", y = "Ratio (allele frequeny / population size)") + theme_bw() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        legend.position = "none", panel.grid.minor = element_blank()) +
  stat_compare_means(method = "t.test", na.rm = T, label.y = 0.15, 
                     comparisons = list(c("Positive sites", "Negative sites")), size = 4)
dev.off()
pd <- pe.re.svms$PEs %>% group_by(V7) %>% 
  summarise(snvs.count = length(V5),
            pos.count = sum(V14 > 0),
            pos.ratio = sum(V14 > 0)/length(V14)) %>% as.data.frame() %>% 
  mutate(group = case_when(pos.ratio >= 0.50 ~ "Positive sites",
                           pos.ratio <= 0.25 ~ "Negative sites")) %>% na.omit()
table(pd$group)
pd <- pe.re.alle$PEs_in_CHB %>% 
  filter(V16 > 50) %>% group_by(V7) %>% 
  summarise(value = mean(V16/V14)) %>% as.data.frame() %>% 
  mutate(group = case_when(V7 %in% subset(pd, group == "Positive sites")$V7 ~ "Positive sites",
                           V7 %in% subset(pd, group == "Negative sites")$V7 ~ "Negative sites")) %>% na.omit() %>% 
  mutate(group = factor(group, levels = c("Positive sites", "Negative sites")))
pdf(file.path(outdir, "Evo_selection_allele_frequency_in_Chinese_CHB_population_in_PEs_positive_and_negative_sites.pdf"), width = 2.5, height = 3.5)
do.call("rbind", replicate(1, pd, simplify = FALSE)) %>% 
  ggplot(aes(x = group, y = value)) +
  geom_boxplot(aes(fill = group), width  = 0.6) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  labs(x = "", y = "Ratio (allele frequeny / population size)") + theme_bw() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        legend.position = "none", panel.grid.minor = element_blank()) +
  stat_compare_means(method = "t.test", na.rm = T, label.y = 0.15, 
                     comparisons = list(c("Positive sites", "Negative sites")), size = 4)
dev.off()
# plotting for CHS
pd <- pe.re.svms$REs_H3K9me3 %>% group_by(V7) %>% 
  summarise(snvs.count = length(V5),
            pos.count = sum(V14 > 0),
            pos.ratio = sum(V14 > 0)/length(V14)) %>% as.data.frame() %>% 
  mutate(group = case_when(pos.ratio >= 0.50 ~ "Positive sites",
                           pos.ratio <= 0.25 ~ "Negative sites")) %>% na.omit()
table(pd$group)
pd <- pe.re.alle$REs_H3K9me3_in_CHS %>% 
  filter(V16 > 60) %>% group_by(V7) %>% 
  summarise(value = mean(V16/V14)) %>% as.data.frame() %>% 
  mutate(group = case_when(V7 %in% subset(pd, group == "Positive sites")$V7 ~ "Positive sites",
                           V7 %in% subset(pd, group == "Negative sites")$V7 ~ "Negative sites")) %>% na.omit() %>% 
  mutate(group = factor(group, levels = c("Positive sites", "Negative sites")))
pdf(file.path(outdir, "Evo_selection_allele_frequency_in_Chinese_CHS_population_in_REs_H3K9me3_positive_and_negative_sites.pdf"), width = 2.5, height = 3.5)
do.call("rbind", replicate(1, pd, simplify = FALSE)) %>% 
  ggplot(aes(x = group, y = value)) +
  geom_boxplot(aes(fill = group), width  = 0.6) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  labs(x = "", y = "Ratio (allele frequeny / population size)") + theme_bw() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        legend.position = "none", panel.grid.minor = element_blank()) +
  stat_compare_means(method = "t.test", na.rm = T, label.y = 0.15, 
                     comparisons = list(c("Positive sites", "Negative sites")), size = 4)
dev.off()
pd <- pe.re.svms$REs_H3K27me3 %>% group_by(V7) %>% 
  summarise(snvs.count = length(V5),
            pos.count = sum(V14 > 0),
            pos.ratio = sum(V14 > 0)/length(V14)) %>% as.data.frame() %>% 
  mutate(group = case_when(pos.ratio >= 0.50 ~ "Positive sites",
                           pos.ratio <= 0.25 ~ "Negative sites")) %>% na.omit()
table(pd$group)
pd <- pe.re.alle$REs_H3K27me3_in_CHS %>% 
  filter(V16 > 60) %>% group_by(V7) %>% 
  summarise(value = mean(V16/V14)) %>% as.data.frame() %>% 
  mutate(group = case_when(V7 %in% subset(pd, group == "Positive sites")$V7 ~ "Positive sites",
                           V7 %in% subset(pd, group == "Negative sites")$V7 ~ "Negative sites")) %>% na.omit() %>% 
  mutate(group = factor(group, levels = c("Positive sites", "Negative sites")))
pdf(file.path(outdir, "Evo_selection_allele_frequency_in_Chinese_CHS_population_in_REs_H3K27me3_positive_and_negative_sites.pdf"), width = 2.5, height = 3.5)
do.call("rbind", replicate(1, pd, simplify = FALSE)) %>% 
  ggplot(aes(x = group, y = value)) +
  geom_boxplot(aes(fill = group), width  = 0.6) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  labs(x = "", y = "Ratio (allele frequeny / population size)") + theme_bw() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        legend.position = "none", panel.grid.minor = element_blank()) +
  stat_compare_means(method = "t.test", na.rm = T, label.y = 0.15, 
                     comparisons = list(c("Positive sites", "Negative sites")), size = 4)
dev.off()
pd <- pe.re.svms$PEs %>% group_by(V7) %>% 
  summarise(snvs.count = length(V5),
            pos.count = sum(V14 > 0),
            pos.ratio = sum(V14 > 0)/length(V14)) %>% as.data.frame() %>% 
  mutate(group = case_when(pos.ratio >= 0.50 ~ "Positive sites",
                           pos.ratio <= 0.25 ~ "Negative sites")) %>% na.omit()
table(pd$group)
pd <- pe.re.alle$PEs_in_CHS %>% 
  filter(V16 > 60) %>% group_by(V7) %>% 
  summarise(value = mean(V16/V14)) %>% as.data.frame() %>% 
  mutate(group = case_when(V7 %in% subset(pd, group == "Positive sites")$V7 ~ "Positive sites",
                           V7 %in% subset(pd, group == "Negative sites")$V7 ~ "Negative sites")) %>% na.omit() %>% 
  mutate(group = factor(group, levels = c("Positive sites", "Negative sites")))
pdf(file.path(outdir, "Evo_selection_allele_frequency_in_Chinese_CHS_population_in_PEs_positive_and_negative_sites.pdf"), width = 2.5, height = 3.5)
do.call("rbind", replicate(1, pd, simplify = FALSE)) %>% 
  ggplot(aes(x = group, y = value)) +
  geom_boxplot(aes(fill = group), width  = 0.6) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  labs(x = "", y = "Ratio (allele frequeny / population size)") + theme_bw() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        legend.position = "none", panel.grid.minor = element_blank()) +
  stat_compare_means(method = "t.test", na.rm = T, label.y = 0.15, 
                     comparisons = list(c("Positive sites", "Negative sites")), size = 4)
dev.off()


### >>> 5. SNP density + Recombination rate
# - define positive and negative sites
pd <- pe.re.svms$REs_H3K9me3 %>% 
  mutate(pos.id = paste(V1, ":", V2, "-", V3, sep = "")) %>% 
  group_by(pos.id) %>% 
  summarise(pos.ratio = sum(V14 > 0)/length(V14)) %>% 
  as.data.frame() %>% 
  mutate(group = case_when(pos.ratio >= 0.35 ~ "Positive sites",
                           pos.ratio <= 0.3 ~ "Negative sites")) %>% 
  na.omit()
# - load data
resk9.sva.snpden <- ParseSNPdensity("/home/yhw/bioinfo/project-mine/Embryo.93/evo/revised_8c_enhancer/snpdensity/REs_K9_sva")
resk9.sva.rate <- ParseRate("/home/yhw/bioinfo/project-mine/Embryo.93/evo/revised_8c_enhancer/snpdensity/REs_K9_sva")
resk27.sva.snpden <- ParseSNPdensity("/home/yhw/bioinfo/project-mine/Embryo.93/evo/revised_8c_enhancer/snpdensity/REs_K27_sva")
resk27.sva.rate <- ParseRate("/home/yhw/bioinfo/project-mine/Embryo.93/evo/revised_8c_enhancer/snpdensity/REs_K27_sva")
pes.sva.snpden <- ParseSNPdensity("/home/yhw/bioinfo/project-mine/Embryo.93/evo/revised_8c_enhancer/snpdensity/PEs_sva")
pes.sva.rate <- ParseRate("/home/yhw/bioinfo/project-mine/Embryo.93/evo/revised_8c_enhancer/snpdensity/PEs_sva")
# - plotting for REs H3K9me3
pd.pos <- resk9.sva.snpden$region$region.id[resk9.sva.snpden$region$pos.id %in% subset(pd, group == "Positive sites")$pos.id]
pd.neg <- resk9.sva.snpden$region$region.id[resk9.sva.snpden$region$pos.id %in% subset(pd, group == "Negative sites")$pos.id]
resk9.sva.snpden$value %>% 
  rownames_to_column(var = "pos.id") %>% 
  mutate(Group = case_when(pos.id %in% pd.pos & bins_21 <= 0.02 ~ "Positive site",
                           pos.id %in% pd.neg ~ "Negative site")) %>% 
  na.omit() %>% 
  gather(key = "bins", value = "ratio", -Group, -pos.id) %>% 
  group_by(Group) %>% arrange(Group) %>% 
  mutate(bins = gsub("bins_", "", bins)) %>% 
  mutate(bins = as.numeric(bins)) -> pd.snp
table(pd.snp$Group, pd.snp$bins)
resk9.sva.rate$value %>% 
  rownames_to_column(var = "pos.id") %>% 
  mutate(Group = case_when(pos.id %in% pd.pos ~ "Positive site",
                           pos.id %in% pd.neg ~ "Negative site")) %>% 
  na.omit() %>% 
  gather(key = "bins", value = "ratio", -Group, -pos.id) %>% 
  group_by(Group) %>% arrange(Group) %>% 
  mutate(bins = gsub("bins_", "", bins)) %>% 
  mutate(bins = as.numeric(bins)) -> pd.rate
pd.snp$Data <- "SNP density"
pd.rate$Data <- "Recombination rate"
rbind(pd.snp, pd.rate) -> pd
pdf(file.path(outdir, "SNP_density_and_Recombination_rate_in_positive_and_negative_sites_around_REs_K9_SVAs_10kb.pdf"), width = 4.5, height = 3.5)
pd %>% 
  ggplot(aes(x = bins, y = ratio)) +
  stat_smooth(method = "loess", span = 0.05, se = T, alpha = 0.3, aes(color = Group, fill = Group), lwd = 1.25, level = 0.25) +
  scale_fill_manual(values = c("Positive site" = "#df233d", "Negative site" = "#3093dc")) + 
  scale_color_manual(values = c("Positive site" = "#df233d", "Negative site" = "#3093dc")) + 
  geom_vline(xintercept = 21, colour = "#000000", linetype = "longdash", size = 0.75, alpha = 0.75) +
  labs(x = "SVAs (± 10Kb)") + ylab(NULL) + theme_bw() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        legend.position = "top", panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "plain", size = 12, colour = "black", angle = 0),
        strip.placement = "outside") +
  facet_wrap(. ~ Data, scales = 'free', nrow = 2, strip.position = "left",
             labeller = as_labeller(c(`SNP density` = "SNPs density", `Recombination rate` = "Recombination rate")))
dev.off()
pdf(file.path(outdir, "SNP_density_and_Recombination_rate_in_positive_and_negative_sites_around_REs_K9_SVAs_5kb.pdf"), width = 4.5, height = 3.5)
pd %>% 
  filter(bins %in% 11:31) %>% 
  ggplot(aes(x = bins, y = ratio)) +
  stat_smooth(method = "loess", span = 0.05, se = T, alpha = 0.3, aes(color = Group, fill = Group), lwd = 1.25, level = 0.25) +
  scale_fill_manual(values = c("Positive site" = "#df233d", "Negative site" = "#3093dc")) + 
  scale_color_manual(values = c("Positive site" = "#df233d", "Negative site" = "#3093dc")) + 
  geom_vline(xintercept = 21, colour = "#000000", linetype = "longdash", size = 0.75, alpha = 0.75) +
  labs(x = "SVAs (± 5Kb)") + ylab(NULL) + theme_bw() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        legend.position = "top", panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "plain", size = 12, colour = "black", angle = 0),
        strip.placement = "outside") +
  facet_wrap(. ~ Data, scales = 'free', nrow = 2, strip.position = "left",
             labeller = as_labeller(c(`SNP density` = "SNPs density", `Recombination rate` = "Recombination rate")))
dev.off()
# - plotting for REs H3K27me3
pd <- pe.re.svms$REs_H3K27me3 %>% 
  mutate(pos.id = paste(V1, ":", V2, "-", V3, sep = "")) %>% 
  group_by(pos.id) %>% 
  summarise(pos.ratio = sum(V14 > 0)/length(V14)) %>% 
  as.data.frame() %>% 
  mutate(group = case_when(pos.ratio >= 0.35 ~ "Positive sites",
                           pos.ratio <= 0.3 ~ "Negative sites")) %>% 
  na.omit()
pd.pos <- resk27.sva.snpden$region$region.id[resk27.sva.snpden$region$pos.id %in% subset(pd, group == "Positive sites")$pos.id]
pd.neg <- resk27.sva.snpden$region$region.id[resk27.sva.snpden$region$pos.id %in% subset(pd, group == "Negative sites")$pos.id]
resk27.sva.snpden$value %>% 
  rownames_to_column(var = "pos.id") %>% 
  mutate(Group = case_when(pos.id %in% pd.pos & bins_21 <= 0.02 ~ "Positive site",
                           pos.id %in% pd.neg ~ "Negative site")) %>% 
  na.omit() %>% 
  gather(key = "bins", value = "ratio", -Group, -pos.id) %>% 
  group_by(Group) %>% arrange(Group) %>% 
  mutate(bins = gsub("bins_", "", bins)) %>% 
  mutate(bins = as.numeric(bins)) -> pd.snp
table(pd.snp$Group, pd.snp$bins)
resk27.sva.rate$value %>% 
  rownames_to_column(var = "pos.id") %>% 
  mutate(Group = case_when(pos.id %in% pd.pos ~ "Positive site",
                           pos.id %in% pd.neg ~ "Negative site")) %>% 
  na.omit() %>% 
  gather(key = "bins", value = "ratio", -Group, -pos.id) %>% 
  group_by(Group) %>% arrange(Group) %>% 
  mutate(bins = gsub("bins_", "", bins)) %>% 
  mutate(bins = as.numeric(bins)) -> pd.rate
pd.snp$Data <- "SNP density"
pd.rate$Data <- "Recombination rate"
rbind(pd.snp, pd.rate) -> pd
pdf(file.path(outdir, "SNP_density_and_Recombination_rate_in_positive_and_negative_sites_around_REs_K27_SVAs_10kb.pdf"), width = 4.5, height = 3.5)
pd %>% 
  ggplot(aes(x = bins, y = ratio)) +
  stat_smooth(method = "loess", span = 0.05, se = T, alpha = 0.3, aes(color = Group, fill = Group), lwd = 1.25, level = 0.25) +
  scale_fill_manual(values = c("Positive site" = "#df233d", "Negative site" = "#3093dc")) + 
  scale_color_manual(values = c("Positive site" = "#df233d", "Negative site" = "#3093dc")) + 
  geom_vline(xintercept = 21, colour = "#000000", linetype = "longdash", size = 0.75, alpha = 0.75) +
  labs(x = "SVAs (± 10Kb)") + ylab(NULL) + theme_bw() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        legend.position = "top", panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "plain", size = 12, colour = "black", angle = 0),
        strip.placement = "outside") +
  facet_wrap(. ~ Data, scales = 'free', nrow = 2, strip.position = "left",
             labeller = as_labeller(c(`SNP density` = "SNPs density", `Recombination rate` = "Recombination rate")))
dev.off()
pdf(file.path(outdir, "SNP_density_and_Recombination_rate_in_positive_and_negative_sites_around_REs_K27_SVAs_5kb.pdf"), width = 4.5, height = 3.5)
pd %>% 
  filter(bins %in% 11:31) %>% 
  ggplot(aes(x = bins, y = ratio)) +
  stat_smooth(method = "loess", span = 0.05, se = T, alpha = 0.3, aes(color = Group, fill = Group), lwd = 1.25, level = 0.25) +
  scale_fill_manual(values = c("Positive site" = "#df233d", "Negative site" = "#3093dc")) + 
  scale_color_manual(values = c("Positive site" = "#df233d", "Negative site" = "#3093dc")) + 
  geom_vline(xintercept = 21, colour = "#000000", linetype = "longdash", size = 0.75, alpha = 0.75) +
  labs(x = "SVAs (± 5Kb)") + ylab(NULL) + theme_bw() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        legend.position = "top", panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "plain", size = 12, colour = "black", angle = 0),
        strip.placement = "outside") +
  facet_wrap(. ~ Data, scales = 'free', nrow = 2, strip.position = "left",
             labeller = as_labeller(c(`SNP density` = "SNPs density", `Recombination rate` = "Recombination rate")))
dev.off()
# - plotting for REs H3K9me3
pd <- pe.re.svms$PEs %>% 
  mutate(pos.id = paste(V1, ":", V2, "-", V3, sep = "")) %>% 
  group_by(pos.id) %>% 
  summarise(pos.ratio = sum(V14 > 0)/length(V14)) %>% 
  as.data.frame() %>% 
  mutate(group = case_when(pos.ratio >= 0.35 ~ "Positive sites",
                           pos.ratio <= 0.3 ~ "Negative sites")) %>% 
  na.omit()
pd.pos <- pes.sva.snpden$region$region.id[pes.sva.snpden$region$pos.id %in% subset(pd, group == "Positive sites")$pos.id]
pd.neg <- pes.sva.snpden$region$region.id[pes.sva.snpden$region$pos.id %in% subset(pd, group == "Negative sites")$pos.id]
pes.sva.snpden$value %>% 
  rownames_to_column(var = "pos.id") %>% 
  mutate(Group = case_when(pos.id %in% pd.pos & bins_21 <= 0.02 ~ "Positive site",
                           pos.id %in% pd.neg ~ "Negative site")) %>% 
  na.omit() %>% 
  gather(key = "bins", value = "ratio", -Group, -pos.id) %>% 
  group_by(Group) %>% arrange(Group) %>% 
  mutate(bins = gsub("bins_", "", bins)) %>% 
  mutate(bins = as.numeric(bins)) -> pd.snp
table(pd.snp$Group, pd.snp$bins)
pes.sva.rate$value %>% 
  rownames_to_column(var = "pos.id") %>% 
  mutate(Group = case_when(pos.id %in% pd.pos ~ "Positive site",
                           pos.id %in% pd.neg ~ "Negative site")) %>% 
  na.omit() %>% 
  gather(key = "bins", value = "ratio", -Group, -pos.id) %>% 
  group_by(Group) %>% arrange(Group) %>% 
  mutate(bins = gsub("bins_", "", bins)) %>% 
  mutate(bins = as.numeric(bins)) -> pd.rate
pd.snp$Data <- "SNP density"
pd.rate$Data <- "Recombination rate"
rbind(pd.snp, pd.rate) -> pd
pdf(file.path(outdir, "SNP_density_and_Recombination_rate_in_positive_and_negative_sites_around_PEs_SVAs_10kb.pdf"), width = 4.5, height = 3.5)
pd %>% 
  ggplot(aes(x = bins, y = ratio)) +
  stat_smooth(method = "loess", span = 0.05, se = T, alpha = 0.3, aes(color = Group, fill = Group), lwd = 1.25, level = 0.25) +
  scale_fill_manual(values = c("Positive site" = "#df233d", "Negative site" = "#3093dc")) + 
  scale_color_manual(values = c("Positive site" = "#df233d", "Negative site" = "#3093dc")) + 
  geom_vline(xintercept = 21, colour = "#000000", linetype = "longdash", size = 0.75, alpha = 0.75) +
  labs(x = "SVAs (± 10Kb)") + ylab(NULL) + theme_bw() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        legend.position = "top", panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "plain", size = 12, colour = "black", angle = 0),
        strip.placement = "outside") +
  facet_wrap(. ~ Data, scales = 'free', nrow = 2, strip.position = "left",
             labeller = as_labeller(c(`SNP density` = "SNPs density", `Recombination rate` = "Recombination rate")))
dev.off()
pdf(file.path(outdir, "SNP_density_and_Recombination_rate_in_positive_and_negative_sites_around_PEs_SVAs_5kb.pdf"), width = 4.5, height = 3.5)
pd %>% 
  filter(bins %in% 11:31) %>% 
  ggplot(aes(x = bins, y = ratio)) +
  stat_smooth(method = "loess", span = 0.05, se = T, alpha = 0.3, aes(color = Group, fill = Group), lwd = 1.25, level = 0.25) +
  scale_fill_manual(values = c("Positive site" = "#df233d", "Negative site" = "#3093dc")) + 
  scale_color_manual(values = c("Positive site" = "#df233d", "Negative site" = "#3093dc")) + 
  geom_vline(xintercept = 21, colour = "#000000", linetype = "longdash", size = 0.75, alpha = 0.75) +
  labs(x = "SVAs (± 5Kb)") + ylab(NULL) + theme_bw() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        legend.position = "top", panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "plain", size = 12, colour = "black", angle = 0),
        strip.placement = "outside") +
  facet_wrap(. ~ Data, scales = 'free', nrow = 2, strip.position = "left",
             labeller = as_labeller(c(`SNP density` = "SNPs density", `Recombination rate` = "Recombination rate")))
dev.off()


### >>> 6. Variance of gene expression level
# - output regions for GREAT (REs H3K9me3 SVA positive and negative)
pd <- pe.re.svms$REs_H3K9me3 %>% 
  mutate(pos.id = paste(V1, ":", V2, "-", V3, sep = "")) %>% 
  group_by(pos.id) %>% 
  summarise(pos.ratio = sum(V14 > 0)/length(V14)) %>% 
  as.data.frame() %>% 
  mutate(group = case_when(pos.ratio >= 0.50 ~ "Positive sites",
                           pos.ratio <= 0.25 ~ "Negative sites")) %>% 
  na.omit() %>% 
  separate(pos.id, c("chr", "start", "end"))
write.table(subset(pd, group == "Positive sites")[,1:3], file.path(outdir, "REs_H3K9me3_SVA_Positive_sites.bed"), 
            sep = "\t", col.names = F, row.names = F, quote = F)
write.table(subset(pd, group == "Negative sites")[,1:3], file.path(outdir, "REs_H3K9me3_SVA_Negative_sites.bed"), 
            sep = "\t", col.names = F, row.names = F, quote = F)
# - output regions for GREAT (REs H3K27me3 SVA positive and negative)
pd <- pe.re.svms$REs_H3K27me3 %>% 
  mutate(pos.id = paste(V1, ":", V2, "-", V3, sep = "")) %>% 
  group_by(pos.id) %>% 
  summarise(pos.ratio = sum(V14 > 0)/length(V14)) %>% 
  as.data.frame() %>% 
  mutate(group = case_when(pos.ratio >= 0.50 ~ "Positive sites",
                           pos.ratio <= 0.25 ~ "Negative sites")) %>% 
  na.omit() %>% 
  separate(pos.id, c("chr", "start", "end"))
write.table(subset(pd, group == "Positive sites")[,1:3], file.path(outdir, "REs_H3K27me3_SVA_Positive_sites.bed"), 
            sep = "\t", col.names = F, row.names = F, quote = F)
write.table(subset(pd, group == "Negative sites")[,1:3], file.path(outdir, "REs_H3K27me3_SVA_Negative_sites.bed"), 
            sep = "\t", col.names = F, row.names = F, quote = F)
# - output regions for GREAT (PEs SVA positive and negative)
pd <- pe.re.svms$PEs %>% 
  mutate(pos.id = paste(V1, ":", V2, "-", V3, sep = "")) %>% 
  group_by(pos.id) %>% 
  summarise(pos.ratio = sum(V14 > 0)/length(V14)) %>% 
  as.data.frame() %>% 
  mutate(group = case_when(pos.ratio >= 0.50 ~ "Positive sites",
                           pos.ratio <= 0.25 ~ "Negative sites")) %>% 
  na.omit() %>% 
  separate(pos.id, c("chr", "start", "end"))
write.table(subset(pd, group == "Positive sites")[,1:3], file.path(outdir, "PEs_SVA_Positive_sites.bed"), 
            sep = "\t", col.names = F, row.names = F, quote = F)
write.table(subset(pd, group == "Negative sites")[,1:3], file.path(outdir, "PEs_SVA_Negative_sites.bed"), 
            sep = "\t", col.names = F, row.names = F, quote = F)
# - load GREAT results
slec.neargene <- list()
for (files in list.files("results/R/Graphs/revised_Fig2/", "*table.txt", full.names = T)) {
  index <- grep(files, list.files("results/R/Graphs/revised_Fig2/", "*table.txt", full.names = T))
  slec.neargene[[index]] <- read.table(files, sep = "\t")
  slec.neargene[[index]]$gene.id <- DbiIDtrans(slec.neargene[[index]]$V1, "SYMBOL", "ENSEMBL", "human")
  slec.neargene[[index]] <- na.omit(slec.neargene[[index]])
  names(slec.neargene)[index] <- gsub("_GREAT_region-gene_table.txt", "", basename(files))
}; rm(files, index)
# - load GSE98150 data
gse98150.ge.count <- read.table("results/R/RawData/featurecount/GSE98150_all_samples_gene_count_id_matrix.txt", header = T, row.names = 1)
colnames(gse98150.ge.count)[-1:-5] <- c(paste("MII.oocyte.", 1:2, sep = ""),
                                        paste("2Cell.", 1:4, sep = ""),
                                        paste("4Cell.", 1:4, sep = ""),
                                        paste("8Cell.", 1:3, sep = ""),
                                        paste("Morula.", 1:2, sep = ""),
                                        paste("ICM.", 1:4, sep = ""),
                                        paste("TE.", 1:4, sep = ""),
                                        paste("E6.5.Epi.", 1:2, sep = ""),
                                        paste("E6.5.Exe.", 1:2, sep = ""))
gse98150.ge.tpm <- CountToTpm(gse98150.ge.count[,-1:-5], gse98150.ge.count$Length)
# - preprocess data
hm.homo.gene <- read.table("results/R/Tables/evolution/human_mouse_homologous_one2one_genes_ensembl.txt", sep = "\t")[,c(1,9)]
hm.homo.gene <- hm.homo.gene[!duplicated(hm.homo.gene$V1),]
p1 <- gse98150.ge.tpm; colnames(p1) <- paste("mm", colnames(p1), sep = "."); p1$gene.id <- rownames(p1)
p2 <- gse101571.ge.tpm; colnames(p2) <- paste("hs", colnames(p2), sep = "."); p2$gene.id <- DbiIDtrans(rownames(p2), "SYMBOL", "ENSEMBL", "human")
hm.homo.gene.tpm <- merge(hm.homo.gene, p2, by.x = "V1", by.y = "gene.id")
hm.homo.gene.tpm <- merge(hm.homo.gene.tpm, p1, by.x = "V9", by.y = "gene.id")
hm.homo.gene.tpm <- hm.homo.gene.tpm[!duplicated(hm.homo.gene.tpm$V1),]
rownames(hm.homo.gene.tpm) <- hm.homo.gene.tpm$V1
hm.homo.gene.tpm <- hm.homo.gene.tpm[, c(12:13, 18:21)]; rm(p1, p2)
hm.homo.gene.tpm.quan <- as.matrix(hm.homo.gene.tpm) %>% normalize.quantiles() %>% as.data.frame(row.names = rownames(hm.homo.gene.tpm))
colnames(hm.homo.gene.tpm.quan) <- colnames(hm.homo.gene.tpm)
# - plotting REs H3K9me3
pd <- hm.homo.gene.tpm.quan %>% 
  rownames_to_column(var = "gene.id") %>% 
  mutate(hs.8c.mean = rowMeans(hm.homo.gene.tpm.quan[,1:2]),
         mm.2c.mean = rowMeans(hm.homo.gene.tpm.quan[,3:6])) %>% 
  mutate(group = case_when(gene.id %in% slec.neargene$REs_H3K9me3_SVA_Positive_sites$gene.id & hs.8c.mean >= mm.2c.mean ~ "Positive site",
                           gene.id %in% slec.neargene$REs_H3K9me3_SVA_Negative_sites$gene.id & hs.8c.mean <= mm.2c.mean ~ "Negative site")) %>% na.omit()
pdf(file.path(outdir, "Evo_selection_expression_level_of_genes_near_REs_H3K9me3_SVA_positive_and_negative_sites.pdf"), width = 3, height = 3)
pd.title <- paste("N.positive = ", table(pd[,-1:-7]$group)[2], "; N.negative = ", table(pd[,-1:-7]$group)[1], sep = "")
do.call("rbind", replicate(1, pd[,-1:-7], simplify = FALSE)) %>% 
  gather(key = "celltype", value = "TPM", -group) %>% 
  mutate(celltype = gsub("hs.8c.mean", "hs 8C", celltype)) %>% 
  mutate(celltype = gsub("mm.2c.mean", "mm 2C", celltype)) %>%
  ggplot(aes(x = celltype, y = log2(TPM+0.1))) +
  geom_boxplot(aes(fill = celltype), width = 0.6) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  facet_grid(. ~ group) +
  labs(x = "", y = "log2(TPM+0.1)") + theme_bw() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        legend.position = "none", panel.grid.minor = element_blank(),
        strip.background = element_rect(colour = "black", fill = "white", size = 1, linetype = "solid")) +
  stat_compare_means(method = "t.test", na.rm = T, label.y = 10, 
                     comparisons = list(c("hs 8C", "mm 2C")), size = 4) +
  ggtitle(pd.title)
dev.off()
# - plotting REs H3K27me3
pd <- hm.homo.gene.tpm.quan %>% 
  rownames_to_column(var = "gene.id") %>% 
  mutate(hs.8c.mean = rowMeans(hm.homo.gene.tpm.quan[,1:2]),
         mm.2c.mean = rowMeans(hm.homo.gene.tpm.quan[,3:6])) %>% 
  mutate(group = case_when(gene.id %in% slec.neargene$REs_H3K27me3_SVA_Positive_sites$gene.id & hs.8c.mean >= mm.2c.mean ~ "Positive site",
                           gene.id %in% slec.neargene$REs_H3K27me3_SVA_Negative_sites$gene.id & hs.8c.mean <= mm.2c.mean ~ "Negative site")) %>% na.omit()
pdf(file.path(outdir, "Evo_selection_expression_level_of_genes_near_REs_H3K27me3_SVA_positive_and_negative_sites.pdf"), width = 3, height = 3)
pd.title <- paste("N.positive = ", table(pd[,-1:-7]$group)[2], "; N.negative = ", table(pd[,-1:-7]$group)[1], sep = "")
do.call("rbind", replicate(1, pd[,-1:-7], simplify = FALSE)) %>% 
  gather(key = "celltype", value = "TPM", -group) %>% 
  mutate(celltype = gsub("hs.8c.mean", "hs 8C", celltype)) %>% 
  mutate(celltype = gsub("mm.2c.mean", "mm 2C", celltype)) %>%
  ggplot(aes(x = celltype, y = log2(TPM+0.1))) +
  geom_boxplot(aes(fill = celltype), width = 0.6) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  facet_grid(. ~ group) +
  labs(x = "", y = "log2(TPM+0.1)") + theme_bw() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        legend.position = "none", panel.grid.minor = element_blank(),
        strip.background = element_rect(colour = "black", fill = "white", size = 1, linetype = "solid")) +
  stat_compare_means(method = "t.test", na.rm = T, label.y = 10, 
                     comparisons = list(c("hs 8C", "mm 2C")), size = 4) +
  ggtitle(pd.title)
dev.off()
# - plotting PEs
pd <- hm.homo.gene.tpm.quan %>% 
  rownames_to_column(var = "gene.id") %>% 
  mutate(hs.8c.mean = rowMeans(hm.homo.gene.tpm.quan[,1:2]),
         mm.2c.mean = rowMeans(hm.homo.gene.tpm.quan[,3:6])) %>% 
  mutate(group = case_when(gene.id %in% slec.neargene$PEs_SVA_Positive_sites$gene.id ~ "Positive site",
                           gene.id %in% slec.neargene$PEs_SVA_Negative_sites$gene.id ~ "Negative site")) %>% na.omit()
pdf(file.path(outdir, "Evo_selection_expression_level_of_genes_near_PEs_SVA_positive_and_negative_sites.pdf"), width = 3, height = 3)
pd.title <- paste("N.positive = ", table(pd[,-1:-7]$group)[2], "; N.negative = ", table(pd[,-1:-7]$group)[1], sep = "")
do.call("rbind", replicate(1, pd[,-1:-7], simplify = FALSE)) %>% 
  gather(key = "celltype", value = "TPM", -group) %>% 
  mutate(celltype = gsub("hs.8c.mean", "hs 8C", celltype)) %>% 
  mutate(celltype = gsub("mm.2c.mean", "mm 2C", celltype)) %>%
  ggplot(aes(x = celltype, y = log2(TPM+0.1))) +
  geom_boxplot(aes(fill = celltype), width = 0.6) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  facet_grid(. ~ group) +
  labs(x = "", y = "log2(TPM+0.1)") + theme_bw() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 10, angle = 0),
        legend.position = "none", panel.grid.minor = element_blank(),
        strip.background = element_rect(colour = "black", fill = "white", size = 1, linetype = "solid")) +
  stat_compare_means(method = "t.test", na.rm = T, label.y = 10, 
                     comparisons = list(c("hs 8C", "mm 2C")), size = 4) +
  ggtitle(pd.title)
dev.off()



### ===========================
### 9th Part: Hic data analysis
### ===========================

### >>> 1. SVA contacts with major ZGA gene in 2 and 8 cell
# load data
sva.contacts <- list()
for (i in 1:length(list.files("/home/yhw/bioinfo/project-mine/Embryo.93/hic/hs_embryo/results/bedtools/contacts_with_zga", ".txt", full.names = T))) {
  sva.contacts[[i]] <- read.table(list.files("/home/yhw/bioinfo/project-mine/Embryo.93/hic/hs_embryo/results/bedtools/contacts_with_zga", ".txt", full.names = T)[[i]], 
                                  sep = "\t") %>% arrange(V7)
}; rm(i)

# plotting
p1 <- rbind(sva.contacts[[5]], 
            sva.contacts[[9]]) %>%
  separate(V8, c("cell", "enhancer"), ":") %>%
  mutate(len = V2 - V5, V7 = log2(V7),
         cell = gsub("_rep.", "", cell)) %>%
  filter(len >= 400000 & V1 == V4) %>%
  unique() %>%
  ggplot(aes(x = cell, y = V7)) +
  geom_boxplot(aes(fill = cell)) +
  stat_compare_means(method = "t.test", na.rm = T, label.y = 6, comparisons = list(c("8cell", "2cell")), size = 4) +
  theme_bw() +
  scale_fill_brewer() +
  labs(y = NULL) +
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white", size = 1, linetype = "solid"),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 11, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 11, angle = 45, hjust = 1),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 11, angle = 0)) +
  facet_wrap(. ~ enhancer)
p2 <- rbind(sva.contacts[[7]], 
            sva.contacts[[15]][1:nrow(sva.contacts[[7]])+300,]) %>%
  separate(V8, c("cell", "enhancer"), ":") %>%
  mutate(len = V2 - V5, V7 = log2(V7),
         cell = gsub("_rep.", "", cell)) %>%
  filter(len >= 400000 & V1 == V4) %>%
  unique() %>%
  ggplot(aes(x = cell, y = V7)) +
  geom_boxplot(aes(fill = cell)) +
  stat_compare_means(method = "t.test", na.rm = T, label.y = 6, comparisons = list(c("8cell", "2cell")), size = 4) +
  theme_bw() +
  scale_fill_brewer() +
  labs(y = "log2(Normalized count)") +
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white", size = 1, linetype = "solid"),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 11, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 11, angle = 45, hjust = 1),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 11, angle = 0)) +
  facet_wrap(. ~ enhancer)
p3 <- rbind(sva.contacts[[8]], 
            sva.contacts[[12]]) %>%
  separate(V8, c("cell", "enhancer"), ":") %>%
  mutate(len = V2 - V5, V7 = log2(V7),
         cell = gsub("_rep.", "", cell)) %>%
  filter(len >= 400000 & V1 == V4) %>%
  unique() %>%
  ggplot(aes(x = cell, y = V7)) +
  geom_boxplot(aes(fill = cell)) +
  stat_compare_means(method = "t.test", na.rm = T, label.y = 6, comparisons = list(c("8cell", "2cell")), size = 4) +
  theme_bw() +
  scale_fill_brewer() +
  labs(y = NULL) +
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white", size = 1, linetype = "solid"),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 11, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 11, angle = 45, hjust = 1),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 11, angle = 0)) +
  facet_wrap(. ~ enhancer)
pdf(file.path(outdir, "Contact_frequency_of_enhancers_with_ZGA_genes_in_2cell_and_8cell.pdf"), width = 5, height = 3)
plot_grid(p2, p1, p3, nrow = 1)
dev.off()



### =========================================
### 9th Part: 
### =========================================

### >>> . Load GSE101571 data
# - gene
gse101571.ge.count <- read.table("/home/data/publicdata/GSE101571/rnaseq/analysis/results/featurecounts/top10_gene/all_samples_gene_count_name_matrix.txt",
                                 sep = "\t", row.names = 1, header = T)
colnames(gse101571.ge.count)[-1:-5] <- c("GV.oocyte.r1", "GV.oocyte.r2", "MII.oocyte.r1", "MII.oocyte.r2", "2Cell.r1", "2Cell.r2", "2Cell.r3",
                                         "4Cell.r1", "4Cell.r2", "8Cell.r1", "8Cell.r2", "ICM.r1", "ICM.r2")
gse101571.ge.tpm <- CountToTpm(gse101571.ge.count[,-1:-5], gse101571.ge.count$Length)
# - repeat
gse101571.re.count <- read.table("/home/data/publicdata/GSE101571/rnaseq/analysis/results/featurecounts/top10_repeat/all_samples_repeat_count_matrix.txt",
                                 sep = "\t", row.names = 1, header = T)
colnames(gse101571.re.count)[-1:-5] <- c("GV.oocyte.r1", "GV.oocyte.r2", "MII.oocyte.r1", "MII.oocyte.r2", "2Cell.r1", "2Cell.r2", "2Cell.r3",
                                         "4Cell.r1", "4Cell.r2", "8Cell.r1", "8Cell.r2", "ICM.r1", "ICM.r2")
gse101571.re.tpm <- CountToTpm(gse101571.re.count[,-1:-5], gse101571.re.count$Length)
# - 
all(colnames(gse101571.ge.count)[-1:-5]==colnames(gse101571.re.count)[-1:-5])


### >>> . Load regions
res.k9.sva <- read.table("")
res.sva <- read.table("/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer/8cell_reprogrammed_enhancer_derived_by_SVA.bed", sep = "\t")
boxplot(log1p(gse101571.re.tpm[res.sva$V7,]))


### >>> Load geneset
# - ZGA gene list
hs.major.zga <- read.table("/home/yhw/bioinfo/project-mine/Embryo.93/R/Tables/genelist/hs_zga_gene_pc_lncrna_ensembl.bed") %>% 
  separate(V4, c("EnsemblID", "Symbol"), sep = ":")
# - DEGs in "Human embryonic genome activation initiates at the one-cell stage"
pre8c.zga <- read.csv("/home/yhw/bioinfo/project-mine/Embryo.93/R/RawData/zygote_ZGA/1-s2.0-S1934590921004847-mmc2.csv", header = T) %>% 
  filter(logFC > 0 & Gene.Type == "protein_coding")
length(intersect(pre8c.zga$Gene.Symbol, hs.major.zga$Symbol))
# - gene list targeted by SVA-derived enhancers
hs.resk9.gene <- read.table("/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer/genes_interacted_with_sva_c8_RE_k93_enhancers_in_8cell_full.txt", 
                            sep = "\t") %>% 
  separate(col = V16, into = c("Geneid", "Symbol"), sep = ":")
hs.resk27.gene <- read.table("/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer/genes_interacted_with_sva_c8_RE_k273_enhancers_in_8cell_full.txt", 
                             sep = "\t") %>% 
  separate(col = V16, into = c("Geneid", "Symbol"), sep = ":")
hs.pes.gene <- read.table("/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer/genes_interacted_with_sva_c8_PE_enhancers_in_8cell_full.txt", 
                          sep = "\t") %>% 
  separate(col = V16, into = c("Geneid", "Symbol"), sep = ":")
length(unique(intersect(pre8c.zga$Gene.Symbol, intersect(hs.resk9.gene$Symbol, hs.major.zga$Symbol))))
length(unique(intersect(pre8c.zga$Gene.Symbol, intersect(hs.resk27.gene$Symbol, hs.major.zga$Symbol))))
length(unique(intersect(pre8c.zga$Gene.Symbol, intersect(hs.pes.gene$Symbol, hs.major.zga$Symbol))))

