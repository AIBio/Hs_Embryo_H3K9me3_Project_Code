################################################# ~ Good Luck ~ #########################################################
###### >>>>>>                       TITLE: H3K9me3 Project Revised Figure5 analysis                         <<<<<< ######

# - Figure 5 content:
# - 1st part: Library
# - 2nd part: Transcriptional heterogeneity in morula (YHW)
# - 3th part: Open chromatin in morula (YHW)


### =================
### 1st part: Library
### =================
### >>> 1. Library
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("BSgenome.Hsapiens.UCSC.hg38")
library("SingleCellExperiment")
library("Seurat")
library("ggplot2")
library("ChIPseeker")
library("cowplot")
library("tidyr")
library("tibble")
library("tidyverse")
library("stringr")
library("circlize")
library("ComplexHeatmap")
library("amap")
library("karyoploteR")
library("preprocessCore")
library("dendextend")
library("moonBook")
library("webr")
library("ggpubr")
library("ggrepel")
library("VennDiagram")
library("RColorBrewer")
library("FactoMineR")
library("ade4")
library("factoextra")



### =======================================================
### 2nd part: Transcriptional heterogeneity in morula (YHW)
### =======================================================

### >>> 1. Load data
# count matrix
rna.count <- read.table("RawData/scATAC/41467_2018_8205_MOESM9_ESM_supplement_data_6_modified.txt", header = T, sep = "\t")
colnames(rna.count) <- gsub("scCAT_", "", colnames(rna.count))
# gene length
hs.gene.anno <- read.table("RawData/scATAC/Homo_sapiens.GRCh38.97_all_gene_body_embl.bed", sep = "\t") %>%
  mutate(Length = V3-V2)
colnames(hs.gene.anno) <- c("chr", "start", "end", "strand", "ENSEMBL", "SYMBOL", "Type", "Length")
# add annotation information to count matrix
rna.count <- merge(rna.count, hs.gene.anno, by = "ENSEMBL")
rownames(rna.count) <- rna.count$ENSEMBL
rna.count <- rna.count[, c(1:73,80)]
# ICM vs TE DEGs
icm.high <- read.table("RawData/scATAC/hs_te_vs_icm_down_gene_2016.txt", sep = "\t")
te.high <- read.table("RawData/scATAC/hs_te_vs_icm_up_gene_2016.txt", sep = "\t")
# H3K9me3 modifier
c("KDM4A", "KDM4B", "KDM4C", "KDM4D", "KDM4E", "KDM4F", "SETDB1", "SETDB2", "SUV39H1", "SUV39H2", "UHRF1")
hs.k9.factor <- c("ENSG00000066135", "ENSG00000127663", "ENSG00000107077", "ENSG00000186280",
                  "ENSG00000235268", "ENSG00000255855", "ENSG00000143379", "ENSG00000136169",
                  "ENSG00000101945", "ENSG00000152455", "ENSG00000276043")
# GSE36552
tfc.ge.count <- readRDS("RawData/scATAC/gse36552_ge_count.rds")
tfc.ge.count <- tfc.ge.count[, 47:92]
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 11132692 16782216 20631442 21577826 25261317 38518674 
tfc.ge.tpm <- readRDS("RawData/scATAC/gse36552_ge_tpm.rds")
tfc.ge.tpm <- tfc.ge.tpm[, 47:92]
tfc.ge.tpm <- tfc.ge.tpm[rowSums(tfc.ge.tpm) > 5,]
boxplot(log1p(tfc.ge.tpm))
tfc.ge.tpm.quan <- as.matrix(tfc.ge.tpm) %>% normalize.quantiles() %>% as.data.frame(row.names = rownames(tfc.ge.tpm))
colnames(tfc.ge.tpm.quan) <- colnames(tfc.ge.tpm)
boxplot(log1p(tfc.ge.tpm.quan))
tfc.morula.ge.tpm.quan.z <- as.data.frame(t(apply(tfc.ge.tpm.quan[,grep("Morula", colnames(tfc.ge.tpm.quan))], 1, function(x){(x - mean(x))/(sd(x))})), 
                                          row.names = rownames(tfc.ge.tpm.quan))
# 2016 Cell
cell.ge.count <- readRDS("RawData/scATAC/2016Cell_gene_count_e345.rds")
cell.ge.count <- cell.ge.count[rowSums(cell.ge.count)>20,]
cell.ge.count <- cell.ge.count[,colSums(cell.ge.count)>3000000]
table(rownames(cell.ge.count) %in% hs.gene.anno$ENSEMBL)
cell.ge.count$ENSEMBL <- rownames(cell.ge.count)
cell.ge.count <- merge(cell.ge.count, hs.gene.anno, by = "ENSEMBL")
cell.ge.tpm <- CountToTpm(cell.ge.count[,grep("_", colnames(cell.ge.count))], cell.ge.count$Length)
rownames(cell.ge.tpm) <- cell.ge.count$ENSEMBL
cell.ge.tpm.quan <- as.matrix(cell.ge.tpm) %>% normalize.quantiles() %>% as.data.frame(row.names = rownames(cell.ge.tpm))
colnames(cell.ge.tpm.quan) <- colnames(cell.ge.tpm)
cell2016.meta <- readRDS("RawData/scATAC/2016Cell_meta_e345.rds")
cell2016.meta <- cell2016.meta %>% 
  separate("cell_name", c("Stage", "EmbryoID", "EmbryoIndex"), "\\.") %>%
  mutate(Embryo = paste(Stage, EmbryoID, sep = ".")) %>%
  filter(CellId %in% grep("Morula", colnames(cell.ge.count), value = T))
rownames(cell2016.meta) <- cell2016.meta$CellId
# gene list positively or negatively correlated with GATA3
# Gerri, C. et al. Initiation of a conserved trophectoderm program in human, cow and mouse embryos. Nature, doi:10.1038/s41586-020-2759-x (2020).
gata3.cor.gene <- read.table("RawData/scATAC/41586_2020_2759_MOESM3_ESM_supplement_table2.txt", sep = "\t", header = T)
gata3.cor.gene$ENSEMBL <- DbiIDtrans(gata3.cor.gene$Gene, "SYMBOL", "ENSEMBL", "human")


### >>> 2. Quality control
# sequencing depth
barplot(colSums(rna.count[, c(-1,-ncol(rna.count))]))
# filtering sample
keep <- colnames(rna.count)[c(-1,-ncol(rna.count))][colSums(rna.count[, c(-1,-ncol(rna.count))]) >= 3000000]
rna.count <- rna.count[, c("ENSEMBL", keep, "Length")]
# filtering gene
rna.count <- rna.count[rowSums(rna.count[, c(-1,-ncol(rna.count))]) > 10,]
# create metadata
cell.meta <- data.frame(SampleID = colnames(rna.count)[c(-1,-ncol(rna.count))],
                        Stage = str_split_fixed(colnames(rna.count)[c(-1,-ncol(rna.count))], "_", 2)[,1],
                        row.names = colnames(rna.count)[c(-1,-ncol(rna.count))])


### >>> 3. Normalization
rna.tpm <- CountToTpm(rna.count[, c(-1,-ncol(rna.count))], rna.count$Length)
rna.tpm.quan <- as.matrix(rna.tpm) %>% normalize.quantiles() %>% as.data.frame(row.names = rownames(rna.tpm))
colnames(rna.tpm.quan) <- colnames(rna.tpm)
boxplot(log1p(rna.tpm.quan))
# divide samples
rna.morula.tpm <- rna.tpm[,grep("Morula", colnames(rna.tpm))]
rna.morula.tpm <- rna.morula.tpm[rowSums(rna.morula.tpm) > 0,]
boxplot(log1p(rna.morula.tpm))
rna.morula.tpm.quan <- rna.tpm.quan[,grep("Morula", colnames(rna.tpm.quan))]
boxplot(log1p(rna.morula.tpm.quan))
rna.morula.tpm.quan <- rna.morula.tpm.quan[rowSums(rna.morula.tpm.quan) > 0,]
rna.morula.tpm.quan.z <- as.data.frame(t(apply(rna.morula.tpm.quan, 1, function(x){(x - mean(x))/(sd(x))})), row.names = rownames(rna.morula.tpm.quan))


### >>> 4. Highly variable genes
rownames(HVGfinderPv(rna.morula.tpm.quan, 0.05))
rownames(HVGfinderPv(tfc.ge.tpm.quan[,grep("Morula", colnames(tfc.ge.tpm.quan))], 0.05))
rownames(HVGfinderPv(cell.ge.tpm.quan[,grep("Morula", colnames(cell.ge.tpm.quan))], 0.01))
venn.diagram(x = list(d1 = rownames(HVGfinderPv(rna.morula.tpm.quan, 0.05)), 
                      d2 = rownames(HVGfinderPv(tfc.ge.tpm.quan[,grep("Morula", colnames(tfc.ge.tpm.quan))], 0.05)), 
                      d3 = rownames(HVGfinderPv(cell.ge.tpm.quan[,grep("Morula", colnames(cell.ge.tpm.quan))], 0.01))),
             output = TRUE, imagetype = "png", filename = "Graphs/test.png",
             height = 2000, width = 2000, resolution = 500,
             cat.cex = 1, cat.default.pos = "outer",
             cat.fontfamily = "sans", cat.col = c("#000000", "#000000", "#000000"),
             fill = c(alpha("#ed5113", 1), alpha("#dc2056", 1), alpha("#085fbb", 1)),
             lty = 1, lwd = 2.5, col = c("#ffffff", "#ffffff", "#ffffff"),
             cex = 1, label.col = "#000000", margin = 0.1)


### >>> 5. Seurat pipeline
# single ATAC+RNA
rna.seurat <- CreateSeuratObject(counts = rna.count[,c(-1,-ncol(rna.count))], meta.data = cell.meta)
rna.seurat <- NormalizeData(rna.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
rna.seurat <- FindVariableFeatures(rna.seurat, selection.method = "vst", nfeatures = 1000)
rna.seurat <- ScaleData(rna.seurat, features = rownames(rna.seurat))
rna.seurat <- RunPCA(rna.seurat, features = VariableFeatures(object = rna.seurat))
rna.seurat <- JackStraw(rna.seurat, num.replicate = 100)
rna.seurat <- ScoreJackStraw(rna.seurat, dims = 1:20)
ElbowPlot(rna.seurat)
rna.seurat <- FindNeighbors(rna.seurat, dims = 1:8)
rna.seurat <- FindClusters(rna.seurat, resolution = 1.5)
rna.seurat <- RunUMAP(rna.seurat, dims = 1:8)
DimPlot(rna.seurat, reduction = "umap", label = TRUE, group.by = "ident", label.size = 4, pt.size = 2)
UMAPPlot(rna.seurat, reduction = "umap", label = TRUE, group.by = "Stage", label.size = 4, pt.size = 2)
FeaturePlot(rna.seurat, features = hs.gene.anno[match(c("NANOG", "SOX2", "POU5F1", "CDX2", "GATA2", "GATA3"), hs.gene.anno$SYMBOL),5], 
            cols = c("#c0bda5", "#cd0e0f"), reduction = "umap", label = FALSE, label.size = 4, pt.size = 1)
# 2016 cell
cell2016.seurat <- CreateSeuratObject(counts = cell.ge.count[,cell2016.meta$CellId], meta.data = cell2016.meta)
cell2016.seurat <- NormalizeData(cell2016.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
cell2016.seurat <- FindVariableFeatures(cell2016.seurat, selection.method = "vst", nfeatures = 1000)
cell2016.seurat <- ScaleData(cell2016.seurat, features = rownames(cell2016.seurat))
cell2016.seurat <- RunPCA(cell2016.seurat, features = VariableFeatures(object = cell2016.seurat))
cell2016.seurat <- JackStraw(cell2016.seurat, num.replicate = 100)
cell2016.seurat <- ScoreJackStraw(cell2016.seurat, dims = 1:20)
ElbowPlot(cell2016.seurat)
cell2016.seurat <- FindNeighbors(cell2016.seurat, dims = 1:8)
cell2016.seurat <- FindClusters(cell2016.seurat, resolution = 1.5)
cell2016.seurat <- RunUMAP(cell2016.seurat, dims = 1:8)
DimPlot(cell2016.seurat, reduction = "umap", label = TRUE, group.by = "ident", label.size = 4, pt.size = 2)
UMAPPlot(cell2016.seurat, reduction = "umap", label = TRUE, group.by = "Embryo", label.size = 4, pt.size = 2)
FeaturePlot(cell2016.seurat, features = hs.gene.anno[match(c("NANOG", "SOX2", "POU5F1", "CDX2", "GATA2", "GATA3", "PARD3"), hs.gene.anno$SYMBOL),5], 
            cols = c("#c0bda5", "#cd0e0f"), reduction = "umap", label = FALSE, label.size = 4, pt.size = 1)


### >>> 6. Define outer and inner cell in morula
# scATAC-seq + scRNA-seq (Liu, L. et al.2019)
rna.morula.tpm.quan[hs.gene.anno[match(c("GATA3"), hs.gene.anno$SYMBOL),5],] %>%
  rownames_to_column(var = "ENSEMBL") %>%
  gather(key = "cell", value = "TPM", -ENSEMBL) %>%
  filter(TPM < 0.01) -> rna.morula.inner
rna.morula.tpm.quan[hs.gene.anno[match(c("GATA3"), hs.gene.anno$SYMBOL),5],] %>%
  rownames_to_column(var = "ENSEMBL") %>%
  gather(key = "cell", value = "TPM", -ENSEMBL) %>%
  filter(TPM > 0.8) -> rna.morula.outer

# plotting
pd <- rna.morula.tpm.quan[,c(rna.morula.inner$cell, rna.morula.outer$cell)]
pd <- pd[rowSums(pd)>0,]
pd <- na.omit(pd[na.omit(subset(gata3.cor.gene, Spearman >= 0.25 | Spearman <= -0.25)$ENSEMBL),])
pd$SYMBOL <- DbiIDtrans(rownames(pd), "ENSEMBL", "SYMBOL", "human")
rownames(pd) <- pd$SYMBOL
pd.cor <- corr.test(t(pd[,grep("Morula", colnames(pd))]), method = "spearman", adjust = "fdr")
pd <- pd.cor$r[c("GATA3"),] %>% as.data.frame() %>%
  rownames_to_column(var = "gene") %>% na.omit()
colnames(pd) <- c("gene", "cor")
pdf("Graphs/ICM_and_TE/scATAC/inner_outer_cell_gene_correlation_with_GATA3.pdf", height = 3, width = 5)
pd %>%
  mutate(fill = case_when(cor > 0 ~ "pos", cor < 0 ~ "neg",)) %>%
  filter(abs(cor) > 0.25) %>%
  filter(gene %in% c(intersect(subset(gata3.cor.gene, Spearman < 0)$Gene, subset(pd, cor <= 0)$gene),
                     intersect(subset(gata3.cor.gene, Spearman > 0)$Gene, subset(pd, cor > 0)$gene))) %>%
  ggplot(aes(x = fct_reorder(gene, cor, .desc = T), y = cor)) +
  geom_bar(aes(fill = fill), stat = "identity", color = "#ecf0f1") +
  scale_fill_manual(values = c(pos = "#AB0419", neg = "#0B488A")) +
  theme_bw() + labs(x = "Gene", y = "Correlation (vs GATA3)") +
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 11, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 11, angle = 45, hjust = 1),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 11, angle = 0))
dev.off(); rm(pd, pd.cor)
# Yan, L. et al. 2013
tfc.ge.tpm[hs.gene.anno[match(c("GATA3"), hs.gene.anno$SYMBOL),5], grep("Morula", colnames(tfc.ge.tpm))] %>%
  rownames_to_column(var = "gene") %>%
  gather(key = "cell", value = "TPM", -gene) %>%
  ggplot(aes(x = fct_reorder(cell, TPM), y = TPM)) +
  geom_bar(stat = "identity", fill = "#3498db", width = 1) +
  facet_wrap(. ~ gene) +
  theme_classic()
# Petropoulos, S. et al. 2016
cell.ge.tpm.quan[hs.gene.anno[match(c("GATA3"), hs.gene.anno$SYMBOL),5], grep("Morula", colnames(cell.ge.tpm.quan))] %>%
  rownames_to_column(var = "gene") %>%
  gather(key = "cell", value = "TPM", -gene) %>%
  ggplot(aes(x = fct_reorder(cell, TPM), y = TPM)) +
  geom_bar(stat = "identity", fill = "#3498db", width = 1) +
  facet_wrap(. ~ gene) +
  theme_classic()

### >>> 7. Find differential genes
library("limma")
rna.morula.tpm
group <- factor(epi.morula.filter.meta$SampleGroup)
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
rownames(design) <- colnames(rna.morula.tpm[,epi.morula.filter.meta$SampleID])
contr.matrix <- makeContrasts(contrasts = paste("Inner", "Outer", sep = "-"), levels = design)
fit <- lmFit(log2(rna.morula.tpm[,epi.morula.filter.meta$SampleID]+1), design)
fit <- contrasts.fit(fit, contr.matrix)
fit <- eBayes(fit, trend = TRUE)
rna.deg <- topTable(fit, sort.by = "P", n = Inf); rm(fit, group, design, contr.matrix)
rna.deg <- rna.deg %>% rownames_to_column(var = "ENSEMBL") %>%
  mutate(SYMBOL = DbiIDtrans(ENSEMBL, "ENSEMBL", "SYMBOL", "human"))
rna.deg$log10Pvalue <- -log10(rna.deg$P.Value)
rna.deg <- na.omit(rna.deg)
subset(rna.deg, ENSEMBL %in% hs.k9.factor)


### >>> 8. GO
rna.deg.inner.go <- FullSet.GO("human", na.omit(subset(rna.deg, P.Value <= 0.05 & logFC >= log2(2))$SYMBOL), "morula_inner_up_gene",
                               "SYMBOL", "Tables/GO/Morula_inner_up_gene")
rna.deg.outer.go <- FullSet.GO("human", na.omit(subset(rna.deg, P.Value <= 0.05 & logFC <= -log2(2))$SYMBOL), "morula_outer_up_gene",
                               "SYMBOL", "Tables/GO/Morula_outer_up_gene")
pdf("Graphs/ICM_and_TE/scATAC/inner_outer_cell_scATAC_differential_inner_up_gene_GO.pdf", height = 4, width = 5)
rbind(rna.deg.inner.go[[1]][c("GO:0008380", "GO:0000398", "GO:0031123", "GO:0043487", "GO:0019827", "GO:0098722"),]) %>% 
  mutate(Log10Pvalue = -log10(pvalue)) %>%
  ggplot(aes(x = Log10Pvalue, y = fct_reorder(Description, Log10Pvalue))) + 
  labs(x = "-Log10.pvalue", y = "GO biological process") +
  geom_bar(stat = "identity", fill = "#3498db", color = "#ecf0f1") + 
  theme_classic() +
  geom_vline(xintercept = c(-log10(0.05)), color = "#FF1300", linetype = "longdash", size = 1, alpha = 0.8) + 
  theme(axis.title.x = element_text(face="plain", colour = "#1D1D1E", size = 10, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#1D1D1E", size = 10, angle = 90, hjust = 1), 
        axis.text.x  = element_text(face="plain", colour = "#1D1D1E", size = 10, angle = 0, hjust = 1),
        axis.text.y  = element_text(face="plain", colour = "#1D1D1E", size = 10, angle = 0, hjust = 1))
dev.off()
pdf("Graphs/ICM_and_TE/scATAC/inner_outer_cell_scATAC_differential_outer_up_gene_GO.pdf", height = 4, width = 8)
rbind(rna.deg.outer.go[[1]][c("GO:0019080", "GO:0019083", "GO:0006613", "GO:0000184", "GO:0006402", "GO:0031935"),]) %>% 
  mutate(Log10Pvalue = -log10(pvalue)) %>%
  ggplot(aes(x = Log10Pvalue, y = fct_reorder(Description, Log10Pvalue))) + 
  labs(x = "-Log10.pvalue", y = "GO biological process") +
  geom_bar(stat = "identity", fill = "#3498db", color = "#ecf0f1") + 
  theme_classic() +
  geom_vline(xintercept = c(-log10(0.05)), color = "#FF1300", linetype = "longdash", size = 1, alpha = 0.8) + 
  theme(axis.title.x = element_text(face="plain", colour = "#1D1D1E", size = 10, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#1D1D1E", size = 10, angle = 90, hjust = 1), 
        axis.text.x  = element_text(face="plain", colour = "#1D1D1E", size = 10, angle = 0, hjust = 1),
        axis.text.y  = element_text(face="plain", colour = "#1D1D1E", size = 10, angle = 0, hjust = 1))
dev.off()

### >>> 9. Plotting
#(MAplot or Volcano)
KAP1.protein <- c("ZNF10","ZNF112","ZNF133","ZNF140","ZNF175","ZNF19","ZNF197","ZNF2",
                  "ZNF224","ZNF25","ZNF256","ZNF274","ZNF30","ZNF320","ZNF333","ZNF33A",
                  "ZNF33B","ZNF350","ZNF354A","ZNF382","ZNF432","ZNF496","ZNF573","ZNF585B",
                  "ZNF616","ZNF689","ZNF708","ZNF720","ZNF721","ZNF732","ZNF738","ZNF74")
viral.gene.expresion <- c("RPL7A","RPL4","RPS15A","POLR2H","RPL24","RPS3","RPLP0","NUP37",
                          "TRIM32","RPS12","RPL19","RPSA","RPL22","RPL12","RPL41","RPS4Y1",
                          "AAAS","RPL15","RPS17","RPL14","EIF3D","SEC13","USF1","RPL10")
rna.splicing <- c("RSR2","PPP2CA","EIF4A3","PAPOLA","SRSF11","CWC15","SCNM1","HNRNPU",
                  "YBX1","ZC3H13","TRA2A","GTF2F2","C1QBP","SRSF2","RBM17","TMBIM6",
                  "SNRPA1","BUD31","PTBP1","SRSF9","SON","DNAJC8","YTHDC1","THOC3","LSM10",
                  "SF1","PRPF38A","PRCC","ZRANB2","SNU13","RBM25","SF3B6","SYF2","LUC7L3",
                  "SCAF11","SNRNP25","KHDRBS1","RBM39","TRA2B","SRSF6","SRSF10","HNRNPH1",
                  "IK","FRG1","SNRNP27","LSM1","SRRM1","HNRNPM","RP9","DDX20","CIR1","NUP98",
                  "POLR2F","PQBP1","U2AF1","ALYREF","PPIG","IWS1","HNRNPA1L2","RRAGC","UPF3B","SNRPB2")
pd <- rna.deg %>% 
  mutate(group = case_when((logFC >= log2(2) & P.Value <= 0.05) ~ "Up.regulated", 
                           (logFC <= -log2(2) & P.Value <= 0.05) ~ "Down.regulated",
                           (abs(logFC) < log2(2) | P.Value > 0.05) ~ "Non.changed")) %>% 
  mutate(label = case_when(group == "Up.regulated" & logFC >= log2(6) & P.Value <= 0.000025 ~ "Up.label",
                           group == "Down.regulated" & logFC <= -log2(4) & P.Value <= 0.00025 ~ "down.label"))
pdf("Graphs/ICM_and_TE/scATAC/inner_outer_cell_scATAC_differential_gene_volcano.pdf", height = 4, width = 5)
pd %>%
  ggplot(aes(x = logFC, y = log10Pvalue)) + geom_point(aes(color = group, alpha = group), size = 2) + theme_classic() +
  geom_vline(xintercept = c(-1, 1), colour = "#000000", linetype = "longdash", size = 1, alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), colour = "#000000", linetype = "longdash", size = 1, alpha = 0.5) +
  scale_color_manual(values = c(Up.regulated = "#EA2027", Down.regulated = "#0652DD", Non.changed = "#000000")) +
  scale_alpha_manual(values = c(Up.regulated = 0.5, Down.regulated = 0.5, Non.changed = 0.25)) +
  theme(legend.position = "none") + labs(x = expression(log[2]("Mean FPKM")), y = expression(log[2]("Fold of Change"))) +
  theme(axis.title.x = element_text(face="plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face="plain", colour = "#000000", size = 10, angle = 0),
        axis.text.y  = element_text(face="plain", colour = "#000000", size = 10, angle = 0)) +
  geom_label_repel(aes(label = SYMBOL),
                   data = rbind(na.omit(pd), subset(pd, SYMBOL %in% c("GATA3", "NANOG"))),
                   box.padding = 0.2, point.padding = 0.2, segment.color = "#130f40",
                   label.size = 0.2, fill = "white", size = 3)
dev.off; rm(pd)
# heatmap (TPM)
pdf("Graphs/ICM_and_TE/scATAC/inner_outer_cell_scATAC_differential_peaks_heatmap_TPM.pdf", height = 5, width = 5)
pd <- rbind(rna.morula.tpm[subset(rna.deg, logFC >= log2(2) & P.Value <= 0.05)$ENSEMBL,epi.morula.filter.meta$SampleID],
            rna.morula.tpm[subset(rna.deg, logFC <= -log2(2) & P.Value <= 0.05)$ENSEMBL,epi.morula.filter.meta$SampleID])
Heatmap(pd,
        col = colorRamp2(c(0, 10, 20, 30, 40), c("#481b6d", "#414487", "#1e9b8a", "#72cf56", "#f7e620")), name = "TPM",
        cluster_rows = T, show_row_dend = T, row_dend_side = "left", row_dend_width = unit(1, "cm"), 
        clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
        show_column_names = T, column_names_rot = 45, column_names_side = "bottom",
        cluster_columns = T, show_column_dend = T, column_dend_height = unit(1, "cm"), column_dend_side = "top",
        clustering_distance_columns = "euclidean", clustering_method_columns = "complete",
        show_row_names = F, row_names_rot = 0, row_names_side = "right",
        width = unit(5, "cm"), height = unit(7, "cm"),
        column_split = epi.morula.filter.meta$SampleGroup, column_gap = unit(2, "mm"), border = T,
        row_split = c(rep("Up", nrow(subset(rna.deg, logFC >= log2(2) & P.Value <= 0.05))),
                      rep("Down", nrow(subset(rna.deg, logFC <= -log2(2) & P.Value <= 0.05)))), row_gap = unit(2, "mm"))
dev.off(); rm(pd)

### >>> Overlap with ICMvsTE differential expression genes
dea.icm.vs.te.ge <- list()
for (i in 1:length(list.files("Tables/DEGs/ICM_vs_TE/", pattern = ".txt"))) {
  dea.icm.vs.te.ge[[i]] <- read.table(list.files("Tables/DEGs/ICM_vs_TE/", pattern = ".txt", full.names = T)[i], sep = "\t", header = T)
  #dea.icm.vs.te.ge[[i]] <- dea.icm.vs.te.ge[[i]][,grep("id|ave|var|med", colnames(dea.icm.vs.te.ge[[i]]))]
  colnames(dea.icm.vs.te.ge[[i]]) <- gsub("X", "", colnames(dea.icm.vs.te.ge[[i]]))
  dea.icm.vs.te.ge[[i]]$symbol <- DbiIDtrans(dea.icm.vs.te.ge[[i]]$id, "ENSEMBL", "SYMBOL", "human")
  dea.icm.vs.te.ge[[i]] <- na.omit(dea.icm.vs.te.ge[[i]])
}; rm(i)
names(dea.icm.vs.te.ge) <- c("ICM.spe", "TE.spe")
length(intersect(dea.icm.vs.te.ge$ICM.spe$symbol, na.omit(subset(rna.deg, P.Value <= 0.05 & logFC >= log2(2))$SYMBOL)))
length(intersect(dea.icm.vs.te.ge$TE.spe$symbol, na.omit(subset(rna.deg, P.Value <= 0.05 & logFC <= -log2(2))$SYMBOL)))



### ========================================
### 3th part: Open chromatin in morula (YHW)
### ========================================

### >>> 1. Load data
# count matrix
epi.count <- read.table("RawData/scATAC/41467_2018_8205_MOESM8_ESM_supplement_data_5.txt", header = T, sep = "\t")
colnames(epi.count) <- gsub("scCAT_", "", colnames(epi.count))
rownames(epi.count) <- epi.count$Pos.ID


### >>> 2. Filtering
# sequencing depth
barplot(colSums(epi.count[, colnames(epi.count)%in%c(rna.morula.inner$cell, rna.morula.outer$cell)]))
# filtering sample
keep <- colnames(epi.count)[-1][colSums(epi.count[, -1]) < 200000]
epi.count <- epi.count[, c("Pos.ID", keep)]
# filtering gene
epi.count <- epi.count[rowSums(epi.count[, -1]) > 10,]


### >>> 3. Normalization
epi.len <- as.numeric(str_split_fixed(epi.count$Pos.ID, "_", 3)[,3])-as.numeric(str_split_fixed(epi.count$Pos.ID, "_", 3)[,2])
epi.tpm <- CountToTpm(epi.count[,-1], epi.len)
epi.tpm.quan <- as.matrix(epi.tpm) %>% normalize.quantiles() %>% as.data.frame(row.names = rownames(epi.tpm))
colnames(epi.tpm.quan) <- colnames(epi.tpm)
boxplot(log1p(epi.tpm[,colnames(epi.tpm) %in% c(rna.morula.inner$cell, rna.morula.outer$cell)]))
boxplot(log1p(epi.tpm.quan[,colnames(epi.tpm.quan) %in% c(rna.morula.inner$cell, rna.morula.outer$cell)]))


### >>> 4. Correlation analysis
pd <- epi.tpm.quan[,colnames(epi.tpm.quan) %in% c(rna.morula.inner$cell, rna.morula.outer$cell)]
pd <- pd[,!colnames(pd) %in% c("Morula_520", "Morula_528", "Morula_545")]
pd <- pd[rowSums(pd)>10,]
pdf("Graphs/ICM_and_TE/scATAC/inner_outer_cell_scATAC_clustering.pdf", height = 5, width = 5)
dend <- pd %>% na.omit %>%
  t %>% scale %>% hcluster(method = "pearson", link = "complete") %>% as.dendrogram
dend %>%
  set("branches_k_color", value = c("#E3BF03", "#437ECE"), k = 2) %>% set("branches_lwd", 3) %>%
  set("labels_cex", 1.5) %>% set("labels_col", value = c("#E3BF03", "#437ECE"), k = 2) %>% set("hang_leaves", 0.3) %>%
  set("leaves_cex", 1.5) %>% set("leaves_pch", 19) %>% set("leaves_col", "#E9405D") %>% plot(horiz = TRUE, axes = T) 
dend %>% rect.dendrogram(k = 2, horiz = TRUE, border = 8, lty = 5, lwd = 2)
dev.off(); rm(pd, dend)


### >>> 5. Filter samples and peaks
# inner cells
pd.inner <- epi.tpm.quan[,c("Morula_501", "Morula_502", "Morula_503", "Morula_507")]
pd.inner <- pd.inner[apply(pd.inner, 1, median) > 5,]
nrow(pd.inner)
# outer cells
pd.outer <- epi.tpm.quan[,c("Morula_525", "Morula_526", "Morula_534", "Morula_539")]
pd.outer <- pd.outer[apply(pd.outer, 1, median) > 5,]
nrow(pd.outer)
# merge
epi.morula.filter <- epi.tpm.quan[unique(c(rownames(pd.inner), rownames(pd.outer))), 
                                  c("Morula_501", "Morula_502", "Morula_503", "Morula_507",
                                    "Morula_525", "Morula_526", "Morula_534", "Morula_539")]
write.table(rownames(epi.morula.filter),
            "Tables/scATAC/Morula_all_open_chromatin_peaks.txt", col.names = F, row.names = F, quote = F)
epi.morula.filter.meta <- data.frame(SampleID = colnames(epi.morula.filter),
                                     SampleGroup = c(rep("Inner", 4), rep("Outer", 4)),
                                     row.names = colnames(epi.morula.filter))


### >>> 6. ICM and TE specific gained H3K9me3 peaks
# load peaks
icm.te.spe.k9 <- list()
for (i in seq(1, length(list.files("Tables/scATAC", pattern = "*in_morula.bed")))) {
  icm.te.spe.k9[[i]] <- read.table(list.files("Tables/scATAC", pattern = "*in_morula.bed", full.names = T)[i], sep = "\t")
  icm.te.spe.k9[[i]]$pos.id <- paste(icm.te.spe.k9[[i]]$V5, icm.te.spe.k9[[i]]$V6, icm.te.spe.k9[[i]]$V7, sep = "_")
}; rm(i)
names(icm.te.spe.k9) <- c("ICM.specific", "TE.specific")


boxplot(log2(epi.morula.filter[c(icm.te.spe.k9$ICM.specific$pos.id),]+0.1))
boxplot(log2(epi.morula.filter[c(icm.te.spe.k9$TE.specific$pos.id),]+0.1))

# ATAC signals
pdf("Graphs/ICM_and_TE/specific_pks/filtered/ICM_TE_specific_gained_H3K9me3_peaks_in_morula_ATAC_signal.pdf", height = 3, width = 6)
pd <- epi.morula.filter %>%
  rownames_to_column(var = "pos.id") %>% 
  mutate(Inner = apply(epi.morula.filter[,3:4], 1, mean),
         Outer = apply(epi.morula.filter[,5:7], 1, mean)) %>%
  mutate(Inner = log2(Inner+0.1), Outer = log2(Outer+0.1))
p1 <- pd[,c("pos.id", "Inner", "Outer")] %>% gather(key = "Cell", value = "Log2TPM", -pos.id) %>% 
  mutate(type = case_when(pos.id %in% icm.te.spe.k9$ICM.specific$pos.id ~ "ICM specific",
                          pos.id %in% icm.te.spe.k9$TE.specific$pos.id ~ "TE specific")) %>% na.omit() %>% 
  filter(type == "ICM specific") %>% 
  ggplot(aes(x = Cell, y = Log2TPM)) +
  geom_boxplot(aes(fill = Cell)) +
  scale_fill_brewer() + theme_bw() +
  stat_compare_means(method = "t.test", na.rm = T, label.y = 11, comparisons = list(c("Inner", "Outer")), size = 4) +
  facet_wrap(. ~ type) +
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_rect(colour = "black", fill = "white", size = 1, linetype = "solid"),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 11, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 11, angle = 45, hjust = 1),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 11, angle = 0))
pd <- epi.morula.filter %>%
  rownames_to_column(var = "pos.id") %>% 
  mutate(Inner = apply(epi.morula.filter[,1:2], 1, mean),
         Outer = apply(epi.morula.filter[,c(6,8)], 1, mean)) %>%
  mutate(Inner = log2(Inner+0.1), Outer = log2(Outer+0.1))
p2 <- pd[,c("pos.id", "Inner", "Outer")] %>% gather(key = "Cell", value = "Log2TPM", -pos.id) %>% 
  mutate(type = case_when(pos.id %in% icm.te.spe.k9$ICM.specific$pos.id ~ "ICM specific",
                          pos.id %in% icm.te.spe.k9$TE.specific$pos.id ~ "TE specific")) %>% na.omit() %>% 
  filter(type == "TE specific") %>% 
  ggplot(aes(x = Cell, y = Log2TPM)) +
  geom_boxplot(aes(fill = Cell)) +
  scale_fill_brewer() + theme_bw() +
  stat_compare_means(method = "t.test", na.rm = T, label.y = 11, comparisons = list(c("Inner", "Outer")), size = 4) +
  facet_wrap(. ~ type) +
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_rect(colour = "black", fill = "white", size = 1, linetype = "solid"),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 11, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 11, angle = 45, hjust = 1),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 11, angle = 0))
plot_grid(p1, p2, ncol = 2)
dev.off(); rm(p1, p2, pd)


### >>> 7. Find differential peaks
library("limma")
group <- factor(epi.morula.filter.meta$SampleGroup)
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
rownames(design) <- colnames(epi.morula.filter)
contr.matrix <- makeContrasts(contrasts = paste("Inner", "Outer", sep = "-"), levels = design)
fit <- lmFit(log2(epi.morula.filter+1), design)
fit <- contrasts.fit(fit, contr.matrix)
fit <- eBayes(fit, trend = TRUE)
res <- topTable(fit, sort.by = "P", n = Inf); rm(fit, group, design, contr.matrix)
# output tables (pvalue 0.15: )
write.table(rownames(subset(res, logFC >= log2(1.5))),
            "Tables/scATAC/Morula_inner_cell_specific_open_chromatin_fc1.5.txt", col.names = F, row.names = F, quote = F)
write.table(rownames(subset(res, logFC <= -log2(1.5))),
            "Tables/scATAC/Morula_outer_cell_specific_open_chromatin_fc1.5.txt", col.names = F, row.names = F, quote = F)
# output tables (pvalue 0.15: )
write.table(rownames(subset(res, logFC >= log2(1.5) & P.Value <= 0.15)),
            "Tables/scATAC/Morula_inner_cell_specific_open_chromatin_fc1.5_p0.15.txt", col.names = F, row.names = F, quote = F)
write.table(rownames(subset(res, logFC <= -log2(1.5) & P.Value <= 0.15)),
            "Tables/scATAC/Morula_outer_cell_specific_open_chromatin_fc1.5_p0.15.txt", col.names = F, row.names = F, quote = F)
# output tables (pvalue 0.1: )
write.table(rownames(subset(res, logFC >= log2(1.5) & P.Value <= 0.1)),
            "Tables/scATAC/Morula_inner_cell_specific_open_chromatin_fc1.5_p0.1.txt", col.names = F, row.names = F, quote = F)
write.table(rownames(subset(res, logFC <= -log2(1.5) & P.Value <= 0.1)),
            "Tables/scATAC/Morula_outer_cell_specific_open_chromatin_fc1.5_p0.1.txt", col.names = F, row.names = F, quote = F)
# output tables (pvalue 0.05: )
write.table(rownames(subset(res, logFC >= log2(1.5) & P.Value <= 0.05)),
            "Tables/scATAC/Morula_inner_cell_specific_open_chromatin_fc1.5_p0.05.txt", col.names = F, row.names = F, quote = F)
write.table(rownames(subset(res, logFC <= -log2(1.5) & P.Value <= 0.05)),
            "Tables/scATAC/Morula_outer_cell_specific_open_chromatin_fc1.5_p0.05.txt", col.names = F, row.names = F, quote = F)


### >>> 8. Annotate peaks
library("ChIPseeker")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb.hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene 
morula.spe.pks <- list.files("Tables/scATAC", pattern = ".bed", full.names = T)
morula.spe.pks <- as.list(morula.spe.pks)
names(morula.spe.pks) <- c("Inner", "Outer")
pks.anno <- lapply(morula.spe.pks, annotatePeak, TxDb = txdb.hg38, tssRegion = c(-3000, 3000), verbose = FALSE, annoDb = "org.Hs.eg.db")
pdf("Graphs/Annotation_peak_rawplot.pdf", width = 5, height = 5)
p1 <- plotAnnoBar(pks.anno, title = "Distribution of Peakss in Genome")
p2 <- plotDistToTSS(pks.anno, title = "Distribution of Peaks Relative to TSS")
plot_grid(p1, p2, ncol = 1)
dev.off(); rm(p1, p2)


### >>> 9. Plotting
# MA plot
pdf("Graphs/ICM_and_TE/scATAC/inner_outer_cell_scATAC_differential_peaks_MAplot.pdf", height = 3, width = 3)
na.omit(res) %>% mutate(group = case_when((logFC >= log2(2) & P.Value <= 0.1) ~ "Up.regulated", 
                                          (logFC <= -log2(2) & P.Value <= 0.1) ~ "Down.regulated",
                                          (abs(logFC) < log2(2) | P.Value > 0.1) ~ "Non.changed")) %>%
  ggplot(aes(x = AveExpr, y = logFC)) + geom_point(aes(color = group, alpha = group), size = 1) + theme_classic() +
  geom_hline(yintercept = 1, colour = "#000000", linetype = "longdash", size = 1, alpha = 0.5) +
  geom_hline(yintercept = -1, colour = "#000000", linetype = "longdash", size = 1, alpha = 0.5) +
  scale_color_manual(values = c(Up.regulated = "#EA2027", Down.regulated = "#0652DD", Non.changed = "#000000")) +
  scale_alpha_manual(values = c(Up.regulated = 0.5, Down.regulated = 0.5, Non.changed = 0.25)) +
  theme(legend.position = "none") + labs(x = expression(log[2]("Mean FPKM")), y = expression(log[2]("Fold of Change"))) +
  theme(axis.title.x = element_text(face="plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face="plain", colour = "#000000", size = 10, angle = 0),
        axis.text.y  = element_text(face="plain", colour = "#000000", size = 10, angle = 0)) +
  xlim(c(1.5, 10))
dev.off()
# heatmap (TPM)
pdf("Graphs/ICM_and_TE/scATAC/inner_outer_cell_scATAC_differential_peaks_heatmap_TPM.pdf", height = 5, width = 5)
pd <- rbind(epi.morula.filter[rownames(subset(res, logFC >= log2(0.5) & P.Value <= 0.05)),],
            epi.morula.filter[rownames(subset(res, logFC <= -log2(0.5) & P.Value <= 0.05)),])
Heatmap(pd,
        col = colorRamp2(c(0, 10, 20, 30, 40), c("#481b6d", "#414487", "#1e9b8a", "#72cf56", "#f7e620")), name = "TPM",
        cluster_rows = T, show_row_dend = T, row_dend_side = "left", row_dend_width = unit(1, "cm"), 
        clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
        show_column_names = T, column_names_rot = 45, column_names_side = "bottom",
        cluster_columns = T, show_column_dend = T, column_dend_height = unit(1, "cm"), column_dend_side = "top",
        clustering_distance_columns = "euclidean", clustering_method_columns = "complete",
        show_row_names = F, row_names_rot = 0, row_names_side = "right",
        width = unit(5, "cm"), height = unit(7, "cm"),
        column_split = epi.morula.filter.meta$SampleGroup, column_gap = unit(2, "mm"), border = T,
        row_split = c(rep("Up", nrow(subset(res, logFC >= log2(0.5) & P.Value <= 0.05))),
                      rep("Down", nrow(subset(res, logFC <= -log2(0.5) & P.Value <= 0.05)))), row_gap = unit(2, "mm"))
dev.off(); rm(pd)
# heatmap (zscore)
pdf("Graphs/ICM_and_TE/scATAC/inner_outer_cell_scATAC_differential_peaks_heatmap_zscore.pdf", height = 5, width = 5)
pd <- rbind(epi.morula.filter[rownames(subset(res, logFC >= log2(0.5) & P.Value <= 0.05)),],
            epi.morula.filter[rownames(subset(res, logFC <= -log2(0.5) & P.Value <= 0.05)),])
pd <- as.data.frame(t(apply(pd, 1, function(x){(x - mean(x))/(sd(x))})), 
                    row.names = rownames(pd))
Heatmap(pd,
        col = colorRamp2(c(-1, 0, 1), c("#20bfdc", "#ffffff", "#dc2056")), name = "Z-Score",
        cluster_rows = T, show_row_dend = T, row_dend_side = "left", row_dend_width = unit(1, "cm"), 
        clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
        show_column_names = T, column_names_rot = 45, column_names_side = "bottom",
        cluster_columns = T, show_column_dend = T, column_dend_height = unit(1, "cm"), column_dend_side = "top",
        clustering_distance_columns = "euclidean", clustering_method_columns = "complete",
        show_row_names = F, row_names_rot = 0, row_names_side = "right",
        width = unit(5, "cm"), height = unit(7, "cm"),
        column_split = epi.morula.filter.meta$SampleGroup, column_gap = unit(2, "mm"), border = T,
        row_split = c(rep("Up", nrow(subset(res, logFC >= log2(0.5) & P.Value <= 0.05))),
                      rep("Down", nrow(subset(res, logFC <= -log2(0.5) & P.Value <= 0.05)))), row_gap = unit(2, "mm"))
dev.off(); rm(pd)
