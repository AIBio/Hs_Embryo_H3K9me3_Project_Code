################################################# ~ Good Luck ~ #########################################################
###### >>>>>>                       TITLE: H3K9me3 Project Revised Figure3 analysis                         <<<<<< ######

# - Figure 3 content:
# - 1st part: Global setting
# - 2nd part: Library
# - 3rd part: Load and preprocessing various data (YHW)
# - 4th part: Integrated with Tang et.al data (YHW)
# - 5th part: SVA expression level (YHW)
# - 6th part: Seurat pipeline (YHW)
# - 7th part: Tracjetory analysis (YHW)



### ========================
### 1st part: Global setting
### ========================
setwd("/home/yhw/bioinfo/project-mine/Embryo.93")
#save.image("R/CodeData/revised_Fig3.RData")
#load("R/CodeData/revised_Fig3.RData")



### =================
### 2nd part: Library
### =================
.libPaths("/home/cmq/software/anaconda3/envs/rs.p37.r4111/lib/R/library")
library(devtools)
library(BiocManager)
### >>> 1. CRAN packages
cran.pks <- c("forcats", "dplyr", "tidyr", "tidyverse", "stringr", "ggplot2", "ggrepel", "circlize", "cowplot", "ggpubr", 
              "preprocessCore", "FactoMineR", "factoextra", "amap", "dendextend", "biclust", "Seurat")
for (pks in cran.pks) {
  if (!require(pks, character.only = TRUE)) {
    install.packages(pks)
    library(pks, character.only = T)
  } else {
    library(pks, character.only = T)
  }
}; rm(pks)
### >>> 2. Bioconductor packages
bioc.pks <- c("limma", "sva", "ComplexHeatmap", "monocle", "preprocessCore")
for (pks in bioc.pks) {
  if (!require(pks, character.only = TRUE)) {
    BiocManager::install(pks)
    library(pks, character.only = T)
  } else {
    library(pks, character.only = T)
  }
}; rm(pks)
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
FullSet.GO <- function(species, genelist, basename, genetype, outdir){
  library("ggplot2"); library("forcats"); library("dplyr"); library("stringr")
  dir.create(outdir, showWarnings = FALSE, recursive = T)
  if (species == "human"){
    go.data <- clusterProfiler.GO(species, genelist, genetype)
  } else if (species == "mouse"){
    go.data <- clusterProfiler.GO(species, genelist, genetype)
  } else if (species == "rat"){
    go.data <- clusterProfiler.GO(species, genelist, genetype)
  }
  for (i in seq(1, length(go.data), 1)){
    if (!is.null(go.data[[i]])){
      plot.data <- as.data.frame(go.data[[i]]) %>%
        mutate(GeneRatio.new = sapply(GeneRatio, function(x) eval(parse(text=as.character(x)))),
               Description = as.factor(Description), `-Log10(pvalue)` = -log10(pvalue))
      pdf.name <- c(paste(basename, "_GO_BP.pdf", sep = ""), paste(basename, "_GO_CC.pdf", sep = ""),
                    paste(basename, "_GO_MF.pdf", sep = ""), paste(basename, "_KEGG_pathway.pdf", sep = ""),
                    paste(basename, "_Reactome_pathway.pdf", sep = ""))
      title.name <- c("GO Enrichment Analysis: Biological Process", "GO Enrichment Analysis: Cellular Component",
                      "GO Enrichment Analysis: Molecular Function", "Enrichment Analysis: KEGG Pathways",
                      "Enrichment Analysis: Reactome Pathways")
      write.table(plot.data, file.path(outdir, paste(str_split_fixed(pdf.name[i], "\\.", 2)[, 1], ".txt", sep = "")), sep = "\t", row.names = F, quote = F)
      pdf(file.path(outdir, pdf.name[i]), height = 6, width = 10)
      print(plot.data[1:20, ] %>% mutate(Description = fct_reorder(Description, GeneRatio.new)) %>%
              ggplot(aes(x = GeneRatio.new, y = Description)) + 
              geom_point(aes(size = Count, color = `-Log10(pvalue)`), alpha = 1, shape = 16) + 
              scale_colour_gradient(low = "#0652DD", high = "#EA2027") + theme_classic() +
              ggtitle(label = title.name[i]) + ylab("GO Term Discription") + xlab("Gene Ratio") +
              theme(plot.title = element_text(size = 18)) +
              labs(colour = "-log10(Pvalue)", size = "Gene Counts"))
      dev.off()
    } else {
      print("Find a NULL results")
    }
  }
  return(go.data)
}
clusterProfiler.GO <- function(species, genelist, genetype){
  library("GO.db"); library("ReactomePA"); library("DOSE"); library("org.Hs.eg.db")
  library("org.Mm.eg.db"); library("org.Rn.eg.db"); library("clusterProfiler"); library("AnnotationDbi")
  if (species == "human"){
    go.database <- "org.Hs.eg.db"
    org.short <- "hsa"
    entrezid <- mapIds(x = org.Hs.eg.db, keys = unique(as.character(genelist)),
                       column = "ENTREZID", keytype = genetype, multiVals = "first")
  } else if (species == "mouse"){
    go.database <- "org.Mm.eg.db"
    org.short <- "mmu"
    entrezid <- mapIds(x = org.Mm.eg.db, keys = unique(as.character(genelist)),
                       column = "ENTREZID", keytype = genetype, multiVals = "first")
  } else if (species == "rat"){
    go.database <-  "org.Rn.eg.db"
    org.short <- "rno"
    entrezid <- mapIds(x = org.Rn.eg.db, keys = unique(as.character(genelist)),
                       column = "ENTREZID", keytype = genetype, multiVals = "first")
  } else {
    print("You may need to define a new species dataset in the function")
  }
  all.results <- list()
  for (index in seq(1, 3, 1)){
    GO.names <- c("BP", "CC", "MF")
    bcf <- enrichGO(gene = unique(as.character(genelist)), keyType = genetype,
                    OrgDb = go.database, ont = GO.names[index], pvalueCutoff = 1)
    if (!is.null(bcf)){
      bcf.res <- as.data.frame(bcf@result)
      all.results[[index]] <- bcf.res[order(bcf.res$pvalue), ]
    } else {
      all.results[[index]] <- NULL
    }
  }
  kegg <- enrichKEGG(unique(as.character(na.omit(entrezid))), organism = org.short, pvalueCutoff = 1, minGSSize = 1)
  if (!is.null(kegg)){
    kegg.res <- as.data.frame(kegg@result)
    all.results[[4]] <- kegg.res[order(kegg.res$pvalue), ]
  } else {
    all.results[[4]] <- NULL
  }
  reactome <- enrichPathway(unique(as.character(na.omit(entrezid))), organism = species, pvalueCutoff = 1, minGSSize = 1)
  if (!is.null(reactome)){
    reactome.res <- as.data.frame(reactome@result)
    all.results[[5]] <- reactome.res[order(reactome.res$pvalue), ]
  } else {
    all.results[[5]] <- NULL
  }
  # Disease Ontology (DO) Semantic and Enrichment analysis
  #do <- enrichDO(unique(as.character(na.omit(entrezid))), ont = "DO", pvalueCutoff = 0.05)
  #do_results <- as.data.frame(do@result[1:10, ])
  #do_results$Positive_log10pvalue <- -log10(do_results$pvalue)
  names(all.results) <- c("GO.BP", "GO.CC", "GO.MF", "KEGG", "REACTOME")
  return(all.results)
}
HVGfinderTopn <- function(expr.mat, n){
  expr.mat.means <- rowMeans(expr.mat)
  expr.mat.vars <- apply(expr.mat, 1, var)
  expr.mat.cv2 <- expr.mat.vars/expr.mat.means^2
  library("statmod")
  minMeanForFit <- unname(quantile(expr.mat.means[which(expr.mat.cv2 > 0.3)], 0.95))
  useForFit <- expr.mat.means >= minMeanForFit
  fit <- glmgam.fit(cbind(a0 = 1, altilde = 1/expr.mat.means[useForFit]), expr.mat.cv2[useForFit])
  a0 <- unname(fit$coefficients["a0"])
  a1 <- unname(fit$coefficients["altilde"])
  df <- ncol(expr.mat) - 1
  afit <- a1/expr.mat.means+a0
  varFitRatio <- expr.mat.vars/(afit*expr.mat.means^2)
  varorder <- order(varFitRatio, decreasing = T)
  oed <- expr.mat[varorder, ]
  expr.mat.hvg <- oed[1:n, ]
  return(rownames(expr.mat.hvg))
}
edgeR.scRNA <- function(count, meta, g1, g2, lfc, sig, dir){
  # - count: count matrix of gene or repeat;
  # - meta: sample annotation
  # - comparison: g1 vs g2;
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
  contrast <- factor(c(rep(1, length(g2.cols)), rep(2, length(g1.cols))))
  dge <- DGEList(counts = tar.count, samples = tar.meta, group = contrast)
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
  names(degs.all.res) <- c("all", "sig", "up", "down")
  write.csv(tt$table,  paste(dir, "/edgeR_t1_", g1, "_Vs_", g2, "_all_gene_results.csv", sep = ""), row.names = F)
  write.csv(degs.all,  paste(dir, "/edgeR_t2_", g1, "_Vs_", g2, "_all_sig_results.csv", sep = ""), row.names = F)
  write.csv(degs.up,   paste(dir, "/edgeR_t3_", g1, "_Vs_", g2, "_up_sig_gene_lfc", lfc, "_sig", sig, ".csv", sep = ""), row.names = F)
  write.csv(degs.down, paste(dir, "/edgeR_t4_", g1, "_Vs_", g2, "_down_sig_gene_lfc", lfc, "_sig", sig, ".csv", sep = ""), row.names = F)
  return(degs.all.res)
}
edgeR.Ge.ET.LRT <- function(count, meta, g1, g2, lfc, sig, dir, spe, nor, with.rep){
  library("edgeR")
  library("dplyr")
  library("ggplot2")
  library("cowplot")
  dir.create(dir, showWarnings = F, recursive = T)
  con.index <- grep(paste("^", g2, "$", sep = ""), as.character(meta$SampleGroup))
  exp.index <- grep(paste("^", g1, "$", sep = ""), as.character(meta$SampleGroup))
  count <- count[, c(con.index, exp.index)]
  meta <- meta[c(con.index, exp.index), ]
  nor <- nor[, c(con.index, exp.index)]
  nor$ENSEMBL <- rownames(nor)
  if (isTRUE(with.rep)){
    contrast <- factor(c(rep(1, length(con.index)), rep(2, length(exp.index))))
    deg <- DGEList(counts = count, group = contrast)
  } else if (isFALSE(with.rep)){
    deg <- DGEList(counts = count, samples = meta)
    deg$samples$group <- meta$SampleGroup
  }
  keep <- filterByExpr(deg, min.count = 10)
  deg <- deg[keep, , keep.lib.sizes = FALSE]
  deg <- calcNormFactors(deg, method = "TMM")
  lcpm <- as.data.frame(cpm(deg, log = F))
  lcpm$ENSEMBL <- rownames(lcpm)
  if (isTRUE(with.rep)){
    design <- model.matrix(~contrast)
    deg <- estimateDisp(deg, design)
    deg.fit <- glmFit(deg, design)
    deg.test <- glmLRT(deg.fit, coef = 2)
  } else if (isFALSE(with.rep)){
    design <- model.matrix(~meta$SampleGroup)
    deg <- estimateGLMCommonDisp(deg, method = "deviance", robust = TRUE, subset = NULL)
    deg.test <- exactTest(deg, pair = c(g2, g1), rejection.region = "doubletail")
  }
  test.res <- topTags(deg.test, sort.by = "logFC", adjust.method = "fdr", n = nrow(deg))
  test.res$table$ENSEMBL <- rownames(test.res$table)
  test.res$table$SYMBOL <- DbiIDtrans(rownames(test.res$table), "ENSEMBL", "SYMBOL", spe)
  compare.info <- data.frame(Name = c("Comparison", "Testing for DE genes", "Method to adjust p.value"),
                             info = c(paste(g1, "vs", g2, sep = " "), test.res$test, test.res$adjust.method))
  res <- merge(test.res$table, lcpm, by = "ENSEMBL")
  res <- merge(res, nor, by = "ENSEMBL")
  colnames(res) <- gsub("\\.x", "\\.TMM", colnames(res))
  colnames(res) <- gsub("\\.y", "\\.TPM", colnames(res))
  res <- res %>% mutate(Log2MeanTMM = log2(rowMeans(res[, grep("TMM", colnames(res))])),
                        Log2MeanTPM = log2(rowMeans(res[, grep("TPM", colnames(res))])),
                        color = case_when((logFC >= lfc & FDR <= sig) ~ "Up.regulated", 
                                          (logFC <= -lfc & FDR <= sig) ~ "Down.regulated",
                                          (abs(logFC) < lfc | FDR > sig) ~ "Non.changed"),
                        trans = case_when((logFC >= lfc & FDR <= sig) ~ 1, 
                                          (logFC <= -lfc & FDR <= sig) ~ 1,
                                          (abs(logFC) < lfc | FDR > sig) ~ 0.5))
  pdf(paste(dir, "/edgeR_p1_", g1, "_Vs_", g2, "_all_gene_MAplot_rawplot.pdf", sep = ""), height = 5, width = 15)
  for (value in c("Log2MeanTMM", "Log2MeanTPM", "logCPM")){
    p <- ggplot(res, aes_string(x = value, y = "logFC")) + geom_point(aes(color = color, alpha = trans), size = 2) + theme_classic() +
      geom_hline(yintercept = lfc, colour = "#000000", linetype = "longdash", size = 1, alpha = 0.5) +
      geom_hline(yintercept = -lfc, colour = "#000000", linetype = "longdash", size = 1, alpha = 0.5) +
      scale_color_manual(name = "Color", values = c(Up.regulated = "#EA2027", Down.regulated = "#0652DD", Non.changed = "#000000")) +
      theme(legend.position = "none") + labs(x = expression(log[2]("Normalized Count")), y = expression(log[2]("Fold of Change"))) +
      theme(axis.title.x = element_text(face="plain", colour = "#000000", size = 12, angle = 0),
            axis.title.y = element_text(face="plain", colour = "#000000", size = 12, angle = 90),
            axis.text.x  = element_text(face="plain", colour = "#000000", size = 10, angle = 0),
            axis.text.y  = element_text(face="plain", colour = "#000000", size = 10, angle = 0))
    assign(paste("p", grep(value, c("Log2MeanTMM", "Log2MeanTPM", "logCPM")), sep = ""), p)
  }
  print(plot_grid(p1, p2, p3, ncol = 3, align = "h", labels = c("Log2MeanTMM", "Log2MeanTPM", "Log2MeanCPM")))
  dev.off()
  all.sig <- subset(res, PValue <= sig & abs(logFC) >= lfc)
  up.sig <- subset(res, PValue <= sig & logFC >= lfc)
  down.sig <- subset(res, PValue <= sig & logFC <= -(lfc))
  write.csv(compare.info, paste(dir, "/edgeR_t1_", g1, "_Vs_", g2, "_comparision_info.csv", sep = ""),  row.names = F)
  write.csv(lcpm, paste(dir, "/edgeR_t2_", g1, "_Vs_", g2, "_normalized_count_matrix.csv", sep = ""), row.names = F)
  write.csv(res, paste(dir, "/edgeR_t3_", g1, "_Vs_", g2, "_all_gene_results.csv", sep = ""), row.names = F)
  write.csv(all.sig, paste(dir, "/edgeR_t4_", g1, "_Vs_", g2, "_all_sig_gene_results.csv", sep = ""), row.names = F)
  write.csv(up.sig, paste(dir, "/edgeR_t5_", g1, "_Vs_", g2, "_up_fc", lfc, "_sig", sig, ".csv", sep = ""), row.names = F)
  write.csv(down.sig, paste(dir, "/edgeR_t6_", g1, "_Vs_", g2, "_down_fc", lfc, "_sig", sig, ".csv", sep = ""), row.names = F)
  output <- list(compare.info, lcpm, res, all.sig, up.sig, down.sig)
  names(output) <- c("compare.info", "lcpm", "res", "all.sig", "up.sig", "down.sig")
  return(output)
}


### ===================================================
### 3rd part: Load and preprocessing various data (YHW)
### ===================================================

### >>> 1. Setting output dir
outdir <- file.path(getwd(), "R/Graphs/revised_Fig3")
if (! dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}


### >>> 2. Load count data (RNA-seq)
# - gene
abembryo.ge.count <- read.table("R/RawData/featurecount/Abnormal_embryo_all_samples_gene_count_id_matrix.txt", header = T, sep = "\t")
colnames(abembryo.ge.count)[-1:-6] <- gsub("gene.", "embryo", str_split_fixed(colnames(abembryo.ge.count)[-1:-6], "_", 8)[,7])
# - repeat (gene + repeat)
abembryo.re.count <- read.table("R/RawData/featurecount/Abnormal_embryo_all_samples_repeat_count_matrix.txt", header = T, sep = "\t")
colnames(abembryo.re.count)[-1:-6] <- gsub("repeat.", "embryo", str_split_fixed(colnames(abembryo.re.count)[-1:-6], "_", 8)[,7])
abembryo.re.tpm <- CountToTpm(abembryo.re.count[,-1:-6], abembryo.re.count$Length)
abembryo.re.tpm$Repeatid <- abembryo.re.count$Geneid
# - Tang et al. 2013 (gene)
tang.ge <- readRDS("R/RawData/featurecount/GSE36552_sce_obj_from_CX.rds")
tang.ge.count <- as.data.frame(counts(tang.ge))
tang.ge.count$Geneid <- rownames(tang.ge.count)
# - GSE101571 (gene + repeat)
gse101571.ge.count <- read.table("R/RawData/featurecount/GSE101571_all_samples_gene_count_id_matrix.txt", header = T, row.names = 1)
colnames(gse101571.ge.count)[-1:-5] <- c("GV.oocyte.1", "GV.oocyte.2", "MII.oocyte.1", "MII.oocyte.2", "2Cell.1", "2Cell.2", "2Cell.3",
                                         "4Cell.1", "4Cell.2", "8Cell.1", "8Cell.2", "ICM.1", "ICM.2")
gse101571.ge.tpm <- CountToTpm(gse101571.ge.count[,-1:-5], gse101571.ge.count$Length)
gse101571.re.count <- read.table("R/RawData/featurecount/GSE101571_all_samples_repeat_count_matrix.txt", header = T, row.names = 1)
colnames(gse101571.re.count)[-1:-5] <- c("GV.oocyte.1", "GV.oocyte.2", "MII.oocyte.1", "MII.oocyte.2", "2Cell.1", "2Cell.2", "2Cell.3",
                                         "4Cell.1", "4Cell.2", "8Cell.1", "8Cell.2", "ICM.1", "ICM.2")
gse101571.re.tpm <- CountToTpm(gse101571.re.count[,-1:-5], gse101571.re.count$Length)
gse101571.re.count$Geneid <- rownames(gse101571.re.count)


### >>> 3. Load gene list
# - ZGA gene list
hs.major.zga <- read.table("R/Tables/genelist/hs_zga_gene_pc_lncrna_ensembl.bed") %>% separate(V4, c("EnsemblID", "Symbol"), sep = ":")
# - gene list targeted by SVA-derived enhancers
hs.resk9.gene <- read.table("/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer/genes_interacted_with_c8_RE_k93_enhancers_in_8cell_full.txt", 
                            sep = "\t") %>% 
  separate(col = V16, into = c("Geneid", "Symbol"), sep = ":")
hs.resk27.gene <- read.table("/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer/genes_interacted_with_c8_RE_k273_enhancers_in_8cell_full.txt", 
                             sep = "\t") %>% 
  separate(col = V16, into = c("Geneid", "Symbol"), sep = ":")
hs.pes.gene <- read.table("/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer/genes_interacted_with_c8_PE_enhancers_in_8cell_full.txt", 
                          sep = "\t") %>% 
  separate(col = V16, into = c("Geneid", "Symbol"), sep = ":")
# - enhancers (SVAs from PEs and REs)
hs.pes.sva <- read.table("/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer/8cell_primed_enhancer_derived_by_SVA.bed", 
                         sep = "\t")
hs.res.k27.sva <- read.table("/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer/8cell_reprogrammed_enhancer_H3K27me3_derived_by_SVA.bed", 
                             sep = "\t")
hs.res.k9.sva <- read.table("/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer/8cell_reprogrammed_enhancer_H3K9me3_derived_by_SVA.bed", 
                            sep = "\t")
# - enhancers-gene pairs (SVAs-derived enhancers)
hs.sva.gene <- list()
for (i in 1:4) {
  files <- list.files("/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer", "_repeat.txt", full.names = T)
  hs.sva.gene[[i]] <- read.table(files[i], sep = "\t") %>% separate(V7, c("ENSEMBL", "SYMBOL"), sep = ":")
  colnames(hs.sva.gene[[i]])[13] <- "Repeatid"
}; rm(i, files)
names(hs.sva.gene) <- gsub("_enhancers_in_8cell_full_anno_by_repeat.txt", "",
                           gsub("genes_interacted_with_sva_", "", 
                                list.files("/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer", "_repeat.txt")))
# - APOPTOSIS gene (from KEGG)
hs.apoptosis <- read.table("R/Tables/cited_data/KEGG_APOPTOSIS.txt")
hs.apoptosis$Geneid <- DbiIDtrans(hs.apoptosis$V1, "SYMBOL", "ENSEMBL", "human")
# - cell cycle genes
cell.cycle.genes <- read.csv("R/Tables/cited_data/human_cell_cycle.csv")
s.genes <- subset(cell.cycle.genes, phase == "S")$geneID
g2m.genes <- subset(cell.cycle.genes, phase == "G2/M")$geneID
# - K9-related modifier
#ENSG00000066135 KDM4A
#ENSG00000127663 KDM4B
#ENSG00000107077 KDM4C
#ENSG00000186280 KDM4D
#ENSG00000235268 KDM4E
#ENSG00000255855 KDM4F
#ENSG00000143379 SETDB1
#ENSG00000136169 SETDB2
#ENSG00000101945 SUV39H1
#ENSG00000152455 SUV39H2
#ENSG00000276043 UHRF1
#ENSG00000094916 HP1α
#ENSG00000108468 HP1β
hs.k9.factor <- c("ENSG00000066135", "ENSG00000127663", "ENSG00000107077", "ENSG00000186280", "ENSG00000235268", "ENSG00000255855", 
                  "ENSG00000143379", "ENSG00000136169", "ENSG00000101945", "ENSG00000152455", 
                  "ENSG00000276043", "ENSG00000094916", "ENSG00000108468")
names(hs.k9.factor) <- c("KDM4A", "KDM4B", "KDM4C", "KDM4D", "KDM4E", "KDM4F",
                         "SETDB1", "SETDB2", "SUV39H1", "SUV39H2", "UHRF1", "HP1α", "HP1β")
#EZH1 ENSG00000108799
#EZH2 ENSG00000106462
#KDM6A ENSG00000147050
#KDM6B ENSG00000132510
#KDM7A ENSG00000006459
#PHF8/KDM7B ENSG00000172943
hs.k27.factor <- c("ENSG00000108799", "ENSG00000106462", "ENSG00000147050", 
                   "ENSG00000132510", "ENSG00000006459", "ENSG00000172943")
#EZH1 ENSG00000108799
#EZH2 ENSG00000106462
#EED ENSG00000074266
#RBBP4/RbAp48 ENSG00000162521
#SUZ12 ENSG00000178691
#JARID2 ENSG00000008083
#AEBP2 ENSG00000139154
hs.prc.complex <- c("ENSG00000108799", "ENSG00000106462", "ENSG00000074266", "ENSG00000162521", 
                    "ENSG00000178691", "ENSG00000008083", "ENSG00000139154")
#CREBBP/KAT3A ENSG00000005339
#EP300/KAT3B ENSG00000100393
hs.k27ac.factor <- c("ENSG00000005339", "ENSG00000100393")

### >>> 3. Filtering cells with abnormal karyotype
abembryo.ge.count <- abembryo.ge.count[,-grep("embryo8", colnames(abembryo.ge.count))]
abembryo.re.count <- abembryo.re.count[,-grep("embryo8", colnames(abembryo.re.count))]
abembryo.re.tpm <- abembryo.re.tpm[,-grep("embryo8", colnames(abembryo.re.tpm))]



### ===============================================
### 4th part: Integrated with Tang et.al data (YHW)
### ===============================================

# - Object:
# - 1. integrated abnormal scRNA-seq data with Tang et.al scRNA-seq data;
# - 2. PCA analysis to infer cell state;
# - 3. Perform differential expression analysis of gene and repeats;
# - 4. scRNA-seq pipeline to infer cell cycle;
# - 5. scRNA-seq pipeline to perform trajectory analysis;

### >>> 1. Setting output dir
outdir <- file.path(getwd(), "R/Graphs/revised_Fig3/integrat_with_tang")
if (! dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}


### >>> 2. Merge scRNA-seq data (abnormal + Tang et al)
# prepare count data and metadata
abembryo.tang <- merge(abembryo.ge.count, tang.ge.count, by = "Geneid")
rownames(abembryo.tang) <- abembryo.tang$Geneid
abembryo.tang.tpm <- CountToTpm(abembryo.tang[,-1:-6], abembryo.tang$Length)
abembryo.tang <- abembryo.tang[,c(-1:-6, -52:-55)]
abembryo.tang <- abembryo.tang[,-grep("ESC", colnames(abembryo.tang))]
abembryo.tang.meta <- data.frame(cell_name = colnames(abembryo.tang),
                                 cell_type = c(rep("unknown", 11), str_split_fixed(colnames(abembryo.tang)[-1:-11], " ", 2)[,1]), 
                                 batch = c(rep("A", 11), rep("B", ncol(abembryo.tang)-11)),
                                 row.names = colnames(abembryo.tang))
abembryo.tang.meta$CellType <- gsub("-", "", abembryo.tang.meta$cell_type)
abembryo.tang.tpm <- abembryo.tang[,colnames(abembryo.tang)]


### >>> 3. Normalization and remove batch effect
# quantitle normalization
abembryo.tang.tpm.quan <- as.data.frame(normalize.quantiles(as.matrix(abembryo.tang.tpm)), row.names = rownames(abembryo.tang.tpm))
colnames(abembryo.tang.tpm.quan) <- colnames(abembryo.tang.tpm)
# batch correction and PCA
sva.modcombat <- model.matrix(~1, data = abembryo.tang.meta)
abembryo.tang.tpm.sva <- ComBat(dat = as.matrix(abembryo.tang.tpm[rowSums(abembryo.tang.tpm)>=30,]), batch = abembryo.tang.meta$batch,
                                mod = sva.modcombat, par.prior = TRUE, prior.plots = FALSE)
abembryo.tang.tpm.sva <- abembryo.tang.tpm.sva[rowSums(abembryo.tang.tpm.sva)>=30,]


### >>> 4. PCA analysis
plot.data.pca <- PCA(t(abembryo.tang.tpm.sva), scale.unit = TRUE, ncp = 5, graph = FALSE)
pdf(file.path(outdir, "PCA_plot_abnormal_human_embryo_integrated_with_Tang_scRNAseq_remove_hESCs.pdf"), height = 4, width = 6)
fviz_pca_ind(plot.data.pca, axes = c(1, 2),
             geom.ind = c("point"),
             pointshape = 16, pointsize = 3, alpha.ind = 0.7,
             col.ind = c(rep("abEmbryo", 11), str_split_fixed(colnames(abembryo.tang.tpm.sva)[-1:-11], " ", 2)[,1]),
             repel = T, label = "all", labelsize = 3, legend.title = "Groups") + 
  ggtitle("2D PCA-plot") + theme_classic()
dev.off(); rm(plot.data.pca)
plot.data.pca <- PCA(t(abembryo.tang.tpm.sva[HVGfinderTopn(abembryo.tang.tpm.sva, 2000),]), scale.unit = TRUE, ncp = 5, graph = FALSE)
pdf(file.path(outdir, "PCA_plot_abnormal_human_embryo_integrated_with_Tang_scRNAseq_remove_hESCs_top2000_HVGs.pdf"), height = 4, width = 6)
fviz_pca_ind(plot.data.pca, axes = c(1, 2),
             geom.ind = c("point"),
             pointshape = 16, pointsize = 3, alpha.ind = 0.7,
             col.ind = c(rep("abEmbryo", 13), str_split_fixed(colnames(abembryo.tang.tpm.sva)[-1:-13], " ", 2)[,1]),
             repel = F, label = "all", labelsize = 3, legend.title = "Groups") + 
  ggtitle("2D PCA-plot") + theme_classic()
dev.off(); rm(plot.data.pca)


### >>> 5. Differential expression analysis: abembryo vs 8cell
# gene
abembryo.tang.dea.vs.8c <- edgeR.scRNA(abembryo.tang, abembryo.tang.meta, "unknown", "8cell", 1, 0.05, "R/Tables/DEGs/dea_unknown_vs_tang8cell")
write.table(subset(pd, group == "Up.regulated")[,1:7], "Graphs/reprogrammed/abnormal_embryo/DEGs_abembryo_vs_8cell_volcano_up_gene.txt", sep = "\t", col.names = T, row.names = F)
write.table(subset(pd, group == "Down.regulated")[,1:7], "Graphs/reprogrammed/abnormal_embryo/DEGs_abembryo_vs_8cell_volcano_down_gene.txt", sep = "\t", col.names = T, row.names = F)
# repeat
abembryo.re.count <- merge(abembryo.re.count, gse101571.re.count[,-1:-5], by = "Geneid")
rownames(abembryo.re.count) <- abembryo.re.count$Geneid
abembryo.re.tpm <- CountToTpm(abembryo.re.count[,-1:-6], abembryo.re.count$Length)
abembryo.meta <- data.frame(SampleID = colnames(abembryo.re.tpm),
                            CellType = c(rep("unknown", 13), str_split_fixed(colnames(abembryo.re.tpm)[-1:-13], "\\.", 2)[,1]),
                            row.names = colnames(abembryo.re.tpm))

### >>> 6. Plotting
# venn diagram
venn.diagram(x = list(ZGA = hs.major.zga$V1,
                      REs = hs.res.zga.gene$Geneid,
                      PEs = hs.res.zga.gene$Geneid,
                      abembryo.down = subset(pd, group == "Down.regulated")$Geneid,
                      abembryo.up = subset(pd, group == "Up.regulated")$Geneid),
             output = TRUE, imagetype="png", filename = "Graphs/reprogrammed/abnormal_embryo/Venn_diagram.png",
             height = 2000, width = 2000, resolution = 400,
             cat.cex = 1, cat.default.pos = "outer",
             cat.fontfamily = "sans", cat.col = c("#000000", "#000000", "#000000", "#000000", "#000000"),
             fill = c(alpha("#ed5113", 1), alpha("#34e1a5", 1), alpha("#dc2056", 1), alpha("#085fbb", 1), alpha("#34ade1", 1)),
             lty = 1, lwd = 2.5, col = c("#ffffff", "#ffffff", "#ffffff", "#ffffff", "#ffffff"),
             cex = 1, label.col = "#000000", margin = 0.1)

# volcano plot
pd.tpm <- abembryo.ge.tpm
pd.tpm$Geneid <- rownames(pd.tpm)
pd <- abembryo.tang.dea.vs.8c$all %>% rownames_to_column(var = "Geneid")
pd$Symbol <- DbiIDtrans(genelist = pd$Geneid, intype = "ENSEMBL", outtype = "SYMBOL", "human")
pd <- merge(pd, pd.tpm, by = "Geneid"); rm(pd.tpm)
pd <- pd %>% mutate(abembryo = rowMeans(pd[,grep("embryo", colnames(pd))]),
                    `8cell` = rowMeans(pd[,grep("8-cell", colnames(pd))]))
pd <- pd %>% na.omit() %>%
  mutate(log10Pvalue = -log10(PValue)) %>% 
  mutate(group = case_when((logFC >= log2(2) & PValue <= 0.05) ~ "Up.regulated", 
                           (logFC <= -log2(2) & PValue <= 0.05) ~ "Down.regulated",
                           (abs(logFC) < log2(2) | PValue > 0.05) ~ "Non.changed")) %>% 
  mutate(discard = case_when(group == "Down.regulated" & `8cell` >= 10 ~ "keep", 
                             group == "Up.regulated" & abembryo >= 20 ~ "keep",
                             group == "Non.changed" & abembryo >= 5 ~ "keep")) %>% na.omit
pdf(file.path(outdir, "Volcano_plot_to_show_DEGs_between_abnormal_embryo_vs_8cell.pdf"), height = 4, width = 5)
pd %>%
  ggplot(aes(x = logFC, y = log10Pvalue)) + geom_point(aes(color = group, alpha = group), size = 2) + theme_classic() +
  geom_vline(xintercept = c(-1, 1), colour = "#000000", linetype = "longdash", size = 1, alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), colour = "#000000", linetype = "longdash", size = 1, alpha = 0.5) +
  scale_color_manual(values = c(Up.regulated = "#EA2027", Down.regulated = "#0652DD", Non.changed = "#000000")) +
  scale_alpha_manual(values = c(Up.regulated = 0.5, Down.regulated = 0.5, Non.changed = 0.25)) +
  theme(legend.position = "none") + labs(x = expression(log[2]("Fold of Change")), y = "-log2(Pvalue)") +
  theme(axis.title.x = element_text(face="plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face="plain", colour = "#000000", size = 10, angle = 0),
        axis.text.y  = element_text(face="plain", colour = "#000000", size = 10, angle = 0)) +
  geom_label_repel(aes(label = Symbol),
                   data = subset(pd, Symbol %in% c("ZSCAN4", "ZSCAN5B", "DUXA", "YY1", "CTCF", "KLF4", "KLF5", "SOX2", "KDM4D",
                                                   "TUBB4B", "KHDC3L", "WAPL", "TPX2",
                                                   "TENT4A", "TENT4B", "EIF3A", 
                                                   "HNRNPD", "SNORD70", "SNORD11", "SRSF1", "SRSF2", "SRSF3")),
                   box.padding = 0.2, point.padding = 0.2, segment.color = "#130f40",
                   label.size = 0.2, fill = "white", size = 3, max.overlaps = 50)
dev.off(); rm(pd)
# - GO analysis based on filtered DEGs
abembryo.tang.dea.vs.8c.up <- FullSet.GO(species = "human", genelist = subset(pd, group == "Up.regulated")$Symbol, 
                                         basename = "deg_unknown_vs_tang8cell_up", genetype = "SYMBOL", outdir = "R/Tables/GO/dea_unknown_vs_tang8cell")
abembryo.tang.dea.vs.8c.down <- FullSet.GO(species = "human", genelist = subset(pd, group == "Down.regulated")$Symbol, 
                                           basename = "deg_unknown_vs_tang8cell_down", genetype = "SYMBOL", outdir = "R/Tables/GO/dea_unknown_vs_tang8cell")
# - visualization
# GO:0022613                                                          ribonucleoprotein complex biogenesis 8.115156e-11
# GO:0033044                                                         regulation of chromosome organization 5.005769e-10
# GO:0098813                                                                nuclear chromosome segregation 9.438331e-09
# GO:1901990                                             regulation of mitotic cell cycle phase transition 1.775292e-08
# GO:0000070                                                          mitotic sister chromatid segregation 1.276758e-07
# GO:0000280                                                                              nuclear division 1.553370e-07
# GO:0008380                                                                                  RNA splicing 1.689539e-07
# GO:0000082                                                         G1/S transition of mitotic cell cycle 1.733995e-07
# GO:0006260                                                                               DNA replication 2.445850e-07
# GO:0042254                                                                           ribosome biogenesis 1.358467e-06
pdf(file.path(outdir, "Bar_plot_to_show_GO_term_down_genes_between_abnormal_embryo_vs_8cell.pdf"), height = 4, width = 8)
abembryo.tang.dea.vs.8c.down$GO.BP[c("GO:0022613", "GO:0033044", "GO:0098813", "GO:1901990", "GO:0000070", 
                                     "GO:0000280", "GO:0008380", "GO:0000082", "GO:0006260", "GO:0042254",
                                     "GO:0016570"),] %>% 
  mutate(Log10Pvalue = -log10(pvalue)) %>%
  ggplot(aes(x = Log10Pvalue, y = fct_reorder(Description, Log10Pvalue))) + 
  labs(x = "-Log10.pvalue", y = "GO biological process") +
  geom_bar(stat = "identity", fill = "#3498db", color = "#ecf0f1", width = 0.75) + 
  theme_classic() +
  geom_vline(xintercept = c(-log10(0.05)), color = "#FF1300", linetype = "longdash", size = 1, alpha = 0.8) + 
  theme(axis.title.x = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 90, hjust = 1), 
        axis.text.x  = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 0, hjust = 1),
        axis.text.y  = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 0, hjust = 1))
dev.off()
pdf(file.path(outdir, "Bar_plot_to_show_GO_term_up_genes_between_abnormal_embryo_vs_8cell.pdf"), height = 4, width = 8)
abembryo.tang.dea.vs.8c.up$GO.BP[c("GO:0022613", "GO:0033044", "GO:0098813", "GO:1901990", "GO:0000070", 
                                   "GO:0000280", "GO:0008380", "GO:0000082", "GO:0006260", "GO:0042254"),] %>% 
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

### >>> 6. Gene expression level
# - K9-regulated factors
abembryo.tang.tpm[hs.k9.factor,] %>% 
  rownames_to_column(var = "gene.id") %>% 
  gather(key = "cell.id", value = "tpm", -gene.id) %>% 
  mutate(cell.type = str_split_fixed(cell.id, " ", 2)[,1]) %>% 
  mutate(cell.type = gsub("embryo...", "abEmbryo", cell.type)) %>% 
  mutate(cell.type = gsub("-cell", "Cell", cell.type)) %>% 
  mutate(cell.type = factor(cell.type, c("Oocyte", "Zygote", "2Cell", "4Cell", "8Cell", "abEmbryo", "Morula", "LB", "hESC"))) %>% 
  filter(cell.type != "hESC") -> pd
pd$symbol <- DbiIDtrans(pd$gene.id, "ENSEMBL", "SYMBOL", "human")
pdf(file.path(outdir, "Box_plot_to_show_K9_factors_expr_during_embryo_development.pdf"), height = 4, width = 20)
do.call("rbind", replicate(3, pd, simplify = FALSE)) %>% 
  ggplot(aes(x = cell.type, y = log2(tpm+0.1))) +
  geom_boxplot(aes(fill = cell.type), width = 0.75, lwd = 0.5) +
  facet_wrap(. ~ symbol, ncol = 10) +
  scale_fill_brewer("Blues") +
  theme_bw() +
  labs(x = "Group", y = "Log2(TPM+0.1)") +
  stat_compare_means(method = "t.test", na.rm = T, label.y = 8, comparisons = list(c("8Cell", "abEmbryo")), size = 4) +
  theme(axis.title.x = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 90), 
        axis.text.x  = element_blank(), panel.grid = element_blank(),
        axis.text.y  = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 0),
        strip.background = element_rect(colour = "black", fill = "white", size = 0.75, linetype = "solid"),
        strip.text.x = element_text(face = "plain", size = 12, colour = "black", angle = 0),
        strip.text.y = element_text(face = "plain", size = 12, colour = "black", angle = 0))
dev.off()
# - Apoptosis-related genes
abembryo.tang.tpm[hs.apoptosis$Geneid,] %>% 
  rownames_to_column(var = "gene.id") %>% 
  gather(key = "cell.id", value = "tpm", -gene.id) %>% 
  mutate(cell.type = str_split_fixed(cell.id, " ", 2)[,1]) %>% 
  mutate(cell.type = gsub("embryo...", "abEmbryo", cell.type)) %>% 
  mutate(cell.type = gsub("-cell", "Cell", cell.type)) %>% 
  mutate(cell.type = factor(cell.type, c("Oocyte", "Zygote", "2Cell", "4Cell", "8Cell", "abEmbryo", "Morula", "LB", "hESC"))) %>% 
  filter(cell.type != "hESC") -> pd
pd$symbol <- DbiIDtrans(pd$gene.id, "ENSEMBL", "SYMBOL", "human")
pdf(file.path(outdir, "Box_plot_to_show_apoptosis_genes_expr_during_embryo_development.pdf"), height = 17, width = 20)
do.call("rbind", replicate(1, pd, simplify = FALSE)) %>% 
  ggplot(aes(x = cell.type, y = log2(tpm+0.1))) +
  geom_boxplot(aes(fill = cell.type), width = 0.75, lwd = 0.5) +
  facet_wrap(. ~ symbol) +
  scale_fill_brewer("Blues") +
  theme_bw() +
  labs(x = "Group", y = "Log2(TPM+0.1)") +
  stat_compare_means(method = "t.test", na.rm = T, label.y = 8, comparisons = list(c("8Cell", "abEmbryo")), size = 4) +
  theme(axis.title.x = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 90), 
        axis.text.x  = element_blank(), panel.grid = element_blank(),
        axis.text.y  = element_text(face="plain", colour = "#1D1D1E", size = 12, angle = 0),
        strip.background = element_rect(colour = "black", fill = "white", size = 0.75, linetype = "solid"),
        strip.text.x = element_text(face = "plain", size = 12, colour = "black", angle = 0),
        strip.text.y = element_text(face = "plain", size = 12, colour = "black", angle = 0))
dev.off(); rm(pd)



### ====================================
### 5th part: SVA expression level (YHW)
### ====================================

### >>> 1. Setting output dir
outdir <- file.path(getwd(), "R/Graphs/revised_Fig3/integrat_with_xie")
if (! dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}


### >>> 2. DEA
# - create merged count data
if (all(abembryo.ge.count$Geneid == rownames(gse101571.ge.count))) {
  inte.xie.ge.count <- cbind(abembryo.ge.count[,-1:-6], gse101571.ge.count[,-1:-7])
  inte.xie.ge.tpm <- CountToTpm(inte.xie.ge.count, gse101571.ge.count$Length)
}
if (all(abembryo.re.count$Geneid == rownames(gse101571.re.count))) {
  inte.xie.re.count <- cbind(abembryo.re.count[,-1:-6], gse101571.re.count[,-1:-7])
  inte.xie.re.tpm <- CountToTpm(inte.xie.re.count, gse101571.re.count$Length)
}
# - create metadata
inte.xie.meta <- data.frame(SampleID = colnames(inte.xie.ge.count),
                            SampleGroup = c(rep("Unknown", 11), gsub("\\.(1|2|3)", "", colnames(inte.xie.ge.count)[-1:-11])),
                            CellType = c(rep("Unknown", 11), gsub("\\.(1|2|3)", "", colnames(inte.xie.ge.count)[-1:-11])),
                            row.names = colnames(inte.xie.ge.count))
# - differential expression analysis
xie.deg.ab.vs.4c <- edgeR.Ge.ET.LRT(count = inte.xie.ge.count, meta = inte.xie.meta, g1 = "Unknown", g2 = "4Cell", lfc = 1, sig = 0.05, 
                                    dir = "R/Tables/DEGs/dea_abembryo_vs_Xie4Cell", spe = "human", nor = inte.xie.ge.tpm, with.rep = T)
xie.deg.ab.vs.8c <- edgeR.Ge.ET.LRT(count = inte.xie.ge.count, meta = inte.xie.meta, g1 = "Unknown", g2 = "8Cell", lfc = 1, sig = 0.05, 
                                    dir = "R/Tables/DEGs/dea_abembryo_vs_Xie8Cell", spe = "human", nor = inte.xie.ge.tpm, with.rep = T)
length(intersect(xie.deg.ab.vs.8c$down.sig$ENSEMBL, hs.major.zga$EnsemblID))
length(intersect(xie.deg.ab.vs.8c$down.sig$SYMBOL, hs.major.zga$Symbol))
length(intersect(xie.deg.ab.vs.8c$up.sig$ENSEMBL, hs.major.zga$EnsemblID))
venn.diagram(x = list(c8.up = rownames(abembryo.xie.dea.8c.vs.4c[[3]]),
                      c8.down = rownames(abembryo.xie.dea.8c.vs.4c[[4]]),
                      abembryo.down = rownames(abembryo.xie.dea.vs.8c[[3]]),
                      abembryo.up = rownames(abembryo.xie.dea.vs.8c[[4]])),
             output = TRUE, imagetype="png", filename = "Graphs/reprogrammed/abnormal_embryo/Venn_diagram_DERs_between_abembryo_8cell.png",
             height = 2000, width = 2000, resolution = 500,
             cat.cex = 1, cat.default.pos = "outer",
             cat.fontfamily = "sans", cat.col = c("#000000", "#000000", "#000000", "#000000"),
             fill = c(alpha("#ed5113", 1), alpha("#34e1a5", 1), alpha("#dc2056", 1), alpha("#085fbb", 1)),
             lty = 1, lwd = 2.5, col = c("#ffffff", "#ffffff", "#ffffff", "#ffffff"),
             cex = 1, label.col = "#000000", margin = 0.1)


### >>> 2. Correlation between SVAs and SVAs-targeted genes
# - load repeat data
pd1 <- abembryo.re.tpm
pd2 <- gse101571.re.tpm[,-1:-2]
pd2$Repeatid <- rownames(pd2)
# - merge data
pd <- merge(pd1, pd2, by = "Repeatid"); rm(pd1, pd2)
rownames(pd) <- pd$Repeatid
pd <- pd[,-1]
pd <- pd[rowSums(pd)>=10,]
pd.meta <- data.frame(SampleID = colnames(pd),
                      SampleGroup = str_split_fixed(colnames(pd), "\\.", 2)[,1],
                      Batch = c(rep("1", 11), rep("2", 11)),
                      row.names = colnames(pd))
# - remove batch effect
sva.modcombat <- model.matrix(~1, data = pd.meta)
pd.sva <- ComBat(dat = as.matrix(pd), batch = pd.meta$Batch,
                 mod = sva.modcombat, par.prior = TRUE, prior.plots = FALSE)
pd.sva[pd.sva<0] <- 0
pd.sva <- pd.sva %>% 
  as.data.frame() %>% 
  mutate(Oocyte = rowMeans(pd.sva[,grep("oocyte", colnames(pd.sva))]),
         `2 Cell` = rowMeans(pd.sva[,grep("2Cell", colnames(pd.sva))]),
         `4 Cell` = rowMeans(pd.sva[,grep("4Cell", colnames(pd.sva))]),
         `8 Cell` = rowMeans(pd.sva[,grep("8Cell", colnames(pd.sva))]),
         ICM = rowMeans(pd.sva[,grep("ICM", colnames(pd.sva))]),
         abEmbryo = rowMeans(pd.sva[,grep("embryo", colnames(pd.sva))]))
pd.sva <- pd.sva[,23:28]
# - quantitle normalization
length(grep("SVA", rownames(pd)))
pd.quan <- as.data.frame(normalize.quantiles(as.matrix(pd)), row.names = rownames(pd))
colnames(pd.quan) <- colnames(pd)
pd.quan <- pd.quan %>% 
  mutate(Oocyte = rowMeans(pd.quan[,grep("oocyte", colnames(pd.quan))]),
         `2 Cell` = rowMeans(pd.quan[,grep("2Cell", colnames(pd.quan))]),
         `4 Cell` = rowMeans(pd.quan[,grep("4Cell", colnames(pd.quan))]),
         `8 Cell` = rowMeans(pd.quan[,grep("8Cell", colnames(pd.quan))]),
         ICM = rowMeans(pd.quan[,grep("ICM", colnames(pd.quan))]),
         abEmbryo = rowMeans(pd.quan[,grep("embryo", colnames(pd.quan))]))
pd.quan <- pd.quan[,23:28]
# - prepare gene expression
pd.gene <- CountToTpm(abembryo.ge.count[,-1:-6], abembryo.ge.count$Length)
rownames(pd.gene) <- abembryo.ge.count$Geneid
pd.gene <- gse101571.ge.tpm %>% 
  mutate(Oocyte = rowMeans(gse101571.ge.tpm[,grep("MII.oocyte", colnames(gse101571.ge.tpm))]),
         `2 Cell` = rowMeans(gse101571.ge.tpm[,grep("2Cell", colnames(gse101571.ge.tpm))]),
         `4 Cell` = rowMeans(gse101571.ge.tpm[,grep("4Cell", colnames(gse101571.ge.tpm))]),
         `8 Cell` = rowMeans(gse101571.ge.tpm[,grep("8Cell", colnames(gse101571.ge.tpm))]),
         ICM = rowMeans(gse101571.ge.tpm[,grep("ICM", colnames(gse101571.ge.tpm))]),
         abEmbryo = rowMeans(pd.gene[,grep("embryo", colnames(pd.gene))]))
pd.gene <- pd.gene[,14:19]
pd.gene <- as.data.frame(t(apply(pd.gene, 1, function(x){(x - mean(x))/(sd(x))})), row.names = rownames(pd.gene)) %>% na.omit()
pd.gene$ENSEMBL <- rownames(pd.gene)
# - prepare K9-REs
pd.k9.sva <- pd.sva[unique(hs.res.k9.sva$V7),] %>% 
  na.omit() %>% 
  rownames_to_column(var = "Repeatid") %>% 
  mutate(Repeatid = gsub("\\..", "", Repeatid)) %>% 
  filter(`4 Cell` <= 1)
pd.k9.sva <- merge(pd.k9.sva, hs.sva.gene$c8_RE_k93[,c(7:8,13)], by = "Repeatid")
pd.k9.sva <- merge(pd.k9.sva, pd.gene, by = "ENSEMBL") %>% filter(abEmbryo.x < 2 & abEmbryo.y < 1)
# - prepare K27-REs
pd.k27.sva <- pd.sva[unique(hs.res.k27.sva$V7),] %>% 
  na.omit() %>% 
  rownames_to_column(var = "Repeatid") %>% 
  mutate(Repeatid = gsub("\\..", "", Repeatid)) %>% 
  filter(`4 Cell` <= 1)
pd.k27.sva <- merge(pd.k27.sva, hs.sva.gene$c8_RE_k273[,c(7:8,13)], by = "Repeatid")
pd.k27.sva <- merge(pd.k27.sva, pd.gene, by = "ENSEMBL") %>% filter(abEmbryo.x < 2 & abEmbryo.y < 1)
# - prepare PEs
pd.pes.sva <- pd.sva[unique(hs.pes.sva$V7),] %>% 
  na.omit() %>% 
  rownames_to_column(var = "Repeatid") %>% 
  mutate(Repeatid = gsub("\\..", "", Repeatid)) %>% 
  filter(`4 Cell` <= 1)
pd.pes.sva <- merge(pd.pes.sva, hs.sva.gene$c8_PE[,c(7:8,13)], by = "Repeatid")
pd.pes.sva <- merge(pd.pes.sva, pd.gene, by = "ENSEMBL") %>% filter(abEmbryo.x > 1.5)
# - plot heatmap
pd <- rbind(pd.k9.sva, pd.k27.sva, pd.pes.sva)
row.anno.ha <- rowAnnotation(`Type` = c(rep("REs.K9", nrow(pd.k9.sva)), rep("REs.K27", nrow(pd.k27.sva)), rep("PEs", nrow(pd.pes.sva))),
                             col = list(`Type` = c("REs.K9" = "#f53e3e", "REs.K27" = "#f5f56d", "PEs" = "#2893e2")))
pdf(file.path(outdir, "Heatmap_to_show_SVAs_and_targeted_genes_expression.pdf"), height = 5, width = 5)
p <- Heatmap(as.matrix(pd[,c("Oocyte.x", "2 Cell.x", "4 Cell.x", "8 Cell.x", "abEmbryo.x", "ICM.x")]), name = "TPM",
             col = colorRamp2(seq(0, 4, 0.5), c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B")),
             cluster_rows = T, show_row_dend = T, row_dend_side = "left", row_dend_width = unit(1, "cm"), 
             clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
             cluster_columns = F, show_column_dend = T, column_dend_height = unit(1, "cm"), column_dend_side = "top",
             clustering_distance_columns = "euclidean", clustering_method_columns = "complete",
             show_column_names = T, column_names_rot = 45, column_names_side = "bottom", column_title = "SVAs (TPM)",
             show_row_names = F, row_names_rot = 0, row_names_side = "right", use_raster = T,
             left_annotation = row.anno.ha, border = T, border_gp = gpar(col = "#000000", lwd = 1.5),
             row_split = c(rep("REs.K9", nrow(pd.k9.sva)), rep("REs.K27", nrow(pd.k27.sva)), rep("PEs", nrow(pd.pes.sva))))
p + Heatmap(as.matrix(pd[,c("8 Cell.y", "abEmbryo.y", "ICM.y")]), name = "Z-score",
            col = colorRamp2(c(-1, 0, 1), c("#20bfdc", "#ffffff", "#dc2056")),
            cluster_rows = T, show_row_dend = T, row_dend_side = "left", row_dend_width = unit(1, "cm"), 
            clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
            cluster_columns = F, show_column_dend = T, column_dend_height = unit(1, "cm"), column_dend_side = "top",
            clustering_distance_columns = "euclidean", clustering_method_columns = "complete",
            show_column_names = T, column_names_rot = 45, column_names_side = "bottom", column_title = "Gene (Z-score)",
            show_row_names = F, row_names_rot = 0, row_names_side = "right", use_raster = T, 
            border = T, border_gp = gpar(col = "#000000", lwd = 1.5),
            row_split = c(rep("REs.K9", nrow(pd.k9.sva)), rep("REs.K27", nrow(pd.k27.sva)), rep("PEs", nrow(pd.pes.sva))))
dev.off(); rm(p, row.anno.ha, pd.meta, pd.k27.sva, pd.k9.sva, pd.pes.sva, pd.sva, pd.quan)



### ===============================
### 6th part: Seurat pipeline (YHW)
### ===============================

### >>> 1. Setting output dir
outdir <- file.path(getwd(), "R/Graphs/revised_Fig3/integrat_with_tang")
if (! dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}

### >>> 2. Running seurat pipeline
sr.ob <- CreateSeuratObject(counts = abembryo.tang, meta.data = abembryo.tang.meta)
sr.ob <- NormalizeData(sr.ob) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
sr.ob <- CellCycleScoring(sr.ob, s.features=s.genes, g2m.features=g2m.genes)
sr.ob <- ScaleData(sr.ob, verbose = FALSE, vars.to.regress=c('S.Score','G2M.Score')) %>% 
  RunPCA(npcs = 30, verbose = FALSE) %>% 
  JackStraw(num.replicate = 100, dims = 30) %>% 
  ScoreJackStraw(dims = 1:30)
ElbowPlot(sr.ob, ndims = 50)
dim.n <- 20
sr.ob <- RunUMAP(sr.ob, reduction = "pca", dims = 1:dim.n) %>% RunTSNE(reduction = "pca", dims = 1:dim.n) %>% 
  FindNeighbors(reduction = "pca", dims = 1:dim.n) %>% FindClusters(resolution = 0.5)
sr.harmony <- RunHarmony(sr.ob, group.by.vars = "batch") %>% FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% 
  FindClusters(resolution = 0.5) %>% RunTSNE(reduction = "harmony", dims = 1:dim.n)
sr.harmony <-  ScaleData(sr.harmony, features = rownames(sr.harmony), verbose = FALSE)
print(DimPlot(sr.harmony, reduction = "tsne", label = TRUE, repel = TRUE, group.by = c("cell_type", "ident")))

### >>> 3. Gene expression level
pd.gene <- subset(hs.apoptosis, V1 %in% c("CASP10", "CASP3", "CASP6", "CASP7", "CASP8", "CASP9", "TP53", "TNF", "NFKB1", "BAX", "BCL2"))
pdf(file.path(outdir, "Violin_plot_integrated_with_tang_to_check_apoptosis_gene_expression.pdf"), height = 5, width = 6.5)
VlnPlot(sr.harmony, features = pd.gene$Geneid, stack = T, group.by = "cell_type", flip = T)
dev.off()
pdf(file.path(outdir, "Violin_plot_integrated_with_tang_to_check_K9_earser_expression.pdf"), height = 5, width = 6.5)
VlnPlot(sr.harmony, features = hs.k9.factor[1:6], stack = T, group.by = "cell_type", flip = T)
dev.off()
pdf(file.path(outdir, "Violin_plot_integrated_with_tang_to_check_K9_writer_expression.pdf"), height = 5.5, width = 6.5)
VlnPlot(sr.harmony, features = hs.k9.factor[-1:-6], stack = T, group.by = "cell_type", flip = T)
dev.off()
pdf(file.path(outdir, "Violin_plot_integrated_with_tang_to_check_K27_modifier_expression.pdf"), height = 5, width = 6.5)
VlnPlot(sr.harmony, features = hs.k27.factor, stack = T, group.by = "cell_type", flip = T)
dev.off()
pdf(file.path(outdir, "Violin_plot_integrated_with_tang_to_check_K27_PRC_complex_expression.pdf"), height = 6, width = 6.5)
VlnPlot(sr.harmony, features = hs.prc.complex, stack = T, group.by = "cell_type", flip = T)
dev.off()
pdf(file.path(outdir, "Violin_plot_integrated_with_tang_to_check_K27ac_modifier_expression.pdf"), height = 3, width = 6.5)
VlnPlot(sr.harmony, features = hs.k27ac.factor, stack = T, group.by = "cell_type", flip = T)
dev.off()
pdf(file.path(outdir, "PieDonut_plot_to_show_Cell_cycle_and_cell_type_integrated_with_tang.pdf"), height = 5, width = 5)
PieDonut(sr.harmony@meta.data,
         r0 = getOption("PieDonut.r0", 0), r1 = getOption("PieDonut.r1", 0.7), r2 = getOption("PieDonut.r2", 0.9),
         aes(pies = cell_type, donuts = Phase), addDonutLabel = TRUE, showRatioDonut = F,
         showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.001), color = "#ffffff",
         ratioByGroup = F, labelposition = 1, pieLabelSize = 2, donutLabelSize = 2)
dev.off()



### ===================================
### 7th part: Tracjetory analysis (YHW)
### ===================================

### >>> 1. Create a directory to save Monocle results
outdir <- file.path(getwd(), "R/Graphs/revised_Fig3/monocle_outputs")
if (! dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}

### >>> 2. Create a new CellDataSet object as follows
expr.matrix <- GetAssayData(sr.harmony, assay= "RNA", slot = "counts")
sample.meta <- sr.harmony@meta.data
gene.anno <- data.frame(gene_short_name = DbiIDtrans(row.names(sr.harmony), "ENSEMBL", "SYMBOL", "human"),
                        row.names = row.names(sr.harmony),
                        gene_id = row.names(sr.harmony))

pd <- new("AnnotatedDataFrame", data = sample.meta)
fd <- new("AnnotatedDataFrame", data = gene.anno)
mncle.sub <- newCellDataSet(as.matrix(expr.matrix), phenoData = pd, featureData = fd, expressionFamily=negbinomial.size())

### >>> 3. Estimate size factors and dispersions
mncle.sub <- estimateSizeFactors(mncle.sub)
mncle.sub <- estimateDispersions(mncle.sub)

### >>> 5. Filtering low-quality cells
mncle.sub <- detectGenes(mncle.sub, min_expr = 0.1)
expressed.genes <- row.names(subset(fData(mncle.sub), num_cells_expressed >= 8))
#expressed.genes <- VariableFeatures(sr.harmony)


### >>> 6. Choose genes that define a cell's progress
diff.test.res <- differentialGeneTest(mncle.sub[expressed.genes,], fullModelFormulaStr = "~cell_type", cores = 30)
ordering.genes <- row.names(subset(diff.test.res, qval <= 0.005))
#ordering.genes <- expressed.genes
mncle.sub <- setOrderingFilter(mncle.sub, ordering.genes)
plot_ordering_genes(mncle.sub)

### >>> 7. Reduce data dimensionality
mncle.sub <- reduceDimension(mncle.sub, max_components = 2, method = 'DDRTree')

### >>> 8. Order cells along the trajectory
mncle.sub <- orderCells(mncle.sub)
# visualize the trajectory in the reduced dimensional space
pdf(file.path(outdir, "Plot_cell_trajectory_color_by_cell_type.pdf"), height = 4.5, width = 4)
plot_cell_trajectory(mncle.sub, color_by = "cell_type", cell_size = 3, cell_name_size = 2)
dev.off()

### >>> 9. Re-order the cell by manually setting root state
plot_cell_trajectory(mncle.sub, color_by = "State", cell_size = 2, cell_name_size = 2)
mncle.sub.order <- orderCells(mncle.sub, root_state = 2)
pdf(file.path(outdir, "Plot_cell_trajectory_color_by_Pseudotime.pdf"), height = 4.5, width = 4)
plot_cell_trajectory(mncle.sub.order, color_by = "Pseudotime", cell_size = 3) + scale_color_viridis_c()
dev.off()

### >>> 10. Finding Genes that Change as a Function of Pseudotime
# perform differential expression analysis
diff.test.pseu <- differentialGeneTest(mncle.sub.order[expressed.genes,], fullModelFormulaStr = "~sm.ns(Pseudotime)", cores  = 20)
diff.test.pseu <- diff.test.pseu[order(diff.test.pseu$qval),]

# plot marker genes (in scatter and heatmap plots)
pd.gene <- c("DUX4", "ZSCAN4", "NFYA", "LEUTX", "OTX2", "NANOG", "POU5F1", "SOX2")
to.be.tested <- row.names(subset(fData(mncle.sub.order), gene_short_name %in% pd.gene))
cds_subset <- mncle.sub.order[to.be.tested,]
pdf(file.path(outdir, "Plot_marker_genes_expr_along_Pseudotime_scatter.pdf"), height = 1*length(pd.gene), width = 4)
plot_genes_in_pseudotime(cds_subset, color_by = "cell_type", trend_formula = "~ sm.ns(Pseudotime, df=3)")
dev.off()

### >>> 12. Delete useless variables
rm(pd, fd, pd.gene, to.be.tested, cds_subset)
