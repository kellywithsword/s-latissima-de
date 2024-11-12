## INITIALIZATION ##
# Remove everything from memory
rm(list=ls())
# Set working directory
workdir <- "/Users/kellydeweese/Documents/GradSchool/NuzhdinData/S_latissima/lindell_RNAseq/DE"
be_so_dir <- paste(workdir, "/just_BE_SO", sep="")
setwd(workdir)


## PREPARE DATA FOR DESEQ2 ##
# Load in sample conditions and IDs
samples <- read.table("deseq_samples_file.txt", header=FALSE)
# samples <- read.table("deseq_samples_file_just_BE_SO.txt", header=FALSE)
names(samples) <- c("lifestage", "run")
rownames(samples) <- samples$run
# Reorder samples data frame by run column sorted alphabetically
samples <- samples[order(samples$run),]
files <- list.files(workdir, pattern=".genes.results")
# files <- list.files(be_so_dir, pattern=".genes.results")
names(files) <- unlist(lapply(strsplit(files, "\\."),
                              function(l) l[[1]]))
names(files) <- samples$run
samples$lifestage <- as.factor(samples$lifestage)
# Save DESeq2-formatted samples table as text file
write.table(samples, file="samples_DESeq2.txt", quote=FALSE)
# write.table(samples, file="just_BE_SO_samples_DESeq2.txt", quote=FALSE)
# Use tximport to extract quantification data relevant to DESeq2
library("tximport")
library("readr")
txi <- tximport(files, type="rsem",
                txIn=FALSE, txOut=FALSE)
zero_length <- (apply(txi$length, 1, min) == 0)
txi$length <- txi$length[!zero_length,]
txi$abundance <- txi$abundance[!zero_length,]
txi$counts <- txi$counts[!zero_length,]
# Save txi to disc
saveRDS(txi, file="txi_all.rds")
# saveRDS(txi, file="txi_just_BE_SO.rds")


## RUN DESEQ2 ##
txi <- readRDS("txi_all.rds")
library("DESeq2")
samples <- read.table("samples_DESeq2.txt", stringsAsFactors=TRUE)
input_dds <- DESeqDataSetFromTximport(txi,
                                colData = samples,
                                design = ~ lifestage)
dds <- DESeq(input_dds)
res <- results(dds)
# Save DESEq2 objects to disc
saveRDS(input_dds, file="rounded_counts_deseq_object_all.rds")
saveRDS(dds, file="dds_all.rds")
saveRDS(res, file="res_all.rds")


## SHRINKING ##
resApe <- lfcShrink(dds, coef=2, type="apeglm")
resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")
# Save each shrunken result
saveRDS(resApe, file="res_all_ape.rds")
saveRDS(resNorm, file="res_all_norm.rds")
saveRDS(resAsh, file="res_all_ash.rds")


## IMPORT DESEQ2 OBJECTS ##
library("DESeq2")
samples <- read.table("samples_DESeq2.txt")
# samples <- read.table("just_BE_SO_samples_DESeq2.txt")
# DESeq2 object (DESeq2 not yet run)
input_dds <- readRDS("rounded_counts_deseq_object_all.rds")
# DESeq2 object after running DESeq2 on rounded counts
dds <- readRDS("dds_all.rds")
# results of running DESeq2 on rounded counts
res <- readRDS("res_all.rds")
# shrunken results
resApe <- readRDS("res_all_ape.rds")
resNorm <- readRDS("res_all_norm.rds")
resAsh <- readRDS("res_all_ash.rds")


## VISUALIZATION ##
# MA plots
# Raw data
png("all_MA_plot.png", height=900, width=900, res=110)
# png("just_BE_SO_MA_plot.png", height=900, width=900, res=110)
DESeq2::plotMA(res, main="Sorus vs. non-reproductive tissue")
dev.off()
# Shrunken LFCs
png("all_MA_plot_Ape.png", height=900, width=900, res=110)
# png("just_BE_SO_MA_plot_Ape.png", height=900, width=900, res=110)
DESeq2::plotMA(resApe, 
               main="Sorus vs. non-reproductive tissue\nshrink = approx. posterior est.")
dev.off()
png("all_MA_plot_Norm.png", height=900, width=900, res=110)
# png("just_BE_SO_MA_plot_Norm.png", height=900, width=900, res=110)
DESeq2::plotMA(resNorm, 
               main="Sorus vs. non-reproductive tissue\nshrink = normal dist.")
dev.off()
png("all_MA_plot_Ash.png", height=900, width=900, res=110)
# png("just_BE_SO_MA_plot_Ash.png", height=900, width=900, res=110)
DESeq2::plotMA(resAsh, 
               main="Sorus vs. non-reproductive tissue\nshrink = adaptive")
dev.off()
# Volcano plot
library(EnhancedVolcano)
volcano_label <-  gsub(".*\\|", "", rownames(res))
# volcano_label <-  gsub(".*\\#_", "", volcano_label)
png("all_volcano_plot.png", height=900, width=900, res=110)
# png("just_BE_SO_volcano_plot.png", height=900, width=900, res=110)
EnhancedVolcano(res,
                lab = volcano_label,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Sorus vs. non-reproductive tissue')
dev.off()


## EXTRACT GENE LISTS ##
# Remove NAs from analysis
no_na_res <- res[!is.na(res$log2FoldChange),]
no_na_res <- res[!is.na(res$padj),]
no_na_dds <- dds[!is.na(res$log2FoldChange),]
no_na_dds <- dds[!is.na(res$padj),]
# Plot normalized counts of genes of interest
png("glycosyltransferase_gene_counts.png", height=900, width=900, res=110)
plotCounts(no_na_dds, gene=which.min(no_na_res$padj), intgroup="lifestage")
dev.off()
# Inspect low adjusted P-values & high LFC values
rownames(no_na_res[-log10(no_na_res$padj) > 5 & no_na_res$log2FoldChange > 0,])


