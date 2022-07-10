setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#load libs
source("functions_HD.R")
library(dplyr)
library(DESeq2)
library(stringr)
library(data.table)
library(SummarizedExperiment)
library(sva)
#load data
df <- read.csv("../mRNA_counts_geneNames.csv", row.names = 1)
pheno <- df[,c(1:4)]
df <- df[,-c(1:4)]
df <- df %>% mutate_all(as.integer)
df <- as.data.frame(t(df))
df <- df[rowSums(df < 10) <= 56 , , drop = FALSE]
#remove 2m data
df <- df %>% dplyr::select(-contains("2m"))
pheno <- pheno[- grep("2m", pheno$age),]
pheno$genotype <- sub(pheno$genotype, pattern = "Q20", replacement = "WT")
pheno$genotype <- str_replace_all(pheno$genotype, c("Q111" = "HD", 
                                                    "Q140" = "HD",
                                                    "Q175" = "HD",
                                                    "Q80" = "HD",
                                                    "Q92" = "HD"))
df <- as.matrix(df)
pheno <- as.matrix(pheno)
#create multiple DE contrast objects
DE.gender <- DESeqDataSetFromMatrix(countData = df,
                                            colData = pheno,
                                            design = ~ Sex)
#DE.gender$Sex <- relevel(DE.gender$Sex, "male")
#DE.gender$Sex
DE.genotype <- DESeqDataSetFromMatrix(countData = df,
                                    colData = pheno,
                                    design = ~ genotype)
#keep WT as condition 1 as it will be denominator
DE.genotype$genotype <- relevel(DE.genotype$genotype, "WT" )
DE.both <- DESeqDataSetFromMatrix(countData = df,
                                      colData = pheno,
                                      design = ~ genotype + Sex)
DE.both$genotype <- relevel(DE.both$genotype, "WT" )
as.data.frame(colData(DE.both))
#run DE
# gender
DE.gender.deseq <- DESeq(DE.gender)
resultsNames(DE.gender.deseq)
DE.gender.res <- results(DE.gender.deseq, contrast = c("Sex", "male", "female"))
DE.gender.table <- as.data.frame(DE.gender.res)
# genotype
DE.genotype.deseq <- DESeq(DE.genotype)
resultsNames(DE.genotype.deseq)
DE.genotype.res <- results(DE.genotype.deseq, contrast = c("genotype", "WT", "HD"))
DE.genotype.table <- as.data.frame(DE.genotype.res)
#both
DE.both.deseq <- DESeq(DE.both)
resultsNames(DE.both.deseq)
DE.both.res <- results(DE.both.deseq, contrast = c("genotype", "WT", "HD"))
DE.both.table <- as.data.frame(DE.both.res)
#ranks
sum(DE.gender.table$padj < 0.05)
sum(DE.gender.table$pvalue < 0.05)
sum(DE.genotype.table$padj < 0.05)
sum(DE.genotype.table$pvalue < 0.05)
sum(DE.both.table$padj < 0.05)
sum(DE.both.table$pvalue < 0.05)
################################################################################
#comBAT
batch <- pheno[,1]
batch <- sub(batch, pattern = "female", replacement = 1)
batch <- sub(batch, pattern = "male", replacement = 2)
batch <- as.numeric(batch)
df.adjusted <- ComBat_seq(df, batch=batch)
#re-perform DE
#create multiple DE contrast objects
DE.gender <- DESeqDataSetFromMatrix(countData = df.adjusted,
                                    colData = pheno,
                                    design = ~ Sex)
DE.genotype <- DESeqDataSetFromMatrix(countData = df.adjusted,
                                      colData = pheno,
                                      design = ~ genotype)
DE.genotype$genotype <- relevel(DE.genotype$genotype, "WT")
DE.both <- DESeqDataSetFromMatrix(countData = df.adjusted,
                                  colData = pheno,
                                  design = ~ genotype + Sex)
DE.both$Genotype <- relevel(DE.both$genotype, "WT")
as.data.frame(colData(DE.both))
#run DE
#gender
DE.gender.deseq <- DESeq(DE.gender)
resultsNames(DE.gender.deseq)
DE.gender.res <- results(DE.gender.deseq, contrast = c("Sex", "male", "female"))
DE.gender.table <- as.data.frame(DE.gender.res)
# genotype
DE.genotype.deseq <- DESeq(DE.genotype)
resultsNames(DE.genotype.deseq)
DE.genotype.res <- results(DE.genotype.deseq, contrast = c("genotype", "WT", "HD"))
DE.genotype.table <- as.data.frame(DE.genotype.res)
#both
DE.both.deseq <- DESeq(DE.both)
resultsNames(DE.both.deseq)
DE.both.res <- results(DE.both.deseq, contrast = c("genotype", "WT", "HD"))
DE.both.table <- as.data.frame(DE.both.res)
#ranks
sum(DE.gender.table$padj < 0.05)
sum(DE.gender.table$pvalue < 0.05)
sum(DE.genotype.table$padj < 0.05)
sum(DE.genotype.table$pvalue < 0.05)
sum(DE.both.table$padj < 0.05)
sum(DE.both.table$pvalue < 0.05)
#keep genes with a log2fc value > 0.2 or < 0.2
variantgenesup <- subset(DE.both.table, log2FoldChange > 0.2)
variantgenesdwn <- subset(DE.both.table, log2FoldChange < -0.2)
variantgenes <- rbind(variantgenesup, variantgenesdwn)
#read out normalised DE values
normalized_counts <- counts(DE.both.deseq, normalized=TRUE)
normalized_counts <- normalized_counts[which(rownames(normalized_counts) %in% rownames(variantgenes)),]
normalized_counts <- round(normalized_counts, 0)
normalized_counts <- t(normalized_counts)
normalized_counts <- as.data.frame(normalized_counts)
Samples <- rownames(normalized_counts)
Samples <- sub(Samples, pattern = "female_", replacement = "")
Samples <- sub(Samples, pattern = "male_", replacement = "")
Samples <- sub(Samples, pattern = "_10m...", replacement = "")
Samples <- stringr::str_replace_all(Samples, "Q111", "HD")
Samples <- stringr::str_replace_all(Samples, "Q175", "HD")
Samples <- stringr::str_replace_all(Samples, "Q140", "HD")
Samples <- stringr::str_replace_all(Samples, "Q92", "HD")
Samples <- stringr::str_replace_all(Samples, "Q80", "HD")
Samples <- stringr::str_replace_all(Samples, "Q20", "WT")
normalized_counts <- cbind(normalized_counts, Samples)
write.csv(normalized_counts, "../10m_mRNA.csv")
################################################################################
#extract 2m data
df.2m <- read.csv("../mRNA_counts_geneNames.csv", row.names = 1)
df.10m <- read.csv("../10m_mRNA.csv", row.names = 1)
df.10m <- as.data.frame(t(df.10m))
pheno.2m <- df.2m[,c(1:4)]
df.2m <- df.2m[,-c(1:4)]
df.2m <- as.data.frame(t(df.2m))
df.2m <- df.2m %>% mutate_all(as.integer)
df.2m <- df.2m %>% select(-contains("10m"))
pheno.2m <- pheno.2m[- grep("10m", pheno.2m$age),]
pheno.2m <- as.data.frame(pheno.2m)
df.2m <- df.2m[which(rownames(df.2m) %in% rownames(df.10m)),]
df.2m <- as.matrix(df.2m)
#create pheno object
pheno.2m$genotype <- sub(pheno.2m$genotype, pattern = "Q20", 
                         replacement = "WT")
pheno.2m$genotype <- str_replace_all(pheno.2m$genotype, c("Q111" = "HD", 
                                                    "Q140" = "HD",
                                                    "Q175" = "HD",
                                                    "Q80" = "HD",
                                                    "Q92" = "HD"))
pheno.2m <- as.matrix(pheno.2m)
batch.2m <- pheno.2m[,1]
batch.2m <- sub(batch.2m, pattern = "female", replacement = 1)
batch.2m <- sub(batch.2m, pattern = "male", replacement = 2)
batch.2m <- as.numeric(batch.2m)
df.adjusted.2m <- ComBat_seq(df.2m, batch=batch.2m)
#creatr DE objects
DE.gender.2m <- DESeqDataSetFromMatrix(countData = df.adjusted.2m,
                                     colData = pheno.2m,
                                     design = ~ Sex)
DE.both.2m <- DESeqDataSetFromMatrix(countData = df.adjusted.2m,
                                  colData = pheno.2m,
                                  design = ~ genotype + Sex)
DE.both.2m$genotype <- relevel(DE.both.2m$genotype, "WT" )
DE.both.2m$Sex <- relevel(DE.both.2m$Sex, "male")
as.data.frame(colData(DE.both.2m))
#DE to check combat worked
DE.gender.deseq.2m <- DESeq(DE.gender.2m)
resultsNames(DE.gender.deseq.2m)
DE.gender.res.2m <- results(DE.gender.deseq.2m, contrast = c("Sex", "male", "female"))
DE.gender.table.2m <- as.data.frame(DE.gender.res.2m)
DE.both.deseq.2m <- DESeq(DE.both.2m)
resultsNames(DE.both.deseq.2m)
DE.both.res.2m <- results(DE.both.deseq.2m, contrast = c("Sex", "male", "female"))
DE.both.table.2m <- as.data.frame(DE.both.res.2m)
sum(DE.gender.table.2m$padj < 0.05)
sum(DE.gender.table.2m$pvalue < 0.05)
sum(DE.both.table.2m$padj < 0.05)
sum(DE.both.table.2m$pvalue < 0.05)
DE.both.res.2m <- results(DE.both.deseq.2m, contrast = c("genotype", "WT", "HD"))
DE.both.table.2m <- as.data.frame(DE.both.res.2m)
sum(DE.both.table.2m$padj < 0.05)
sum(DE.both.table.2m$pvalue < 0.05)
#save data and only export genes (+samples) from 10m ananlysis
DE.both.deseq.2m <- DESeq(DE.both.2m)
normalized_counts.2m <- counts(DE.both.deseq.2m, normalized=TRUE)
normalized_counts.2m <- normalized_counts.2m[which(rownames(normalized_counts.2m) %in% rownames(df.10m)),]
normalized_counts.2m <- round(normalized_counts.2m, 0)
normalized_counts.2m <- t(normalized_counts.2m)
Samples <- rownames(normalized_counts.2m)
Samples <- sub(Samples, pattern = "female_", replacement = "")
Samples <- sub(Samples, pattern = "male_", replacement = "")
Samples <- sub(Samples, pattern = "_2m...", replacement = "")
Samples <- stringr::str_replace_all(Samples, "Q111", "HD")
Samples <- stringr::str_replace_all(Samples, "Q175", "HD")
Samples <- stringr::str_replace_all(Samples, "Q140", "HD")
Samples <- stringr::str_replace_all(Samples, "Q92", "HD")
Samples <- stringr::str_replace_all(Samples, "Q80", "HD")
Samples <- stringr::str_replace_all(Samples, "Q20", "WT")
normalized_counts.2m <- cbind(normalized_counts.2m, Samples)
write.csv(normalized_counts.2m, "../2m_mRNA.csv")
