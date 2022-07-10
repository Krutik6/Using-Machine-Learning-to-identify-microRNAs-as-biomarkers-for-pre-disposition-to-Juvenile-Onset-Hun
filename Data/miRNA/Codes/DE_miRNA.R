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
df <- read.csv("../miRNA_counts_outliersRemoves.csv", row.names = 1)
pheno <- df[,c(1:4)]
df <- df[,-c(1:5)]
df <- as.data.frame(t(df))
df <- df[rowSums(df < 10) <= 56 , , drop = FALSE]
#remove 2m data
df <- df %>% select(-contains("2m"))
pheno <- pheno[- grep("2m", pheno$AGE),]
pheno$Genotype <- sub(pheno$Genotype, pattern = "Q20", replacement = "WT")
pheno$Genotype <- str_replace_all(pheno$Genotype, c("Q111" = "HD", 
                                                    "Q140" = "HD",
                                                    "Q175" = "HD",
                                                    "Q80" = "HD",
                                                    "Q92" = "HD"))
df <- as.matrix(df)
pheno <- as.matrix(pheno)
#create multiple DE contrast objects
DE.gender <- DESeqDataSetFromMatrix(countData = df,
                                            colData = pheno,
                                            design = ~ sex)
DE.gender$sex <- relevel(DE.gender$sex, "male")
DE.genotype <- DESeqDataSetFromMatrix(countData = df,
                                    colData = pheno,
                                    design = ~ Genotype)
DE.genotype$Genotype <- relevel(DE.genotype$Genotype, "WT" )
DE.both <- DESeqDataSetFromMatrix(countData = df,
                                      colData = pheno,
                                      design = ~ Genotype + sex)
DE.both$Genotype <- relevel(DE.both$Genotype, "WT" )
DE.both$sex <- relevel(DE.both$sex, "male")
as.data.frame( colData(DE.both))
#run DE
# gender
DE.gender.deseq <- DESeq(DE.gender)
resultsNames(DE.gender.deseq)
DE.gender.res <- results(DE.gender.deseq, contrast = c("sex", "male", "female"))
DE.gender.table <- as.data.frame(DE.gender.res)
# genotype
DE.genotype.deseq <- DESeq(DE.genotype)
resultsNames(DE.genotype.deseq)
DE.genotype.res <- results(DE.genotype.deseq, contrast = c("Genotype", "HD", "WT"))
DE.genotype.table <- as.data.frame(DE.genotype.res)
#both
DE.both.deseq <- DESeq(DE.both)
resultsNames(DE.both.deseq)
DE.both.res <- results(DE.both.deseq, contrast = c("Genotype", "HD", "WT"))
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
batch <- pheno[,4]
batch <- sub(batch, pattern = "female", replacement = 1)
batch <- sub(batch, pattern = "male", replacement = 2)
batch <- as.numeric(batch)
df.adjusted <- ComBat_seq(df, batch=batch)
#re-perform DE
#create multiple DE contrast objects
DE.gender <- DESeqDataSetFromMatrix(countData = df.adjusted,
                                    colData = pheno,
                                    design = ~ sex)
DE.gender$sex <- relevel(DE.gender$sex, "male")
DE.genotype <- DESeqDataSetFromMatrix(countData = df.adjusted,
                                      colData = pheno,
                                      design = ~ Genotype)
DE.genotype$Genotype <- relevel(DE.genotype$Genotype, "WT" )
DE.both <- DESeqDataSetFromMatrix(countData = df.adjusted,
                                  colData = pheno,
                                  design = ~ Genotype + sex)
DE.both$Genotype <- relevel(DE.both$Genotype, "WT" )
DE.both$sex <- relevel(DE.both$sex, "male")
as.data.frame( colData(DE.both))
#run DE
# gender
DE.gender.deseq <- DESeq(DE.gender)
resultsNames(DE.gender.deseq)
DE.gender.res <- results(DE.gender.deseq, contrast = c("sex", "male", "female"))
DE.gender.table <- as.data.frame(DE.gender.res)
# genotype
DE.genotype.deseq <- DESeq(DE.genotype)
resultsNames(DE.genotype.deseq)
DE.genotype.res <- results(DE.genotype.deseq, contrast = c("Genotype", "HD", "WT"))
DE.genotype.table <- as.data.frame(DE.genotype.res)
#both
DE.both.deseq <- DESeq(DE.both)
resultsNames(DE.both.deseq)
DE.both.res <- results(DE.both.deseq, contrast = c("Genotype", "HD", "WT"))
DE.both.table <- as.data.frame(DE.both.res)
#ranks
sum(DE.gender.table$padj < 0.05)
sum(DE.gender.table$pvalue < 0.05)
sum(DE.genotype.table$padj < 0.05)
sum(DE.genotype.table$pvalue < 0.05)
sum(DE.both.table$padj < 0.05)
sum(DE.both.table$pvalue < 0.05)
#read out normalised DE values
normalized_counts <- counts(DE.both.deseq, normalized=TRUE)
normalized_counts <- round(normalized_counts, 0)
normalized_counts <- t(normalized_counts)
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
write.csv(normalized_counts, "../10m_miRNA.csv")
#extract 2m data and perform DE using this data
df.2m <- read.csv("../miRNA_counts_outliersRemoves.csv", row.names = 1)
pheno.2m <- df.2m[,c(1:4)]
df.2m <- df.2m[,-c(1:5)]
df.2m <- as.data.frame(t(df.2m))
df.2m <- df.2m[rowSums(df.2m < 10) <= 56 , , drop = FALSE]
df.2m <- df.2m %>% select(-contains("10m"))
df.2m <- as.matrix(df.2m)
pheno.2m <- pheno.2m[- grep("10m", pheno.2m$AGE),]
pheno.2m$Genotype <- sub(pheno.2m$Genotype, pattern = "Q20",
                         replacement = "WT")
pheno.2m$Genotype <- str_replace_all(pheno.2m$Genotype, c("Q111" = "HD", 
                                                    "Q140" = "HD",
                                                    "Q175" = "HD",
                                                    "Q80" = "HD",
                                                    "Q92" = "HD"))
pheno.2m <- as.matrix(pheno.2m)
batch.2m <- pheno.2m[,4]
batch.2m <- sub(batch.2m, pattern = "female", replacement = 1)
batch.2m <- sub(batch.2m, pattern = "male", replacement = 2)
batch.2m <- as.numeric(batch.2m)
df.adjusted.2m <- ComBat_seq(df.2m, batch=batch.2m)
DE.gender.2m <- DESeqDataSetFromMatrix(countData = df.adjusted.2m,
                                     colData = pheno.2m,
                                     design = ~ sex)
DE.both.2m <- DESeqDataSetFromMatrix(countData = df.adjusted.2m,
                                  colData = pheno.2m,
                                  design = ~ Genotype + sex)
#DE to check combat worked
DE.gender.deseq.2m <- DESeq(DE.gender.2m)
resultsNames(DE.gender.deseq.2m)
DE.gender.res.2m <- results(DE.gender.deseq.2m, contrast = c("sex", "male", "female"))
DE.gender.table.2m <- as.data.frame(DE.gender.res.2m)
DE.both.deseq.2m <- DESeq(DE.both.2m)
resultsNames(DE.both.deseq.2m)
DE.both.res.2m <- results(DE.both.deseq.2m, contrast = c("sex", "male", "female"))
DE.both.table.2m <- as.data.frame(DE.both.res.2m)
sum(DE.gender.table.2m$padj < 0.05)
sum(DE.gender.table.2m$pvalue < 0.05)
sum(DE.both.table.2m$padj < 0.05)
sum(DE.both.table.2m$pvalue < 0.05)
#save genes found in the 10m analysis only!
df.10m <- read.csv("../10m_miRNA.csv", row.names = 1)
genes.10m <- colnames(df.10m)
genes.2m <- rownames(df.2m)
#should only be one difference - samples
diff <- genes.10m[which(genes.10m %in% genes.2m == FALSE)]
DE.both.deseq.2m <- DESeq(DE.both.2m)
normalized_counts.2m <- counts(DE.both.deseq.2m, normalized=TRUE)
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
write.csv(normalized_counts.2m, "../2m_miRNA.csv")
