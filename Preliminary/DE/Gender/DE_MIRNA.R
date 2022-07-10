library(DESeq2)
library(limma)
library(dplyr)
library(tidyverse) 
setwd("/home/b7053098/Documents/HD/TimiRGeN/counts")
miRNA <- read.csv("miRNA_counts.csv", row.names = 1)
colnames(miRNA)
Pheno <- read.csv("miRNA_Pheno.csv", row.names = 1)
#remove outliers
miRNA <- miRNA[,-c(63, 96, 127, 145, 157)]
Pheno <- Pheno[-c(63, 96, 127, 145, 157),]
#remove gender, convert >Q20 -> HD and others to WT
Genderonly <- function(X){
    HD <- sub(X$Name, pattern = "Q20", replacement = "WT")
    HD <- sub(HD, pattern = "Q111", replacement = "HD")
    HD <- sub(HD, pattern = "Q140", replacement = "HD")
    HD <- sub(HD, pattern = "Q175", replacement = "HD")
    HD <- sub(HD, pattern = "Q80", replacement = "HD")
    HD <- sub(HD, pattern = "Q92", replacement = "HD")
    return(HD)
}
gender <- Genderonly(Pheno)
# If a third of the samples (37) have a rowSum of more than 50
X <- miRNA[!rowSums(miRNA < 50) >= 37, , drop = FALSE]
m <- mapply(X, FUN=as.integer)
rownames(m) <- rownames(X)
#create Conditions
Samples <- colnames(m)
Conditions <- gender
colData <- cbind(Samples, Conditions)
rownames(colData) <- colnames(m)
#make DESeq object
dds <- DESeqDataSetFromMatrix(countData = m, colData = colData, 
                              design = ~Conditions)
# keep <- rowSums(counts(dds)) >= 1000
# keep
# dds <- dds[keep,]

plotMDS(dds@assays@data@listData$counts, col = as.numeric(dds$Conditions))
boxplot(dds@assays@data@listData$counts, col = as.numeric(dds$Conditions))

### normalised counts for ML
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="normalized_miRNA_counts.txt", sep="\t",
            quote=F, col.names=NA)

###
dds$Conditions <- factor(dds$Conditions, levels = unique(dds$Conditions))
dds$Conditions
dds <- DESeq(dds)
resultsNames(dds)

getDE <- function(numC, denC){
    res <- results(dds, contrast= c("Conditions", numC, denC))
    res_B <- suppressMessages(as.data.frame(lfcShrink(dds=dds, 
                                                      contrast=c("Conditions",
                                                                 numC,
                                                                 denC), 
                                                      res=res,
                                                      type = 'ashr')))
    return(res_B)
}
HD2 <- getDE(numC = 'male_HD_2m', denC = 'female_HD_2m')
hist(HD2$padj)
WT2 <- getDE(numC = 'male_WT_2m', denC = 'female_WT_2m')
hist(WT2$padj)
HD6 <- getDE(numC = 'male_HD_6m', denC = 'female_HD_6m')
hist(HD6$padj)
WT6 <- getDE(numC = 'male_WT_6m', denC = 'female_WT_6m')
hist(WT6$padj)
HD10 <- getDE(numC = 'male_HD_10m', denC = 'female_HD_10m')
hist(HD10$padj)
WT10 <- getDE(numC = 'male_WT_10m', denC = 'female_WT_10m')
hist(WT10$padj)
################################################################################
PHD2 <- HD2[which(HD2$padj < 0.05),]
PWT2 <- WT2[which(WT2$padj < 0.05),]
PHD6 <- HD6[which(HD6$padj < 0.05),]
PWT6 <- WT6[which(WT6$padj < 0.05),]
PHD10 <- HD10[which(HD10$padj < 0.05),]
PWT10 <- WT10[which(WT10$padj < 0.05),]

U <- unique(c(rownames(PHD2), rownames(PWT2), rownames(PHD6), 
              rownames(PWT6), rownames(PHD10), rownames(PWT10)))

UHD2 <- HD2[which(rownames(HD2) %in% U == TRUE),]
UWT2 <- WT2[which(rownames(WT2) %in% U == TRUE),]
UHD6 <- HD6[which(rownames(HD6) %in% U == TRUE),]
UWT6 <- WT6[which(rownames(WT6) %in% U == TRUE),]
UHD10 <- HD10[which(rownames(HD10) %in% U == TRUE),]
UWT10 <- WT10[which(rownames(WT10) %in% U == TRUE),]

DF <- data.frame(row.names = rownames(UHD2), 
                 UHD2$log2FoldChange, UHD2$padj,
                 UWT2$log2FoldChange, UWT2$padj,
                 UHD6$log2FoldChange, UHD6$padj,
                 UWT6$log2FoldChange, UWT6$padj,
                 UHD10$log2FoldChange, UHD10$padj,
                 UWT10$log2FoldChange, UWT10$padj)

colnames(DF) <- sub(colnames(DF), pattern = "U", replacement = "")

setwd("/home/b7053098/Documents/HD/TimiRGeN/Gender//")
write.csv(DF, "miRNA_gender.csv")
