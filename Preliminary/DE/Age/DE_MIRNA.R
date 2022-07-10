library(DESeq2)
library(limma)
library(dplyr)
library(tidyverse) 
library(factoextra)
setwd("/home/b7053098/Documents/HD/TimiRGeN/counts")
miRNA <- read.csv("miRNA_counts.csv", row.names = 1)
colnames(miRNA)
Pheno <- read.csv("miRNA_Pheno.csv", row.names = 1)
#remove outliers
miRNA <- miRNA[,-c(63, 96, 127, 145, 157)]
Pheno <- Pheno[-c(63, 96, 127, 145, 157),]
#remove gender, convert >Q20 -> HD and others to WT
ageonly <- function(X){
    HD <- sub(X$Name, pattern = "fe", replacement = "")
    HD <- sub(HD, pattern = "male_", replacement = "")
    HD <- sub(HD, pattern = "Q20", replacement = "WT")
    HD <- sub(HD, pattern = "Q111", replacement = "HD")
    HD <- sub(HD, pattern = "Q140", replacement = "HD")
    HD <- sub(HD, pattern = "Q175", replacement = "HD")
    HD <- sub(HD, pattern = "Q80", replacement = "HD")
    HD <- sub(HD, pattern = "Q92", replacement = "HD")
    return(HD)
}
age <- ageonly(Pheno)
# If a third of the samples (37) have a rowSum of more than 50
X <- miRNA[!rowSums(miRNA < 50) >= 37, , drop = FALSE]
m <- mapply(X, FUN=as.integer)
rownames(m) <- rownames(X)
#create Conditions
Samples <- colnames(m)
Conditions <- age
colData <- cbind(Samples, Conditions)
rownames(colData) <- colnames(m)
#make DESeq object
dds <- DESeqDataSetFromMatrix(countData = m, colData = colData, 
                              design = ~Conditions)
keep <- rowSums(counts(dds)) >= 1000
keep
dds <- dds[keep,]

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
M2 <- getDE(numC = 'HD_2m', denC = 'WT_2m')
hist(M2$padj)
M6 <- getDE(numC = 'HD_6m', denC = 'WT_6m')
hist(M6$padj)
M10 <- getDE(numC = 'HD_10m', denC = 'WT_10m')
hist(M10$padj)
################################################################################
P2 <- M2[which(M2$padj < 0.05),]
P6 <- M6[which(M6$padj < 0.000001),]
P10 <- M10[which(M10$padj < 0.05),]

U <- unique(c(rownames(P2), rownames(P6), rownames(P10)))

U2 <- M2[which(rownames(M2) %in% U == TRUE),]
U6 <- M6[which(rownames(M6) %in% U == TRUE),]
U10 <- M10[which(rownames(M10) %in% U == TRUE),]

colnames(U2) <- sub(colnames(U2), pattern = "log", replacement = "M2.log")
colnames(U2) <- sub(colnames(U2), pattern = "padj", replacement = "M2.padj")
colnames(U6) <- sub(colnames(U6), pattern = "log", replacement = "M6.log")
colnames(U6) <- sub(colnames(U6), pattern = "padj", replacement = "M6.padj")
colnames(U10) <- sub(colnames(U10), pattern = "log", replacement = "M10.log")
colnames(U10) <- sub(colnames(U10), pattern = "padj", replacement = "M10.padj")

DF <- data.frame(row.names = rownames(U2),U2$M2.log2FoldChange, U2$M2.padj,
                 U6$M6.log2FoldChange, U6$M6.padj, U10$M10.log2FoldChange, 
                 U10$M10.padj)

setwd("/home/b7053098/Documents/HD/TimiRGeN/Age/")
write.csv(DF, "miRNA_age.csv")