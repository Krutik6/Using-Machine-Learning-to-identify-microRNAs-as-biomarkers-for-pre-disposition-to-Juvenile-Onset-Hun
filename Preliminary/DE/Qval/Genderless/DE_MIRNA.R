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
Qs <- function(X){
    HD <- sub(X$Name, pattern = "fe", replacement = "")
    HD <- sub(HD, pattern = "male_", replacement = "")
    return(HD)
}
Qs <- Qs(Pheno)
# If a third of the samples (37) have a rowSum of more than 50
X <- miRNA[!rowSums(miRNA < 50) >= 37, , drop = FALSE]
m <- mapply(X, FUN=as.integer)
rownames(m) <- rownames(X)
#create Conditions
Samples <- colnames(m)
Conditions <- Qs
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

x <- unique(dds$Conditions)
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
Q202M <- getDE(numC = 'Q20_2m', denC = 'WT_2m')
Q802M <- getDE(numC = 'Q80_2m', denC = 'WT_2m')
Q922M <- getDE(numC = 'Q92_2m', denC = 'WT_2m')
Q1112M <- getDE(numC = 'Q111_2m', denC = 'WT_2m')
Q1402M <- getDE(numC = 'Q140_2m', denC = 'WT_2m')
Q1752M <- getDE(numC = 'Q175_2m', denC = 'WT_2m')

Q206m <- getDE(numC = 'Q20_6m', denC = 'WT_6m')
Q806m <- getDE(numC = 'Q80_6m', denC = 'WT_6m')
Q926m <- getDE(numC = 'Q92_6m', denC = 'WT_6m')
Q1116m <- getDE(numC = 'Q111_6m', denC = 'WT_6m')
Q1406m <- getDE(numC = 'Q140_6m', denC = 'WT_6m')
Q1756m <- getDE(numC = 'Q175_6m', denC = 'WT_6m')

Q2010m <- getDE(numC = 'Q20_10m', denC = 'WT_10m')
Q8010m <- getDE(numC = 'Q80_10m', denC = 'WT_10m')
Q9210m <- getDE(numC = 'Q92_10m', denC = 'WT_10m')
Q11110m <- getDE(numC = 'Q111_10m', denC = 'WT_10m')
Q14010m <- getDE(numC = 'Q140_10m', denC = 'WT_10m')
Q17510m <- getDE(numC = 'Q175_10m', denC = 'WT_10m')
################################################################################
Q202M_P <- Q202M[which(Q202M$padj < 0.05),]
Q802M_P <- Q802M[which(Q802M$padj < 0.05),]
Q922M_P <- Q922M[which(Q922M$padj < 0.05),]
Q1112M_P <- Q922M[which(Q1112M$padj < 0.05),]
Q1402M_P <- Q922M[which(Q1402M$padj < 0.05),]
Q1752M_P <- Q922M[which(Q1752M$padj < 0.05),]

Q206m_P <- Q206m[which(Q206m$padj < 0.000001),]
Q806m_P <- Q806m[which(Q806m$padj < 0.000001),]
Q926m_P <- Q926m[which(Q926m$padj < 0.000001),]
Q1116m_P <- Q926m[which(Q1116m$padj < 0.000001),]
Q1406m_P <- Q926m[which(Q1406m$padj < 0.000001),]
Q1756m_P <- Q926m[which(Q1756m$padj < 0.000001),]

Q2010m_P <- Q2010m[which(Q2010m$padj < 0.05),]
Q8010m_P <- Q8010m[which(Q8010m$padj < 0.05),]
Q9210m_P <- Q9210m[which(Q9210m$padj < 0.05),]
Q11110m_P <- Q9210m[which(Q11110m$padj < 0.05),]
Q14010m_P <- Q9210m[which(Q14010m$padj < 0.05),]
Q17510m_P <- Q9210m[which(Q17510m$padj < 0.05),]

U <- unique(c(rownames(Q202M_P), rownames(Q802M_P), rownames(Q922M_P), 
              rownames(Q1112M_P), rownames(Q1402M_P), rownames(Q1752M_P),
              rownames(Q206m_P), rownames(Q806m_P), rownames(Q926m_P), 
              rownames(Q1116m_P), rownames(Q1406m_P), rownames(Q1756m_P),
              rownames(Q2010m_P), rownames(Q8010m_P), rownames(Q9210m_P), 
              rownames(Q11110m_P), rownames(Q14010m_P), rownames(Q17510m_P)))

Q202M <- Q202M[which(rownames(Q202M) %in% U == TRUE),]
Q802M <- Q202M[which(rownames(Q802M) %in% U == TRUE),]
Q922M <- Q202M[which(rownames(Q922M) %in% U == TRUE),]
Q1112M <- Q202M[which(rownames(Q1112M) %in% U == TRUE),]
Q1402M <- Q202M[which(rownames(Q1402M) %in% U == TRUE),]
Q1752M <- Q202M[which(rownames(Q1752M) %in% U == TRUE),]

Q206m <- Q206m[which(rownames(Q206m) %in% U == TRUE),]
Q806m <- Q206m[which(rownames(Q806m) %in% U == TRUE),]
Q926m <- Q206m[which(rownames(Q926m) %in% U == TRUE),]
Q1116m <- Q206m[which(rownames(Q1116m) %in% U == TRUE),]
Q1406m <- Q206m[which(rownames(Q1406m) %in% U == TRUE),]
Q1756m <- Q206m[which(rownames(Q1756m) %in% U == TRUE),]

Q2010m <- Q2010m[which(rownames(Q2010m) %in% U == TRUE),]
Q8010m <- Q2010m[which(rownames(Q8010m) %in% U == TRUE),]
Q9210m <- Q2010m[which(rownames(Q9210m) %in% U == TRUE),]
Q11110m <- Q2010m[which(rownames(Q11110m) %in% U == TRUE),]
Q14010m <- Q2010m[which(rownames(Q14010m) %in% U == TRUE),]
Q17510m <- Q2010m[which(rownames(Q17510m) %in% U == TRUE),]


DF <- data.frame(row.names = rownames(Q202M), 
                 Q202M$log2FoldChange, Q202M$padj,
                 Q802M$log2FoldChange, Q802M$padj,
                 Q922M$log2FoldChange, Q922M$padj,
                 Q1112M$log2FoldChange, Q1112M$padj,
                 Q1402M$log2FoldChange, Q1402M$padj,
                 Q1752M$log2FoldChange, Q1752M$padj,
                 Q206m$log2FoldChange, Q206m$padj,
                 Q806m$log2FoldChange, Q806m$padj,
                 Q926m$log2FoldChange, Q926m$padj,
                 Q1116m$log2FoldChange, Q1116m$padj,
                 Q1406m$log2FoldChange, Q1406m$padj,
                 Q1756m$log2FoldChange, Q1756m$padj,
                 Q2010m$log2FoldChange, Q2010m$padj,
                 Q8010m$log2FoldChange, Q8010m$padj,
                 Q9210m$log2FoldChange, Q9210m$padj,
                 Q11110m$log2FoldChange, Q11110m$padj,
                 Q14010m$log2FoldChange, Q14010m$padj,
                 Q17510m$log2FoldChange, Q17510m$padj)


setwd("/home/b7053098/Documents/HD/TimiRGeN/Qval/Genderless/")
write.csv(DF, "miRNA_qval_genderless.csv")
