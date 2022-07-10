library(DESeq2)
library(limma)
library(dplyr)
library(tidyverse) 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
mRNA <- read.csv("../../counts/mRNA_counts.csv", row.names = 1)
colnames(mRNA)
Pheno <- read.csv("../../counts/mRNA_Pheno.csv", row.names = 1)
#remove outliers
mRNA <- mRNA[,-c(63, 96, 127, 145, 157)]
Pheno <- Pheno[-c(63, 96, 127, 145, 157),]
#remove gender, convert >Q20 -> HD and others to WT
Qs <- Pheno$Name
# If a third of the samples (37) have a rowSum of more than 50
X <- mRNA[!rowSums(mRNA < 50) >= 37, , drop = FALSE]
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
write.table(normalized_counts, file="normalized_mRNA_counts.txt", sep="\t",
            quote=F, col.names=NA)

###
dds$Conditions <- factor(dds$Conditions, levels = unique(dds$Conditions))
dds$Conditions
dds <- DESeq(dds)
resultsNames(dds)

x <- unique(dds$Conditions)
x


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
MQ202M <- getDE(numC = 'male_Q20_2m', denC = 'male_WT_2m')
MQ802M <- getDE(numC = 'male_Q80_2m', denC = 'male_WT_2m')
MQ922M <- getDE(numC = 'male_Q92_2m', denC = 'male_WT_2m')
MQ1112M <- getDE(numC = 'male_Q111_2m', denC = 'male_WT_2m')
MQ1402M <- getDE(numC = 'male_Q140_2m', denC = 'male_WT_2m')
MQ1752M <- getDE(numC = 'male_Q175_2m', denC = 'male_WT_2m')
FQ202M <- getDE(numC = 'female_Q20_2m', denC = 'female_WT_2m')
FQ802M <- getDE(numC = 'female_Q80_2m', denC = 'female_WT_2m')
FQ922M <- getDE(numC = 'female_Q92_2m', denC = 'female_WT_2m')
FQ1112M <- getDE(numC = 'female_Q111_2m', denC = 'female_WT_2m')
FQ1402M <- getDE(numC = 'female_Q140_2m', denC = 'female_WT_2m')
FQ1752M <- getDE(numC = 'female_Q175_2m', denC = 'female_WT_2m')

MQ206m <- getDE(numC = 'male_Q20_6m', denC = 'male_WT_6m')
MQ806m <- getDE(numC = 'male_Q80_6m', denC = 'male_WT_6m')
MQ926m <- getDE(numC = 'male_Q92_6m', denC = 'male_WT_6m')
MQ1116m <- getDE(numC = 'male_Q111_6m', denC = 'male_WT_6m')
MQ1406m <- getDE(numC = 'male_Q140_6m', denC = 'male_WT_6m')
MQ1756m <- getDE(numC = 'male_Q175_6m', denC = 'male_WT_6m')
FQ206m <- getDE(numC = 'female_Q20_6m', denC = 'female_WT_6m')
FQ806m <- getDE(numC = 'female_Q80_6m', denC = 'female_WT_6m')
FQ926m <- getDE(numC = 'female_Q92_6m', denC = 'female_WT_6m')
FQ1116m <- getDE(numC = 'female_Q111_6m', denC = 'female_WT_6m')
FQ1406m <- getDE(numC = 'female_Q140_6m', denC = 'female_WT_6m')
FQ1756m <- getDE(numC = 'female_Q175_6m', denC = 'female_WT_6m')

MQ2010m <- getDE(numC = 'male_Q20_10m', denC = 'male_WT_10m')
MQ8010m <- getDE(numC = 'male_Q80_10m', denC = 'male_WT_10m')
MQ9210m <- getDE(numC = 'male_Q92_10m', denC = 'male_WT_10m')
MQ11110m <- getDE(numC = 'male_Q111_10m', denC = 'male_WT_10m')
MQ14010m <- getDE(numC = 'male_Q140_10m', denC = 'male_WT_10m')
MQ17510m <- getDE(numC = 'male_Q175_10m', denC = 'male_WT_10m')
FQ2010m <- getDE(numC = 'female_Q20_10m', denC = 'female_WT_10m')
FQ8010m <- getDE(numC = 'female_Q80_10m', denC = 'female_WT_10m')
FQ9210m <- getDE(numC = 'female_Q92_10m', denC = 'female_WT_10m')
FQ11110m <- getDE(numC = 'female_Q111_10m', denC = 'female_WT_10m')
FQ14010m <- getDE(numC = 'female_Q140_10m', denC = 'female_WT_10m')
FQ17510m <- getDE(numC = 'female_Q175_10m', denC = 'female_WT_10m')
################################################################################
MQ202M_P <- MQ202M[which(MQ202M$padj < 0.05),]
MQ802M_P <- MQ802M[which(MQ802M$padj < 0.05),]
MQ922M_P <- MQ922M[which(MQ922M$padj < 0.05),]
MQ1112M_P <- MQ922M[which(MQ1112M$padj < 0.05),]
MQ1402M_P <- MQ922M[which(MQ1402M$padj < 0.05),]
MQ1752M_P <- MQ922M[which(MQ1752M$padj < 0.05),]

MQ206m_P <- MQ206m[which(MQ206m$padj < 0.05),]
MQ806m_P <- MQ806m[which(MQ806m$padj < 0.05),]
MQ926m_P <- MQ926m[which(MQ926m$padj < 0.05),]
MQ1116m_P <- MQ926m[which(MQ1116m$padj < 0.05),]
MQ1406m_P <- MQ926m[which(MQ1406m$padj < 0.05),]
MQ1756m_P <- MQ926m[which(MQ1756m$padj < 0.05),]

MQ2010m_P <- MQ2010m[which(MQ2010m$padj < 0.05),]
MQ8010m_P <- MQ8010m[which(MQ8010m$padj < 0.05),]
MQ9210m_P <- MQ9210m[which(MQ9210m$padj < 0.05),]
MQ11110m_P <- MQ9210m[which(MQ11110m$padj < 0.05),]
MQ14010m_P <- MQ9210m[which(MQ14010m$padj < 0.05),]
MQ17510m_P <- MQ9210m[which(MQ17510m$padj < 0.05),]

FQ202M_P <- FQ202M[which(FQ202M$padj < 0.05),]
FQ802M_P <- FQ802M[which(FQ802M$padj < 0.05),]
FQ922M_P <- FQ922M[which(FQ922M$padj < 0.05),]
FQ1112M_P <- FQ922M[which(FQ1112M$padj < 0.05),]
FQ1402M_P <- FQ922M[which(FQ1402M$padj < 0.05),]
FQ1752M_P <- FQ922M[which(FQ1752M$padj < 0.05),]

FQ206m_P <- FQ206m[which(FQ206m$padj < 0.05),]
FQ806m_P <- FQ806m[which(FQ806m$padj < 0.05),]
FQ926m_P <- FQ926m[which(FQ926m$padj < 0.05),]
FQ1116m_P <- FQ926m[which(FQ1116m$padj < 0.05),]
FQ1406m_P <- FQ926m[which(FQ1406m$padj < 0.05),]
FQ1756m_P <- FQ926m[which(FQ1756m$padj < 0.05),]

FQ2010m_P <- FQ2010m[which(FQ2010m$padj < 0.05),]
FQ8010m_P <- FQ8010m[which(FQ8010m$padj < 0.05),]
FQ9210m_P <- FQ9210m[which(FQ9210m$padj < 0.05),]
FQ11110m_P <- FQ9210m[which(FQ11110m$padj < 0.05),]
FQ14010m_P <- FQ9210m[which(FQ14010m$padj < 0.05),]
FQ17510m_P <- FQ9210m[which(FQ17510m$padj < 0.05),]

U <- unique(c(rownames(MQ202M_P), rownames(MQ802M_P), rownames(MQ922M_P), 
              rownames(MQ1112M_P), rownames(MQ1402M_P), rownames(MQ1752M_P),
              rownames(MQ206m_P), rownames(MQ806m_P), rownames(MQ926m_P), 
              rownames(MQ1116m_P), rownames(MQ1406m_P), rownames(MQ1756m_P),
              rownames(MQ2010m_P), rownames(MQ8010m_P), rownames(MQ9210m_P), 
              rownames(MQ11110m_P), rownames(MQ14010m_P), rownames(MQ17510m_P),
              rownames(FQ202M_P), rownames(FQ802M_P), rownames(FQ922M_P), 
              rownames(FQ1112M_P), rownames(FQ1402M_P), rownames(FQ1752M_P),
              rownames(FQ206m_P), rownames(FQ806m_P), rownames(FQ926m_P), 
              rownames(FQ1116m_P), rownames(FQ1406m_P), rownames(FQ1756m_P),
              rownames(FQ2010m_P), rownames(FQ8010m_P), rownames(FQ9210m_P), 
              rownames(FQ11110m_P), rownames(FQ14010m_P), rownames(FQ17510m_P)))

MQ202M <- MQ202M[which(rownames(MQ202M) %in% U == TRUE),]
MQ802M <- MQ802M[which(rownames(MQ802M) %in% U == TRUE),]
MQ922M <- MQ922M[which(rownames(MQ922M) %in% U == TRUE),]
MQ1112M <- MQ1112M[which(rownames(MQ1112M) %in% U == TRUE),]
MQ1402M <- MQ1402M[which(rownames(MQ1402M) %in% U == TRUE),]
MQ1752M <- MQ1752M[which(rownames(MQ1752M) %in% U == TRUE),]

MQ206m <- MQ206m[which(rownames(MQ206m) %in% U == TRUE),]
MQ806m <- MQ806m[which(rownames(MQ806m) %in% U == TRUE),]
MQ926m <- MQ926m[which(rownames(MQ926m) %in% U == TRUE),]
MQ1116m <- MQ1116m[which(rownames(MQ1116m) %in% U == TRUE),]
MQ1406m <- MQ1406m[which(rownames(MQ1406m) %in% U == TRUE),]
MQ1756m <- MQ1756m[which(rownames(MQ1756m) %in% U == TRUE),]

MQ2010m <- MQ2010m[which(rownames(MQ2010m) %in% U == TRUE),]
MQ8010m <- MQ8010m[which(rownames(MQ8010m) %in% U == TRUE),]
MQ9210m <- MQ9210m[which(rownames(MQ9210m) %in% U == TRUE),]
MQ11110m <- MQ11110m[which(rownames(MQ11110m) %in% U == TRUE),]
MQ14010m <- MQ14010m[which(rownames(MQ14010m) %in% U == TRUE),]
MQ17510m <- MQ17510m[which(rownames(MQ17510m) %in% U == TRUE),]

FQ202M <- FQ202M[which(rownames(FQ202M) %in% U == TRUE),]
FQ802M <- FQ802M[which(rownames(FQ802M) %in% U == TRUE),]
FQ922M <- FQ922M[which(rownames(FQ922M) %in% U == TRUE),]
FQ1112M <- FQ1112M[which(rownames(FQ1112M) %in% U == TRUE),]
FQ1402M <- FQ1402M[which(rownames(FQ1402M) %in% U == TRUE),]
FQ1752M <- FQ1752M[which(rownames(FQ1752M) %in% U == TRUE),]

FQ206m <- FQ206m[which(rownames(FQ206m) %in% U == TRUE),]
FQ806m <- FQ806m[which(rownames(FQ806m) %in% U == TRUE),]
FQ926m <- FQ926m[which(rownames(FQ926m) %in% U == TRUE),]
FQ1116m <- FQ1116m[which(rownames(FQ1116m) %in% U == TRUE),]
FQ1406m <- FQ1406m[which(rownames(FQ1406m) %in% U == TRUE),]
FQ1756m <- FQ1756m[which(rownames(FQ1756m) %in% U == TRUE),]

FQ2010m <- FQ2010m[which(rownames(FQ2010m) %in% U == TRUE),]
FQ8010m <- FQ8010m[which(rownames(FQ8010m) %in% U == TRUE),]
FQ9210m <- FQ9210m[which(rownames(FQ9210m) %in% U == TRUE),]
FQ11110m <- FQ11110m[which(rownames(FQ11110m) %in% U == TRUE),]
FQ14010m <- FQ14010m[which(rownames(FQ14010m) %in% U == TRUE),]
FQ17510m <- FQ17510m[which(rownames(FQ17510m) %in% U == TRUE),]


DF <- data.frame(row.names = rownames(MQ202M), 
                 MQ202M$log2FoldChange, MQ202M$padj,
                 MQ802M$log2FoldChange, MQ802M$padj,
                 MQ922M$log2FoldChange, MQ922M$padj,
                 MQ1112M$log2FoldChange, MQ1112M$padj,
                 MQ1402M$log2FoldChange, MQ1402M$padj,
                 MQ1752M$log2FoldChange, MQ1752M$padj,
                 MQ206m$log2FoldChange, MQ206m$padj,
                 MQ806m$log2FoldChange, MQ806m$padj,
                 MQ926m$log2FoldChange, MQ926m$padj,
                 MQ1116m$log2FoldChange, MQ1116m$padj,
                 MQ1406m$log2FoldChange, MQ1406m$padj,
                 MQ1756m$log2FoldChange, MQ1756m$padj,
                 MQ2010m$log2FoldChange, MQ2010m$padj,
                 MQ8010m$log2FoldChange, MQ8010m$padj,
                 MQ9210m$log2FoldChange, MQ9210m$padj,
                 MQ11110m$log2FoldChange, MQ11110m$padj,
                 MQ14010m$log2FoldChange, MQ14010m$padj,
                 MQ17510m$log2FoldChange, MQ17510m$padj,
                 FQ202M$log2FoldChange, FQ202M$padj,
                 FQ802M$log2FoldChange, FQ802M$padj,
                 FQ922M$log2FoldChange, FQ922M$padj,
                 FQ1112M$log2FoldChange, FQ1112M$padj,
                 FQ1402M$log2FoldChange, FQ1402M$padj,
                 FQ1752M$log2FoldChange, FQ1752M$padj,
                 FQ206m$log2FoldChange, FQ206m$padj,
                 FQ806m$log2FoldChange, FQ806m$padj,
                 FQ926m$log2FoldChange, FQ926m$padj,
                 FQ1116m$log2FoldChange, FQ1116m$padj,
                 FQ1406m$log2FoldChange, FQ1406m$padj,
                 FQ1756m$log2FoldChange, FQ1756m$padj,
                 FQ2010m$log2FoldChange, FQ2010m$padj,
                 FQ8010m$log2FoldChange, FQ8010m$padj,
                 FQ9210m$log2FoldChange, FQ9210m$padj,
                 FQ11110m$log2FoldChange, FQ11110m$padj,
                 FQ14010m$log2FoldChange, FQ14010m$padj,
                 FQ17510m$log2FoldChange, FQ17510m$padj)

write.csv(DF, "mRNA_qval_gender.csv")
