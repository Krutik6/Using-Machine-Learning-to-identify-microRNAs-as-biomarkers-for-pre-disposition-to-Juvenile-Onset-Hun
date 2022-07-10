library(DESeq2)
library(limma)
setwd("~/Documents/HD/Data/miRNA")
D_data <- read.csv("~/Documents/HD/Data/miRNA/miRNA_counts.csv", row.names = 1)
D_data <- D_data[,-c(63, 157)]
colnames(D_data)
Pheno <- read.csv("Pheno.csv", row.names = 1)
Pheno <- Pheno[-c(63, 157),]
#Create colData
Samples <- colnames(D_data)
Conditions <- Pheno$Name
Gender <- Pheno$sex
colData <- cbind(Samples, Conditions, Gender)
rownames(colData) <- colnames(D_data)
# Make DE object
dds <- DESeqDataSetFromMatrix(countData = D_data, colData = colData, 
                              design = ~Gender)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
#Qcheck
plotMDS(dds@assays@data$counts, col = as.numeric(dds$Conditions))
boxplot(dds@assays@data$counts, col = as.numeric(dds$Conditions))
#Add batch 1
B <- factor(c(rep(c(rep(c("A","B"),each=8), rep("A", 8)), each = 6), 
            rep(c("A", "C", "A"), each = 8)))
B <- B[-c(63, 157)]
dds$batch <- B
# Add Conditions and check
dds$Conditions <- factor(dds$Conditions, levels = c(unique(Pheno$Name)))
dds$Conditions
setwd("DE/Batch/")
#### code to plot from dds data
# Dispersion estimate plot
pdf(file="Dispersion_estimates_for_deseq_data.pdf")
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
plotDispEsts(dds)
dev.off()
##### variance stabilising transformation of data
vst_data = DESeq2::vst(dds, blind = TRUE, nsub = 500)
#### check colnames
colnames(vst_data)
##### PCA plot of vsd coloured by batch
library(ggplot2)
pdf(file="PCA_plot_vsd_by_batch.pdf")
plot_data = plotPCA(vst_data, intgroup="batch", returnData=TRUE)
percentVar = round(100 * attr(plot_data, "percentVar"))
ggplot(plot_data, aes(PC1, PC2, color=batch, label = name)) +
    geom_point(size=1.5) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()
#### rename data and create a data frame copy of expression data 
Fleischer <- vst_data
Fleischer_df <- as.data.frame(assay(vst_data))
#### calculate principal components
pca <- prcomp(assay(vst_data), scale.=TRUE)
tpca <- t(pca$x)
rownames(tpca) <- colnames(Fleischer_df)
##################### remove batch