setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(DESeq2)
library(sva)
library(stringr)
library(plyr)
library(limma)
#prepare data
df <- read.table("GSE108396_PREDICTHD_smallRNAseq_CSF.RPM.tsv", row.names = 1, header = TRUE, 
                 fill=TRUE)
names <- read.table("names.txt")
names <- sub(names, pattern = " [miRNA-seq]", replacement = "", fixed = TRUE)
colnames(df) <- names
#check for outliers and normalisations
boxplot(df)
plotMDS(df)
remove <- c("CSF-HD-NS00061044-batch4-rep2",
            "CSF-HD-NS00025003-batch3-rep1")
df <- df[ , -which(names(df) %in% remove)]
plotMDS(df)
# no clear batch effect
names <- read.table("names.txt")
names <- sub(names, pattern = " [miRNA-seq]", replacement = "", fixed = TRUE)
names <- as.data.frame(names)
names <- names[which(names$names %in% remove == FALSE),]
x <- str_match(names, "-\\s*(.*?)\\s*-")
x <- gsub(x, pattern = "-", replacement = "")
colData <- cbind(names, x[,1])
rownames(colData) <- colData[,1]
colnames(colData) <- c("Samples", "Conditions")
#prepare data for DE
df <- round(df, 0)
DEobj <- DESeqDataSetFromMatrix(countData = df,
                                colData = colData,
                                design = ~ Conditions)
#perform DE
DEobj.deseq <- DESeq(DEobj)
resultsNames(DEobj.deseq)
DEobj.res <- results(DEobj.deseq, contrast = c("Conditions", "HD", "control"))
DEobj.table <- as.data.frame(DEobj.res)
DEobj.table <- as.data.frame(DEobj.table)
write.csv(DEobj.table, "CSF_dong_DE.csv")
