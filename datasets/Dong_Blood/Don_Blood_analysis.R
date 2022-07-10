setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(sva)
library(stringr)
library(plyr)
library(limma)
library(dplyr)
df <- read.table("GSE108395_table.txt", row.names = 1, header = TRUE, 
                 fill=TRUE)

names <- read.table("GSE108395_names.txt")
colnames(names) <- names
comparenames <- cbind(colnames(df), colnames(names))
colnames(df) <- names
#check normaised and outliers
boxplot(df)
plotMDS(df)
remove <- c("blood-control-NS00026404-batch1-rep1", 
           "blood-control-NS00026404-batch1-rep2",
           "blood-HD-NS00026563-batch1-rep1",
           "blood-control-NS00026037-batch1-rep1",
           "blood-HD-NS00026311-batch1-rep1",
           "blood-HD-NS00026555-batch1-rep1")
df <- df[ , -which(names(df) %in% remove)]
boxplot(df)
plotMDS(df)
#correct for batch 1 and batch 2
names <- as.data.frame(t(names))
names$condition <- str_match(names$V1, "-\\s*(.*?)\\s*-")
names$condition <- gsub(names$condition, pattern = "-", replacement = "")
#names$bool <- names$condition
#names$bool <- revalue(names$bool, c("HD"=1))
#names$bool <- revalue(names$bool, c("control"=0))
#names$batch <- gsub(names$V1, pattern = ".....$", replacement = "")
#names$batch <- str_sub(names$batch, -1)
# df.adjusted <- ComBat(dat = df, batch = names$batch)  - not clear batch effect from MDS plots
#perform DE
#create design matrix
names <- names[which(names$V1 %in% remove == FALSE),]
conditions <- names$condition[,1]
x <- as.data.frame(conditions)
x$seq <- ""
des <- x %>% 
  mutate(seq = ifelse(conditions == "HD", 1, seq))
des <- des %>% 
  mutate(seq = ifelse(conditions == "control",0, seq))
des <- as.data.frame(cbind(rep(1, 227), des$seq))
rownames(des) <- colnames(df)
colnames(des) <- c("Samples", "Conditions")
des$Samples <- as.numeric(des$Samples)
des$Conditions <- as.numeric(des$Conditions)
des2 <- as.matrix(des)

df <- df[grepl('hsa', rownames(df)),]
df <- as.matrix(df)

fit <- lmFit(df, des2)
fit <- eBayes(fit)
result <- topTable(fit, adjust = "BH", coef = "Conditions", n=Inf)  
write.csv(result, "blood_Dong_DE_limma.csv")

#try DESeq2 due to large disparity in log2fc values seen (high variance and poisson models
# did not fit well)
#create colData object
names <- read.table("GSE108395_names.txt")
names <- as.data.frame(t(names))
x <- str_match(names$V1, "-\\s*(.*?)\\s*-")
x <- gsub(x, pattern = "-", replacement = "")
names <- cbind(names, x[,1])
names <- names[which(names$V1 %in% remove == FALSE),]
rownames(names) <- names$V1
colnames(names) <- c("Samples", "Conditions")
colData <- as.matrix(names)
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
write.csv(DEobj.table, "blood_dong_DE_Deseq2.csv")

