#libraries
library(ggfortify)
library(ggplot2)
library(ggbiplot)
library(dplyr)
source("functions_HD.R")
#set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
df <- read.csv("../miRNA_counts.csv", row.names = 1)
pheno <- read.csv("../miRNA_Pheno.csv")
#remove 6 month data
df <- df %>% select(-contains("6m"))
pheno <- pheno[- grep("6m", pheno$AGE),]
#remove miRNAs with only 0s
df.no0 <- df[rowSums(df == 0) <= 100 , , drop = FALSE]
#add phenotypic data
df <- as.data.frame(t(df))
df <- cbind(pheno[,3-5], t(df.no0))
rownames(df) <- colnames(df.no0)
df$X <- rownames(df)
df$X <- sub(df$X, pattern = "female", replacement = "F")
df$X <- sub(df$X, pattern = "male", replacement = "M")
#perform PCA
pr <- prcomp(x = df[,-c(1,2,3,4,5)], scale. = TRUE, center = TRUE)
makepca(probj = pr, col = df$AGE, df = df)
#makepca(probj = pr, col = df$Genotype, df = df)
#makepca(probj = pr, col = df$sex, df = df)
#remove outliers
threshold <- mean(pr$sdev) * 6
pr$sdev < threshold
prdat <- pr$x
#save file
write.csv(df, "../miRNA_counts_outliersRemoves.csv")
