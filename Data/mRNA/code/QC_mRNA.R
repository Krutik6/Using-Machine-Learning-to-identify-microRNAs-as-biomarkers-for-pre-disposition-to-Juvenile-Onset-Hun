#libraries
library(ggfortify)
library(ggplot2)
library(ggbiplot)
library(dplyr)
source("functions_HD.R")
#set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
df <- read.csv("../mRNA_counts.csv", row.names = 1)
pheno <- read.csv("../mRNA_Pheno.csv", row.names = 1)
pheno <- pheno[,c(1,2,3,4,5)]
#remove 6 month data
df <- df %>% select(-contains("6m"))
pheno <- pheno[- grep("6m", pheno$age),]
#remove miRNAs with only 0s
df.no0 <- df[rowSums(df == 0) <= 100 , , drop = FALSE]
#add phenotypic data
df <- as.data.frame(t(df))
df <- cbind(pheno[,2:5], t(df.no0))
rownames(df) <- colnames(df.no0)
rownames(df) <- sub(rownames(df), pattern = ".*?_", replacement = "")
df$Name <- sub(df$Name, pattern = "female", replacement = "F")
df$Name <- sub(df$Name, pattern = "male", replacement = "M")
#perform PCA
pr <- prcomp(x = df[,-c(1,2,3,4)], scale. = TRUE, center = TRUE)
makepca(probj = pr, col = df$age, df = df)
#remove outliers
threshold <- mean(pr$sdev) * 6
pr$sdev < threshold
rot <- pr$rotation
prdat <- pr$x
#systematically remove each value as the potential outlier
# names <- rownames(df)
# allThreshs <- list()
# allBins <- list()
# for (i in seq_along(names)) {
#   k <- names[i]
#   x <- df[which(rownames(df) != k),]
#   pr <- prcomp(x = x[,-c(1,2,3,4)], scale. = TRUE, center = TRUE)
#   threshold <- mean(pr$sdev) * 6
#   allThreshs[[i]] <- threshold
#   names(allThreshs)[i] <- k
#   allBins[[i]] <- pr$sdev < threshold
#   names(allBins)[i] <- k
# }
#remove outlier with highest variance

makepca(probj = pr, col = df$age, df = df)

df <- df[which(rownames(df) != "male_Q92_2m_3R"),]
pr <- prcomp(x = df[,-c(1,2,3,4)], scale. = TRUE, center = TRUE)
threshold <- mean(pr$sdev) * 6
pr$sdev < threshold
prdat <- pr$x
makepca(probj = pr, col = df$age, df = df)

#save file
write.csv(df, "../mRNA_counts_outliersRemoves.csv")
