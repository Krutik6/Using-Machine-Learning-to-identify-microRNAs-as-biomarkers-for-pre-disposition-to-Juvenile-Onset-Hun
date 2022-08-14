setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source(file = "plothelper.R")

library(dplyr)
library(reshape2)
library(stringr)
library(ggplot2)
library(wesanderson)
library(plotly)
################################################################################
#prep for fig 1A
# df <- read.csv("recursiveelimdata.csv")
# df <- melt(df)
# df$variable <- as.character(df$variable)
# df <- data.frame(df,do.call(rbind,strsplit(df$variable,".",fixed = TRUE)))
# colnames(df) <- c("classifier", "whole", "accuracy", "type", "age", "genes")
# df$pairs <- paste0(df$classifier,"_", df$age, "_", df$genes)
# df$accuracy <- round(df$accuracy, digits = 2)
# df$age <- sub(df$age, pattern = "10", replacement = "10M Samples")
# df$age <- sub(df$age, pattern = "2", replacement = "2M Samples")
# df$age <- sub(df$age, pattern = "AGE", replacement = "Predisposition")
# df$classifier <- sub(df$classifier, pattern = "Gaussian_Process",
#                      replacement = "Gaussian_NB")
# df$classifier <- sub(df$classifier, pattern = "Polynomial_SVM",
#                      replacement = "Poly_SVM")


#################################################################################
#prep for fig 1B
df <- read.csv("kbestresults_mat.csv")
names(df)[1] <- "classifier"
df$Age <- sub(df$Age, pattern = "10", replacement = "10M Samples")
df$Age <- sub(df$Age, pattern = "2", replacement = "2M Samples")
df$Age <- sub(df$Age, pattern = "AGE", replacement = "Predisposition")
df$classifier <- sub(df$classifier, pattern = "RandomForest", replacement = "Random Forest")
df$classifier <- sub(df$classifier, pattern = "ExtraTrees", replacement = "Extra Trees")

miRNA <- df[which(df$Genes == "miRNA"),]
mRNA <- df[which(df$Genes == "mRNA"),]

class(df$precision)
################################################################################

