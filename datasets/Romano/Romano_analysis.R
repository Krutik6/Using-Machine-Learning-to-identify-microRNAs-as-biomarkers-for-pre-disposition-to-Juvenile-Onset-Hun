setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(limma)
df <- read.table("GSE167630_table.txt", fill=TRUE, row.names = 1,
                 header = TRUE)
newcolnames <- c("HD1", "HD2", "HD3", "HD4",
                 "HD5", "HD6", "HD7", "HD8",
                 "HD9", "C1", "C2", "C3",
                 "C4", "C5", "C6", "C7", "C8",
                 "S1", "S2", "S3", "S4", "S5")

comparenames <- cbind(colnames(df), newcolnames)
colnames(df) <- newcolnames
#keep controls and HD
df <- df[,1:17]
#keep hsa samples
df <- df[grepl('hsa', rownames(df)),]
df <- df[- grep("hp_", rownames(df)),]
df <- df[- grep("v11_", rownames(df)),]
rownames(df) <- sub(rownames(df), pattern = "_st", replacement = "")
#check data has been normalised
boxplot(df) #looks normalised
plotMDS(df) # no outliers
#perform DE
conditions <- c(rep(0, 9), rep(1, 8))
samples <- rep(1, 17)
des <- cbind(samples, conditions)
rownames(des) <- colnames(df)
des <- as.matrix(des)
fit <- lmFit(df, des)
fit <- eBayes(fit)
result <- topTable(fit, adjust = "BH", coef = "conditions", n=Inf)  
write.csv(result, "blood_Romano_DE.csv")
