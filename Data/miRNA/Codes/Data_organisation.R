setwd("~/Documents/HD/Data/miRNA")

X <- read.table("SraRunTable_miRNAs.txt", sep = ',', header = TRUE)

useful_names <- X[,c(1, 2, 15, 26)]

useful_names$Genotype <- sub(useful_names$Genotype, pattern = "Wild Type", 
                             replacement = "WT")

useful_names$Genotype <- sub(useful_names$Genotype, pattern = "Het (",
                             replacement = "", fixed = TRUE)

useful_names$Genotype <- sub(useful_names$Genotype, pattern = ") Knock-In",
                             replacement = "", fixed = TRUE)

useful_names$AGE <- sub(useful_names$AGE, pattern = " month",
                             replacement = "m", fixed = TRUE)

useful_names$Name <- paste0(useful_names$sex, '_', useful_names$Genotype, 
                            '_', useful_names$AGE)

write.csv(useful_names, "Pheno.csv")

Col <- stats::ave(as.character(useful_names$Name), useful_names$Name,
                      FUN=function(x) if (length(x)>1) paste0(x[1],
                                                              '_',
                                                              seq_along(x),
                                                              'R') else x[1])
useful_names$Repeated <- Col

useful_names$filename <- paste0(useful_names$Run, '_', useful_names$Repeated, '.csv')

setwd("Counts/")

Files <- list.files()

useful_names <- cbind(useful_names, Files)

# for (i in seq_along(useful_names$Files)) {
#     file.rename(from = useful_names$Files[i], to = useful_names$filename[i])
# }

################################################################################
getwd()
Files <- list.files()
myfiles = lapply(Files, read.delim)
names(myfiles) <- useful_names$Repeated
counts <- lapply(myfiles, function(x) x[c(2)])
X <- do.call(cbind,counts)
colnames(X) <- names(counts)
head(X)
Genes <- myfiles[[1]][1]
Data <- cbind(Genes, X)
rownames(Data) <- Data$X.miRNA
Ordered_data <- Data[order(Data$X.miRNA),]

colnames(Data)[2]
O_data <- Data[with(Data, order(X.miRNA, -female_Q111_2m_1R)), ]
D_data <- O_data[! duplicated(O_data$X.miRNA),]
rownames(D_data) <- D_data$X.miRNA
D_data$X.miRNA <- NULL
write.csv(D_data, "miRNA_counts.csv")
################################################################################
