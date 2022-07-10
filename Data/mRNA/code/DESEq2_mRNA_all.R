library(DESeq2)
library(limma)

setwd("~/Documents/HD/Data/mRNA")

Data <- read.csv("org_counts.csv", row.names = 1)
Pheno <- read.csv("mRNA_Phenotype.csv", row.names = 1)
Data <- Data[,-c(63, 96, 127, 145, 157)]
Pheno <- Pheno[-c(63, 96, 127, 145, 157),]

m <- mapply(Data, FUN=as.integer)
rownames(m) <- rownames(Data)

Samples <- colnames(Data)
Conditions <- Pheno$Name
Gender <- Pheno$Sex
colData <- cbind(Samples, Conditions, Gender)
rownames(colData) <- colnames(Data)

dds <- DESeqDataSetFromMatrix(countData = m, colData = colData, 
                              design = ~Conditions)

keep <- rowSums(counts(dds)) >= 50
dds <- dds[keep,]

plotMDS(dds@assays@data@listData$counts, col = as.numeric(dds$Conditions))
boxplot(dds@assays@data@listData$counts, col = as.numeric(dds$Conditions))
dds$Conditions
dds$Conditions <- factor(dds$Conditions, levels = unique(dds$Conditions))

### normalised counts for ML
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="normalized_mRNA_Nooutliers.txt", sep="\t", quote=F, col.names=NA)
###
dds <- DESeq(dds)
resultsNames(dds)

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

F_WT_6M <- getDE(numC = 'female_WT_6m', denC = 'female_WT_2m')
F_WT_10M <- getDE(numC = 'female_WT_10m', denC = 'female_WT_2m')

F_Q20_6M <- getDE(numC = 'female_Q20_6m', denC = 'female_Q20_2m')
F_Q20_10M <- getDE(numC = 'female_Q20_10m', denC = 'female_Q20_2m')

F_Q80_6M <- getDE(numC = 'female_Q80_6m', denC = 'female_Q80_2m')
F_Q80_10M <- getDE(numC = 'female_Q80_10m', denC = 'female_Q80_2m')

F_Q111_6M <- getDE(numC = 'female_Q111_6m', denC = 'female_Q111_2m')
F_Q111_10M <- getDE(numC = 'female_Q111_10m', denC = 'female_Q111_2m')

F_Q140_6M <- getDE(numC = 'female_Q140_6m', denC = 'female_Q140_2m')
F_Q140_10M <- getDE(numC = 'female_Q140_10m', denC = 'female_Q140_2m')

F_Q175_6M <- getDE(numC = 'female_Q175_6m', denC = 'female_Q175_2m')
F_Q175_10M <- getDE(numC = 'female_Q175_10m', denC = 'female_Q175_2m')
#
M_WT_6M <- getDE(numC = 'male_WT_6m', denC = 'male_WT_2m')
M_WT_10M <- getDE(numC = 'male_WT_10m', denC = 'male_WT_2m')

M_Q20_6M <- getDE(numC = 'male_Q20_6m', denC = 'male_Q20_2m')
M_Q20_10M <- getDE(numC = 'male_Q20_10m', denC = 'male_Q20_2m')

M_Q80_6M <- getDE(numC = 'male_Q80_6m', denC = 'male_Q80_2m')
M_Q80_10M <- getDE(numC = 'male_Q80_10m', denC = 'male_Q80_2m')

M_Q111_6M <- getDE(numC = 'male_Q111_6m', denC = 'male_Q111_2m')
M_Q111_10M <- getDE(numC = 'male_Q111_10m', denC = 'male_Q111_2m')

M_Q140_6M <- getDE(numC = 'male_Q140_6m', denC = 'male_Q140_2m')
M_Q140_10M <- getDE(numC = 'male_Q140_10m', denC = 'male_Q140_2m')

M_Q175_6M <- getDE(numC = 'male_Q175_6m', denC = 'male_Q175_2m')
M_Q175_10M <- getDE(numC = 'male_Q175_10m', denC = 'male_Q175_2m')
#
MF_WT_2M <- getDE(numC = 'male_WT_2m', denC = 'female_WT_2m')
MF_WT_6M <- getDE(numC = 'male_WT_6m', denC = 'female_WT_6m')
MF_WT_10M <- getDE(numC = 'male_WT_10m', denC = 'female_WT_10m')

MF_Q20_2M <- getDE(numC = 'male_Q20_2m', denC = 'female_Q20_2m')
MF_Q20_6M <- getDE(numC = 'male_Q20_6m', denC = 'female_Q20_6m')
MF_Q20_10M <- getDE(numC = 'male_Q20_10m', denC = 'female_Q20_10m')

MF_Q80_2M <- getDE(numC = 'male_Q80_2m', denC = 'female_Q80_2m')
MF_Q80_6M <- getDE(numC = 'male_Q80_6m', denC = 'female_Q80_6m')
MF_Q80_10M <- getDE(numC = 'male_Q80_10m', denC = 'female_Q80_10m')

MF_Q92_2M <- getDE(numC = 'male_Q92_2m', denC = 'female_Q92_2m')
MF_Q92_6M <- getDE(numC = 'male_Q92_6m', denC = 'female_Q92_6m')
MF_Q92_10M <- getDE(numC = 'male_Q92_10m', denC = 'female_Q92_10m')

MF_Q111_2M <- getDE(numC = 'male_Q111_2m', denC = 'female_Q111_2m')
MF_Q111_6M <- getDE(numC = 'male_Q111_6m', denC = 'female_Q111_6m')
MF_Q111_10M <- getDE(numC = 'male_Q111_10m', denC = 'female_Q111_10m')

MF_Q140_2M <- getDE(numC = 'male_Q140_2m', denC = 'female_Q140_2m')
MF_Q140_6M <- getDE(numC = 'male_Q140_6m', denC = 'female_Q140_6m')
MF_Q140_10M <- getDE(numC = 'male_Q140_10m', denC = 'female_Q140_10m')

MF_Q175_2M <- getDE(numC = 'male_Q175_2m', denC = 'female_Q175_2m')
MF_Q175_6M <- getDE(numC = 'male_Q175_6m', denC = 'female_Q175_6m')
MF_Q175_10M <- getDE(numC = 'male_Q175_10m', denC = 'female_Q175_10m')

setwd("~/Documents/HD/Data/mRNA/DE/MF")
write.csv(MF_WT_2M, "MF_WT_2M.csv")
write.csv(MF_WT_6M, "MF_WT_6M.csv")
write.csv(MF_WT_10M, "MF_WT_10M.csv")
write.csv(MF_Q20_2M, "MF_Q20_2M.csv")
write.csv(MF_Q20_6M, "MF_Q20_6M.csv")
write.csv(MF_Q20_10M, "MF_Q20_10M.csv")
write.csv(MF_Q80_2M, "MF_Q80_2M.csv")
write.csv(MF_Q80_6M, "MF_Q80_6M.csv")
write.csv(MF_Q80_10M, "MF_Q80_10M.csv")
write.csv(MF_Q92_2M, "MF_Q92_2M.csv")
write.csv(MF_Q92_6M, "MF_Q92_6M.csv")
write.csv(MF_Q92_10M, "MF_Q92_10M.csv")
write.csv(MF_Q111_2M, "MF_Q111_2M.csv")
write.csv(MF_Q111_6M, "MF_Q111_6M.csv")
write.csv(MF_Q111_10M, "MF_Q111_10M.csv")
write.csv(MF_Q140_2M, "MF_Q140_2M.csv")
write.csv(MF_Q140_6M, "MF_Q140_6M.csv")
write.csv(MF_Q140_10M, "MF_Q140_10M.csv")
write.csv(MF_Q175_2M, "MF_Q175_2M.csv")
write.csv(MF_Q175_6M, "MF_Q175_6M.csv")
write.csv(MF_Q175_10M, "MF_Q175_10M.csv")

#
F_Q20_WT_2M <- getDE(numC = 'female_Q20_2m', denC = 'female_WT_2m')
F_Q20_WT_6M <- getDE(numC = 'female_Q20_6m', denC = 'female_WT_6m')
F_Q20_WT_10M <- getDE(numC = 'female_Q20_10m', denC = 'female_WT_10m')

F_Q80_WT_2M <- getDE(numC = 'female_Q80_2m', denC = 'female_WT_2m')
F_Q80_WT_6M <- getDE(numC = 'female_Q80_6m', denC = 'female_WT_6m')
F_Q80_WT_10M <- getDE(numC = 'female_Q80_10m', denC = 'female_WT_10m')

F_Q92_WT_2M <- getDE(numC = 'female_Q92_2m', denC = 'female_WT_2m')
F_Q92_WT_6M <- getDE(numC = 'female_Q92_6m', denC = 'female_WT_6m')
F_Q92_WT_10M <- getDE(numC = 'female_Q92_10m', denC = 'female_WT_10m')

F_Q111_WT_2M <- getDE(numC = 'female_Q111_2m', denC = 'female_WT_2m')
F_Q111_WT_6M <- getDE(numC = 'female_Q111_6m', denC = 'female_WT_6m')
F_Q111_WT_10M <- getDE(numC = 'female_Q111_10m', denC = 'female_WT_10m')

F_Q140_WT_2M <- getDE(numC = 'female_Q140_2m', denC = 'female_WT_2m')
F_Q140_WT_6M <- getDE(numC = 'female_Q140_6m', denC = 'female_WT_6m')
F_Q140_WT_10M <- getDE(numC = 'female_Q140_10m', denC = 'female_WT_10m')

F_Q175_WT_2M <- getDE(numC = 'female_Q175_2m', denC = 'female_WT_2m')
F_Q175_WT_6M <- getDE(numC = 'female_Q175_6m', denC = 'female_WT_6m')
F_Q175_WT_10M <- getDE(numC = 'female_Q175_10m', denC = 'female_WT_10m')

setwd("~/Documents/HD/Data/mRNA/DE/WT/Female")
write.csv(F_Q20_WT_2M, "F_Q20_WT_2M.csv")
write.csv(F_Q20_WT_6M, "F_Q20_WT_6M.csv")
write.csv(F_Q20_WT_10M, "F_Q20_WT_10M.csv")
write.csv(F_Q80_WT_2M, "F_Q80_WT_2M.csv")
write.csv(F_Q80_WT_6M, "F_Q80_WT_6M.csv")
write.csv(F_Q80_WT_10M, "F_Q80_WT_10M.csv")
write.csv(F_Q92_WT_2M, "F_Q92_WT_2M.csv")
write.csv(F_Q92_WT_6M, "F_Q92_WT_6M.csv")
write.csv(F_Q92_WT_10M, "F_Q92_WT_10M.csv")
write.csv(F_Q111_WT_2M, "F_Q111_WT_2M.csv")
write.csv(F_Q111_WT_6M, "F_Q111_WT_6M.csv")
write.csv(F_Q111_WT_10M, "F_Q111_WT_10M.csv")
write.csv(F_Q140_WT_2M, "F_Q140_WT_2M.csv")
write.csv(F_Q140_WT_6M, "F_Q140_WT_6M.csv")
write.csv(F_Q140_WT_10M, "F_Q140_WT_10M.csv")
write.csv(F_Q175_WT_2M, "F_Q175_WT_2M.csv")
write.csv(F_Q175_WT_6M, "F_Q175_WT_6M.csv")
write.csv(F_Q175_WT_10M, "F_Q175_WT_10M.csv")

#
M_Q20_WT_2M <- getDE(numC = 'male_Q20_2m', denC = 'male_WT_2m')
M_Q20_WT_6M <- getDE(numC = 'male_Q20_6m', denC = 'male_WT_6m')
M_Q20_WT_10M <- getDE(numC = 'male_Q20_10m', denC = 'male_WT_10m')

M_Q80_WT_2M <- getDE(numC = 'male_Q80_2m', denC = 'male_WT_2m')
M_Q80_WT_6M <- getDE(numC = 'male_Q80_6m', denC = 'male_WT_6m')
M_Q80_WT_10M <- getDE(numC = 'male_Q80_10m', denC = 'male_WT_10m')

M_Q92_WT_2M <- getDE(numC = 'male_Q92_2m', denC = 'male_WT_2m')
M_Q92_WT_6M <- getDE(numC = 'male_Q92_6m', denC = 'male_WT_6m')
M_Q92_WT_10M <- getDE(numC = 'male_Q92_10m', denC = 'male_WT_10m')

M_Q111_WT_2M <- getDE(numC = 'male_Q111_2m', denC = 'male_WT_2m')
M_Q111_WT_6M <- getDE(numC = 'male_Q111_6m', denC = 'male_WT_6m')
M_Q111_WT_10M <- getDE(numC = 'male_Q111_10m', denC = 'male_WT_10m')

M_Q140_WT_2M <- getDE(numC = 'male_Q140_2m', denC = 'male_WT_2m')
M_Q140_WT_6M <- getDE(numC = 'male_Q140_6m', denC = 'male_WT_6m')
M_Q140_WT_10M <- getDE(numC = 'male_Q140_10m', denC = 'male_WT_10m')

M_Q175_WT_2M <- getDE(numC = 'male_Q175_2m', denC = 'male_WT_2m')
M_Q175_WT_6M <- getDE(numC = 'male_Q175_6m', denC = 'male_WT_6m')
M_Q175_WT_10M <- getDE(numC = 'male_Q175_10m', denC = 'male_WT_10m')


setwd("~/Documents/HD/Data/mRNA/DE/WT/Male")
write.csv(M_Q20_WT_2M, "M_Q20_WT_2M.csv")
write.csv(M_Q20_WT_6M, "M_Q20_WT_6M.csv")
write.csv(M_Q20_WT_10M, "M_Q20_WT_10M.csv")
write.csv(M_Q80_WT_2M, "M_Q80_WT_2M.csv")
write.csv(M_Q80_WT_6M, "M_Q80_WT_6M.csv")
write.csv(M_Q80_WT_10M, "M_Q80_WT_10M.csv")
write.csv(M_Q92_WT_2M, "M_Q92_WT_2M.csv")
write.csv(M_Q92_WT_6M, "M_Q92_WT_6M.csv")
write.csv(M_Q92_WT_10M, "M_Q92_WT_10M.csv")
write.csv(M_Q111_WT_2M, "M_Q111_WT_2M.csv")
write.csv(M_Q111_WT_6M, "M_Q111_WT_6M.csv")
write.csv(M_Q111_WT_10M, "M_Q111_WT_10M.csv")
write.csv(M_Q140_WT_2M, "M_Q140_WT_2M.csv")
write.csv(M_Q140_WT_6M, "M_Q140_WT_6M.csv")
write.csv(M_Q140_WT_10M, "M_Q140_WT_10M.csv")
write.csv(M_Q175_WT_2M, "M_Q175_WT_2M.csv")
write.csv(M_Q175_WT_6M, "M_Q175_WT_6M.csv")
write.csv(M_Q175_WT_10M, "M_Q175_WT_10M.csv")
