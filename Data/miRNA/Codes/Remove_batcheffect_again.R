library(DESeq2)
library(limma)
setwd("~/Documents/HD/Data/miRNA")
D_data <- read.csv("~/Documents/HD/Data/miRNA/miRNA_counts.csv", row.names = 1)
D_data <- D_data[,-c(63, 127, 157)]
colnames(D_data)
Pheno <- read.csv("Pheno.csv", row.names = 1)
Pheno <- Pheno[-c(63, 127, 157),]
#Create colData
Samples <- colnames(D_data)
Conditions <- sub(Pheno$Name, pattern = "female_", replacement = "")
Conditions <- sub(Conditions, pattern = "male_", replacement = "")
Repeats <- rep(c(1,2,3,4), 42)
Repeats <- Repeats[-c(63, 127, 157)]
Batch <- c(rep(c(rep(c("A","B"), each = 8), rep("A", 8)), 6), 
           rep("A", 24))
Batch <- Batch[-c(63, 127, 157)]
colData <- cbind(Conditions, Repeats, Batch)
rownames(colData) <- colnames(D_data)
# Make DE object
dds <- DESeqDataSetFromMatrix(countData = D_data, colData = colData, 
                              design = ~Conditions + Conditions:Batch + Conditions:Repeats)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
# quality control
plotMDS(dds@assays@data@listData$counts, col = as.numeric(dds$Conditions))
boxplot(dds@assays@data@listData$counts, col = as.numeric(dds$Conditions))
#Add batch 1

dds$Batch <- Batch
# Remove Batch1

# set condition
dds$Condition <- factor(dds$Condition, levels = unique(dds$Condition))
dds$Condition
dds <- DESeq(dds)
resultsNames(dds)

# function
getDE <- function(numC, denC){
    res <- results(dds, contrast= c("Conditions", numC, denC))
    res_B <- suppressMessages(as.data.frame(lfcShrink(dds=dds, 
                                                      contrast=c("Condition",
                                                                 numC,
                                                                 denC), 
                                                      res=res,
                                                      type = 'ashr')))
    return(res_B)
}

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