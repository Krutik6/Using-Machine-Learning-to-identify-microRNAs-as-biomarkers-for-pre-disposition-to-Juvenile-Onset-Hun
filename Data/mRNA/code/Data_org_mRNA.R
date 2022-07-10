setwd("~/Documents/HD/Data/mRNA/counts")

SRA <- read.table("../SraRunTable.txt", sep = '\t', header = TRUE)
SRA_pheno <- SRA[,c(10, 13:15)]

SRA_pheno$genotype <- sub(SRA_pheno$genotype, pattern = "Wild Type", 
                             replacement = "WT")

SRA_pheno$genotype <- sub(SRA_pheno$genotype, pattern = "Het (",
                             replacement = "", fixed = TRUE)

SRA_pheno$genotype <- sub(SRA_pheno$genotype, pattern = ") Knock-In",
                             replacement = "", fixed = TRUE)

SRA_pheno$age <- sub(SRA_pheno$age, pattern = " month",
                        replacement = "m", fixed = TRUE)

SRA_pheno$Name <- paste0(SRA_pheno$Sex, '_', SRA_pheno$genotype, 
                            '_', SRA_pheno$age)

Col <- stats::ave(as.character(SRA_pheno$Name), SRA_pheno$Name,
                  FUN=function(x) if (length(x)>1) paste0(x[1],
                                                          '_',
                                                          seq_along(x),
                                                          'R') else x[1])
SRA_pheno$Repeated <- Col

SRA_pheno$filename <- paste0(SRA_pheno$Run, '_', SRA_pheno$Repeated)

Files <- list.files()

SRA_pheno <- cbind(SRA_pheno, Files)

# for (i in seq_along(SRA_pheno$Files)) {
#     file.rename(from = SRA_pheno$Files[i], to = SRA_pheno$filename[i])
# }

library(tximport)
library(tximeta)
files <- file.path(getwd(), SRA_pheno$Files, "quant.sf")
names(files) <- paste0(SRA_pheno$Files)
files
# need to summarise to Gene level
coldata <- data.frame(files, names = "SRR",condition = SRA_pheno$Name, stringsAsFactors=FALSE)
se <- tximeta(coldata, type="salmon")
# create tx2gene
a <- se@rowRanges$tx_id
b <- se@rowRanges$gene_id
tx2gene <- as.data.frame(cbind(a, b))
# tximport
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
dim(txi.salmon$counts)
mRNA <- txi.salmon$counts
write.csv(mRNA, "../mRNA_counts.csv")
write.csv(SRA_pheno, "../Phenotype.csv")
