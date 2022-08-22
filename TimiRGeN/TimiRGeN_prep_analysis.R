#perform TimiRGeN analysis on 2M and 10M data per genotype 
#looking for gender based differences per genotype
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#load libs
library(TimiRGeN)
library(org.Mm.eg.db)
library(dplyr)
library(clusterProfiler)
#prep miRNA data
miRNA <- read.csv("miRNA_qval_gender.csv", row.names = 1, header = TRUE)
miRNA.6m <- dplyr::select(miRNA, -contains("6M")) #remove 6 month
miRNA.6m <- miRNA.6m[complete.cases(miRNA.6m), ] # remove NAs
colnames(miRNA.6m) <- sub(colnames(miRNA.6m), pattern = "M.", replacement=".")
colnames(miRNA.6m) <- sub(colnames(miRNA.6m), pattern = "m.", replacement=".")
miRNA.F <- dplyr::select(miRNA.6m, contains("FQ")) #female
colnames(miRNA.F) <- sub(colnames(miRNA.F), pattern = "FQ", replacement="QM")
miRNA.M <- dplyr::select(miRNA.6m, -contains("FQ")) #male
colnames(miRNA.M) <- sub(colnames(miRNA.M), pattern = "^.", replacement="")
colnames(miRNA.M) <- sub(colnames(miRNA.M), pattern = "M", replacement="")
colnames(miRNA.M) <- sub(colnames(miRNA.M), pattern = "^", replacement="QM")
#prep mRNA data
mRNA <- read.csv("mRNA_qval_gender.csv", row.names = 1, header = TRUE)
mRNA.6m <- dplyr::select(mRNA, -contains("6M")) #remove 6 month
mRNA.6m <- mRNA.6m[complete.cases(mRNA.6m), ] # remove NAs
colnames(mRNA.6m) <- sub(colnames(mRNA.6m), pattern = "M.", replacement=".")
colnames(mRNA.6m) <- sub(colnames(mRNA.6m), pattern = "m.", replacement=".")
mRNA.F <- dplyr::select(mRNA.6m, contains("FQ")) #female
colnames(mRNA.F) <- sub(colnames(mRNA.F), pattern = "FQ", replacement="QM")
mRNA.M <- dplyr::select(mRNA.6m, -contains("FQ")) #male
colnames(mRNA.M) <- sub(colnames(mRNA.M), pattern = "^.", replacement="")
colnames(mRNA.M) <- sub(colnames(mRNA.M), pattern = "M", replacement="")
colnames(mRNA.M) <- sub(colnames(mRNA.M), pattern = "^", replacement="QM")
#begin TimiRGeN for female samples
F.MAE <- startObject(miRNA.F, mRNA.F)
F.MAE <- getIdsMir(MAE = F.MAE, assay(F.MAE, 1), orgDB = org.Mm.eg.db, 
                   miRPrefix = "mmu")
F.MAE <- getIdsMrna(MAE = F.MAE, assay(F.MAE, 2), mirror = "useast",
                    species = "mmusculus", orgDB = org.Mm.eg.db)
F.MAE <- combineGenes(MAE = F.MAE, miR_data = assay(F.MAE, 1), assay(F.MAE, 2))
F.MAE <- genesList(MAE = F.MAE, method = "c", genetic_data = assay(F.MAE, 9),
                   timeString = "QM")
F.MAE <- significantVals(MAE = F.MAE, method = "c",
                       geneList = metadata(F.MAE)[[1]], maxVal = 0.05, 
                       stringVal = "padj")
F.MAE <- addIds(MAE = F.MAE, method = "c", filtered_genelist = metadata(F.MAE)[[2]],
                miR_IDs = assay(F.MAE, 3), mRNA_IDs = assay(F.MAE, 7))
F.MAE <- eNames(MAE = F.MAE, method = "c", gene_IDs = metadata(F.MAE)[[3]])
#check which samples will lead to no identified significant genes - and remove
x <- metadata(F.MAE)[[4]]
x <- x[lapply(x,length)>0]
F.MAE2 <- MultiAssayExperiment()
F.MAE2 <- dloadGmt(MAE = F.MAE2, species = "Mus musculus")
# F.MAE2 <- enrichWiki(MAE = F.MAE2, method = "c", 
#                      ID_list = x, orgDB = org.Mm.eg.db, 
#                      path_gene = assay(F.MAE2, 1), path_name = assay(F.MAE2, 2), 
#                      ID = "ENTREZID", universe = assay(F.MAE2, 1)[[2]])
#required additional step due to a bug - will add to tool!
bg <- as.list(F.MAE2[[1]][2])
bg <- unlist(bg)
names(bg) <- NULL
lst2 <- lapply(x, function(x) {
  enricher(x, TERM2GENE = F.MAE2[[1]], TERM2NAME = F.MAE2[[2]], 
           universe = bg)
})
lst2 <- lst2[lapply(lst2,length)>0]
sigwiki <- lapply(lst2, function(x) {
  setReadable(x, org.Mm.eg.db, keyType = "ENTREZID")
})
names(sigwiki) <- paste0(names(sigwiki), "_wikipathways", sep = "")
metadata(F.MAE2)[["sigwiki"]] <- sigwiki
#plot out results
setwd("f/")
savePlots(largeList = metadata(F.MAE2)[[1]], maxInt = 4, fileType = "jpeg")
##################################################################################
#male samples
M.MAE <- startObject(miRNA.M, mRNA.M)
M.MAE <- getIdsMir(MAE = M.MAE, assay(M.MAE, 1), orgDB = org.Mm.eg.db, 
                   miRPrefix = "mmu")
M.MAE <- getIdsMrna(MAE = M.MAE, assay(M.MAE, 2), mirror = "useast",
                    species = "mmusculus", orgDB = org.Mm.eg.db)
M.MAE <- combineGenes(MAE = M.MAE, miR_data = assay(M.MAE, 1), assay(M.MAE, 2))
M.MAE <- genesList(MAE = M.MAE, method = "c", genetic_data = assay(M.MAE, 9),
                   timeString = "QM")
M.MAE <- significantVals(MAE = M.MAE, method = "c",
                         geneList = metadata(M.MAE)[[1]], maxVal = 0.05, 
                         stringVal = "padj")
M.MAE <- addIds(MAE = M.MAE, method = "c", filtered_genelist = metadata(M.MAE)[[2]],
                miR_IDs = assay(M.MAE, 3), mRNA_IDs = assay(M.MAE, 7))
M.MAE <- eNames(MAE = M.MAE, method = "c", gene_IDs = metadata(M.MAE)[[3]])
x <- metadata(M.MAE)[[4]]
x <- x[lapply(x,length)>0]
M.MAE2 <- MultiAssayExperiment()
M.MAE2 <- dloadGmt(MAE = M.MAE2, species = "Mus musculus")
bg <- as.list(M.MAE2[[1]][2])
bg <- unlist(bg)
names(bg) <- NULL
lst2 <- lapply(x, function(x) {
  enricher(x, TERM2GENE = M.MAE2[[1]], TERM2NAME = M.MAE2[[2]], 
           universe = bg)
})
lst2 <- lst2[lapply(lst2,length)>0]
sigwiki <- lapply(lst2, function(x) {
  setReadable(x, org.Mm.eg.db, keyType = "ENTREZID")
})
names(sigwiki) <- paste0(names(sigwiki), "_wikipathways", sep = "")
metadata(M.MAE2)[["sigwiki"]] <- sigwiki
#plot out results
setwd("../m/")
savePlots(largeList = metadata(M.MAE2)[[1]], maxInt = 3, fileType = "jpeg")
