setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#use biomaRt to extract gene names
library(biomaRt)
#start biomart
df <- read.csv("../mRNA_counts_outliersRemoves.csv", row.names = 1)
pheno <- df[,c(1:4)]
df <- df[,-c(1:4)]
df <- as.data.frame(t(df))
df$names <- rownames(df)
#create mart and extract data
mart <- biomaRt::useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", 
                    host = "www.ensembl.org")
a <- listAttributes(mart)
f <- listFilters(mart)
x <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"), 
           filters = "ensembl_gene_id", mart = mart, values = df$names, 
           uniqueRows = TRUE)
x <- x[!duplicated(x$mgi_symbol),]
x <- x[!duplicated(x$ensembl_gene_id),]
#make sure correct and non duplicated data gets attached in the right place
m <- merge(x, df, by.x = "ensembl_gene_id", by.y = "names", all = TRUE)
m[is.na(m)] <- ""
m$mgi_symbol <- ifelse(m$mgi_symbol == "", m$ensembl_gene_id, m$mgi_symbol)
rownames(m) <- m$mgi_symbol
m$ensembl_gene_id <- m$mgi_symbol <- NULL
#write out data
df <- as.data.frame(t(m))
df <- cbind(pheno, df)
write.csv(df, "../mRNA_counts_geneNames.csv", quote = FALSE)
