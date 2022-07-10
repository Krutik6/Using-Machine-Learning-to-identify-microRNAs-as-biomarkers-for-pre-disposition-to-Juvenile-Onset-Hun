setwd("~/Documents/HD/Data/mRNA/")
counts <- read.csv("mRNA_counts.csv", row.names = 1)
colnames(counts) <- gsub(colnames(counts), pattern = '^...........', replacement = '')
head(counts)

library(biomaRt)
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
G_list <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                filters = "ensembl_gene_id", values = rownames(counts), 
                mart = ensembl)

counts$ensembl_gene_id <- rownames(counts)
named_counts <- merge(G_list, counts, 'ensembl_gene_id') 
head(named_counts)

ordered_counts <- named_counts[order(named_counts$external_gene_name),]
no_dup_counts <- ordered_counts[! duplicated(ordered_counts$external_gene_name),]
rownames(no_dup_counts) <- no_dup_counts$external_gene_name
no_dup_counts$ensembl_gene_id <- no_dup_counts$external_gene_name <- NULL

write.csv(no_dup_counts, "org_counts.csv")
