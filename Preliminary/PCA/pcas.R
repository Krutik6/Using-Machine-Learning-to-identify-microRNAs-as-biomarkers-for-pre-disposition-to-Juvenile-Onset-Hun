#Create PCA plots for all data to show homogeneiry
#begin with all data (2m, 6m 10m) * miRNA + mRNA
library(pca3d)
library(dplyr)
library(ggplot2)
library(plotly)
library(stats)
data(iris)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# ALL
#mRNA - all
#df <- read.csv("../DE/Qval/Gender/normalized_mRNA_counts.txt", sep = "\t", 
#                  row.names = 1)
#miRNA - all
#df <- read.csv("../DE/Qval/Gender/normalized_miRNA_counts.txt", sep = "\t", 
#               row.names = 1)
df.pca <- prcomp(t(df), scale. = FALSE)
prop <- summary(df.pca)
dat <- prop$importance[,1:2]
dat1 <- round(dat[2,1], 2) * 100
dat2 <- round(dat[2,2], 2) * 100
names <- colnames(df)
Gender <- sub(names, pattern = "_.*$", replacement = "")
names <- sub(names, pattern = "female_", replacement = "")
names <- sub(names, pattern = "male_", replacement = "")
names <- sub(names, pattern = ".{3}$", replacement = "")
names <- sub(names, pattern = "Q20", replacement = "WT")
names <- sub(names, pattern = "Q80", replacement = "JOHD")
names <- sub(names, pattern = "Q92", replacement = "JOHD")
names <- sub(names, pattern = "Q111", replacement = "JOHD")
names <- sub(names, pattern = "Q140", replacement = "JOHD")
names <- sub(names, pattern = "Q175", replacement = "JOHD")
Genotype <- sub(names, pattern = "\\_.*", replacement = "")
Age <- sub(names, pattern = ".*_", replacement = "")
x <- as.data.frame(df.pca$x)

p <- ggplot(data = x, aes(x = PC1, y = PC2, colour=Age,
                          shape=Genotype, fill=Gender))+
  geom_point(show.legend = TRUE, size=4.5, alpha=0.6)+
  scale_colour_manual(values=c("Red", "Blue", "Orange"))+
  scale_fill_manual(values = c("Black", "Light Gray"))+
  theme_classic()+
  labs(title="PCA - All mRNAs",
       x=paste("PC1 (",dat1,"%)", sep = ""),
       y=paste("PC2 (",dat2,"%)", sep = ""),
       )+
  theme(plot.title=element_text(size=18, face="bold",hjust = 0.7, vjust=0.2,
                                margin = margin(t = 0, r = 10, b = 10, l = 0)),
        axis.text.x=element_text(size=14, face="bold"),
        axis.text.y=element_text(size=14, face="bold"),
        axis.title.x=element_text(size=16, face="bold",
                                  margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y=element_text(size=16, face="bold",
                                  margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.text=element_text(size=14, face="bold"),
        legend.title =element_text(size=14, face="bold"))
ggplotly(p) %>%
  layout(legend = list(orientation = "h", x =-0.1, y =-0.2),
         title = list(y=0.98, x=0.48),
         margin = list(l = 50, t = 50, r= 50, b=100),
         plot_bgcolor='white')

# For ML
# young miRNA
#df <- read.csv("../../Data/miRNA/2m_miRNA.csv", row.names = 1)
# young mRNA
#df <- read.csv("../../Data/mRNA/2m_mRNA.csv", row.names = 1)
# aged miRNA
df <- read.csv("../../Data/miRNA/10m_miRNA_outlier.csv", row.names = 1)
# aged mRNA
#df <- read.csv("../../Data/mRNA/10m_mRNA_outlier.csv", row.names = 1)
df$Samples <- NULL
df.pca <- prcomp(df, scale. = FALSE)
prop <- summary(df.pca)
dat <- prop$importance[,1:2]
dat1 <- round(dat[2,1], 2) * 100
dat2 <- round(dat[2,2], 2) * 100

names <- rownames(df)
Gender <- sub(names, pattern = "_.*$", replacement = "")
names <- sub(names, pattern = "female_", replacement = "")
names <- sub(names, pattern = "male_", replacement = "")
names <- sub(names, pattern = ".{3}$", replacement = "")
names <- sub(names, pattern = "Q20", replacement = "WT")
names <- sub(names, pattern = "Q80", replacement = "JOHD")
names <- sub(names, pattern = "Q92", replacement = "JOHD")
names <- sub(names, pattern = "Q111", replacement = "JOHD")
names <- sub(names, pattern = "Q140", replacement = "JOHD")
names <- sub(names, pattern = "Q175", replacement = "JOHD")
Genotype <- sub(names, pattern = "\\_.*", replacement = "")
Age <- sub(names, pattern = ".*_", replacement = "")
x <- as.data.frame(df.pca$x)

p <- ggplot(data = x, aes(x = PC1, y = PC2, colour=Age,
                          shape=Genotype, fill=Gender))+
  geom_point(show.legend = TRUE, size=4.5, alpha=0.6)+
  scale_colour_manual(values=c("Red"))+
  scale_fill_manual(values = c("Black", "Light Gray"))+
  theme_classic()+
  labs(title="PCA - 10m miRNA Samples after Processing",
       x=paste("PC1 (",dat1,"%)", sep = ""),
       y=paste("PC2 (",dat2,"%)", sep = ""),
  )+
  theme(plot.title=element_text(size=18, face="bold",hjust = 0.7, vjust=0.2,
                                margin = margin(t = 0, r = 10, b = 10, l = 0)),
        axis.text.x=element_text(size=14, face="bold"),
        axis.text.y=element_text(size=14, face="bold"),
        axis.title.x=element_text(size=16, face="bold",
                                  margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y=element_text(size=16, face="bold",
                                  margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.text=element_text(size=14, face="bold"),
        legend.title =element_text(size=14, face="bold"))
ggplotly(p) %>%
  layout(legend = list(orientation = "h", x =-0.1, y =-0.3),
         title = list(y=0.98, x=0.44),
         margin = list(l = 50, t = 50, r= 50, b=100),
         plot_bgcolor='white')
