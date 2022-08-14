#plot scatter diarams which contrast True and Predicted labels
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../Current/Info/besty/")
#import libraries
library(reshape2)
library(ggplot2)
library(plotly)
library(dplyr)
library(tidyr)
library(ggfortify)
library(wesanderson)
library(stringr)
#prepare data for plotting
files <- list.files()
flist <- lapply(files, read.csv)
######################################################################
#retrieve all the gene names for table 1
names <- lapply(flist, function(x) colnames(x))
names <- rapply(names,str_replace,pattern="X",replacement="",
                how='replace')
names <- rapply(names,str_replace,pattern="Predicted.Labels",replacement="",
                how='replace')
names <- rapply(names,str_replace,pattern="True.Labels",replacement="",
                how='replace')
names <- lapply(names,function(x) x[order(x, decreasing = FALSE)])
names <- lapply(names,function(x) {x <- unlist(noquote(x))
return(x)})
listnames <- lapply(names, function(x) paste(x, collapse = ", "))
######################################################################
flist <- lapply(flist, function(x) transform(x, age=""))
flist <- lapply(flist, function(x) transform(x, genes=""))
geneslist <- rep(c("miRNA", "mRNA"), each=3)
agelist <- rep(c("10m", "AGE", "2m"), 2)

for (i in seq_along(geneslist)) {
  flist[[i]]$age <- agelist[i]
  flist[[i]]$genes <- geneslist[i]
}

flist <- lapply(flist, function(x) transform(x, correct=""))
flist <- lapply(flist, function(x){
x <- x %>% mutate(across(True.Labels:Predicted.Labels, ~as.character(.x)), 
              across(correct, ~ ifelse(True.Labels != Predicted.Labels,
                                       "FALSE", as.character(.x))), 
              across(correct, ~ ifelse(True.Labels == Predicted.Labels,
                                       "TRUE", as.character(.x))), 
              across(True.Labels:Predicted.Labels, ~as.factor(.x)))
return(x)}
)
df <- flist[[5]]
minus <- ncol(df)-5
df2 <- df[,2:minus] 
df2 <- as.matrix(df2)
pallete = c('red', 'blue')
pca_res <- prcomp(df2, scale. = FALSE)
p <- autoplot(pca_res, data=df, colour='True.Labels', size=5, shape='correct')
p <- p + 
  theme_classic()+
  scale_colour_manual(values=pallete)+
  labs(title="PCA - mRNA Predisposed Samples")+
  theme(plot.title=element_text(size=18, face="bold",hjust = 0.5, vjust=0.2,
                                margin = margin(t = 0, r = 10, b = 10, l = 0)),
        axis.text.x=element_text(size=14, face="bold"),
        axis.text.y=element_text(size=14, face="bold"),
        axis.title.x=element_text(size=16, face="bold",
                                  margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y=element_text(size=16, face="bold",
                                  margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.text=element_text(size=16, face="bold"),
        legend.title =element_text(size=16, face="bold"),
        strip.text.x = element_text(size=18, colour = "red"),
        strip.text.y = element_text(size=18, colour = "red"))
 ggplotly(p) %>%
  layout(legend = list(orientation = "h", x =0, y =-0.4, title="x"),
         title = list(y=0.95, x=0.08, xref="plot"),
         margin = list(l = 50, t = 80, r= 20, b=50),
         plot_bgcolor='White')


  
