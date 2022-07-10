#functions
#biplot (PCA)
makepca <- function(probj, col, df){
title=paste0("Hippocampus miRNA data")
ggbiplot(probj, labels = rownames(df), groups=col, labels.size=2, 
         var.axes = F, obs.scale = 0) +
  labs(title =  title)+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8), 
        legend.text = element_text(size = 10))+
  theme_classic()+
  theme(plot.margin=unit(c(0,0,0,0), 'cm'))
}


#DE
getDE <- function(dds, factor, numC, denC) {
  res <- results(dds, contrast= c(factor, numC, denC))
  res_B <- suppressMessages(as.data.frame(lfcShrink(dds=dds, 
                                                    contrast=c("Conditions",
                                                               numC,
                                                               denC), 
                                                    res=res,
                                                    type = 'ashr')))
  return(res_B)
}
