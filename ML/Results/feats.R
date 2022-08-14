setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
x <- read.csv("bestFeatures.csv")

a <- paste0(x[,3], sep="")
b <- noquote(a[1:80])
c <- b[order(b)]
c
