setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

miRNAnames <- read.csv("predisp.txt")
miRNAs <- as.data.frame(colnames(miRNAnames))
#look into MGI for homologues
miRNAs$humannames <- ""
colnames(miRNAs) <- c("mousenames", "humannames")
miRNAs$mousenames <- gsub(miRNAs$mousenames, pattern = ".", replacement = "-", 
                          fixed = TRUE)
miRNAs$mousenames[41]
miRNAs$humannames <- c("hsa-let-7a-2-3p", "hsa-let-7i-3p", "hsa-let-7i-5p",
                       "NA", "hsa-miR-1247-5p", "hsa-miR-130b-5p",
                       "hsa-134-5p", "hsa-miR-135b-5p", "hsa-miR-153-5p",
                       "hsa-miR-154-5p", "hsa-miR-181a-1-3p", "hsa-miR-181a-2-3p",
                       "hsa-miR-181a-5p", "NA", "hsa-miR-212-3p", "hsa-miR-212-5p",
                       "hsa-miR-221-3p", "hsa-miR-223-3p", "NA", 
                       "NA", "NA", "NA", "NA", "NA", "NA", "hsa-miR-370-5p",
                       "hsa-miR-378b", "hsa-miR-382-5p", "hsa-miR-410-3p",
                       "hsa-miR-433-5p", "hsa-miR-488-5p", "NA", "NA",
                       "hsa-miR-543-3p", "NA", "NA", "hsa-miR-665-3p", "NA",
                       "NA", "hsa-miR-92b-5p", "hsa-miR-99b-5p")

blood_romano <- read.csv("../Romano/blood_Romano_DE.csv")
blood_dong_limma <- read.csv("../Dong_Blood/blood_Dong_DE_limma.csv")
blood_dong_deseq2 <- read.csv("../Dong_Blood/blood_dong_DE_Deseq2.csv")
csf_dong <- read.csv("../Dong_CSF/CSF_dong_DE.csv")



blood_romano.int <- blood_romano[which(blood_romano$X %in% miRNAs$humannames == TRUE),]
csf_dong.int <- csf_dong[which(csf_dong$X %in% miRNAs$humannames == TRUE),]
blood_dong_limma.int <- blood_dong_limma[which(blood_dong_limma$X %in% miRNAs$humannames == TRUE),]
blood_dong_deseq2.int <- blood_dong_deseq2[which(blood_dong_deseq2$X %in% miRNAs$humannames == TRUE),]
#DESeq2 for blood_dong much more intuitive 