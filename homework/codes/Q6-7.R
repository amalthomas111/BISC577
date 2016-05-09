#### Question 6 ####
#source("https://bioconductor.org/biocLite.R")
#biocLite()
# AnnotationHub
#biocLite("AnnotationHub")
# rtracklayer (might be required)
#biocLite("rtracklayer")
# Reference genome
#biocLite("BSgenome.Mmusculus.UCSC.mm10")
#install.packages(c('ROCR','e1071'))
#library(ROCR)
#library(e1071
library(AnnotationHub)
library(rtracklayer)
library(BSgenome.Mmusculus.UCSC.mm10)
library(DNAshapeR)

## Information retreival
ah = AnnotationHub()
## Genome sequence retreival
chipseqdata = ah[['AH28451']]
seqlevelsStyle(chipseqdata) <- "UCSC"
#getFasta(chipseqdata, BSgenome = Mmusculus, width = 150, filename = "AH28451.fa")
#getFasta(chipseqdata, BSgenome = Mmusculus, width = 300, filename = "AH28451.fa")
getFasta(chipseqdata, BSgenome = Mmusculus, width = 400, filename = "AH28451.fa")
chipseqdata

#### Question 7 ####
#plots
pred = getShape(filename = "AH28451.fa", shapeType = 'All', parse = TRUE)
# Generate ensemble plots
png(file="MGW.png")
plotShape(pred$MGW)
title("Minor Groove Width")
dev.off()
png(file="PropellorTwist.png")
plotShape(pred$ProT)
title("Propellor Twist")
dev.off()
png(file="Roll.png")
plotShape(pred$Roll)
title("Roll")
dev.off()
png(file="HelixTwist.png")
plotShape(pred$HelT)
title("Helix Twist")
dev.off()
