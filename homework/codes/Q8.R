### QUESTION 8 ####
#loading libraries
library(AnnotationHub)
library(rtracklayer)
library(BSgenome.Mmusculus.UCSC.mm10)
library(caret)
library(ROCR)
library(pROC)
library(DNAshapeR)

## Data retreive
seqLength = 30
sampleSize = 1000

# Bound (ChIP-seq)
ah = AnnotationHub()
ctcfPeaks = ah[["AH28451"]]
seqlevelsStyle(ctcfPeaks) = "UCSC"
#save sequences to fasta file:bound.fa
getFasta( GR = ctcfPeaks, BSgenome = Mmusculus, width = seqLength, filename = "bound.fa" )

# Unbound (random regions w/o overlapping)
chrName = names(Mmusculus)[1:22]
chrLength = seqlengths(Mmusculus)[1:22]
#create bed file with genome size
d = data.frame(chrName,chrLength)
write(t(d),file="genome.bed",ncolumns = 2,sep="\t")

#create a bed file for ctcf cordinates
start_end = data.frame(ctcfPeaks@ranges)[,1:2]
d=data.frame(ctcfPeaks@seqnames,start_end)
write(t(d),file="ctcf_cordinates.bed",ncolumns = 3,sep="\t")

###### Using BEDTools to generate random 1000 sequences of length 30 nt#########
##sort genome.bed using 
#sort -k 1,1 -k2,2n genome.bed >sorted_genome.bed
##sort ctcf_cordinates.bed using 
#sortBed -i ctcf_cordinates.bed > sorted_ctcf_cordinates.bed
##find complement of ctcf cordinates by
#bedtools complement -i sorted_ctcf_cordinates.bed -g sorted_genome.bed > complement_cordinates.bed
##Get random 1000 cordinates with size 30 using
#bedtools random -g complement_cordinates.bed -n 1000 -l 30 > random_cordinates.bed
##Get fasta file for the random-cordinates using:
##Downlaoded mm10 file from UCSC converted twoBit file to fasta using
#twoBitToFa mm10.2bit mm10.fa
##Extracted the fasta sequences using bedtools getFastaFromBed
#bedtools getfasta -fi mm10.fa -bed random_cordinates.bed -fo unbound.fa
###################################################################################

## Merge bound and unbound data
# Combine two datasets and generate one file
boundFasta = readBStringSet("bound.fa")
boundFasta = sample(boundFasta, sampleSize)
unboundFasta = readBStringSet("unbound.fa")
names(unboundFasta) = paste0( names(unboundFasta), "_unbound")
names(unboundFasta) = paste0( names(unboundFasta), "_unbound")
writeXStringSet( c(boundFasta, unboundFasta), "ctcf.fa" )

## DNAshapeR prediction
shapePred = getShape("ctcf.fa")
## Encode feature vectors
featureType_1 = c("1-mer")
featureType_2 = c("1-mer", "1-shape")
featureVector_1 = encodeSeqShape("ctcf.fa", shapePred, featureType_1)
featureVector_2 = encodeSeqShape("ctcf.fa", shapePred, featureType_2)

## Linear Regression model 
## Data preparation
bound_status = c(rep("Y", sampleSize), rep("N", sampleSize))
df_1 = data.frame(isBound = bound_status, featureVector_1)
df_2 = data.frame(isBound = bound_status, featureVector_2)
# Set parameters for Caret
trainControl = trainControl(method = "cv", number = 2, savePredictions = TRUE,
                             classProbs = TRUE)
# Perform prediction
model_1 = train(isBound~ ., data = df_1, trControl = trainControl,
               method = "glm", family = binomial, metric ="ROC")
# Perform prediction
model_2 = train(isBound~ ., data = df_2, trControl = trainControl,
               method = "glm", family = binomial, metric ="ROC")

model_1
model_2
## Plot AUROC
png(file = "AUROC.png")
prediction1 = prediction( model_1$pred$Y, model_1$pred$obs )
performance1 = performance( prediction1, "tpr", "fpr" )
plot(performance1,col="red")
par(new=TRUE)
prediction2 = prediction( model_2$pred$Y, model_2$pred$obs )
performance2 = performance( prediction2, "tpr", "fpr" )
plot(performance2,col="blue")
abline(0,1,lty=2)
legend(0.4,0.4, lty = c(1,1), cex=0.8, c("1-mer", "1-mer+shape"), col=c("red", "Blue"),bty='n')
title("AUROC")
dev.off()

## Caluculate AUC
auc1 = performance(prediction1, "auc")
auc1 = unlist(slot(auc1, "y.values"))
auc1

auc2 = performance(prediction2, "auc")
auc2 = unlist(slot(auc2, "y.values"))
auc2
