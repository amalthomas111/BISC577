#loading the libraries
library(DNAshapeR)
library(caret)

fasta_mad="/home/master/hw-data/MAD.txt.fa"
fasta_myc="/home/master/hw-data/Myc.txt.fa"
fasta_max="/home/master/hw-data/Max.txt.fa"

phred_mad=getShape(fasta_mad)
phred_myc=getShape(fasta_myc)
phred_max=getShape(fasta_max)

featureType_1 = c("1-mer")
featureType_2 = c("1-mer", "1-shape")