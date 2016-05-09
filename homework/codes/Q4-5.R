#loading the libraries
library(DNAshapeR)
library(caret)

fasta_mad ="/home/master/hw-data/Mad.txt.fa"
fasta_myc ="/home/master/hw-data/Myc.txt.fa"
fasta_max ="/home/master/hw-data/Max.txt.fa"

pred_mad=getShape(fasta_mad)
pred_myc=getShape(fasta_myc)
pred_max=getShape(fasta_max)

## Encode feature vectors
featureType_1 = c("1-mer")
featureType_2 = c("1-mer", "1-shape")

featureVector_mad = encodeSeqShape(fasta_mad, pred_mad, featureType_1)
featureVector_myc = encodeSeqShape(fasta_myc, pred_myc, featureType_1)
featureVector_max = encodeSeqShape(fasta_max, pred_max, featureType_1)

featureVector_mad1 = encodeSeqShape(fasta_mad, pred_mad, featureType_2)
featureVector_myc1 = encodeSeqShape(fasta_myc, pred_myc, featureType_2)
featureVector_max1 = encodeSeqShape(fasta_max, pred_max, featureType_2)

## Build MLR model by using Caret
# Data preparation

exp_mad = "/home/master/hw-data/Mad.txt"
exp_mad_data = read.table(exp_mad)
exp_myc = "/home/master/hw-data/Myc.txt"
exp_myc_data = read.table(exp_myc)
exp_max = "/home/master/hw-data/Max.txt"
exp_max_data = read.table(exp_max)

df_mad = data.frame(affinity=exp_mad_data$V2, featureVector_mad)
df_mad1 = data.frame(affinity=exp_mad_data$V2, featureVector_mad1)

df_myc = data.frame(affinity=exp_myc_data$V2, featureVector_myc)
df_myc1 = data.frame(affinity=exp_myc_data$V2, featureVector_myc1)

df_max = data.frame(affinity=exp_max_data$V2, featureVector_max)
df_max1 = data.frame(affinity=exp_max_data$V2, featureVector_max1)

# Arguments setting for Caret
trainControl = trainControl(method = "cv", number = 10, savePredictions = TRUE)

# Prediction with L2-regularized
model_mad = train(affinity~., data = df_mad, trControl=trainControl, 
                method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
model_mad1 = train(affinity~., data = df_mad1, trControl=trainControl, 
                   method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
model_myc = train(affinity~., data = df_myc, trControl=trainControl, 
                   method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
model_myc1 = train(affinity~., data = df_myc1, trControl=trainControl, 
                   method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
model_max = train(affinity~., data = df_max, trControl=trainControl, 
                   method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
model_max1 = train(affinity~., data = df_max1, trControl=trainControl, 
                 method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))

#Coefficient of determination
model_mad$results
model_mad1$results
#rs_mad = mean(model_mad$results$Rsquared,na.rm=T)
#rs_mad1 = mean(model_mad1$results$Rsquared,na.rm=T)
rs_mad = max(model_mad$results$Rsquared,na.rm=T)
rs_mad1 = max(model_mad1$results$Rsquared,na.rm=T)

model_myc$results
model_myc1$results
rs_myc = max(model_myc$results$Rsquared,na.rm=T)
rs_myc1 = max(model_myc1$results$Rsquared,na.rm=T)

model_max$results
model_max1$results
rs_max = max(model_max$results$Rsquared,na.rm=T)
rs_max1 = max(model_max1$results$Rsquared,na.rm=T)

print("For mad 1-mer")
print(rs_mad)
print("For mad 1-mer+shape")
print(rs_mad1)
print("For myc 1-mer")
print(rs_myc)
print("For myc 1-mer+shape")
print(rs_myc1)
print("For max 1-mer")
print(rs_max)
print("For max 1-mer+shape")
print(rs_max1)

###Plotting the R^2 values
x = c(rs_mad,rs_myc,rs_max)
y = c(rs_mad1,rs_myc1,rs_max1)

#plot(x,y,main="Comparison of models",xlab="R^2 sequence model",ylab="R^2 sequence + shape model")#,xlim=c(0.70,0.74), ylim=c(0.70,0.74))
png(file="R^2 sequence model.png")
plot(x[1],y[1],col = "red", pch = 15, xlim=c(0.70,0.90), ylim=c(0.70,0.90),
     main="Comparison of models",xlab="R^2 sequence model",ylab="R^2 sequence + shape model")
points(x[2],y[2], pch= 16, col="blue")
points(x[3],y[3],pch =17, col="green")
#abline(lm(y~x), col="black",lty=2)
abline(0,1,lty=2)
legend(0.873,0.875,pch=c(15,16,17),c("Mad","Myc","Max"),col=c('red','blue','green'),bty='n')
dev.off()

### Question 5 b ####
#Mann–Whitney–Wilcoxon (MWW) or Wilcoxon rank-sum test
seq = c(rs_mad,rs_myc,rs_max)
seq_shape = c(rs_mad1,rs_myc1,rs_max1)

test_data = data.frame(seq,seq_shape)
test_data1 = t(test_data)
test_data1
#Null hypothesis the sequence and sequenc + structure model follow same distribution
wilcox.test(seq,seq_shape)
