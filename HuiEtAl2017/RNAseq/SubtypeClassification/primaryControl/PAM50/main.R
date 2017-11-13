library("beeswarm")
library("beeswarm")
library("RColorBrewer")
library("gplots")
library("caret")
library("xlsx")
library(leaflet)

prep_data <- function(data,pam50.genes) {

  temp <- row.names(data)
  temp[temp=="NUF2"]<-"CDCA1"
  temp[temp=="NDC80"]<-"KNTC2"
  row.names(data) <- temp
  data <- data[pam50.genes,]

  #mark any normal control samples as normal-like

  sample_classes <- colnames(data)
  sample_classes[grep("-11",sample_classes)] <- "normal"
  sample_classes[sample_classes!="normal"] <- "tumor"

  temp <- unlist(
    lapply(
      colnames(data),
      function(x) return(
          paste(strsplit(x,"-")[[1]][1:3],collapse="-")
        )
      )
  )
  temp[sample_classes=="normal"] <- "normal"
  colnames(data) <- temp

  subtypes <- read.csv("/Users/charlesmurphy/Desktop/Research/data/papers/2012_comprehensiveMolecularPortraitsOfHumanBreastCancer/subtype.csv",header=T,check.names=F,colClasses="character")
  subtypes <- subtypes[!is.na(subtypes[,5]),]

  data <- data[,((colnames(data) %in% subtypes[,1]) | (colnames(data)=="normal"))]

  class_labels <- c()
  for (i in colnames(data)) {
    if (grepl("normal",i)) {
      class_labels <- c(class_labels,"Normal-like")
    } else {
      if (length(as.character(subtypes[subtypes[,1]==i,5]))==0) {
        print(i)
      }
      class_labels <- c(class_labels,as.character(subtypes[subtypes[,1]==i,5]))
    }
  }
  return(list(data,class_labels))
}

class_centroids <- function(data,class_labels,shrinkage=1,classes) {

  overall_centroid <- apply(data,1,mean)
  class_centroids <- matrix(0,ncol=length(classes),nrow=nrow(data))
  class_shrunken_centroids <- matrix(0,ncol=length(classes),nrow=nrow(data))
  colnames(class_centroids) <- classes
  colnames(class_shrunken_centroids) <- classes

  s0 <- apply(data,1,median)
  n <- ncol(data)
  k = length(classes)
  prior_class_probability <- c()

  #compute class centroids

  for (i in classes) {
    data.class <- data[,class_labels==i]
    class_centroids[,i] <- apply(data.class,1,mean)
    prior_class_probability <- c(prior_class_probability,ncol(data.class)/n)
  }

  #compute variance

  si = matrix(0,ncol=1,nrow=nrow(data))
  for (i in classes) {
    data.class <- data[,class_labels==i]
    temp <- apply((data.class - class_centroids[,i])**2,1,sum)
    si = si + temp
  }
  si = (1/(n-k)) * si

  #compute shrinkage

  for (i in classes) {
    data.class <- data[,class_labels==i]

    nk <- ncol(data.class)
    mk <- sqrt((1/nk) + (1/n))
    dk <- (class_centroids[,i] - overall_centroid) / (mk * (si + s0))

    temp <- abs(dk) - shrinkage
    temp[temp<0]<-0

    dk.prime <- sign(dk) * temp

    class_shrunken_centroids[,i] <- overall_centroid + mk * (si + s0) * dk.prime
  }

  return(list(class_shrunken_centroids,si,s0,prior_class_probability))
}

classify <- function(test_vector,class_shrunken_centroids,si,s0,prior_class_probability,classes) {
  discriminant_score <- c()

  for (i in 1:length(prior_class_probability)) {
    discriminant_score <- c(discriminant_score,sum((test_vector - class_shrunken_centroids[,i])**2 / (si + s0)**2) - 2 * log(prior_class_probability[i]))
  }

  class_probability <- exp(-0.5 * discriminant_score) / sum(exp(-0.5*discriminant_score))

  return(classes[class_probability==max(class_probability)])
}

classify_all_prob <- function(test_vector,class_shrunken_centroids,si,s0,prior_class_probability,classes) {
  discriminant_score <- c()

  for (i in 1:length(prior_class_probability)) {
    discriminant_score <- c(discriminant_score,sum((test_vector - class_shrunken_centroids[,i])**2 / (si + s0)**2) - 2 * log(prior_class_probability[i]))
  }

  class_probability <- exp(-0.5 * discriminant_score) / sum(exp(-0.5*discriminant_score))

  return(class_probability)
}


cross_validation <- function(data,class_labels,classes) {
  accuracy <- c()
  folds <- createFolds(1:ncol(data),10)
  for (test_set in folds) {
    model_set <- setdiff(1:ncol(data),test_set)
    r <- class_centroids(data=data[,model_set],class_labels=class_labels[model_set],shrinkage=1,classes)
    class_shrunken_centroids <- r[[1]]
    si <- r[[2]]
    s0 <- r[[3]]
    prior_class_probability <- r[[4]]
    class_predictions <- as.character(unlist(apply(data[,test_set],2,function(x) return(classify(x,class_shrunken_centroids,si,s0,prior_class_probability,classes)))))
    class_true <- class_labels[test_set]

    class_predictions <- factor(class_predictions,levels=classes)
    class_true <- factor(class_true,levels=classes)

    r <- confusionMatrix(class_predictions,class_true)
    accuracy <- c(accuracy,as.numeric(r$overall["Accuracy"]))
  }
  return(accuracy)
}

classify_mouse <- function(human_data,human_class_labels,classes, mouse_data) {

  r <- class_centroids(data=human_data,class_labels=human_class_labels,shrinkage=1,classes)
  class_shrunken_centroids <- r[[1]]
  si <- r[[2]]
  s0 <- r[[3]]
  prior_class_probability <- r[[4]]


  r <- unlist(apply(mouse_data,2,function(x) return(classify(x,class_shrunken_centroids,si,s0,prior_class_probability,classes))))

  r2 <- apply(mouse_data,2,function(x) return(classify_all_prob(x,class_shrunken_centroids,si,s0,prior_class_probability,classes)))

  return(list(r,r2))
}

pam50.genes <- read.table("/Users/charlesmurphy/Desktop/Research/data/papers/PAM50/PAM50/bioclassifier_R/pam50_centroids.txt",sep="\t",row.names=1,header=T)
pam50.genes <- row.names(pam50.genes)
classes <- c("Basal-like","Luminal A","Luminal B","HER2-enriched","Normal-like")

load("/Users/charlesmurphy/Desktop/Research/data/GDAC/BRCA/gdac.broadinstitute.org_BRCA.Merge_rnaseq__illuminahiseq_rnaseq__unc_edu__Level_3__gene_expression__data.Level_3.2016012800.0.0/process_data/vst.RData")
data.vst <- data.vst[,grep("TCGA-E2-A15A",colnames(data.vst),invert=T)]
data.vst <- data.vst[,grep("TCGA-E2-A15K",colnames(data.vst),invert=T)]
r <- prep_data(data.vst,pam50.genes)
data.vst <- r[[1]]
class_labels <- r[[2]]

accuracy.vst <- cross_validation(data.vst,class_labels,classes)


load('/Users/charlesmurphy/Desktop/Research/data/GDAC/BRCA/gdac.broadinstitute.org_BRCA.Merge_rnaseq__illuminahiseq_rnaseq__unc_edu__Level_3__gene_expression__data.Level_3.2016012800.0.0/process_data/rpkm.RData')
rpkm <- rpkm[,grep("TCGA-E2-A15A",colnames(rpkm),invert=T)]
rpkm <- rpkm[,grep("TCGA-E2-A15K",colnames(rpkm),invert=T)]
r <- prep_data(rpkm,pam50.genes)
rpkm <- r[[1]]
class_labels <- r[[2]]

accuracy.rpkm <- cross_validation(rpkm,class_labels,classes)
accuracy.logrpkm <- cross_validation(log2(rpkm+0.1),class_labels,classes)

rpkm.rank <- apply(rpkm,2,function(x) return(rank(x)/max(rank(x))))
accuracy.rank <- cross_validation(rpkm.rank,class_labels,classes)

rpkm.median <- log2(rpkm+0.1) - apply(log2(rpkm+0.1),1,median)
accuracy.log2median <- cross_validation(rpkm.median,class_labels,classes)

#plot accuracies

r <- list()
r[["RPKM"]] <- accuracy.rpkm
r[["log2(RPKM+0.1)"]] <- accuracy.logrpkm
r[["Rank normalized"]] <- accuracy.rank
r[["Median-centered\nlog2(RPKM+0.1)"]] <- accuracy.log2median
r[["VST"]] <- accuracy.vst

pdf("Classifier-performace.pdf")
par(mar=c(8.1,4.1,4.1,2.1))
beeswarm(r,lty=1,ylim=c(0,1),ylab="Accuracy",las=2,cex=0.5,main="Comparison of classification accuracy\nby normalization method")
boxplot(r, add=TRUE,outline=FALSE,names=rep("",length(r)),pars=list(boxwex=0.5),boxlty=0,whisklty=0,staplelty=0,yaxt="n")
dev.off()

#classify mouse

fpkm <- read.csv("/Users/charlesmurphy/Desktop/Research/0914_hui/results/RNAseq/Cufflinks/genes-humanSymbols-fpkm.csv",row.names=1,header=T,check.names=F)
temp <- row.names(fpkm)
temp[temp=="NUF2"]<-"CDCA1"
temp[temp=="NDC80"]<-"KNTC2"
temp[temp=="ORC6"]<-"ORC6L"
row.names(fpkm) <- temp
fpkm <- fpkm[pam50.genes,]

fpkm.rank <- apply(fpkm,2,function(x) return(rank(x)/max(rank(x))))
r <- classify_mouse(rpkm.rank,class_labels,classes,fpkm.rank)
class_predictions <- r[[1]]
all_class_prob <- r[[2]]
row.names(all_class_prob) <- classes

samples <- read.xlsx("/Users/charlesmurphy/Desktop/Research/0914_hui/data/samples.xlsx","samples")
samples <- samples[(samples['byl']==0) & (samples['bkm']==0) & (samples['bmn']==0),]
samples <- samples[samples$sequencing_type=="RNAseq",]
row.names(samples) <- as.character(samples$sample_ID)

common <- intersect(as.character(samples$sample_ID),colnames(all_class_prob))
print(length(common))
all_class_prob <- all_class_prob[,common]
all_class_prob <- all_class_prob[,order(all_class_prob[1,],decreasing=T)]
samples <- samples[colnames(all_class_prob),]

write.csv(t(all_class_prob),"per-sample-probabilities.csv",quote=F)

mouse.markers <- read.csv("/Users/charlesmurphy/Desktop/Research/0914_hui/results/RNAseq/SubtypeClassification/Markers/mouse-marker-status-prediction.csv",header=T,row.names=1,colClasses="character")
mouse.markers <- mouse.markers[colnames(all_class_prob),]
mouse.markers[mouse.markers=="Positive"]<-"black"
mouse.markers[mouse.markers=="Negative"]<-"white"

sample_type = as.character(samples$tissue)
sample_type[sample_type=="primary"]<-"#999999"
sample_type[sample_type=="implant"]<-"#999999"
sample_type[sample_type=="normal breast"]<-"green"

N <- ncol(all_class_prob)

png('Mouse-classification-PAM50.png',width=1400,height=220)
par(mar=c(6,7,2,3),cex=0.9)
pp <- barplot(
  all_class_prob,
  las=2,
  xlab="",
  xlim=c(0,ncol(all_class_prob)+30),
  ylab="Probability of subtype\nassignment",
  names.arg=rep("",ncol(all_class_prob)),
  col=c("red","blue","#00ffff","purple","green"),
  cex.lab=1.5
)
par(xpd=T)
legend(
  "topright",
  c("Basal-like","Luminal A","Luminal B","HER2-enriched","Normal-like"),
  col=c("red","blue","#00ffff","purple","green"),
  pch=15,
  cex=1.4
)
title(xlab="Samples", line=4.5, cex.lab=1.5)

#plot sample type
rect(
  pp-0.5,
  rep(-0.04,N),
  pp+0.5,
  rep(-0.11,N),
  col=sample_type
)

#plot ER/PR/HER2 statuses
rect(
  pp-0.5,
  rep(-0.19,N),
  pp+0.5,
  rep(-0.12,N),
  col=mouse.markers$PR_status
)
rect(
  pp-0.5,
  rep(-0.27,N),
  pp+0.5,
  rep(-0.20,N),
  col=mouse.markers$ER_status
)
rect(
  pp-0.5,
  rep(-0.35,N),
  pp+0.5,
  rep(-0.28,N),
  col=mouse.markers$HER2_status
)

#x_box=127.5 # all
x_box=75.5
y_box=-0.48

text(x_box-6,-0.5,"Sample type",cex=1.4)
text(x_box+1.6,y_box+0.035,"Tumor",cex=1.4,adj=c(0,0.5))
rect(x_box,y_box,x_box+1.2,y_box+0.07,col="#999999")
text(x_box+1.6,y_box-0.12+0.035,"Normal",cex=1.4,adj=c(0,0.5))
rect(x_box,y_box-0.12,x_box+1.2,y_box-0.12+0.07,col="green")

x_box <- x_box + 25

text(x_box-8,-0.5,"PR/ER/HER2 status",cex=1.4)
text(x_box+1.6,y_box+0.035,"Negative",cex=1.4,adj=c(0,0.5))
rect(x_box,y_box,x_box+1.2,y_box+0.07,col="white")
text(x_box+1.6,y_box-0.12+0.035,"Positive",cex=1.4,adj=c(0,0.5))
rect(x_box,y_box-0.12,x_box+1.2,y_box-0.12+0.07,col="black")

text(max(pp)+0.7,-0.075,"Sample type",cex=1,adj=c(0,0.5))
text(max(pp)+0.7,-0.155,"PR status",cex=1,adj=c(0,0.5))
text(max(pp)+0.7,-0.235,"ER status",cex=1,adj=c(0,0.5))
text(max(pp)+0.7,-0.315,"HER2 status",cex=1,adj=c(0,0.5))
dev.off()
