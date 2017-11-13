library("xlsx")
library("limma")

load_data <- function(data) {

  subtypes <- read.csv("/Users/charlesmurphy/Desktop/Research/data/GDAC/BRCA/gdac.broadinstitute.org_BRCA.Merge_Clinical.Level_1.2016012800.0.0/markers.csv",header=T,check.names=F,colClasses="character")
  colnames(subtypes) <- c("samples","ER_status","PR_status","HER2_status")
  row.names(subtypes) <- subtypes[,"samples"]

  common <- intersect(row.names(subtypes),colnames(data))
  subtypes <- subtypes[common,]
  data <- data[,common]

  return(list(data,subtypes))
}

load('/Users/charlesmurphy/Desktop/Research/data/GDAC/BRCA/gdac.broadinstitute.org_BRCA.Merge_rnaseq__illuminahiseq_rnaseq__unc_edu__Level_3__gene_expression__data.Level_3.2016012800.0.0/process_data/rpkm.RData')
r <- load_data(rpkm)
expression <- r[[1]]
subtypes <- r[[2]]

fpkm <- read.csv("/Users/charlesmurphy/Desktop/Research/0914_hui/results/RNAseq/Cufflinks/genes-humanSymbols-fpkm.csv",row.names=1,header=T,check.names=F)
samples <- read.xlsx("/Users/charlesmurphy/Desktop/Research/0914_hui/data/samples.xlsx","samples")
samples <- samples[(samples['byl']==0) & (samples['bkm']==0) & (samples['bmn']==0),]
samples <- samples[samples$sequencing_type=="RNAseq",]
row.names(samples) <- as.character(samples$sample_ID)
common <- intersect(row.names(samples),colnames(fpkm))
fpkm <- fpkm[,common]
samples <- samples[common,]

common <- intersect(row.names(expression),row.names(fpkm))
fpkm <- fpkm[common,]
expression <- expression[common,]

all_data <- cbind(fpkm,expression)
all_data <- normalizeQuantiles(all_data)
fpkm <- all_data[,1:ncol(fpkm)]
expression <- all_data[,(ncol(fpkm)+1):ncol(all_data)]

#PR model

data <- data.frame(
  PR=as.numeric(expression["PGR",]),
  PR_status=subtypes[,"PR_status"]
)
data <- data[(data[,"PR_status"]=="positive") | (data[,"PR_status"]=="negative"),]
data$PR_status <- factor(data$PR_status,levels=c("negative","positive"))
#data$PR_status <- model.matrix(~PR_status,data=data)[,"PR_statuspositive"]
model.PR <- glm(PR_status~PR,data=data,family=binomial(link='logit'))

#ER model

data <- data.frame(
  ER=as.numeric(expression["ESR1",]),
  ER_status=subtypes[,"ER_status"]
)
data <- data[(data[,"ER_status"]=="positive") | (data[,"ER_status"]=="negative"),]
data$ER_status <- factor(data$ER_status,levels=c("negative","positive"))
model.ER <- glm(ER_status~ER,data=data,family=binomial(link='logit'))

# HER2 model

data <- data.frame(
  HER2=as.numeric(expression["ERBB2",]),
  HER2_status=subtypes[,"HER2_status"]
)
data <- data[!is.na(data$HER2_status),]
data <- data[(data[,"HER2_status"]=="positive") | (data[,"HER2_status"]=="negative"),]
data$HER2_status <- factor(data$HER2_status,levels=c("negative","positive"))
model.HER2 <- glm(HER2_status~HER2,data=data,family=binomial(link='logit'))


data.mouse <- data.frame(
  HER2=as.numeric(fpkm["ERBB2",]),
  PR=as.numeric(fpkm["PGR",]),
  ER=as.numeric(fpkm["ESR1",])
)

tmp <- predict(model.HER2,data.mouse,type='response')
predict.HER2 <- rep("",length(tmp))
predict.HER2[tmp<=0.5] <- "Negative"
predict.HER2[tmp>0.5] <- "Positive"

tmp <- predict(model.PR,data.mouse,type='response')
predict.PR <- rep("",length(tmp))
predict.PR[tmp<=0.5] <- "Negative"
predict.PR[tmp>0.5] <- "Positive"

tmp <- predict(model.ER,data.mouse,type='response')
predict.ER <- rep("",length(tmp))
predict.ER[tmp<=0.5] <- "Negative"
predict.ER[tmp>0.5] <- "Positive"

mouse.markers <- data.frame(
  samples=colnames(fpkm),
  HER2_status=predict.HER2,
  PR_status=predict.PR,
  ER_status=predict.ER
)

write.csv(mouse.markers,"mouse-marker-status-prediction.csv",quote=F,row.names=F)
