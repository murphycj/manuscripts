library(DESeq2)
library(xlsx)

sample_data = read.xlsx("~/Desktop/Research/0816_sam/data/Bakhoum-RNAseq-Samples.xlsx","Bakhoum-RNAseq-Samples.csv")

counts <- read.csv("HTSeq.gene.counts.csv",row.names=1,header=T,check.names=F)
counts <- counts[rowSums(counts)>10,]

sample_data <- sample_data[sample_data[,"Sample.Name"] %in% colnames(counts),]
counts <- counts[,as.character(sample_data[,"Sample.Name"])]
coldata <- data.frame(
  group=sample_data$tissue,
  row.names=colnames(counts)
)
dds <- DESeqDataSetFromMatrix(counts,coldata,~group)
vst <- varianceStabilizingTransformation(dds,blind=T)
write.csv(assay(vst),"HTSeq.gene.vst.csv",quote=F)


sample_data <- sample_data[sample_data[,"Sample.Name"] %in% c("11","12","13","14","15","16","17","18"),]
counts <- counts[,as.character(sample_data[,"Sample.Name"])]
coldata <- data.frame(
  group=sample_data$tissue,
  row.names=colnames(counts)
)
dds <- DESeqDataSetFromMatrix(counts,coldata,~1)
vst <- varianceStabilizingTransformation(dds,blind=T)
write.csv(assay(vst),"HTSeq.gene.vst.cellLines.csv",quote=F)
