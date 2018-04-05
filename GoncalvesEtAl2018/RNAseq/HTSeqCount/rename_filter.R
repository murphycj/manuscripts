library(DESeq2)

sample_data <- read.csv("../../../data/all_samples.csv")
row.names(sample_data) <- sample_data$ID

data <- read.csv("HTSeq.gene.counts.csv",row.names=1,check.names=F,header=T)
colnames(data) <-as.character(sample_data[colnames(data),"sample"])
data <- data[apply(data,1,sum)>2,]

coldata <- data.frame(
  group=as.character(sample_data.tmp$Tissue),
  row.names=colnames(data.tmp)
)

row.names(coldata) <- colnames(data.tmp)
dds <- DESeqDataSetFromMatrix(data.tmp,coldata,~1)
vst <- varianceStabilizingTransformation(dds,blind=T)
vst <- assay(vst)
write.csv(vst,"vst.csv",quote=F)

pdf("skedasticity.pdf")
par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE,cex=1.5)
smoothScatter(
  apply(vst,1,mean),
  apply(vst,1,sd),
  ylab="SD",
  xlab="mean"
)
dev.off()
