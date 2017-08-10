library(DESeq2)
library(xlsx)
library(RColorBrewer)
library(ggplot2)

perform_pcas <- function(prefix,data,samples) {

	re <- prcomp(t(data))
  pcs=re$sdev**2
	pcs = pcs/sum(pcs) * 100

	write.csv(re$x,paste(prefix,"PC.values.csv",sep=""),quote=F)

  for (i in 1:4) {
    for (j in 1:4) {
      if (i >= j) {
        next
      }
      data <- data.frame(
        Name=row.names(samples),
        PC1=re$x[,i],
        PC2=re$x[,j],
        tissue=samples$tissue,
        cell_line=samples$cell_line_origin
      )

      pdf(paste(prefix,"PC",i,"PC",j,".pdf",sep=""))
      print(
        ggplot(data,aes(y=PC2,x=PC1,color=cell_line,shape=tissue,label=Name)) +
        geom_point() +
        xlab(paste("PC",i," (",round(pcs[1],2),"%)",sep="")) +
        ylab(paste("PC",j," (",round(pcs[2],2),"%)",sep="")) +
        geom_text(aes(label=Name),hjust=0, vjust=0)
      )
      dev.off()
    }
  }
}

samples = read.xlsx("~/Desktop/Research/0816_sam/data/Bakhoum-RNAseq-Samples.xlsx","Bakhoum-RNAseq-Samples.csv")
row.names(samples) <- samples[,"Sample.Name"]
data <- read.csv("../HTSeqCount/HTSeq.gene.vst.csv",header=T,row.names=1,check.names=F)
samples <- samples[colnames(data),]

perform_pcas("./PCAall/",data,samples)
data.sub <- data[order(apply(data,1,sd),decreasing=T),]
perform_pcas("./PCAtop1000Var/",data.sub[1:1000,],samples)
perform_pcas("./PCAtop500Var/",data.sub[1:500,],samples)
data.sub <- data[order(apply(data,1,median),decreasing=T),]
perform_pcas("./PCAtop1000median/",data.sub[1:1000,],samples)
perform_pcas("./PCAtop500median/",data.sub[1:500,],samples)
