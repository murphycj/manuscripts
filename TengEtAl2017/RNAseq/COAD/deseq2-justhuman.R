library(DESeq2)
library("BiocParallel")
register(MulticoreParam(26))

name="COAD_normal_vs_tumor"

counts <- read.csv("COAD.counts.csv",row.names=1,header=T,check.names=F)

samples <- colnames(counts)
samples[grep("-11",samples)] <- "normal"
samples[samples!="normal"] <- "tumor"

coldata <- data.frame(tumor=samples, row.names=colnames(counts))
countTable <- DESeqDataSetFromMatrix(countData=counts,colData=coldata,design=~tumor)
result <- DESeq(countTable,parallel=T)
res <- results(result,parallel=T)


dd = res[with(res,order(padj)),]

write.csv(dd,paste("./",name,"/",name,"_results.csv",sep=""))
pdf(paste("./",name,"/",name,"-pvalue-hist.pdf",sep=""))
hist(res$pvalue, breaks=20, col="grey")
dev.off()

pdf(paste("./",name,"/",name,"-padj-hist.pdf",sep=""))
hist(res$padj, breaks=20, col="grey")
dev.off()

pdf(paste("./",name,"/",name,"-MA-plot.pdf",sep=""))
plotMA(result)
dev.off()
