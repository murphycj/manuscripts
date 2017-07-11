library(pvclust)
library(snow)
library(dendextend)
library("beeswarm")
library("gplots")

vst.mouse <- read.csv("../../HTSeqCount/HTSeq.humanSymbols.vst.csv",row.names=1,header=T)

rsem <- read.csv("../COAD.gene.rsem.csv",row.names=1,header=T,check.names=F)
rsem <- rsem[apply(rsem,1,mean)>1,]

de.genes <- read.csv("COAD_normal_vs_tumor_results.csv",row.names=1,header=T)
de.genes <- de.genes[!is.na(de.genes$padj),]
de.genes <- de.genes[(de.genes$padj<=0.01) & (abs(de.genes$log2FoldChange)>=3) & (row.names(de.genes) %in% row.names(rsem)),]

de.genes.up <- row.names(de.genes[de.genes$log2FoldChange>0,])
print(length(de.genes.up))
de.genes.down <- row.names(de.genes[de.genes$log2FoldChange<0,])
print(length(de.genes.down))

sink("COAD_gene_signatures.gmt")
cat("COAD_UP\tNULL")
for (i in de.genes.up) {
  cat("\t")
  cat(i)
}
cat("\n")

cat("COAD_DOWN\tNULL")
for (i in de.genes.down) {
  cat("\t")
  cat(i)
}
cat("\n")
sink()


common <- intersect(row.names(de.genes),row.names(vst.mouse))


png("de-genes-mouse_coadNormal_vs_tumor.png",width=1000,height=1000)
heatmap.2(
  as.matrix(vst.mouse[common,]),
  trace="none",
  keysize=1.0,
  scale="row",
  margins=c(8,12),
  labRow="",
  col=colorpanel(75,"blue","yellow"),
	distfun=function(x) as.dist(1-cor(t(x),method="spearman", use = "complete.obs")),
  key.title="Gene expression",
  main=""
)
dev.off()
