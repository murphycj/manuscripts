library(DESeq2)

samples <- read.csv("../../../data/samples.csv")

counts <- read.csv("HTSeq.geneSymbols.counts.csv",row.names=1,header=T)
counts <- counts[rowSums(counts)>1,]
counts <- counts[,samples$sample]

coldata <- data.frame(
  group=samples$group,
  row.names=colnames(counts)
)

dds <- DESeqDataSetFromMatrix(counts,coldata,~group)

vst <- varianceStabilizingTransformation(dds,blind=T)

write.csv(assay(vst),"HTSeq.geneSymbols.counts.vst.csv",quote=F)




counts <- read.csv("HTSeq.gene.counts.csv",row.names=1,header=T)
counts <- counts[rowSums(counts)>1,]
counts <- counts[,samples$sample]

coldata <- data.frame(
  group=samples$group,
  row.names=colnames(counts)
)

dds <- DESeqDataSetFromMatrix(counts,coldata,~group)

vst <- varianceStabilizingTransformation(dds,blind=T)

write.csv(assay(vst),"HTSeq.gene.counts.vst.csv",quote=F)
