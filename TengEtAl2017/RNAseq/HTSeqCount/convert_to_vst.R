library(DESeq2)

samples <- read.csv("../../../data/samples.csv")

counts <- read.csv("HTSeq.gene.counts.csv",row.names=1,header=T)
counts <- counts[rowSums(counts)>1,]
counts <- counts[,samples$Sample]

coldata <- data.frame(
  group=samples$Group.name,
  row.names=colnames(counts)
)

dds <- DESeqDataSetFromMatrix(counts,coldata,~group)

vst <- varianceStabilizingTransformation(dds,blind=T)

write.csv(assay(vst),"HTSeq.gene.counts.vst.csv",quote=F)



counts <- read.csv("HTSeq.geneSymbol.counts.csv",row.names=1,header=T)
counts <- counts[rowSums(counts)>1,]
counts <- counts[,samples$Sample]

coldata <- data.frame(
  group=samples$Group.name,
  row.names=colnames(counts)
)

dds <- DESeqDataSetFromMatrix(counts,coldata,~group)

vst <- varianceStabilizingTransformation(dds,blind=T)

write.csv(assay(vst),"HTSeq.geneSymbol.vst.csv",quote=F)
