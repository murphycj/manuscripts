library(gplots)


plot_progression <- function(PADJ,LOG2FC) {
  fpkm <- as.matrix(read.csv("../Cufflinks/cleanUp/genesSymbols-fpkm.csv",row.names=1,header=T,check.names=F))

  low_vs_high <- read.csv("./CIN_low_vs_CIN_high_symbols/CIN_low_vs_CIN_high_symbols_results.csv",row.names=1)
  low_vs_high <- low_vs_high[!is.na(low_vs_high$padj),]
  low_vs_high <- low_vs_high[(low_vs_high$padj<=PADJ) & (low_vs_high$log2FoldChange>=LOG2FC),]

  low_vs_med <- read.csv("./CIN_low_vs_CIN_medium_symbols/CIN_low_vs_CIN_medium_symbols_results.csv",row.names=1)
  low_vs_med <- low_vs_med[!is.na(low_vs_med$padj),]
  low_vs_med <- low_vs_med[(low_vs_med$padj<=PADJ) & (low_vs_med$log2FoldChange>=LOG2FC),]

  common <- intersect(row.names(low_vs_med),row.names(low_vs_high))
  common <- intersect(common,row.names(fpkm))

  data_up <- cbind(apply(fpkm[common,c("MK","Kb")],1,mean),apply(fpkm[common,c("cont","Ka")],1,mean),fpkm[common,"MKH"])


  pdf("expression_up_regulation.pdf")
  boxplot(log2(data_up+0.01),ylab="log2(FPKM+0.1)",names=c("CIN-low","CIN-medium","CIN-high"),main="Average expression of\ncommon up-regulated genes")
  dev.off()


  fpkm <- as.matrix(read.csv("../Cufflinks/cleanUp/genesSymbols-fpkm.csv",row.names=1,header=T,check.names=F))

  low_vs_high <- read.csv("./CIN_low_vs_CIN_high_symbols/CIN_low_vs_CIN_high_symbols_results.csv",row.names=1)
  low_vs_high <- low_vs_high[!is.na(low_vs_high$padj),]
  low_vs_high <- low_vs_high[(low_vs_high$padj<=PADJ) & (low_vs_high$log2FoldChange<=-LOG2FC),]

  low_vs_med <- read.csv("./CIN_low_vs_CIN_medium_symbols/CIN_low_vs_CIN_medium_symbols_results.csv",row.names=1)
  low_vs_med <- low_vs_med[!is.na(low_vs_med$padj),]
  low_vs_med <- low_vs_med[(low_vs_med$padj<=PADJ) & (low_vs_med$log2FoldChange<=-LOG2FC),]

  common <- intersect(row.names(low_vs_med),row.names(low_vs_high))
  common <- intersect(common,row.names(fpkm))

  data_down <- cbind(apply(fpkm[common,c("MK","Kb")],1,mean),apply(fpkm[common,c("cont","Ka")],1,mean),fpkm[common,"MKH"])


  pdf("expression_down_regulation.pdf")
  boxplot(log2(data_down+0.01),ylab="log2(FPKM+0.1)",names=c("CIN-low","CIN-medium","CIN-high"),main="Average expression of\ncommon down-regulated genes")
  dev.off()

}

plots <- function(PADJ,LOG2FC) {

  fpkm <- as.matrix(read.csv("../Cufflinks/cleanUp/genesSymbols-fpkm.csv",row.names=1,header=T,check.names=F))

  low_vs_high <- read.csv("./CIN_low_vs_CIN_high_symbols/CIN_low_vs_CIN_high_symbols_results.csv",row.names=1)
  low_vs_high <- low_vs_high[!is.na(low_vs_high$padj),]
  low_vs_high <- low_vs_high[(low_vs_high$padj<=PADJ) & (abs(low_vs_high$log2FoldChange)>=LOG2FC),]
  common <- intersect(row.names(low_vs_high),row.names(fpkm))
  low_vs_high <- low_vs_high[common,]
  fpkm <- fpkm[common,]

  low_vs_med <- read.csv("./CIN_low_vs_CIN_high_symbols/CIN_low_vs_CIN_high_symbols_results.csv",row.names=1)
  low_vs_med <- low_vs_med[!is.na(low_vs_med$padj),]
  low_vs_med <- low_vs_med[(low_vs_med$padj<=PADJ) & (abs(low_vs_med$log2FoldChange)>=LOG2FC),]
  common <- intersect(row.names(low_vs_med),row.names(fpkm))
  low_vs_med <- low_vs_med[common,]
  fpkm <- fpkm[common,]

  common <- intersect(row.names(low_vs_med),row.names(low_vs_high))

  png(paste("common-padj",PADJ,"-log2fc",LOG2FC,".png",sep=""),width=1000,height=1000)
  heatmap.2(
    log(fpkm[common,]+0.1),
    trace="none",
    #ColSideColors=colors,
    keysize=1.0,
    scale="row",
    margins=c(8,12),
    cexRow=0.2,
    labRow="",
    cexCol=2.5,
    Colv=F,
    dendrogram="row",
    distfun=function(x) as.dist(1-cor(t(x),method="spearman")),
    col=redgreen(75),
    key.title="Gene expression",
    main=paste("Differentially expressed genes\n(padj<=",PADJ,", log2FC>=",LOG2FC,")",sep="")
  )
  #legend("topright",c("r-spondin","APC-WT","APC-/-"),fill=cols)
  dev.off()


  #cols=rainbow(3)
  #colors = c(rep(cols[1],5),rep(cols[2],6),rep(cols[3],4))
  png(paste("low_vs_high-padj",PADJ,"-log2fc",LOG2FC,".png",sep=""),width=1000,height=1000)
  heatmap.2(
    log(fpkm[row.names(low_vs_high),]+0.1),
    trace="none",
    #ColSideColors=colors,
    keysize=1.0,
    scale="row",
    margins=c(8,12),
    cexRow=0.2,
    labRow="",
    cexCol=1.2,
    Colv=F,
    dendrogram="row",
    distfun=function(x) as.dist(1-cor(t(x),method="spearman")),
    col=redgreen(75),
    key.title="Gene expression",
    main=paste("Differentially expressed genes\n(padj<=",PADJ,", log2FC>=",LOG2FC,")",sep="")
  )
  #legend("topright",c("r-spondin","APC-WT","APC-/-"),fill=cols)
  dev.off()

  low_vs_med <- read.csv("./CIN_low_vs_CIN_high_symbols/CIN_low_vs_CIN_high_symbols_results.csv",row.names=1)
  low_vs_med <- low_vs_med[!is.na(low_vs_med$padj),]
  common <- intersect(row.names(low_vs_med),row.names(fpkm))
  low_vs_med <- low_vs_med[common,]
  fpkm <- fpkm[common,]


  cols=rainbow(3)
  colors = c(rep(cols[1],5),rep(cols[2],6),rep(cols[3],4))
  png(paste("low_vs_med-padj",PADJ,"-log2fc",LOG2FC,".png",sep=""),width=1000,height=1000)
  heatmap.2(
    log(fpkm[row.names(low_vs_med),]+0.1),
    trace="none",
    #ColSideColors=colors,
    keysize=1.0,
    scale="row",
    margins=c(8,12),
    cexRow=0.2,
    labRow="",
    cexCol=1.2,
    Colv=F,
    dendrogram="row",
    distfun=function(x) as.dist(1-cor(t(x),method="spearman")),
    col=redgreen(75),
    key.title="Gene expression",
    main=paste("Differentially expressed genes\n(padj<=",PADJ,", log2FC>=",LOG2FC,")",sep="")
  )
  #legend("topright",c("r-spondin","APC-WT","APC-/-"),fill=cols)
  dev.off()
}

plots(0.05,2)

plot_progression(0.05,2)

low_vs_high <- read.csv("./CIN_low_vs_CIN_high_symbols/CIN_low_vs_CIN_high_symbols_results.csv",row.names=1)
low_vs_high <- low_vs_high[!is.na(low_vs_high$padj),]
low_vs_high <- low_vs_high[low_vs_high$padj<=0.05,]

low_vs_med <- read.csv("./CIN_low_vs_CIN_medium_symbols/CIN_low_vs_CIN_medium_symbols_results.csv",row.names=1)
low_vs_med <- low_vs_med[!is.na(low_vs_med$padj),]
low_vs_med <- low_vs_med[low_vs_med$padj<=0.05,]

data <- list()
data[["CIN-low vs CIN-medium"]] <- row.names(low_vs_med)
data[["CIN-low vs CIN-high"]] <- row.names(low_vs_high)

pdf("venn.pdf")
venn(data)
dev.off()
