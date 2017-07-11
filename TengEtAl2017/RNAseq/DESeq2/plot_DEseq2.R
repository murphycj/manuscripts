library(gplots)

plots <- function(PADJ,LOG2FC) {

  fpkm <- as.matrix(read.csv("../HTSeqCount/HTSeq.geneSymbol.vst.csv",row.names=1,header=T))
  fpkm <- fpkm[,c(
    "EN15d","EN4d","EN_M","EN_F","NAIVE_F",
    "A1","A2","A5","A6","NAIVE","LIPO",
    "A1_plus","A2_plus","A5_plus","A6_plus"
  )]

  apcwt_vs_mut <- read.csv("./APCWT_vs_APCMUT/APCWT_vs_APCMUT_symbols_results.csv",row.names=1)
  common <- intersect(row.names(apcwt_vs_mut),row.names(fpkm))
  apcwt_vs_mut <- apcwt_vs_mut[common,]
  fpkm <- fpkm[common,]


  apcwt_vs_mut <- apcwt_vs_mut[!is.na(apcwt_vs_mut$padj),]
  apcwt_vs_mut <- apcwt_vs_mut[(apcwt_vs_mut$padj<=PADJ) & (abs(apcwt_vs_mut$log2FoldChange)>=LOG2FC),]


  cols=rainbow(3)
  colors = c(rep(cols[1],5),rep(cols[2],6),rep(cols[3],4))
  png(paste("apcwt_vs_mut_Genes-padj",PADJ,"-log2fc",LOG2FC,".png",sep=""),width=500,height=1000)
  heatmap.2(
    log(fpkm[row.names(apcwt_vs_mut),]+0.1),
    trace="none",
    ColSideColors=colors,
    keysize=1.0,
    scale="row",
    margins=c(8,12),
    cexRow=0.2,
    labRow="",
    cexCol=1.2,
    Colv=F,
    dendrogram="row",
    distfun=function(x) as.dist(1-cor(t(x),method="spearman")),
    col=colorpanel(75,"red","white","blue"),
    key.title="Gene expression",
    main=paste("Differentially expressed genes\n(padj<=",PADJ,", log2FC>=",LOG2FC,")",sep="")
  )
  legend("topright",c("r-spondin","APC-WT","APC-/-"),fill=cols)
  dev.off()

  wt_vs_rspo <- read.csv("./WT_vs_Rspondin/WT_vs_Rspondin_results.csv",row.names=1)
  wt_vs_rspo <- wt_vs_rspo[!is.na(wt_vs_rspo$padj),]
  wt_vs_rspo <- wt_vs_rspo[(wt_vs_rspo$padj<=PADJ) & (abs(wt_vs_rspo$log2FoldChange)>=LOG2FC),]


  cols=rainbow(3)
  colors = c(rep(cols[1],5),rep(cols[2],6),rep(cols[3],4))
  png(paste("WT-vs-rspondin-padj",PADJ,"-log2fc",LOG2FC,".png",sep=""),width=500,height=500)
  heatmap.2(
    log(fpkm[row.names(wt_vs_rspo),]+0.1),
    trace="none",
    ColSideColors=colors,
    keysize=1.0,
    scale="row",
    margins=c(8,12),
    cexRow=1.2,
    cexCol=1.2,
    Colv=F,
    dendrogram="row",
    distfun=function(x) as.dist(1-cor(t(x),method="spearman")),
    col=colorpanel(75,"red","white","blue"),
    key.title="Gene expression",
    main=paste("Differentially expressed genes\n(padj<=",PADJ,", log2FC>=",LOG2FC,")",sep="")
  )
  legend("topright",c("r-spondin","APC-WT","APC-/-"),fill=cols,cex=0.8)
  dev.off()
}

plots(0.05,1)
