library(pvclust)
library(snow)
library("beeswarm")
library("RColorBrewer")
library("gplots")
library(AIMS)
library("xlsx")

plot_results <- function(aims.result,fpkm,rpkm.entrez) {

  all.probs <- aims.result$all.probs[["20"]]
  all.probs <- all.probs[,c("Basal","LumA","LumB","Her2","Normal")]
  row.names(all.probs) <- row.names(aims.result$cl)

  samples <- read.xlsx("/Users/charlesmurphy/Desktop/Research/0914_hui/data/samples.xlsx","samples",row.names=1,colClasses="character")
  samples <- samples[samples$sequencing_type=="RNAseq",]
  samples <- samples[as.character(samples$sample_name) %in% row.names(all.probs),]

  all.probs <- all.probs[row.names(all.probs) %in% as.character(samples$sample_name),]
  all.probs <- all.probs[order(as.numeric(all.probs[,1]),decreasing=T),]

  png('Mouse-classification-aim2.png',width=1400,height=250)
  par(mar=c(5,7,4,2),cex=0.9)
  barplot(
    t(all.probs),
    las=2,
    xlab="Samples",
    xlim=c(0,nrow(all.probs)),
    ylab="Probability of subtype\nassignment",
    names.arg=rep("",nrow(all.probs)),
    col=c("red","blue","#00ffff","purple","green"),
    cex.lab=1.5
  )
  par(xpd=T)
  legend(
    par("usr")[2]/2,1.35,
    c("Basal-like","Luminal A","Luminal B","HER2-enriched","Normal-like"),
    col=c("red","blue","#00ffff","purple","green"),
    pch=15,
    cex=1.4,
    ncol=5
  )
  dev.off()

  all.probs <- t(all.probs)

  r <- data.frame(
    sample=as.character(unlist(lapply(colnames(all.probs),function(x) return(rep(x,5))))),
    probability=c(all.probs),
    color=rep(c("red","blue","#00ffff","purple","green"),ncol(all.probs))
    )

  r$sample = factor(r$sample,levels=unique(r$sample))

  png("Mouse-classification-aim.png",width=1400,height=300)
  par(mar=c(5,5,3,3),family="Times")
  beeswarm(
    probability~sample,
    data=r,
    las=2,
    cex=1.2,
    cex.axis=1.5,
    cex.lab=1.5,
    pch=16,
    xlim=c(1,nrow(samples)+30),
    pwcol=as.character(color),
    xlab="Samples",
    ylab="Probability of subtype assignment",
    main="",
    labels=F
    )
  #mtext("Class probability", side=2, line=3, cex=2)
  #mtext("AIMS classification", side=3, line=1, cex=2)

  n=1

  #for (i in levels(r$sample)) {
  #  brca1 <- as.character(samples[samples$sample_name==i,"brca1"])
  #  i = gsub("parental ","P. ",i)
  #  i = gsub("\\*","",i)
  #  i = gsub("KO","",i)
  #  i = gsub("HET","",i)
  #  i = gsub("WT","",i)
  #  i = gsub("WhiteL","W",i)
  #  if (is.na(brca1)) {
  #    mtext(i,side=1,at=n,col="black",las=2,line=2,cex=1.3)
  #  } else if (brca1=="WT") {
  #    mtext(i,side=1,at=n,col="black",las=2,line=2,cex=1.3)
  #  } else {
  #    mtext(i,side=1,at=n,col="red",las=2,line=2,cex=1.3)
  #  }
  #  n=n+1
  #}
  legend(
    "topright",
    legend=c("Basal-like","Luminal A","Luminal B","Her2-enriched","Normal-like"),
    pch=16,
    cex=1.8,
    col=c("red","blue","#00ffff","purple","green")
  )
  dev.off()
}


run_AIMS <- function(fpkm) {

  aims.result <- applyAIMS(as.matrix(fpkm),row.names(fpkm))

  cl <- data.frame(mouse=row.names(aims.result$cl),classification=as.character(aims.result$cl[,"20"]))
  write.csv(cl,"AIMS.mouse.TNBC.classification.csv",quote=F,row.names=F)

  all.prob <- aims.result$all.probs[["20"]]
  row.names(all.prob) <- row.names(aims.result$cl)
  write.csv(all.prob,"AIMS.mouse.TNBC.probabilities.csv",quote=F,row.names=T)

  return(aims.result)

}

samples <- read.xlsx("/Users/charlesmurphy/Desktop/Research/0914_hui/data/samples.xlsx","samples",row.names=1,colClasses="character")
samples <- samples[samples$sequencing_type=="RNAseq",]

genes <- read.csv("genes.csv")[,1]

fpkm <- read.csv(
  "/Users/charlesmurphy/Desktop/Research/0914_hui/results/RNAseq/Cufflinks/rename_filtered/genes-humanEntrez-fpkm.csv",
  row.names=1,
  header=T,
  check.names=F
)
common <- intersect(as.character(samples$sample_name),colnames(fpkm))
fpkm <- fpkm[,common]

load("~/Desktop/Research/data/GDAC/BRCA/gdac.broadinstitute.org_BRCA.Merge_rnaseq__illuminahiseq_rnaseq__unc_edu__Level_3__gene_expression__data.Level_3.2016012800.0.0/process_data/rpkm.entrez.RData")

common <- intersect(as.character(genes),intersect(row.names(fpkm),row.names(rpkm.entrez)))
rpkm.entrez <- rpkm.entrez[common,]
fpkm <- fpkm[common,]

#remove normals

samples <- colnames(rpkm.entrez)
samples[grep("-11",samples)] <- "normal"
samples[samples!="normal"] <- "tumor"
rpkm.entrez <- rpkm.entrez[,samples=="tumor"]

#simplify sample names

colnames(rpkm.entrez) <- unlist(
  lapply(
    colnames(rpkm.entrez),
    function(x) return(
        paste(strsplit(x,"-")[[1]][1:3],collapse="-")
      )
    )
  )

subtypes <- read.csv("/Users/charlesmurphy/Desktop/Research/data/papers/2012_comprehensiveMolecularPortraitsOfHumanBreastCancer/subtype.csv",header=T,check.names=F)
common <- intersect(colnames(rpkm.entrez),as.character(subtypes[,"Complete TCGA ID"]))
rpkm.entrez <- rpkm.entrez[,common]
subtypes <- subtypes[as.character(subtypes[,"Complete TCGA ID"]) %in% common,]
rpkm.entrez <- rpkm.entrez[,as.character(as.character(subtypes[,"Complete TCGA ID"]))]
rpkm.entrez <- rpkm.entrez[,!is.na(subtypes[,"PAM50 mRNA"])]
subtypes <- subtypes[!is.na(subtypes[,"PAM50 mRNA"]),]

aims.result <- run_AIMS(fpkm)

fpkm <- log2(fpkm + 0.1)
rpkm.entrez <- log2(rpkm.entrez + 0.1)

common <- intersect(row.names(fpkm),row.names(rpkm.entrez))
rpkm.entrez <- rpkm.entrez[common,]
fpkm <- fpkm[common,]

plot_results(aims.result,fpkm,rpkm.entrez)
