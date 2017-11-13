library("beeswarm")
library("RColorBrewer")
library("gplots")
library(AIMS)
library("xlsx")
library(leaflet)

plot_results <- function(aims.result,fpkm,samples) {

  all.probs <- aims.result$all.probs[["20"]]
  all.probs <- all.probs[,c("Basal","LumA","LumB","Her2","Normal")]
  row.names(all.probs) <- row.names(aims.result$cl)

  all.probs <- all.probs[row.names(all.probs) %in% as.character(samples$sample_ID),]
  all.probs <- all.probs[order(as.numeric(all.probs[,1]),decreasing=T),]

  samples <- samples[as.character(samples$sample_ID) %in% row.names(all.probs),]
  samples <- samples[row.names(all.probs),]

  mouse.markers <- read.csv("/Users/charlesmurphy/Desktop/Research/0914_hui/results/RNAseq/SubtypeClassification/primaryControl/Markers/mouse-marker-status-prediction.csv",header=T,row.names=1,colClasses="character")
  mouse.markers <- mouse.markers[row.names(all.probs),]
  mouse.markers[mouse.markers=="Positive"]<-"black"
  mouse.markers[mouse.markers=="Negative"]<-"white"

  sample_type = as.character(samples$tissue)
  sample_type[sample_type=="primary"]<-"#999999"
  sample_type[sample_type=="implant"]<-"#999999"
  sample_type[sample_type=="normal breast"]<-"green"

  N <- nrow(all.probs)

  png('Mouse-classification-AIMS.png',width=1400,height=250)
  par(mar=c(6,7,4,3),cex=0.9)
  pp=barplot(
    t(all.probs),
    las=2,
    xlab="",
    xlim=c(0,N+30),
    ylab="Probability of subtype\nassignment",
    names.arg=rep("",N),
    col=c("red","blue","#00ffff","purple","green"),
    cex.lab=1.5
  )
  par(xpd=T)
  legend(
    par("usr")[2]/5,1.35,
    c("Basal-like","Luminal A","Luminal B","HER2-enriched","Normal-like"),
    col=c("red","blue","#00ffff","purple","green"),
    pch=15,
    cex=1.4,
    ncol=5
  )
  title(xlab="Samples", line=4.5, cex.lab=1.5)

  #plot sample type
  rect(
    pp-0.5,
    rep(-0.04,N),
    pp+0.5,
    rep(-0.11,N),
    col=sample_type
  )

  #plot ER/PR/HER2 statuses
  rect(
    pp-0.5,
    rep(-0.19,N),
    pp+0.5,
    rep(-0.12,N),
    col=mouse.markers$PR_status
  )
  rect(
    pp-0.5,
    rep(-0.27,N),
    pp+0.5,
    rep(-0.20,N),
    col=mouse.markers$ER_status
  )
  rect(
    pp-0.5,
    rep(-0.35,N),
    pp+0.5,
    rep(-0.28,N),
    col=mouse.markers$HER2_status
  )

  #x_box=122.5 #all
  x_box=75.5
  y_box=-0.48

  text(x_box-6,-0.5,"Sample type",cex=1.4)
  text(x_box+1.6,y_box+0.035,"Tumor",cex=1.4,adj=c(0,0.5))
  rect(x_box,y_box,x_box+1.2,y_box+0.07,col="#999999")
  text(x_box+1.6,y_box-0.12+0.035,"Normal",cex=1.4,adj=c(0,0.5))
  rect(x_box,y_box-0.12,x_box+1.2,y_box-0.12+0.07,col="green")

  x_box <- x_box + 25

  text(x_box-8,-0.5,"PR/ER/HER2 status",cex=1.4)
  text(x_box+1.6,y_box+0.035,"Negative",cex=1.4,adj=c(0,0.5))
  rect(x_box,y_box,x_box+1.2,y_box+0.07,col="white")
  text(x_box+1.6,y_box-0.12+0.035,"Positive",cex=1.4,adj=c(0,0.5))
  rect(x_box,y_box-0.12,x_box+1.2,y_box-0.12+0.07,col="black")

  text(max(pp)+0.7,-0.075,"Sample type",cex=1,adj=c(0,0.5))
  text(max(pp)+0.7,-0.155,"PR status",cex=1,adj=c(0,0.5))
  text(max(pp)+0.7,-0.235,"ER status",cex=1,adj=c(0,0.5))
  text(max(pp)+0.7,-0.315,"HER2 status",cex=1,adj=c(0,0.5))

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

samples <- read.xlsx("/Users/charlesmurphy/Desktop/Research/0914_hui/data/samples.xlsx","samples",colClasses="character")
samples <- samples[(samples['byl']==0) & (samples['bkm']==0) & (samples['bmn']==0),]
samples <- samples[samples$sequencing_type=="RNAseq",]
row.names(samples) <- as.character(samples$sample_ID)

genes <- read.csv("genes.csv")[,1]

fpkm <- read.csv(
  "/Users/charlesmurphy/Desktop/Research/0914_hui/results/RNAseq/Cufflinks/genes-humanEntrez-fpkm.csv",
  row.names=1,
  header=T,
  check.names=F
)
common <- intersect(as.character(samples$sample_ID),colnames(fpkm))

fpkm <- fpkm[,common]
common <- intersect(as.character(genes),row.names(fpkm))
fpkm <- fpkm[common,]

aims.result <- run_AIMS(fpkm)

plot_results(aims.result,fpkm,samples)
