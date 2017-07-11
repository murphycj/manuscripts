library("beeswarm")
library("xlsx")
source("../lib.R")

chrs_names <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X","Y")
chrs_len <- list()
chrs_len[["1"]] <- 195471971
chrs_len[["2"]] <- 182113224
chrs_len[["3"]] <- 160039680
chrs_len[["4"]] <- 156508116
chrs_len[["5"]] <- 151834684
chrs_len[["6"]] <- 149736546
chrs_len[["7"]] <- 145441459
chrs_len[["8"]] <- 129401213
chrs_len[["9"]] <- 124595110
chrs_len[["10"]] <- 124595110
chrs_len[["11"]] <- 122082543
chrs_len[["12"]] <- 120129022
chrs_len[["13"]] <- 120421639
chrs_len[["14"]] <- 124902244
chrs_len[["15"]] <- 104043685
chrs_len[["16"]] <- 98207768
chrs_len[["17"]] <- 94987271
chrs_len[["18"]] <- 90702639
chrs_len[["19"]] <- 61431566
chrs_len[["X"]] <- 171031299
chrs_len[["Y"]] <- 91744698

tally <- function(cnas) {
  cnas.lines <- c()
  cnas.count <- c()

  n <- 0
  previous = 0
  previous.count <- 0
  for (chrom in chrs_names) {
    tmp <- cnas[cnas[,1]==chrom,]
    cnas.lines <- c(cnas.lines,n)
    cnas.count <- c(cnas.count,previous.count)
    if (nrow(tmp) > 0) {
      for (i in 1:nrow(tmp)) {
        start <- as.numeric(tmp[i,2])
        end <- as.numeric(tmp[i,3])
        count <- as.numeric(tmp[i,4])

        if ((start+n) != (previous)) {
          cnas.lines <- c(cnas.lines,previous)
          cnas.count <- c(cnas.count,0)
          previous.count <- 0
        }

        cnas.lines <- c(cnas.lines,start+n)
        cnas.count <- c(cnas.count,previous.count)

        cnas.lines <- c(cnas.lines,start+n)
        cnas.count <- c(cnas.count,count)

        cnas.lines <- c(cnas.lines,end+n)
        cnas.count <- c(cnas.count,count)

        previous.count <- count
        previous <- end+n
      }
    }

    cnas.lines <- c(cnas.lines,n+chrs_len[[chrom]])
    cnas.count <- c(cnas.count,previous.count)

    n <- n + chrs_len[[chrom]]
    previous.count <- 0
  }
  return(list(cnas.lines,cnas.count))
}

plot_recurrent_cnas <- function(cns,prefix) {

  amp <- data.frame()
  del <- data.frame()

  for (sample in names(cns)) {
    tmp <- cns[[sample]]
    tmp.up <- tmp[(tmp$cn >= 4),]
    amp <- rbind(amp,tmp.up)
    tmp.down <- tmp[(tmp$cn == 0),]
    del <- rbind(del,tmp.down)
  }
  amp <- amp[,c("chromosome","start","end")]
  amp$num <- rep(1,nrow(amp))
  del <- del[,c("chromosome","start","end")]
  del$num <- rep(1,nrow(del))

  write.table(amp,"amplifications.bed",sep="\t",quote=F,row.names=F,col.names=F)
  write.table(del,"deletions.bed",sep="\t",quote=F,row.names=F,col.names=F)

  system("~/Desktop/tools/bedtools2/bin/bedtools sort -i amplifications.bed -faidx names.txt > amplifications.sorted.bed")
  system("~/Desktop/tools/bedtools2/bin/bedtools merge -i amplifications.sorted.bed -c 4 -o sum > amplifications.sorted.merged.bed")
  system("~/Desktop/tools/bedtools2/bin/bedtools sort -i deletions.bed -faidx names.txt > deletions.sorted.bed")
  system("~/Desktop/tools/bedtools2/bin/bedtools merge -i deletions.sorted.bed -c 4 -o sum > deletions.sorted.merged.bed")

  amp <- read.table("amplifications.sorted.merged.bed",sep="\t")
  del <- read.table("deletions.sorted.merged.bed",sep="\t")

  rr <- tally(amp)
  amp.lines <- rr[[1]]
  amp.count <- rr[[2]]
  amp.count <- amp.count/length(cns)

  rr <- tally(del)
  del.lines <- rr[[1]]
  del.count <- rr[[2]]
  del.count <- del.count/length(cns)

  n<-0
  line_positions <- c()
  label_pos <- c()
  for (i in chrs_names) {
    label_pos <- c(label_pos,(n + (chrs_len[[i]]/2)))
    n<-n+chrs_len[[i]]
    line_positions <- c(line_positions,n)
  }

  pdf(paste("recurrentCNAs.",prefix,".pdf",sep=""),width=15,height=5)
  par(mar=c(6,6,4,4))
  plot(amp.lines, amp.count,type="n",ylim=c(-1,1),xaxt='n',yaxt='n',xlab="Chromosomes",xlim=c(0,n),ylab="Frequency")
  lines(amp.lines, amp.count,col="red",lwd=1.5)
  lines(del.lines, del.count*-1,col="blue",lwd=1.5)
  axis(1,labels=chrs_names,at=label_pos)
  axis(2,labels=c(1,0.5,0,0.5,1),at=c(-1,-0.5,0,0.5,1))
  segments(line_positions,-1.5,line_positions,1.5,xlim=c(0,n))
  dev.off()
}

expression_vs_cnv_plots <- function(cnvs,cnv_cn,fpkm,samples) {
  wex_mouse = as.character(samples[(samples$sequencing_type=="WEX"),"mouse"])
  samples_controls <- samples[(samples$control==1) & (samples$mouse %in% wex_mouse) & (samples$sequencing_type=="RNAseq"),]

  rnaseq_common <- c()
  wex_common <- c()

  for (i in row.names(samples_controls)) {
    tmp <- samples_controls[i,]
    tmp <- samples[(samples$mouse==tmp$mouse) & (samples$sequencing_type=="WEX") & (samples$tissue!="tail") & (samples$tissue!="liver"),]

    if (nrow(tmp)>0) {
      rnaseq_common <- c(rnaseq_common, rep(i,nrow(tmp)))
      wex_common <- c(wex_common, row.names(tmp))
    }
  }

  colors <- rep("black",length(rnaseq_common))
  colors[wex_common %in% names(cnv_cn["Met",cnv_cn["Met",]>3])]<-"red"
  pdf("Met-fpkm-cnv.pdf")
  par(mar=c(5,5,2,2),family="Times")
  plot(
    as.numeric(cnvs["Met",wex_common]),
    as.numeric(fpkm["Met",rnaseq_common]),
    pch=16,
    cex=2,
    cex.lab=1.5,
    xlab="log2 copy number",
    ylab="FPKM",
    col=colors,
    cex.axis=1.5
  )
  dev.off()

  colors <- rep("black",length(rnaseq_common))
  colors[wex_common %in% names(cnv_cn["Yap1",cnv_cn["Yap1",]>3])]<-"red"
  pdf("Yap1-fpkm-cnv.pdf")
  par(mar=c(5,5,2,2),family="Times")
  plot(
    as.numeric(cnvs["Yap1",wex_common]),
    as.numeric(fpkm["Yap1",rnaseq_common]),
    pch=16,
    cex=2,
    cex.lab=1.5,
    xlab="log2 copy number",
    ylab="FPKM",
    col=colors,
    cex.axis=1.5
  )
  dev.off()


  colors <- rep("black",length(rnaseq_common))
  colors[wex_common=="HL-WES-09"]<-"red"
  pdf("Pten-1660-fpkm-cnv.pdf")
  par(mar=c(5,5,2,2),family="Times")
  plot(
    as.numeric(cnvs["Pten",wex_common]),
    as.numeric(fpkm["Pten",rnaseq_common]),
    pch=16,
    cex=2,
    cex.lab=1.5,
    xlab="log2 copy number",
    ylab="FPKM",
    col=colors,
    cex.axis=1.5
  )
  dev.off()

  colors <- rep("black",length(rnaseq_common))
  colors[wex_common=="HL-WES-04"]<-"red"
  pdf("Fgfr2-1460-fpkm-cnv.pdf")
  par(mar=c(5,5,2,2),family="Times")
  plot(
    as.numeric(cnvs["Fgfr2",wex_common]),
    as.numeric(fpkm["Fgfr2",rnaseq_common]),
    pch=16,
    cex=2,
    cex.lab=1.5,
    xlab="log2 copy number",
    ylab="FPKM",
    col=colors,
    cex.axis=1.5
  )
  dev.off()

  colors <- rep("black",length(rnaseq_common))
  colors[wex_common=="HL-WES-08"]<-"red"
  pdf("Egfr-1616-fpkm-cnv.pdf")
  par(mar=c(5,5,2,2),family="Times")
  plot(
    as.numeric(cnvs["Egfr",wex_common]),
    as.numeric(fpkm["Egfr",rnaseq_common]),
    pch=16,
    cex=2,
    cex.lab=1.5,
    xlab="log2 copy number",
    ylab="FPKM",
    col=colors,
    cex.axis=1.5
  )
  dev.off()
}

cnvs <- read.csv("../gene-sample-log2.csv",row.names=1,header=T,check.names=F)
cnv_cn <- read.csv("../gene-sample-cn.csv",header=T,row.names=1,check.names=F)

fpkm <- read.csv(
  "/Users/charlesmurphy/Desktop/Research/0914_hui/results/RNAseq/Cufflinks/genes-symbols-fpkm.csv",
  row.names=1,
  check.names=F
)

samples <- read.xlsx("/Users/charlesmurphy/Desktop/Research/0914_hui/data/samples.xlsx","samples",row.names=1)

expression_vs_cnv_plots(cnvs,cnv_cn,fpkm,samples)

cns <- load_cnvkit_data_cns(filterY=T)

plot_recurrent_cnas(cns,"all")

samples <- read.xlsx("/Users/charlesmurphy/Desktop/Research/0914_hui/data/samples.xlsx","samples",row.names=1)
samples <- samples[(samples$sequencing_type=="WEX") & (samples$tissue=="primary"),]

cns.brca1 <- list()
for (sample in row.names(samples[samples$brca1=="flox/flox",])) {
  cns.brca1[[sample]] <- cns[[sample]]
}
plot_recurrent_cnas(cns.brca1,"brcaDEL")

cns.brca1 <- list()
for (sample in row.names(samples[samples$brca1=="WT",])) {
  cns.brca1[[sample]] <- cns[[sample]]
}
plot_recurrent_cnas(cns.brca1,"brcaWT")
