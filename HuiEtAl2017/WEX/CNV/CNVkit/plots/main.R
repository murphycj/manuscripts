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

plot_cna_fpkm <- function(gene,fpkm,cnvs,cnv_cn,amp) {

  samples.rnaseq <- c()
  samples.wex <- c()
  samples.brca1 <- c()
  samples.fpkm <- c()

  cc <- unique(samples[,c("mouse","primary_tumor_index")])
  for (i in 1:nrow(cc)) {
    tmp <- samples[(samples$mouse==cc[i,1]) & (samples$primary_tumor_index==cc[i,2]),]
    tmp.wex <- row.names(tmp[tmp$sequencing_type=="WEX",])

    if (length(tmp.wex)>1) {
      print("more than one WES sample")
    }

    if (("RNAseq" %in% tmp$sequencing_type) & ("WEX" %in% tmp$sequencing_type)) {
      tmp.rnaseq <- row.names(tmp[tmp$sequencing_type=="RNAseq",])

      if (length(tmp.rnaseq)>1) {
        if ("primary" %in% as.character(tmp[tmp$sequencing_type=="RNAseq","tissue"])) {

          tmp.rnaseq <- row.names(tmp[(tmp$sequencing_type=="RNAseq") & (tmp$tissue=="primary"),])

          if (length(tmp.rnaseq)>1) {
            samples.fpkm <- c(samples.fpkm,mean(fpkm[gene,tmp.rnaseq]))
          } else {
            samples.fpkm <- c(samples.fpkm,fpkm[gene,tmp.rnaseq])
          }
        } else {
          tmp.rnaseq <- row.names(tmp[(tmp$sequencing_type=="RNAseq") & (tmp$tissue=="implant"),])
          if (length(tmp.rnaseq)>1) {
            samples.fpkm <- c(samples.fpkm,mean(fpkm[gene,tmp.rnaseq]))
          } else {
            samples.fpkm <- c(samples.fpkm,fpkm[gene,tmp.rnaseq])
          }
        }
      } else {
        samples.fpkm <- c(samples.fpkm,fpkm[gene,tmp.rnaseq])
      }

      samples.rnaseq <- c(samples.rnaseq,paste(tmp.rnaseq,collapse=";"))
      samples.wex <- c(samples.wex,tmp.wex)
      samples.brca1 <- c(samples.brca1,as.character(tmp[tmp.wex,"brca1"]))
    }
  }

  if (length(samples.rnaseq)!=length(samples.wex)) {
    print("BOO")
  }

  data <- data.frame(
    WEX_sample=samples.wex,
    RNAseq_sample=samples.rnaseq,
    BRCA1=samples.brca1,
    FPKM=samples.fpkm,
    log2=as.numeric(as.matrix(cnvs[gene,samples.wex])),
    CN=as.numeric(as.matrix(cnv_cn[gene,samples.wex]))
  )

  write.csv(data, paste(gene,"-fpkm-cnv.csv",sep=""), row.names=F)

  colors <- rep("black",length(samples.wex))
  if (amp) {
    colors[samples.wex %in% colnames(cnv_cn)[as.numeric(as.matrix(cnv_cn[gene,]))>3]]<-"red"
  } else {
    colors[samples.wex %in% colnames(cnv_cn)[as.numeric(as.matrix(cnv_cn[gene,]))==0]]<-"red"
  }

  pdf(paste(gene,"-fpkm-cnv.pdf",sep=""))
  par(mar=c(5,5,2,2),family="Times")
  plot(
    as.numeric(as.matrix(cnvs[gene,samples.wex])),
    samples.fpkm,
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
cnv_cn <- cnv_cn[-c(1),]

fpkm <- read.csv(
  "/Users/charlesmurphy/Desktop/Research/0914_hui/results/RNAseq/Cufflinks/genes-symbols-fpkm.csv",
  row.names=1,
  check.names=F
)

samples <- read.xlsx("/Users/charlesmurphy/Desktop/Research/0914_hui/data/samples.xlsx","samples",row.names=1)
samples <- samples[(samples['byl']==0) & (samples['bkm']==0) & (samples['bmn']==0) & (samples$tissue %in% c("primary","implant")),]

common <- intersect(colnames(fpkm),row.names(samples))
fpkm <- fpkm[,common]
fpkm <- t(scale(t(fpkm)))

plot_cna_fpkm("Met",fpkm,cnvs,cnv_cn,T)
plot_cna_fpkm("Yap1",fpkm,cnvs,cnv_cn,T)
plot_cna_fpkm("Pten",fpkm,cnvs,cnv_cn,F)
plot_cna_fpkm("Fgfr2",fpkm,cnvs,cnv_cn,T)
plot_cna_fpkm("Egfr",fpkm,cnvs,cnv_cn,T)
plot_cna_fpkm("Myc",fpkm,cnvs,cnv_cn,T)

#cns <- load_cnvkit_data_cns(filterY=T)

#plot_recurrent_cnas(cns,"all")

#samples <- read.xlsx("/Users/charlesmurphy/Desktop/Research/0914_hui/data/samples.xlsx","samples",row.names=1)
#samples <- samples[(samples$sequencing_type=="WEX") & (samples$tissue=="primary"),]

#cns.brca1 <- list()
#for (sample in row.names(samples[samples$brca1=="flox/flox",])) {
#  cns.brca1[[sample]] <- cns[[sample]]
#}
#plot_recurrent_cnas(cns.brca1,"brcaDEL")

#cns.brca1 <- list()
#for (sample in row.names(samples[samples$brca1=="WT",])) {
#  cns.brca1[[sample]] <- cns[[sample]]
#}
#plot_recurrent_cnas(cns.brca1,"brcaWT")
