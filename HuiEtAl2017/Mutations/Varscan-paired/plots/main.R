library("beeswarm")
library(gplots)
library(xlsx)

plot_shared_mutations_by_primary <- function() {

  ps_n = list()
  ps_n[["1720"]]<-7
  ps_n[["1513"]]<-11
  ps_n[["1259"]]<-10
  ps_n[["1397"]]<-7
  ps_n[["1512"]]<-10
  ps_n[["1367"]]<-10
  ps_n[["1616"]]<-7
  ps_n[["1460"]]<-6
  ps_n[["1221"]]<-6
  ps_n[["1536"]]<-4
  ps_n[["1461"]]<-4
  ps_n[["1413"]]<-2
  ps_n[["1795"]]<-2
  ps_n[["1660"]]<-2
  ps_n[["1661"]]<-2
  ps_n[["1415"]]<-2
  ps_n[["1614"]]<-2

  data <- read.csv("shared_by_primary.csv",header=T,row.names=1,check.names=F)

  png("shared.byprimary.png",height=2500,width=2000)
  par(mfrow=c(4,5),cex=2.5,mar=c(2.5,2.5,2,2))

  for (i in colnames(data)) {

    temp <- c()
    for (j in row.names(data)) {
      temp <- c(temp,rep(j,data[j,i]))
    }
    temp <- as.numeric(temp)

    hist(
      temp,
      cex.lab=2.5,
      main=i,
      ylab="",
      xlab=""
    )
  }
  dev.off()
}

plot_mutation_types_by_primary_one_figure <- function() {

  samples <- read.xlsx("../../../../data/samples.xlsx","samples",row.names=1)

  variant_types <- read.csv("mutation_types.csv",header=T,row.names=1)

  primaries = as.character(unique(samples$mouse))

  png("mutation.byprimary.png",height=2500,width=2000)
  par(mfrow=c(4,5),cex=1.5)

  for (i in primaries) {
    samples.p = samples[samples$mouse==i,]
    samples.p <- samples.p[row.names(samples.p) %in% row.names(variant_types),]

    if (nrow(samples.p)<2) {
      next
    }

    variant_types.p <- variant_types[row.names(samples.p),]
    variant_types.p <- variant_types.p[order(variant_types.p$total,decreasing=T),]

    barplot(
      t(variant_types.p[,2:6]),
      las=2,
      names.arg=row.names(data),
      col=c("#cc0000","#ffff00","#0000ff","#cc00cc","grey"),
      cex.lab=1.5,
      ylim=c(0,130),
      main=i
    )
  }
  plot.new()
  plot.new()
  plot.new()
  par(xpd=TRUE)
  legend(
    "bottomright",
    c("Missense","Spice site","Frameshift","Truncating","Synonymous"),
    col=c("#cc0000","#ffff00","#0000ff","#cc00cc","grey"),
    pch=15,
    cex=2.3
  )
  dev.off()
}

plot_mutation_types_by_primary <- function() {

  samples <- read.xlsx("../../../../data/samples.xlsx","samples",row.names=1)

  variant_types <- read.csv("mutation_types.csv",header=T,row.names=1)

  primaries = as.character(unique(samples$mouse))

  for (i in primaries) {
    samples.p = samples[samples$mouse==i,]
    samples.p <- samples.p[row.names(samples.p) %in% row.names(variant_types),]

    if (nrow(samples.p)<2) {
      next
    }

    variant_types.p <- variant_types[row.names(samples.p),]
    variant_types.p <- variant_types.p[order(variant_types.p$total,decreasing=T),]

    pdf(paste("mutation.",i,".pdf",sep=""),width=14,height=6)
    par(mar=c(10,6,4,4),cex=0.9)
    barplot(
      t(variant_types.p[,2:6]),
      las=2,
      names.arg=row.names(data),
      col=c("#cc0000","#ffff00","#0000ff","#cc00cc","grey"),
      cex.lab=1.5,
      ylim=c(0,130)
    )
    legend(
      "topright",
      c("Missense","Spice site","Frameshift","Truncating","Synonymous"),
      col=c("#cc0000","#ffff00","#0000ff","#cc00cc","grey"),
      pch=15,
      cex=1.5
    )
    mtext("Count", side=2, line=4,cex=1.5)
    dev.off()
  }
}

plot_mutation_types <- function() {

  variant_types <- read.csv("mutation_types.csv",header=T,row.names=1)

  samples <- read.xlsx("/Users/charlesmurphy/Desktop/Research/0914_hui/data/samples.xlsx","samples")
  samples <- samples[as.character(samples$sample_ID) %in% row.names(variant_types),]

  variant_types <- variant_types[as.character(samples$sample_ID),]

  new_sample_names <- as.character(samples$sample_name)
  new_sample_names <- gsub("parental","P.",new_sample_names)
  row.names(variant_types) <- new_sample_names

  variant_types <- variant_types[order(variant_types$total,decreasing=T),]

  pdf('mutation.types.pdf',width=14,height=6)
  par(mar=c(13,4,4,4),cex=0.5)
  barplot(
    t(variant_types[,2:7]),
    las=2,
    names.arg=row.names(variant_types),
    col=c("#cc0000","#ffff00","#0000ff","black","#cc00cc","grey"),
    cex.lab=1.2
  )
  legend(
    "topright",
    c("Missense","Spice site","Frameshift","Truncating","Indel","Synonymous"),
    col=c("#cc0000","#ffff00","#0000ff","black","#cc00cc","grey"),
    pch=15,
    cex=1.5
  )
  dev.off()

  pdf('mutation.types.nolabel.pdf',width=14,height=6)
  par(mar=c(4,4,4,4),cex=0.9)
  barplot(
    t(variant_types[,2:7]),
    las=2,
    names.arg=rep("",nrow(variant_types)),
    col=c("#cc0000","#ffff00","#0000ff","black","#cc00cc","grey"),
    cex.lab=1.2
  )
  legend(
    "topright",
    c("Missense","Spice site","Frameshift","Truncating","Indel","Synonymous"),
    col=c("#cc0000","#ffff00","#0000ff","black","#cc00cc","grey"),
    pch=15,
    cex=1.5
  )
  dev.off()
}


plot_mutation_types()
#plot_by_genotype()
#plot_by_genotype_indel()
#plot_mutation_types_by_primary()
#plot_mutation_types_by_primary_one_figure()
#plot_shared_mutations_by_primary()
