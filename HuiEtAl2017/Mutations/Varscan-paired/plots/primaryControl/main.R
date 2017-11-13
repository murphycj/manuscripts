library("beeswarm")
library(gplots)
library(xlsx)
library(RColorBrewer)

plot_mutation_types <- function() {

  variant_types <- read.csv("mutation_types.csv",header=T,row.names=1)

  samples <- read.xlsx("/Users/charlesmurphy/Desktop/Research/0914_hui/data/samples.xlsx","samples")
  samples <- samples[as.character(samples$sample_ID) %in% row.names(variant_types),]

  variant_types <- variant_types[as.character(samples$sample_ID),]

  new_sample_names <- as.character(samples$sample_name)
  new_sample_names <- gsub("parental","P.",new_sample_names)
  row.names(variant_types) <- new_sample_names

  variant_types <- variant_types[order(variant_types$total,decreasing=T),]

  colors <- rev(brewer.pal(ncol(variant_types),"Spectral"))
  types <- c("Missense","Synonymous","Indel","Frameshift","Nonsense","Spice site")

  pdf('mutation.types.pdf',width=14,height=6)
  par(mar=c(13,4,4,4),cex=0.5)
  barplot(
    t(variant_types[,2:7]),
    las=2,
    names.arg=row.names(variant_types),
    col=colors,
    cex.lab=1.2
  )
  legend(
    "topright",
    types,
    col=colors,
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
    col=colors,
    cex.lab=1.2
  )
  legend(
    "topright",
    types,
    col=colors,
    pch=15,
    cex=1.5
  )
  dev.off()
}


plot_mutation_types()
