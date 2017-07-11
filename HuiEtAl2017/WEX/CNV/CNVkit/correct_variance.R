library("argparse")

GC_INTERVALS <- seq(0,1,0.01)

median_correction <- function(data,reference) {
  variances <- c()
  for (i in 1:length(GC_INTERVALS)) {
    temp <- reference[(reference$gc>=GC_INTERVALS[i-1]) & (reference$gc<GC_INTERVALS[i]),]
    if (nrow(temp)==0) {
      next
    }

    common <- intersect(row.names(temp),row.names(data))
    if (length(common)==0) {
      next
    }

    temp.data <- data[common,]
    variances <- c(variances,mad(temp.data$log2))
  }
  return(median(variances))
}

new_correction <- function(data,reference) {
  log2s <- data$log2
  log2s <- log2s[order(log2s)]
  quants <- quantile(log2s,c(0.25,0.75))
  data <- data[(data$log2 >= quants[["25%"]]) & (data$log2 <= quants[["75%"]]),]

  variances <- c()
  for (i in 1:length(GC_INTERVALS)) {
    temp <- reference[(reference$gc>=GC_INTERVALS[i-1]) & (reference$gc<GC_INTERVALS[i]),]
    if (nrow(temp)==0) {
      next
    }

    common <- intersect(row.names(temp),row.names(data))
    if (length(common)==0) {
      next
    }

    temp.data <- data[common,]
    variances <- c(variances,mad(temp.data$log2))
  }
  return(median(variances))
}

main <- function(args) {

  data <- read.table(args$cnr,sep="\t",header=T,check.names=F)
  row.names(data) <- paste(data$chromosome,data$start,data$end,sep="-")

  reference <- read.table(args$reference,sep="\t",header=T,check.names=F)
  row.names(reference) <- paste(reference$chromosome,reference$start,reference$end,sep="-")

  common <- intersect(row.names(data),row.names(reference))

  png(paste(args$prefix,"_log2Coverage-GC_uncorrected.png",sep=""))
  smoothScatter(
    as.numeric(reference[common,"gc"]),
    as.numeric(data[common,"log2"]),
    xlab="GC content",
    ylab="log2(tumor/liver)",
    pch=".",
    xlim=c(0,1),
    ylim=c(args$ymin,args$ymax)
  )
  dev.off()

  data.corrected <- data.frame()

  if (args$correction=="median") {
    correction <- median_correction(data,reference)
  } else if (args$correction=="new"){
    correction <- new_correction(data,reference)
  } else {
    correction <- 0.3005
  }

  for (i in 1:length(GC_INTERVALS)) {

    temp <- reference[(reference$gc>=GC_INTERVALS[i-1]) & (reference$gc<GC_INTERVALS[i]),]
    if (nrow(temp)==0) {
      next
    }

    common <- intersect(row.names(temp),row.names(data))
    if (length(common)==0) {
      next
    }

    temp.data <- data[common,]
    temp <- mad(temp.data$log2)
    temp <- correction/temp
    temp.data[,"log2"] <- temp * temp.data[,"log2"]

    data.corrected <- rbind(data.corrected,temp.data)
  }

  common <- intersect(row.names(data),row.names(data.corrected))

  data.corrected <- data.corrected[common,]
  data.corrected <- data.corrected[!is.infinite(data.corrected$log2),]
  data.corrected <- data.corrected[!is.na(data.corrected$log2),]

  png(paste(args$prefix,"_log2Coverage-GC_corrected.png",sep=""))
  smoothScatter(
    as.numeric(reference[common,"gc"]),
    as.numeric(data.corrected[common,"log2"]),
    xlab="GC content",
    ylab="log2(tumor/liver)",
    pch=".",
    xlim=c(0,1),
    ylim=c(args$ymin,args$ymax)
  )
  dev.off()

  write.table(
    data.corrected,
    paste(args$prefix,".corrected.cnr",sep=""),
    row.names=F,
    quote=F,
    sep="\t"
    )
}

parser <- ArgumentParser()

parser$add_argument(
  "-cnr",
  type="character",help="The .cnr file."
)
parser$add_argument(
  "-reference",
  type="character",help="The file."
)
parser$add_argument(
  "-prefix",
  type="character",
  help="The new output prefix."
)
parser$add_argument(
  "-ymin",
  type="integer",
  default=-1,
  help="Min Y"
)
parser$add_argument(
  "-ymax",
  type="integer",
  default=-1,
  help="Max Y"
)
parser$add_argument(
  "-correction",
  type="character",
  default=-1,
  help="Correction method"
)
args <- parser$parse_args()

main(args=args)
