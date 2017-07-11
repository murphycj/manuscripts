library(xlsx)
samples <- read.xlsx("/Users/charlesmurphy/Desktop/Research/0914_hui/data/samples.xlsx","samples",colClasses="character")
samples <- samples[(samples$tissue=="primary") & (samples$sequencing_type=="WEX"),]

initial_brca1_samples <- c(
  "HL-WES-01",
  "HL-WES-02",
  "HL-WES-03",
  "HL-WES-04",
  "HL-WES-05",
  "HL-WES-06",
  "HL-WES-07",
  "HL-WES-08",
  "HL-WES-09",
  "HL-WES-10",
  "HL-WES-11"
  )

load_cnvkit_data_cns <- function(filterY=T) {

  data <- list()

  for (i in samples$sample_ID) {
    if (i %in% initial_brca1_samples) {
      tmp <- read.csv(paste("/Users/charlesmurphy/Desktop/Research/0914_hui/results/WEX/CNV/CNVkit/",i,"/",i,".call.cns",sep=""),sep="\t",header=T)
      if (filterY) {
        tmp <- tmp[tmp$chromosome!="Y",]
      }
      data[[i]] <- tmp
    } else {
      tmp <- read.csv(paste("/Users/charlesmurphy/Desktop/Research/0914_hui/results/WEX/CNV/CNVkit/",i,"/",i,".corrected.call.cns",sep=""),sep="\t",header=T)
      if (filterY) {
        tmp <- tmp[tmp$chromosome!="Y",]
      }
      data[[i]] <- tmp
    }
  }

  return(data)
}

load_cnvkit_data_cnr <- function(filterY=T) {

  data <- list()

  for (i in samples$sample_ID) {
    if (i %in% initial_brca1_samples) {
      tmp <- read.csv(paste("/Users/charlesmurphy/Desktop/Research/0914_hui/results/WEX/CNV/CNVkit/",i,"/",i,".call.cnr",sep=""),sep="\t",header=T)
      if (filterY) {
        tmp <- tmp[tmp$chromosome!="Y",]
      }
      data[[i]] <- tmp
    } else {
      tmp <- read.csv(paste("/Users/charlesmurphy/Desktop/Research/0914_hui/results/WEX/CNV/CNVkit/",i,"/",i,".corrected.call.cnr",sep=""),sep="\t",header=T)
      if (filterY) {
        tmp <- tmp[tmp$chromosome!="Y",]
      }
      data[[i]] <- tmp
    }
  }

  return(data)
}
