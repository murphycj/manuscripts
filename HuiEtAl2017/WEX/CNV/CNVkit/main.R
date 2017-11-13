source("lib.R")
library("xlsx")

mutual_exclusivity <- function() {
  cna_data <- read.csv("gene-sample-cn.csv",header=T,row.names=1,check.names=F)
  cna_data <- cna_data[,((cna_data["Met",]>3) | (cna_data["Yap1",]>3))]
  tmp <- matrix(c(
    length(colnames(cna_data[,(cna_data["Met",]<=3) & (cna_data["Yap1",]<=3)])),
    length(colnames(cna_data[,(cna_data["Met",]<=3) & (cna_data["Yap1",]>3)])),
    length(colnames(cna_data[,(cna_data["Met",]>3) & (cna_data["Yap1",]<=3)])),
    length(colnames(cna_data[,(cna_data["Met",]>3) & (cna_data["Yap1",]>3)]))
    ),nrow=2)
}

summarize_by_gene <- function() {
  cns <- load_cnvkit_data_cns(filterY=T)
  cnr <- load_cnvkit_data_cnr(filterY=T)

  genes <- c()

  for (sample in names(cns)) {
    genes <- c(genes, unlist(strsplit(as.character(cns[[sample]]$gene),",")))
  }
  genes <- unique(genes)

  data <- matrix(NA,nrow=length(genes),ncol=length(cns))
  row.names(data) <- genes
  colnames(data) <- names(cns)

  data_binned <- matrix(NA,nrow=length(genes),ncol=length(cns))
  row.names(data_binned) <- genes
  colnames(data_binned) <- names(cns)

  for (sample in names(cns)) {
    cns_sample <- cns[[sample]]
    cnr_sample <- cnr[[sample]]

    previous_genes <- c()
    previous_cn <- 2
    previous_log2 <- 0

    for (i in 1:nrow(cns_sample)) {
      gg <- unlist(strsplit(as.character(cns_sample[i,"gene"]),","))
      cn <- cns_sample[i,"cn"]
      sample.log2 <- cns_sample[i,"log2"]

      common <- intersect(previous_genes,gg)
      if (length(common) > 0) {
        gg <- gg[!(gg %in% common)]

        tmp <-c(cn,previous_cn)
        cn_to_use <- max(tmp[which(abs(tmp-2)==max(abs(tmp-2)))])
        data_binned[gg,sample] <- cns_sample[i,"cn"]
        data_binned[common,sample] <- cn_to_use

        tmp <-c(sample.log2,previous_log2)
        log2_to_use <- max(tmp[which(abs(tmp-0)==max(abs(tmp-0)))])
        data[gg,sample] <- cns_sample[i,"log2"]
        data[common,sample] <- log2_to_use
      } else {
        data[gg,sample] <- cns_sample[i,"log2"]
        data_binned[gg,sample] <- cns_sample[i,"cn"]
      }

      previous_cn <- cn
      previous_log2 <- sample.log2
      previous_genes <- gg
    }
  }
  samples <- read.xlsx("/Users/charlesmurphy/Desktop/Research/0914_hui/data/samples.xlsx","samples",row.names=1)

  data <- data[order(row.names(data)),]
  data <- data[,order(colnames(data))]
  data <- rbind(as.character(samples[colnames(data),"sample_name"]),data)
  data <- data[row.names(data)!="-",]
  tmp <- row.names(data)
  tmp[tmp==""]<-"sample_names"
  row.names(data) <- tmp
  write.csv(data,"gene-sample-log2.csv",quote=F)

  data_binned <- data_binned[order(row.names(data_binned)),]
  data_binned <- data_binned[,order(colnames(data_binned))]
  data_binned <- rbind(as.character(samples[colnames(data_binned),"sample_name"]),data_binned)
  data_binned <- data_binned[row.names(data_binned)!="-",]
  tmp <- row.names(data_binned)
  tmp[tmp==""]<-"sample_names"
  row.names(data_binned) <- tmp
  write.csv(data_binned,"gene-sample-cn.csv",quote=F)
}

summarize_by_gene()
