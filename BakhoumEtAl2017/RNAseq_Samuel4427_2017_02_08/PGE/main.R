library(circlize)
library(xlsx)

plot_data <- function(prefix,de_genes,pge.results) {

  gene_pos <- read.csv("ensembl89.csv",row.names=1,header=F,colClasses=c("character","character","numeric","numeric","character"))
  gene_pos <- gene_pos[row.names(de_genes),]

  colors <- rep("red",length(de_genes$log2FoldChange))
  colors[de_genes$log2FoldChange<0]<-"blue"

  de_data <- data.frame(
    chr=gene_pos[,1],
    start=gene_pos[,2],
    end=gene_pos[,3],
    value1=de_genes$log2FoldChange,
    color=rep("red",length(de_genes$log2FoldChange)),
    stringsAsFactors=FALSE
  )
  de_data[de_data$value1<0,"color"]<-"blue"
  de_data$chr <- unlist(lapply(de_data$chr,function(x) return(paste("chr",x,sep=""))))
  de_data <- de_data[(de_data$chr!="chrMT") & (de_data$chr!="chrNA"),]


  pdf(paste(prefix,"circos.pdf",sep=""))
  circos.par("gap.degree"=c(rep(1,23),7),"start.degree"=90)
  circos.initializeWithIdeogram(
    cytoband="/Users/charlesmurphy/Desktop/Research/data/UCSC/hg38/cytoBand-canonical.txt"
  )
  circos.genomicTrackPlotRegion(
    de_data[,c("chr","start","end","value1")]
  )
  circos.yaxis(labels.cex = 0.6,sector.index="chr1")
  for (chrom in unique(de_data$chr)) {
    tmp <- de_data[de_data$chr==chrom,]
    circos.genomicPoints(region=tmp[,c("start","end")],value=data.frame(value=tmp$value1),col=tmp$color,sector.index=chrom,pch=".",cex=2.5)
  }
  circos.genomicTrackPlotRegion(
    data.frame(
      chr=pge.results$chrom,
      start=pge.results$start,
      end=pge.results$end,
      stringsAsFactors=FALSE
    ),
    ylim=c(0,max(-log10(pge.results$padj))*1.1)
  )
  circos.yaxis(labels.cex = 0.6,sector.index="chr1")
  for (chrom in unique(pge.results$chrom)) {
    tmp <- pge.results[pge.results$chrom==chrom,]
    colors <- as.character(tmp$cna)
    colors[colors=="amplification"] <- "red"
    colors[colors=="deletion"] <- "blue"
    circos.genomicLines(
      region=tmp[,c("start","end")],
      value=data.frame(value=-log10(tmp$padj)),
      type = "segment",
      col=colors,
      sector.index=chrom,
      lwd = 4
    )
  }
  dev.off()
  circos.clear()
}

comparisons <- c('BoneMet_vs_BrainMet')
padj <- c(0.05,0.1,0.25,1.0)
log2s <- c(0,0.25,0.5)
pvals <- c(0.05)

for (comparison in comparisons) {
  for (pv in pvals) {
    for (pp in padj) {
      for (ll in log2s) {

        de_genes <- read.csv(
          paste("/Users/charlesmurphy/Desktop/Research/0816_sam/results/RNAseq_Samuel4427_2017_02_08/DESeq2_ensembl/",comparison,"/",comparison,"_results.csv",sep=""),
          header=T,
          check.names=F,
          row.names=1,
          colClasses=c("character","numeric","numeric","numeric","numeric","numeric","numeric")
        )
        de_genes <- de_genes[(!is.na(de_genes$padj)) & (de_genes$padj<=pp) & (de_genes$pvalue<=pv) & (abs(de_genes$log2FoldChange)>=ll),]

        if (!file.exists(paste("./",comparison,"/pval",pv,"_padj",pp,"_log2FC",ll,"/",comparison,"_results.xlsx",sep=""))) {
          next
        }

        pge.results <- read.xlsx(
          paste("./",comparison,"/pval",pv,"_padj",pp,"_log2FC",ll,"/",comparison,"_results.xlsx",sep=""),
          "Sheet1",
          colClasses=c("character",rep("numeric",4),"character","numeric","numeric","character")
        )
        pge.results <- pge.results[(pge.results$num_differentially_expressed_genes>3) & (pge.results$padj<=0.01),]
        pge.results$chrom <- unlist(lapply(pge.results$chrom,function(x) return(paste("chr",x,sep=""))))

        plot_data(
          paste("./",comparison,"/pval",pv,"_padj",pp,"_log2FC",ll,"/",sep=""),
          de_genes,
          pge.results
        )
      }
    }
  }
}
