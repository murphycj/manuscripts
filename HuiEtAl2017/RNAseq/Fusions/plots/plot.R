library(gplots)
library(ComplexHeatmap)
library(ggplot2)
options( stringsAsFactors = FALSE)
library(OmicCircos)
library("devtools")
library(xlsx)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

fusions_by_sample <- function(fusions,samples) {

  samples_names <- unique(as.character(fusions[,"sample"]))

  for (s in samples_names) {
    fusions.formated <- c()
    fusions.s <- fusions[fusions["sample"]==s,]
    if (nrow(fusions.s)==0) {
      next
    }
    for (i in 1:nrow(fusions.s)) {
      position1 <- as.character(fusions.s[i,"Fusion_point_for_gene_1.5end_fusion_partner."])
      chr1 <- strsplit(position1,":")[[1]][1]
      pos1 <- as.numeric(strsplit(position1,":")[[1]][2])

      position2 <- as.character(fusions.s[i,"Fusion_point_for_gene_2.3end_fusion_partner."])
      chr2 <- strsplit(position2,":")[[1]][1]
      pos2 <- as.numeric(strsplit(position2,":")[[1]][2])

      gene1 <- fusions.s[i,"Gene_1_symbol.5end_fusion_partner."]
      gene2 <- fusions.s[i,"Gene_2_symbol.3end_fusion_partner."]

      fusions.formated <- rbind(fusions.formated,c(chr1,pos1,gene1,chr2,pos2,gene2))
    }

    labels <- data.frame(rbind(fusions.formated[,1:3],fusions.formated[,4:6]))
    names(labels) <- c("chr","po","Gene")


    pdf(paste("./circos_plot/by_sample/011817.",s,".circos.pdf",sep=""))
    par(mar=c(2,2,2,2))
    plot(c(1,1000),c(1,1000),type="n",axes=F,xlab="",ylab="",main="")
    circos(R=400,xc=500,yc=500,type="chr",cir="mm10",print.chr.lab=T,W=15,scale=F,cex=15)
    circos(R=350,xc=500,yc=500,cir="mm10",W=240,mapping=fusions.formated,type="link",lwd=2)
    circos(R=445,xc=500,yc=500,cir="mm10",W=20,mapping=labels,type="label",side="out",cex=0.5)
    dev.off()
  }

}

plot_heatmap <- function(fusions,samples) {

  row.names(samples) <- samples$sample_name
  samples <- samples[(samples$tissue=="primary") | (samples$control==1),]
  common <- intersect(as.character(samples$sample_name),colnames(fusions))
  samples <- samples[common,]
  fusions <- fusions[,c("Gene-5prime","Gene-3prime","Count",common)]
  samples <- samples[order(samples$brca1),]

  BRCA1 <- as.character(samples$brca1)
  BRCA1[BRCA1=="WT"] <- "#ffcccc"
  BRCA1[BRCA1=="flox/flox"] <- "#ff0000"
  BRCA1[is.na(BRCA1)] <- "#ff0000"
  clab <- cbind(BRCA1)
  colnames(clab)<-c("")
  row.names(fusions) <- paste(fusions["Gene-5prime"][,1],fusions["Gene-3prime"][,1],sep="-")
  fusions <- fusions[,-c(1,2,3)]
  fusions <- fusions[apply(fusions,1,sum)>0,]

  new_sample_names=c()
  for (i in colnames(fusions)) {
    i = gsub("parental ","P. ",i)
    i = gsub("\\*","",i)
    i = gsub("KO","",i)
    i = gsub("HET","",i)
    i = gsub("WhiteL","W",i)
    i = gsub("WT","",i)
    new_sample_names <- c(new_sample_names,i)
  }
  colnames(fusions) <- new_sample_names

  pdf("fusions-global.pdf",width=8,height=12)
  par(mar=c(1,1,1,1),family="Times",cex.main=2.)
  heatmap.3(
    fusions,
    trace="none",
    scale="none",
  	key=F,
    margin=c(3,3),
    cexRow=1,cex.lab=2,
    cexCol=1,
    col= c("white","black"),
  	ColSideColors=clab,
  	Rowv=F,
    Colv=F,
    labRow=NA,
    labCol=NA,xlab="Samples",ylab="Gene fusions",
    dendrogram="none",
    main="",
  	sepcolor="#e6e6e6",
  	colsep=1:ncol(fusions),
  	rowsep=1:nrow(fusions),
  	sepwidth=c(0.0001,0.0001),
    lwid=c(1,6),
    lhei=c(1,10)
  )
  dev.off()


  samples <- read.xlsx("~/Desktop/Research/0914_hui/data/samples.xlsx","samples")
  row.names(samples) <- samples$sample_name
  samples <- samples[samples["sequencing_type"]=="RNAseq",]
  samples <- samples[((samples$tissue=="primary") | (samples$control==1)) & (samples$repeat.==0),]
  #samples <- samples[order(samples$brca1),]

  fusions <- read.csv("../fusions-table-noNormal.noReadthrough.n2.csv",header=T,check.names=F)
  row.names(fusions) <- paste(fusions["Gene-5prime"][,1],fusions["Gene-3prime"][,1],sep="-")

  common <- intersect(as.character(samples$sample_name),colnames(fusions))
  samples <- samples[common,]
  fusions <- fusions[,c("Gene-5prime","Gene-3prime","Count",common)]
  fusions <- fusions[,-c(1,2,3)]
  fusions <- fusions[apply(fusions,1,sum)>0,]
  fusion_count <- apply(fusions,2,sum)
  fusions <- fusions[,names(fusion_count[order(fusion_count,decreasing=T)])]
  fusion_count <- fusion_count[order(fusion_count,decreasing=T)]

  #fusion_count.brca1fl <- fusion_count[samples$brca1=="flox/flox"]
  #fusion_count.brca1wt <- fusion_count[samples$brca1=="WT"]
  #fusions.brca1fl <- fusions[,samples[samples$brca1=="flox/flox","sample_name"]]
  #fusions.brca1wt <- fusions[,samples[samples$brca1=="WT","sample_name"]]

  #fusion_count.brca1fl <- fusion_count.brca1fl[order(fusion_count.brca1fl,decreasing=T)]
  #fusion_count.brca1wt <- fusion_count.brca1wt[order(fusion_count.brca1wt,decreasing=T)]

  #fusions.brca1fl <- fusions.brca1fl[,names(fusion_count.brca1fl)]
  #fusions.brca1wt <- fusions.brca1wt[,names(fusion_count.brca1wt)]

  #fusion_count <- c(fusion_count.brca1fl,fusion_count.brca1wt)
  #fusions <- cbind(fusions.brca1fl,fusions.brca1wt)

  fusions2 <- data.frame()
  tmp <- c()
  for (sample in colnames(fusions)) {
    tmp2 <- setdiff(row.names(fusions[fusions[,sample]==1,]),tmp)
    tmp <- c(tmp,tmp2)
    fusions2 <- rbind(fusions2,fusions[tmp2,])
  }


  ha1 = HeatmapAnnotation(
    barplot = anno_barplot(
      data.frame(fusion_count=fusion_count),
      bar_width = 1,
      gp = gpar(col = NA, fill = "black"),
      border = FALSE, axis = TRUE),
      height = unit(2,"cm")
    )

  tmp <- as.character(samples[names(fusion_count),"brca1"])
  tmp[tmp=="flox/flox"] <- "BRCA1-DEL"
  tmp[tmp=="WT"] <- "BRCA1-WT"
  BRCA1 <- data.frame(
    BRCA1=tmp
    )
  ha2 = HeatmapAnnotation(
    df = BRCA1,
    col=list(BRCA1 = c("BRCA1-DEL" =  "#205247", "BRCA1-WT" = "#72fa4c")),
    name="BRCA1 status",
    annotation_legend_param=list(nrow=1)
  )

  pdf("fusions-global2.pdf",width=10,height=10)
  par(mar=c(2,2,2,2),family="Times",cex.main=2.)
  ht1 <- Heatmap(
    as.matrix(fusions2),
    col=c("white","black"),
    top_annotation=ha1,
    bottom_annotation=ha2,
    cluster_columns=FALSE,
    cluster_rows=FALSE,
    show_row_dend = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE,
    rect_gp = gpar(col = "gray", lwd = 1),
    column_title="Tumors",
    column_title_side="bottom",
    row_title="Gene fusions",
    row_title_side="left",
    name="Gene fusion",
    heatmap_legend_param = list(at = c(0, 1), labels = c("No fusion", "Fusion"),nrow=1,grid_border="black")
  )
  ht_global_opt(
    heatmap_legend_title_gp = gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_labels_gp = gpar(fontsize = 12)
  )
  draw(
    ht1,
    padding = unit(c(2, 20, 2, 2), "mm"),
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
  )
  decorate_annotation(
    "barplot", {
      grid.text("Total\nfusions", -0.055, 0.5,rot=90)
  })
  dev.off()
}

samples <- read.xlsx("~/Desktop/Research/0914_hui/data/samples.xlsx","samples")
samples <- samples[samples["sequencing_type"]=="RNAseq",]
fusions <- read.csv("../fusions-table-noNormal.noReadthrough.n2.csv",header=T,check.names=F)
fusions<-fusions[order(fusions$Count,decreasing=T),]

plot_heatmap(fusions=fusions,samples)
