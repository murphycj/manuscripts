library(gplots)
library(ComplexHeatmap)
library(ggplot2)
options( stringsAsFactors = FALSE)
library(OmicCircos)
library("devtools")
library(xlsx)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")


plot_heatmap <- function(fusions,samples) {


  fusions <- read.csv("../fusions-table-noNormal.noReadthrough.n2.csv",header=T,check.names=F)
  row.names(fusions) <- paste(fusions["Gene-5prime"][,1],fusions["Gene-3prime"][,1],sep="-")

  common <- intersect(as.character(samples$sample_ID),colnames(fusions))
  samples <- samples[common,]
  fusions <- fusions[,c("Gene-5prime","Gene-3prime","Count",common)]
  fusions <- fusions[,-c(1,2,3)]
  fusions <- fusions[apply(fusions,1,sum)>0,]
  fusion_count <- apply(fusions,2,sum)
  fusions <- fusions[,names(fusion_count[order(fusion_count,decreasing=T)])]
  fusion_count <- fusion_count[order(fusion_count,decreasing=T)]

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
      gp = gpar(
        col = NA,
        fill = "black"
      ),
      border = FALSE,
      axis = TRUE,
      row_title_gp = gpar(
        fontsize = 16,
        fontface = "bold"
      )
    ),
    height = unit(2,"cm"),
    annotation_legend_param=list(
      title_gp = gpar(
        fontsize = 16,
        fontface = "bold"
      ),
      labels_gp = gpar(
        fontsize = 16,
        fontface = "bold"
      )
    )
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
    annotation_legend_param=list(
      nrow=1,
      title_gp = gpar(fontsize = 18, fontface = "bold"),
      labels_gp = gpar(fontsize = 18)
    )
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
    heatmap_legend_title_gp = gpar(fontsize = 18, fontface = "bold"),
    heatmap_legend_labels_gp = gpar(fontsize = 18),
    heatmap_row_title_gp = gpar(fontsize = 18, fontface = "bold"),
    heatmap_column_title_gp = gpar(fontsize = 18, fontface = "bold")
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
samples <- samples[(samples['byl']==0) & (samples['bkm']==0) & (samples['bmn']==0),]
samples <- samples[samples["sequencing_type"]=="RNAseq",]
samples <- samples[(samples$tissue=="primary") | (samples$control==1),]
row.names(samples) <- samples$sample_ID

plot_heatmap(fusions=fusions,samples)
