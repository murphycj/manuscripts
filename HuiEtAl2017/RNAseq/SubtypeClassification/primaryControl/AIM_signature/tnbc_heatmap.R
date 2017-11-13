library(ComplexHeatmap)
library(xlsx)
library("limma")

samples <- read.xlsx("/Users/charlesmurphy/Desktop/Research/0914_hui/data/samples.xlsx","samples",row.names=1)
samples <- samples[(samples['byl']==0) & (samples['bkm']==0) & (samples['bmn']==0),]
samples <- samples[samples$sequencing_type=="RNAseq",]
samples <- samples[samples$tissue %in% c("primary","implant"),]
samples <- samples[row.names(samples)!="HL23",]
samples <- samples[!is.na(samples$brca1),]


#load mouse and human exrpession data

genes <- read.csv("genes.csv")[,1]

fpkm <- read.csv(
  "/Users/charlesmurphy/Desktop/Research/0914_hui/results/RNAseq/Cufflinks/genes-humanEntrez-fpkm.csv",
  row.names=1,
  header=T,
  check.names=F
)
common <- intersect(colnames(fpkm),row.names(samples))
samples <- samples[common,]
fpkm <- fpkm[,common]

load("~/Desktop/Research/data/GDAC/BRCA/gdac.broadinstitute.org_BRCA.Merge_rnaseq__illuminahiseq_rnaseq__unc_edu__Level_3__gene_expression__data.Level_3.2016012800.0.0/process_data/rpkm.entrez.RData")

rpkm.entrez <- rpkm.entrez[,grep("-11A",colnames(rpkm.entrez),invert=T)]
rpkm.entrez <- rpkm.entrez[,grep("-11B",colnames(rpkm.entrez),invert=T)]
rpkm.entrez <- rpkm.entrez[,grep("-06A",colnames(rpkm.entrez),invert=T)]
rpkm.entrez <- rpkm.entrez[,colnames(rpkm.entrez)!="TCGA-E2-A15A-06A"]
rpkm.entrez <- rpkm.entrez[,colnames(rpkm.entrez)!="TCGA-E2-A15K-06A"]
colnames(rpkm.entrez) <- unlist(lapply(colnames(rpkm.entrez),function(x) return(paste(strsplit(x,"-")[[1]][1:3],collapse="-"))))


common <- intersect(as.character(genes),intersect(row.names(fpkm),row.names(rpkm.entrez)))
fpkm <- fpkm[common,]
rpkm.entrez <- rpkm.entrez[common,]

all_data <- cbind(fpkm,rpkm.entrez)
all_data <- normalizeQuantiles(all_data)
fpkm <- all_data[,1:ncol(fpkm)]
rpkm.entrez <- all_data[,(ncol(fpkm)+1):ncol(rpkm.entrez)]

write.table(common,"mouse-human-common-genes.csv",quote=F,row.names=F,col.names=F)

#read the subtpye information

subtypes <- read.csv("/Users/charlesmurphy/Desktop/Research/data/papers/2012_comprehensiveMolecularPortraitsOfHumanBreastCancer/subtype.csv",header=T,check.names=F,colClasses="character")
subtypes <- subtypes[!is.na(subtypes[,"PAM50 mRNA"]),]
subtypes <- subtypes[!is.na(subtypes[,"ER Status"]),]
subtypes <- subtypes[!is.na(subtypes[,"PR Status"]),]
subtypes <- subtypes[!is.na(subtypes[,"HER2 Final Status"]),]
subtypes <- subtypes[(subtypes[,"ER Status"]=="Positive") | (subtypes[,"ER Status"]=="Negative"),]
subtypes <- subtypes[(subtypes[,"PR Status"]=="Positive") | (subtypes[,"PR Status"]=="Negative"),]
subtypes <- subtypes[(subtypes[,"HER2 Final Status"]=="Positive") | (subtypes[,"HER2 Final Status"]=="Negative"),]

common <- intersect(colnames(rpkm.entrez),subtypes[,"Complete TCGA ID"])
rpkm.entrez <- rpkm.entrez[,common]
subtypes <- subtypes[subtypes[,"Complete TCGA ID"] %in% common,]
rpkm.entrez <- rpkm.entrez[,subtypes[,"Complete TCGA ID"]]

fpkm <- log2(fpkm + 0.1)
rpkm.entrez <- log2(rpkm.entrez + 0.1)

subtypes_class <- list()
subtypes_class[["TNBC"]] <- as.character(subtypes[subtypes[,"ER Status"]=="Negative" & subtypes[,"PR Status"]=="Negative" & subtypes[,"HER2 Final Status"]=="Negative","Complete TCGA ID"])
subtypes_class[["Basal-like"]] <- setdiff(
                                    as.character(subtypes[subtypes[,"PAM50 mRNA"]=="Basal-like","Complete TCGA ID"]),
                                    subtypes_class[["TNBC"]]
                                 )
subtypes_class[["Luminal A"]] <- setdiff(
                                    as.character(subtypes[subtypes[,"PAM50 mRNA"]=="Luminal A","Complete TCGA ID"]),
                                    subtypes_class[["TNBC"]]
                                 )
subtypes_class[["Luminal B"]] <- setdiff(
                                    as.character(subtypes[subtypes[,"PAM50 mRNA"]=="Luminal B","Complete TCGA ID"]),
                                    subtypes_class[["TNBC"]]
                                 )
subtypes_class[["Normal-like"]] <- setdiff(
                                    as.character(subtypes[subtypes[,"PAM50 mRNA"]=="Normal-like","Complete TCGA ID"]),
                                    subtypes_class[["TNBC"]]
                                 )
subtypes_class[["HER2-enriched"]] <- setdiff(
                                    as.character(subtypes[subtypes[,"PAM50 mRNA"]=="HER2-enriched","Complete TCGA ID"]),
                                    subtypes_class[["TNBC"]]
                                 )


data <- as.matrix(cbind(fpkm,rpkm.entrez))

brca1 <- as.character(samples$brca1)
brca1[brca1=="flox/flox"]<-"BRCA1-DEL"
brca1[brca1=="WT"]<-"BRCA1-WT"

subtype <- colnames(rpkm.entrez)
subtype[subtype %in% subtypes_class[["TNBC"]]]<-"TNBC"
subtype[subtype %in% subtypes_class[["Basal-like"]]]<-"Basal-like"
subtype[subtype %in% subtypes_class[["Luminal A"]]]<-"Luminal A"
subtype[subtype %in% subtypes_class[["Luminal B"]]]<-"Luminal B"
subtype[subtype %in% subtypes_class[["HER2-enriched"]]]<-"HER2-enriched"
subtype[subtype %in% subtypes_class[["Normal-like"]]]<-"Normal-like"
subtype <- c(
  brca1,
  subtype
)

mouse.markers <- read.csv("/Users/charlesmurphy/Desktop/Research/0914_hui/results/RNAseq/SubtypeClassification/primaryControl/Markers/mouse-marker-status-prediction.csv",header=T,row.names=1)


subtype <- data.frame(
  Subtype=subtype,
  PR_status=c(as.character(mouse.markers[colnames(fpkm),"PR_status"]),subtypes[,"PR Status"]),
  ER_status=c(as.character(mouse.markers[colnames(fpkm),"ER_status"]),subtypes[,"ER Status"]),
  HER2_status=c(as.character(mouse.markers[colnames(fpkm),"HER2_status"]),subtypes[,"HER2 Final Status"])
)

ha1 = HeatmapAnnotation(
  df = subtype,
  show_annotation_name=T,
  col = list(
    Subtype = c(
      "TNBC" =  "red",
      "Basal-like" = "black",
      "Luminal A" = "purple",
      "Luminal B" = "orange",
      "HER2-enriched" = "grey",
      "Normal-like" = "blue",
      "BRCA1-DEL" = "#00e600",
      "BRCA1-WT" = "#006600"),
    PR_status=c(
      "NA"="grey",
      "Negative"="white",
      "Positive"="black"
    ),
    ER_status=c(
      "NA"="grey",
      "Negative"="white",
      "Positive"="black"
    ),
    HER2_status=c(
      "NA"="grey",
      "Negative"="white",
      "Positive"="black"
    )
  ),
  annotation_name_gp = gpar(fontsize = 14),
  annotation_legend_param = list(
    Subtype = list(
      title="Subtype",
      title_gp = gpar(fontsize = 14),
      labels_gp = gpar(fontsize = 14)),
    PR_status = list(
      title="PR status",
      title_gp = gpar(fontsize = 14),
      labels_gp = gpar(fontsize = 14)),
    ER_status = list(
      title="ER status",
      title_gp = gpar(fontsize = 14),
      labels_gp = gpar(fontsize = 14)),
    HER2_status = list(
      title="HER2 status",
      title_gp = gpar(fontsize = 14),
      labels_gp = gpar(fontsize = 14)))
)
png("TNBC.png",width=1000,height=1000,family="Times")
Heatmap(
  data,
  show_row_names=F,
  show_column_names=F,
  top_annotation = ha1,
  column_title="Samples",
  column_title_side="bottom",
  column_title_gp = gpar(fontsize = 18),
  row_title="Genes",
  row_title_side="left",
  row_title_gp = gpar(fontsize = 18),
  column_dend_height = unit(25, "mm"),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 14),
    title_gp = gpar(fontsize = 14),
    color_bar = "continuous",
    title="Expression\nz-score",
    title_position="topcenter",
    legend_direction="horizontal"
    )
)
dev.off()
