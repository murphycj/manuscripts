library(gplots)
library(xlsx)

samples <- read.xlsx("/Users/charlesmurphy/Desktop/Research/0914_hui/data/samples.xlsx","samples",row.names=1)
samples <- samples[samples$sequencing_type=="RNAseq",]
samples <- samples[samples$tissue %in% c("primary","implant"),]
samples <- samples[row.names(samples)!="HL23",]
samples <- samples[!is.na(samples$brca1),]
brca1 <- as.character(samples$brca1)

#load mouse and human exrpession data

genes <- read.csv("genes.csv")[,1]

fpkm <- read.csv(
  "/Users/charlesmurphy/Desktop/Research/0914_hui/results/RNAseq/Cufflinks/rename_filtered/genes-humanEntrez-fpkm.csv",
  row.names=1,
  header=T,
  check.names=F
)
common <- intersect(colnames(fpkm),as.character(samples$sample_name))
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

write.table(common,"mouse-human-common-genes.csv",quote=F,row.names=F,colnames=F)

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

colors <- colnames(rpkm.entrez)
colors[colors %in% subtypes_class[["TNBC"]]]<-"red"
colors[colors %in% subtypes_class[["Basal-like"]]]<-"black"
colors[colors %in% subtypes_class[["Luminal A"]]]<-"purple"
colors[colors %in% subtypes_class[["Luminal B"]]]<-"orange"
colors[colors %in% subtypes_class[["HER2-enriched"]]]<-"grey"
colors[colors %in% subtypes_class[["Normal-like"]]]<-"blue"
colors <- c(
  brca1,
  colors
)
colors[colors=="WT"] <- "#00e600"
colors[colors=="flox/flox"] <- "#006600"

data <- as.matrix(cbind(fpkm,rpkm.entrez))

png("TNBC.png",width=1000,height=1000,family="Times")
heatmap.2(
  data,
  trace="none",
  keysize=1.5,
  scale="row",
  labCol=NA,
  labRow=NA,
  margins=c(2,2),
  col= redgreen(75),
  key.title="Gene expression",
  ColSideColors=colors,
  main="",
  lmat=rbind(c(5,0),c(0,4),c(0,1),c(3,2)),
  lhei=c(1,0.75,0.25,4),
  lwid=c(1,4),
  key.par=list(mar=c(1,1,3,1),cex.main=2.,cex.axis=2,cex.lab=2),
  distfun=function(x) as.dist((1-cor(t(x),method="spearman"))),
  hclustfun = function(x) hclust(x,method = 'average')
)
legend(
  0.3,1.,
  c("TNBC","Basal-like", "Luminal A","Luminal B","HER2-enriched","Normal-like","mouse (BRCA1 WT)","mouse (BRCA1 DEL)"),
  pch=16,
  col=c("red","black","purple","orange","grey","blue","#00e600","#006600"),
  title="Subtype",
  cex=1.5,ncol=2
)
dev.off()
