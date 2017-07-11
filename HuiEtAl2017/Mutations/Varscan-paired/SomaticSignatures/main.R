library(lsa)
library(xlsx)
library(gplots)
library("SomaticSignatures")
library(ggplot2)
library(BSgenome.Mmusculus.UCSC.mm10)
library(gplots)

load_coverage <- function() {
  coverage <- read.csv("/Users/charlesmurphy/Desktop/Research/0914_hui/results/Mutations/CountTriNucleotideCoverage/counts.csv",header=T,row.names=1,check.names=F)
  coverage <- coverage[,colnames(coverage)!="HL-test,"]
  samples <- colnames(coverage)
  samples <- gsub("\\[","",samples)
  samples <- gsub("\\]","",samples)
  samples <- gsub(",","",samples)
  colnames(coverage) <- samples

  collapse <- c(
    "ata","atc","atg","att",
    "cta","ctc","ctg","ctt",
    "gta","gtc","gtg","gtt",
    "tta","ttc","ttg","ttt",
    "aca","acc","acg","act",
    "cca","ccc","ccg","cct",
    "gca","gcc","gcg","gct",
    "tca","tcc","tcg","tct"
    )
  collapsed_coverage <- matrix(0,nrow=length(collapse),ncol=ncol(coverage))
  row.names(collapsed_coverage) <- collapse
  colnames(collapsed_coverage) <- colnames(coverage)
  for (i in collapse) {
    complement <- c()
    for (j in strsplit(i,"")[[1]]) {
      if (j=="a") {
        complement <- c(complement,"t")
      } else if (j=="c") {
        complement <- c(complement,"c")
      } else if (j=="g") {
        complement <- c(complement,"g")
      } else {
        complement <- c(complement,"a")
      }
    }
    complement <- paste(complement,collapse="")
    collapsed_coverage[i,] <- as.numeric(coverage[i,] + coverage[complement,])
  }

  return(coverage)
}

compare_to_human <- function(data) {

  sca_vr = VRanges(
      seqnames = data$chr,
      ranges = IRanges(data$pos,data$pos),
      ref = data$ref,
      alt = data$alt,
      sampleNames = data$sample,
      brca1 = data$brca1
  )
  sca_motifs = mutationContext(sca_vr, BSgenome.Mmusculus.UCSC.mm10)
  sca_mm = motifMatrix(sca_motifs,normalize = TRUE)

  # normalize by coverage

  coverage <- load_coverage()

  kmers <- c()
  for (i in row.names(sca_mm)) {
    tmp <- strsplit(i," ")[[1]]
    motif <- tolower(paste(
      subseq(tmp[2],1,1),
      subseq(tmp[1],1,1),
      subseq(tmp[2],3,3),
      sep=""
      ))
    kmers <- c(kmers,motif)
  }

  coverage <- coverage[unique(kmers),]
  coverage <- apply(coverage,2,function(x) return(x/sum(x)))
  coverage <- coverage[kmers,]

  for (sample in colnames(sca_mm)) {
    tmp <- sca_mm[,sample] * coverage[,sample]
    tmp <- tmp / sum(tmp)
    sca_mm[,sample] <- tmp
  }

  # computer normalization factor to compare mouse and human

  mouse_motifs <- read.csv("~/Desktop/Research/0914_hui/data/Mouse_Exome_Design/GRCm38.110624_S0276129.3merCount.csv",header=F,row.names=1)
  human_motifs <- read.csv("~/Desktop/Research/data/refdata/hg19_gatk/hg19.3merCount.csv",header=F,row.names=1)
  motifs <- toupper(row.names(human_motifs))
  mouse_motifs <- as.numeric(mouse_motifs[tolower(motifs),1])
  human_motifs <- as.numeric(human_motifs[tolower(motifs),1])

  mouse_motifs <- mouse_motifs/sum(mouse_motifs)
  human_motifs <- human_motifs/sum(human_motifs)

  names(mouse_motifs) <- motifs
  names(human_motifs) <- motifs

  norms <- human_motifs / mouse_motifs
  sca_mm_norm = normalizeMotifs(sca_mm, norms)

  # average by brca1 genotype

  brca1flfl_samples <- as.character(unique(data[data$brca1=="DEL","sample"]))
  brca1WT_samples <- as.character(unique(data[data$brca1=="WT","sample"]))

  sca_mm_brca1 <- matrix(0,nrow=nrow(sca_mm_norm),ncol=2)
  colnames(sca_mm_brca1) <- c("DEL","WT")
  row.names(sca_mm_brca1) <- row.names(sca_mm_norm)
  sca_mm_brca1[,"DEL"] <- apply(sca_mm_norm[,brca1flfl_samples],1,mean)
  sca_mm_brca1[,"WT"] <- apply(sca_mm_norm[,brca1WT_samples],1,mean)

  # compute distance to human signatures

  humanSignatures <- read.csv("~/Desktop/Research/data/papers/AlexandrovEtAl/signatures.txt",sep="\t",check.names=F)
  row.names(humanSignatures) <- humanSignatures[,"Somatic Mutation Type"]
  humanSignatures <- humanSignatures[,4:ncol(humanSignatures)]

  motifs <- c()
  for (i in row.names(sca_mm_brca1)) {
    tmp <- strsplit(i," ")[[1]]
    mut <- tmp[1]
    mut <- strsplit(mut,"")[[1]]
    context <- tmp[2]
    context <- gsub("\\.",paste("[",mut[1],">",mut[2],"]",sep=""),context)
    motifs <- c(motifs,context)
  }
  row.names(sca_mm_brca1) <- motifs
  sca_mm_brca1 <- sca_mm_brca1[row.names(humanSignatures),]

  humanSignaturesNames <- colnames(humanSignatures)[1:(ncol(humanSignatures)-5)]
  comparison <- cbind(humanSignatures,sca_mm_brca1)

  # cosine similarity

  dd <- cosine(as.matrix(comparison))
  brca1flfl <- dd["DEL",humanSignaturesNames]
  brca1flfl <- brca1flfl[order(brca1flfl,decreasing=T)]

  pdf("./compareHuman/similar_to_brca1flfl_cosine.pdf",width=10,height=5)
  par(mar=c(6,6,3,3))
  barplot(brca1flfl,las=2,ylab="Cosine similarity",ylim=c(0,1))
  dev.off()

  brca1wt <- dd["WT",humanSignaturesNames]
  brca1wt <- brca1wt[order(brca1wt,decreasing=T)]

  pdf("./compareHuman/similar_to_brca1WT_cosine.pdf",width=10,height=5)
  par(mar=c(6,6,3,3))
  barplot(brca1wt,las=2,ylab="Cosine similarity",ylim=c(0,1))
  dev.off()
}

plotSignatureHeatmap <- function(data) {

  #collapse by primary

  sca_vr = VRanges(
      seqnames = data$chr,
      ranges = IRanges(data$pos,data$pos),
      ref = data$ref,
      alt = data$alt,
      sampleNames = data$sample,
      brca1 = data$brca1,
      mouse=as.factor(as.character(data$mouse))
  )
  sca_motifs = mutationContext(sca_vr, BSgenome.Mmusculus.UCSC.mm10)
  sca_mm = motifMatrix(sca_motifs, normalize = F)

  sample_labels <- c()
  samples <- c(
    "P.1221","P.1367","P.1413","P.1415","P.1460","P.1461","P.1512","P.1513","P.1536","P.1601","P.1614","P.1616","P.1660","P.1661","P.1662","P.1795",
    "P.1259","P.1397","P.1531","P.1597","P.1702","P.1720","P.1857","P.1867","P.2005","P.517"
  )
  data_matrix <- matrix(0,nrow=24,ncol=4*length(samples))
  nn <- 1

  for (sample in samples) {
    sample_labels <- c(sample_labels,c(NA,NA,sample,NA))
    sample_tmp <- matrix(0,nrow=24,ncol=4)
    n <- 1
    cc <- c()
    for (mut in c("CA","CG","CT","TA","TC","TG")) {
      sample_tmp[n,] <- sca_mm[paste(mut,c("A.A","C.A","G.A","T.A")),sample]
      n <- n + 1
      sample_tmp[n,] <- sca_mm[paste(mut,c("A.C","C.C","C.C","G.C")),sample]
      n <- n + 1
      sample_tmp[n,] <- sca_mm[paste(mut,c("A.G","C.G","G.G","G.G")),sample]
      n <- n + 1
      sample_tmp[n,] <- sca_mm[paste(mut,c("A.T","C.T","G.T","T.T")),sample]
      n <- n + 1
    }
    #row.names(sample_tmp) <- rep(c("A","C","G","T"),6)
    colnames(sample_tmp) <- c("A","C","G","T")

    data_matrix[,c(nn,nn+1,nn+2,nn+3)] <- sample_tmp
    nn <- nn + 4
  }

  width <- 20
  height <- (width/ncol(data_matrix))*nrow(data_matrix)
  my_palette <- colorRampPalette(c("#ffff99","red"))

  pdf("heatmap_persample.pdf",width=width,height=height*1.45)
  heatmap.2(
    log2(as.matrix(data_matrix)+1),
    trace="none",
    scale="none",
    dendrogram="none",
    col=my_palette,
    Colv=F,
    Rowv=F,
    cexCol=1.5,
    labRow=F,
    lwid=width * c(width*1/10,width*(7/10)),
    lhei=c(0.35*height,0.65*height),
    labCol=sample_labels,
    margins=c(6,8),
    add.expr = c(
      segments(
        4 * 1:(26-1) + 0.5,
        rep(0,26-1),
        4 * 1:(26-1) + 0.5,
        rep(25,26-1)
      ),
      segments(
        rep(0,5),
        4 * 1:5 + 0.5,
        rep(26*4+1,5),
        4 * 1:5 + 0.5,
      ),
      c(
        mtext("C>A",2,at=22,cex=1,las=2,line=1),
        mtext("C>G",2,at=18,cex=1,las=2,line=1),
        mtext("C>T",2,at=14,cex=1,las=2,line=1),
        mtext("T>A",2,at=10,cex=1,las=2,line=1),
        mtext("T>C",2,at=6,cex=1,las=2,line=1),
        mtext("T>G",2,at=2,cex=1,las=2,line=1)
      ),
      lapply(
        seq(1,26*4,4),
          function(x) return(c(
            mtext("A",3,at=x),
            mtext("C",3,at=x+1),
            mtext("G",3,at=x+2),
            mtext("T",3,at=x+3)
          )
        )
      ),
      lapply(
        seq(1,6*4,4),
          function(x) return(c(
            mtext("A",4,at=x+3,las=2,line=0.5),
            mtext("C",4,at=x+2,las=2,line=0.5),
            mtext("G",4,at=x+1,las=2,line=0.5),
            mtext("T",4,at=x,las=2,line=0.5)
          )
        )
      )
    ),
    key.title="log2(# mutations + 1)",
    key.xlab="",
    key.ylab="",
    key.par=list(cex=1)
  )
  text(0.5,0.8,"5' base",cex=2)
  text(0.98,0.5*0.75,"3' base",cex=2,srt=90)
  dev.off()
}

# read the variants file and remove all samples except primary or untreated

data <- read.csv("variant_table.csv",header=T)
sample_data <- read.xlsx("/Users/charlesmurphy/Desktop/Research/0914_hui/data/samples.xlsx","samples",colClasses="character")
data <- data[data$sample %in% as.character(sample_data$sample_ID),]

# compute distance to human signatures

compare_to_human(data)

# collapse by primary tumor

new_data <- data.frame()
primaries <- as.character(unique(data$mouse))

for (pp in primaries) {
  tmp <- data[data$mouse==pp,]
  muts <- paste(tmp$chrom,tmp$pos,sep="-")
  tmp <- tmp[!duplicated(muts),]
  tmp$sample <- rep(pp,nrow(tmp))
  new_data <- rbind(new_data,tmp)
}
row.names(new_data) <- 1:nrow(new_data)
data <- new_data

# plot heatmap of signatures

plotSignatureHeatmap(data)

# plot mutation signature for each BRCA1 genotype

# BRCA1-flox/flox

data_sub <- data[data$brca1=="DEL",]
data_sub$brca1 <- as.character(data_sub$brca1)
sca_vr = VRanges(
    seqnames = data_sub$chr,
    ranges = IRanges(data_sub$pos,data_sub$pos),
    ref = data_sub$ref,
    alt = data_sub$alt,
    sampleNames = data_sub$sample,
    brca1 = data_sub$brca1
)
sca_motifs = mutationContext(sca_vr, BSgenome.Mmusculus.UCSC.mm10)
sca_mm = motifMatrix(sca_motifs, group = "brca1", normalize = TRUE)
pdf("brca1-flfl.pdf",width=10,height=5)
plotMutationSpectrum(sca_motifs,"brca1") +
  ylim(0,0.3) +
  theme_grey(base_size=20) +
  theme(axis.title.y=element_blank(),axis.title.x=element_blank(),axis.text.x=element_text(size=6, angle=90)) +
  theme(legend.position="none") +
  scale_fill_manual(values="red")
dev.off()

# BRCA1-WT

data_sub <- data[data$brca1=="WT",]
data_sub$brca1 <- as.character(data_sub$brca1)
sca_vr = VRanges(
    seqnames = data_sub$chr,
    ranges = IRanges(data_sub$pos,data_sub$pos),
    ref = data_sub$ref,
    alt = data_sub$alt,
    sampleNames = data_sub$sample,
    brca1 = data_sub$brca1
)
sca_motifs = mutationContext(sca_vr, BSgenome.Mmusculus.UCSC.mm10)
sca_mm = motifMatrix(sca_motifs, group = "brca1", normalize = TRUE)
pdf("brca1-WT.pdf",width=10,height=5)
plotMutationSpectrum(sca_motifs,"brca1") +
  ylim(0,0.3) +
  theme_grey(base_size=20) +
  theme(axis.title.y=element_blank(),axis.title.x=element_blank(),axis.text.x=element_text(size=6, angle=90)) +
  theme(legend.position="none") +
  scale_fill_manual(values="red")
dev.off()
