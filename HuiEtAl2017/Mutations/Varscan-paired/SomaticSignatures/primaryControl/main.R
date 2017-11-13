library(lsa)
library(xlsx)
library(gplots)
library("SomaticSignatures")
library(ggplot2)
library(BSgenome.Mmusculus.UCSC.mm10)
library(gplots)

load_coverage <- function(sample_data) {
  coverage <- read.csv("/Users/charlesmurphy/Desktop/Research/0914_hui/results/Mutations/CountTriNucleotideCoverage/counts.csv",header=T,row.names=1,check.names=F)

  samples <- colnames(coverage)
  samples <- gsub("\\[","",samples)
  samples <- gsub("\\]","",samples)
  samples <- gsub(",","",samples)
  colnames(coverage) <- samples
  common <- intersect(as.character(sample_data$sample_ID),colnames(coverage))
  coverage <- coverage[,common]

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

compare_to_human <- function(data, sample_data) {

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

  coverage <- load_coverage(sample_data)

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

# read the variants file and remove all samples except primary or untreated

data <- read.csv("variant_table.csv",header=T)
sample_data <- read.xlsx("/Users/charlesmurphy/Desktop/Research/0914_hui/data/samples.xlsx","samples",colClasses="character")
sample_data <- sample_data[(sample_data['byl']==0) & (sample_data['bkm']==0) & (sample_data['bmn']==0),]
data <- data[data$sample %in% as.character(sample_data$sample_ID),]

# compute distance to human signatures

compare_to_human(data,sample_data)

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
