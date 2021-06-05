library(gplots)
source("/Users/charlesmurphy/Desktop/Research/lib/seqpy/lib/plot_pathways.R")

plot_pathway(
	data="/Users/charlesmurphy/Desktop/Research/0816_sam/results/RNAseq/Cufflinks/genesSymbols-fpkm.csv",
	gmt="/Users/charlesmurphy/Desktop/Research/data/msigdb/c2.all.v5.1.symbols.gmt",
	pathway="TAVAZOIE_METASTASIS",
	samples=c(),
	prefix="TAVAZOIE_METASTASIS",
	main="TAVAZOIE_METASTASIS activity in all samples",
	cexRow=1
)
