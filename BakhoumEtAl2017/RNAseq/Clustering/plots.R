library(calibrate)
library(gplots)
library(genefilter)
library(dendextend)
library("beeswarm")

get_group <- function(samples, sample_data) {

	group <- list()
	results <- c()
	n = 1
	for (i in 1:length(samples)) {
		temp <- sample_data[sample_data$sample==samples[i],]

		temp <- as.character(temp$group)

		if (temp %in% names(group)) {
			results <- c(results,group[[temp]])
		} else {
			group[[temp]] = n
			n = n + 1
			results <- c(results,group[[temp]])
		}
	}

	return(list(results,group))
}

perform_pcas <- function(data,sample_data) {

	re <- prcomp(t(data))

	pcs=re$sdev**2
	pcs = pcs/sum(pcs) * 100

	pdf("PC-percent-variance.pdf")
	par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE,cex=1.5)
	barplot(
		pcs,
		names.arg=unlist(lapply(1:length(pcs),function(x) return(paste("PC",x,sep="")))),
		las=2,
		ylab="Percentage of variance (%)",
		ylim=c(0,100)
	)
	dev.off()

	#save loadings
	write.csv(re$rotation,"PCA-loadings.csv")

	write.csv(re$x,"PCA-PCs.csv")

	cohort <- get_group(colnames(data),sample_data)

	colors <- rainbow(max(cohort[[1]]),start=0.1,end=0.9)[cohort[[1]]]

	legend_colors = rainbow(max(cohort[[1]]),start=0.1,end=0.9)


	#pca, color samples by cohort
	for (i in 1:NUMBER_PCS) {
		for (j in 1: NUMBER_PCS) {
			if (i == j) {
				next
			}
			pdf(paste("./PCA/PC",i,"PC",j,".pdf",sep=""))
			par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE,cex=1.5)
			plot(
				re$x[,i],
				re$x[,j],
				col=colors,
				pch=16,
				ylab=paste("PC",j,sep=""),
				xlab=paste("PC",i,sep=""),
				main="PCA on all samples",
				ylim=c(1.2*min(re$x[,j]),1.2*max(re$x[,j])),
				xlim=c(1.2*min(re$x[,i]),1.2*max(re$x[,i]))
			)
			textxy(re$x[,i],re$x[,j],names(re$x[,i]),cex=0.4)
			legend(x="topright",names(cohort[[2]]),cex=0.85,inset=c(-0.55,0),pch=16,col=legend_colors,title="Group")
			dev.off()
		}
	}
}

clustering <- function(data,sample_data) {

	#cluster on correlation distance

	diss <- 1 - cor(data,use="complete.obs",method="spearman")
	diss.2 <- as.dist(diss)
	dend <- as.dendrogram(hclust(diss.2,method="average"))

	#original labels
	png("./Cluster/cluster-all-samples-correlation-distance.png",width=1200,height=600)
	par(mar=c(11,7,7,7),cex=1.5)
	plot(dend,main="Hierarchical clustering on all genes",ylab="Correlation distance")
	dev.off()

	#clsuter on euclidean distance

	diss <- dist(t(data))
	dend <- as.dendrogram(hclust(diss,method="average"))

	png("./Cluster/cluster-all-samples-euclidean.png",width=1200,height=600)
	par(mar=c(11,7,7,7),cex=1.5)
	plot(dend,main="Hierarchical clustering on all genes",ylab="Euclidean distance")
	dev.off()

	#cluster on top 5000 more variable genes

	vars <- rowVars(data)
	t1 <- sort(vars,index.return=T,decreasing=T)
	data2 <- data[t1$ix,]
	data2 <- data2[1:5000,]

	diss <- dist(t(data2))
	dend <- as.dendrogram(hclust(diss,method="average"))

	png("./Cluster/cluster-all-samples-euclidean-top5000.png",width=1200,height=600)
	par(mar=c(11,7,7,7),cex=1.5)
	plot(dend,main="Hierarchical clustering on all genes (top 5000 genes by variance)",xlab="")
	dev.off()

	#cluster on top 1000 more variable gene

	vars <- rowVars(data)
	t1 <- sort(vars,index.return=T,decreasing=T)
	data2 <- data[t1$ix,]
	data2 <- data2[1:1000,]

	diss <- dist(t(data2))
	dend <- as.dendrogram(hclust(diss,method="average"))

	png("./Cluster/cluster-all-samples-euclidean-top1000.png",width=1200,height=600)
	par(mar=c(11,7,7,7),cex=1.5)
	plot(dend,main="Hierarchical clustering on all genes (top 1000 genes by variance)",xlab="")
	dev.off()

}

NUMBER_PCS = 4
data <- read.csv("../Cufflinks/genes-fpkm.csv",row.names=1,check.names=F)
data[data<=0.1]=0.0
data <- data[rowSums(data)>0,]
data <- log2(data + 0.1)

pdf("skedasticity.pdf")
par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE,cex=1.5)
plot(
	apply(data,1,mean),
	apply(data,1,sd),
	ylab="SD",
	xlab="mean"
)
dev.off()

sample_data <- read.csv("../../../data/samples.csv")

perform_pcas(data=data,sample_data)
clustering(data=data,sample_data)
