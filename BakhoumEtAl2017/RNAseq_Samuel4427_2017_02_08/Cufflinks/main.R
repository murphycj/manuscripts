library(xlsx)

data <- read.csv("genes-fpkm.csv",row.names=1,check.names=F)
data[data<=0.1]=0.0
data <- data[rowSums(data)>0,]

pdf("boxplot-fpkm.pdf",width=15,height=5)
boxplot(as.matrix(log2(data+0.1)),las=2,ylab="log2(FPKM+0.1)")
dev.off()
