library(xlsx)

sample_data <- read.xlsx("../../../data/samples.xlsx","samples")
row.names(sample_data) <- sample_data$sample_ID


data <- read.csv("genes-symbols-fpkm.csv",row.names=1,check.names=F)
data[data<=0.1]=0.0
data <- data[rowSums(data)>0,]
colnames(data) <- as.character(sample_data[colnames(data),"sample_name"])
write.csv(data,"./rename_filtered/genes-symbols-fpkm.csv",quote=F)


data <- read.csv("genes-humanSymbols-fpkm.csv",row.names=1,check.names=F)
data[data<=0.1]=0.0
data <- data[rowSums(data)>0,]
colnames(data) <- as.character(sample_data[colnames(data),"sample_name"])
write.csv(data,"./rename_filtered/genes-humanSymbols-fpkm.csv",quote=F)

data <- read.csv("genes-humanEntrez-fpkm.csv",row.names=1,check.names=F)
data[data<=0.1]=0.0
data <- data[rowSums(data)>0,]
colnames(data) <- as.character(sample_data[colnames(data),"sample_name"])
write.csv(data,"./rename_filtered/genes-humanEntrez-fpkm.csv",quote=F)


data <- read.csv("genes-fpkm.csv",row.names=1,check.names=F)
data[data<=0.1]=0.0
data <- data[rowSums(data)>0,]
colnames(data) <- as.character(sample_data[colnames(data),"sample_name"])
write.csv(data,"./rename_filtered/genes-fpkm.csv",quote=F)
