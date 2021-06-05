data <- read.csv("genesSymbols-fpkm.csv",row.names=1,header=T,check.names=F)
data[data<=0.1]<-0
data <- data[apply(data,1,sum)>0,]
write.csv(data,"./cleanUp/genesSymbols-fpkm.csv")
