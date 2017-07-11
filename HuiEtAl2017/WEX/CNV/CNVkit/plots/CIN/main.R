library("beeswarm")

source("../../lib.R")

chrs <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X")

chrs_len <- list()
chrs_len[["1"]] <- 195471971
chrs_len[["2"]] <- 182113224
chrs_len[["3"]] <- 160039680
chrs_len[["4"]] <- 156508116
chrs_len[["5"]] <- 151834684
chrs_len[["6"]] <- 149736546
chrs_len[["7"]] <- 145441459
chrs_len[["8"]] <- 129401213
chrs_len[["9"]] <- 124595110
chrs_len[["10"]] <- 124595110
chrs_len[["11"]] <- 122082543
chrs_len[["12"]] <- 120129022
chrs_len[["13"]] <- 120421639
chrs_len[["14"]] <- 124902244
chrs_len[["15"]] <- 104043685
chrs_len[["16"]] <- 98207768
chrs_len[["17"]] <- 94987271
chrs_len[["18"]] <- 90702639
chrs_len[["19"]] <- 61431566
chrs_len[["X"]] <- 171031299

data <- load_cnvkit_data_cns()
samples <- samples[(samples$sample_ID %in% names(data)),]

#compare number of breakpoints

num_breakpoints <- list()
num_breakpoints[["BRCA1 DEL"]] <- c()
num_breakpoints[["BRCA1 WT"]] <- c()

sample <- c()

for (i in as.character(samples[samples$brca1=="WT","sample_ID"])) {
  tmp <- data[[i]]
  #tmp <- tmp[tmp$cn!=2,]
  segs = 0
  for (j in chrs) {
    tmp2 <- tmp[tmp[,"chromosome"]==j,]
    if (nrow(tmp2) > 1) {
      segs = segs + (nrow(tmp2)-1)
    }
  }
  num_breakpoints[["BRCA1 WT"]] <- c(num_breakpoints[["BRCA1 WT"]],segs)
  sample <- c(sample,i)
}

for (i in as.character(samples[samples$brca1=="flox/flox","sample_ID"])) {
  tmp <- data[[i]]
  #tmp <- tmp[tmp$cn!=2,]
  segs = 0
  for (j in chrs) {
    tmp2 <- tmp[tmp[,"chromosome"]==j,]
    if (nrow(tmp2) > 1) {
      segs = segs + (nrow(tmp2)-1)
    }
  }
  num_breakpoints[["BRCA1 DEL"]] <- c(num_breakpoints[["BRCA1 DEL"]],segs)
  sample <- c(sample,i)
}

r = wilcox.test(num_breakpoints[["BRCA1 DEL"]],num_breakpoints[["BRCA1 WT"]])
print(r)

pdf("num_breakpoints_BRCA1.pdf")
par(mar=c(3,3,3,3),family="Times",cex=1.5)
beeswarm(num_breakpoints,labels=names(num_breakpoints),cex=0.5,xlab="",ylab="",main="",ylim=c(0,1.1*max(num_breakpoints[["BRCA1 DEL"]])))
boxplot(num_breakpoints, add=TRUE,outline=FALSE,names=rep("",length(num_breakpoints)),pars=list(boxwex=0.5),boxlty=0,whisklty=0,staplelty=0,yaxt="n")
#text(1.5,1.05*max(num_breakpoints[["BRCA1 DEL"]]),paste("p-value = ",round(r$p.value,4),sep=""))
dev.off()

breakpoints <- data.frame(
  samples=sample,
  breakpoints=c(num_breakpoints[["BRCA1 WT"]],num_breakpoints[["BRCA1 DEL"]]),
  brca1=c(
    rep("WT",length(num_breakpoints[["BRCA1 WT"]])),
    rep("flox/flox",length(num_breakpoints[["BRCA1 DEL"]]))
  )
)
write.csv(breakpoints,"breakpoints.csv",quote=F,row.names=F)
