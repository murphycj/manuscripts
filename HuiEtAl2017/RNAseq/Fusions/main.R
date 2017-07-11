library("beeswarm")
library(gplots)
library(xlsx)
options( stringsAsFactors = FALSE)


write_fusion_table <- function(fusions,samples) {
  #
  # if there are fusions with multiple breakpoints, just choose one
  # randomly so the plot looks cleaner
  # then create a table of the fusion events
  #


  fusion.ids <- unique(paste(
    fusions$Gene_1_symbol.5end_fusion_partner.,
    fusions$Gene_2_symbol.3end_fusion_partner.,
    sep="-")
  )
  s <- as.character(unique(fusions[,"sample"]))

  fusions.global <- matrix(0,nrow=length(fusion.ids),ncol=length(s))
  row.names(fusions.global) <- fusion.ids
  colnames(fusions.global) <- s

  for (i in s) {
    temp = fusions[fusions$sample==i,]
    temp.ids <- unique(paste(
      temp$Gene_1_symbol.5end_fusion_partner.,
      temp$Gene_2_symbol.3end_fusion_partner.,
      sep="-")
    )
  	fusions.global[temp.ids,i]<-1
  }

  new_names <- c("Gene-5prime","Gene-3prime","Count",colnames(fusions.global))

  fusions.global <- cbind(rowSums(fusions.global),fusions.global)

  fusions.global <- cbind(
    unlist(lapply(strsplit(fusion.ids,"-"),function(x) return(x[1]))),
    unlist(lapply(strsplit(fusion.ids,"-"),function(x) return(x[2]))),
    fusions.global
  )
  colnames(fusions.global) <- new_names

  return(fusions.global)
}


count_fusions <- function(fusions,samples) {

  fusion.ids <- unique(paste(
    fusions$Gene_1_symbol.5end_fusion_partner.,
    fusions$Gene_2_symbol.3end_fusion_partner.,
    fusions$Fusion_point_for_gene_1.5end_fusion_partner.,
    fusions$Fusion_point_for_gene_2.3end_fusion_partner.,
    sep="-")
  )

  fusions.global <- matrix(nrow=length(fusion.ids),ncol=length(samples$sample_name))
  samples <- as.character(samples$sample_name)
  row.names(fusions.global)<-fusion.ids
  colnames(fusions.global) <- samples

  for (i in samples) {
  	temp = fusions[fusions$sample==i,]
    temp.ids <- unique(paste(
      temp$Gene_1_symbol.5end_fusion_partner.,
      temp$Gene_2_symbol.3end_fusion_partner.,
      temp$Fusion_point_for_gene_1.5end_fusion_partner.,
      temp$Fusion_point_for_gene_2.3end_fusion_partner.,
      sep="-")
    )
  	fusions.global[temp.ids,i]<-1
  	fusions.notsample<-setdiff(fusion.ids,temp.ids)
  	fusions.global[fusions.notsample,i]<-0
  }
  return(fusions.global)
  #fusions.global<-fusions.global[order(rowSums(fusions.global),decreasing=T),]
  #fusions.global <- filter_by_occurance(fusions.global,2,sample_data)
}

collapse_by_primary <- function(fusions,fusion.counts,samples) {

  primaries <- unique(samples$mouse)
  fusions.collapsed <- matrix(nrow=nrow(fusion.counts),ncol=length(primaries))
  row.names(fusions.collapsed) <- row.names(fusion.counts)
  colnames(fusions.collapsed) <- primaries

  for (i in 1:length(primaries)) {
    samples.temp <- as.character(samples[samples$mouse==primaries[i],]$sample_name)
    if (length(samples.temp)>1) {
      fusions.collapsed[,primaries[i]] <- rowSums(fusion.counts[,samples.temp])
    } else {
      fusions.collapsed[,primaries[i]] <- fusion.counts[,samples.temp]
    }
  }

	return(fusions.collapsed)
}

filter_by_occurance_in_primary <- function(fusions,fusions.primary,n) {
  fusions.primary <- fusions.primary > 0
  fusions.primary <- row.names(fusions.primary[rowSums(fusions.primary)<n,])

  fusions.to.keep <- c()

  for (i in 1:nrow(fusions)) {
    temp <- paste(
      fusions[i,]$Gene_1_symbol.5end_fusion_partner.,
      fusions[i,]$Gene_2_symbol.3end_fusion_partner.,
      fusions[i,]$Fusion_point_for_gene_1.5end_fusion_partner.,
      fusions[i,]$Fusion_point_for_gene_2.3end_fusion_partner.,
      sep="-"
    )

    if (temp %in% fusions.primary) {
      fusions.to.keep <- c(fusions.to.keep,i)
    }
  }
  return(fusions[fusions.to.keep,])
}


get_unique_treatments <- function(fusions,fusions.table,sample) {

  fusions.table <- fusions.table[,-c(1,2,3)]
  class(fusions.table) <- "numeric"
  fusions.table.unique <- fusions.table[apply(fusions.table,1,sum)==1,]

  samples.implants <- samples[samples["tissue"]=="implant",]
  samples.implants <- samples.implants[samples.implants[,"sample_name"] %in% colnames(fusions.table.unique),]


}

samples <- read.xlsx("~/Desktop/Research/0914_hui/data/samples.xlsx","samples")
row.names(samples) <- as.character(samples$sample_ID)

fusions <- read.csv("fusions.csv",header=T)
fusions$sample <- samples[as.character(fusions$sample),"sample_name"]
write.csv(fusions,"fusions.rename.csv",row.names=F)

#remove those in normal samples

fusions.normal <- fusions[fusions$sample %in% c("Normal breast 1","Normal breast 2","Normal breast 3"),]
fusions.normal <- unique(paste(fusions.normal$Gene_1_symbol.5end_fusion_partner.,fusions.normal$Gene_2_symbol.3end_fusion_partner.,sep="-"))
fusions.all <- paste(fusions$Gene_1_symbol.5end_fusion_partner.,fusions$Gene_2_symbol.3end_fusion_partner.,sep="-")
fusions <- fusions[!(fusions.all %in% fusions.normal),]

#remove readthoughs

fusions <- fusions[!grepl("readthrough",as.character(fusions[,"Fusion_description"])),]

#remove fusions from repeat samples

repeat_samples <- samples[samples$repeat.==1,]$sample_name
fusions <- fusions[!(fusions$sample %in% repeat_samples),]

write.csv(fusions,"fusions.noNormal.noReadthrough.csv",row.names=F)

fusion.counts <- count_fusions(fusions=fusions,samples=samples)

#remove fusions in more than 5 or 3 primary tumors

fusions.primary <- collapse_by_primary(fusions,fusion.counts,samples)

fusions.filtered <- filter_by_occurance_in_primary(fusions=fusions,fusions.primary=fusions.primary,n=3)
write.csv(fusions.filtered,"fusions-noNormal.noReadthrough.n2.csv",row.names=F)

fusions.table <- write_fusion_table(fusions=fusions.filtered,samples=samples)
write.csv(fusions.table,"fusions-table-noNormal.noReadthrough.n2.csv",row.names=F)

#fusions.filtered <- filter_by_occurance_in_primary(fusions=fusions,fusions.primary=fusions.primary,n=5)
#write.csv(fusions.filtered,"fusions-noNormal.noReadthrough.n5.csv",row.names=F)

#fusions.table <- write_fusion_table(fusions=fusions.filtered,samples=samples)
#write.csv(fusions.table,"fusions-table-noNormal.noReadthrough.n5.csv",row.names=F)
