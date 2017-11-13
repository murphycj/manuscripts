library("IRanges")

for (sd in c("1","2","3","4")) {
  for (sample in list.files(paste("./SD",sd,sep=""),"HL*")) {
    sample_name <- paste(sample,"_40000_",sd,sep="")

    data_file <- paste("./SD",sd,"/",sample,"/",sample,"_result_40000_",sd,"/CNAprofiles/segment.Rdata",sep="")

    load(data_file)

    seg_data <- segment.CNA.object$output[segment.CNA.object$output$ID!="log2.q30diff.bam.vs.log2.q30diff.bam",]
    log2_data <- segment.CNA.object$data

    write.csv(seg_data,paste("./SD",sd,"/copywriter_results/",sample,".seg.csv",sep=""),row.names=F)
    write.csv(log2_data,paste("./SD",sd,"/copywriter_results/",sample,".log2.csv",sep=""),row.names=F)
  }
}
