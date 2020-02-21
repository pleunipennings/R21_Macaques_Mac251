library(Rsamtools)
library(stringr)

source("Rscripts/pileupFreq.R")

#number of sampels to process
#runnames<-list.files("ProcessedData/", pattern= "Run")
bamfiles<-list.files("ProcessedData/",pattern="sort.bam$", recursive = TRUE)
baifiles<-list.files("ProcessedData/",pattern="bai$", recursive = TRUE)

for (i in 1:length(bamfiles)){
#for (i in stockrows){
  t1=Sys.time()
  bam<-paste0("ProcessedData/",bamfiles[i])
  index<-paste0(bam,'.bai')
  bf<-BamFile(bam, index=index)
  p_param <- PileupParam(max_depth=70000)
  result<-pileup(bf, pileupParam = p_param, distinguish_strands=FALSE,ignore_query_Ns=FALSE)
  summary<-pileupFreq(result)
  
  summary$TotalReads<-rowSums(summary[3:6])
  maxr<-max(summary$TotalReads)
  cat("The maximum number of read depth is ", maxr)
  cat("\n")
  write.csv(summary, file=paste0("ProcessedData/CSV/",strsplit(bam,"/")[[1]][2],".csv",collapse=""))
  
  t2=Sys.time()
  print(t2-t1)
}


