SeqDataFiles<-list.files("ProcessedData/",pattern="SeqData.csv$", recursive = TRUE)

nucleotides<-c("a","c","g","t")

for (sdf in SeqDataFiles){
  #sdf=SeqDataFiles[1]
  X<-read.csv(paste0("ProcessedData/",sdf), row.names = 1)
  X$mutFreq<-0
  for (i in 1:nrow(X)){
    mutnucleotides=nucleotides[-which(nucleotides==X$MajNt[i])]
    mutrows=which(names(X)%in%mutnucleotides)
    X$mutFreq[i]<-sum(X[i,mutrows])/X$TotalReads[i]
  }
  write.csv(x = X, file = paste0("ProcessedData/",sdf))
}


