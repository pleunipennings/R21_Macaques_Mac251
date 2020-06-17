SeqDataFiles<-list.files("ProcessedData/",pattern="SeqData.csv$", recursive = TRUE)

nucleotides<-c("a","c","g","t")

for (sdf in SeqDataFiles){
  #sdf=SeqDataFiles[1]
  X<-read.csv(paste0("ProcessedData/",sdf), row.names = 1)
  X$mutFreq<-0
  X$trasitionmutFreq<-0
  for (i in 1:nrow(X)){
    #print(i)
    mutnucleotides=nucleotides[-which(nucleotides==X$MajNt[i])]
    mutcols=which(names(X)%in%mutnucleotides)
    if (!is.na(X$TotalReads[i])){X$mutFreq[i]<-sum(X[i,mutcols])/X$TotalReads[i]}
    transcol=which(names(X)%in%X$transition.maj[i])
    if (!is.na(X$TotalReads[i])){X$trasitionmutFreq[i]<-X[i,transcol]/X$TotalReads[i]}
    if (is.na(X$TotalReads[i])){X$mutFreq[i]<-NA; X$trasitionmutFreq[i]<-NA}
    }
  write.csv(x = X, file = paste0("ProcessedData/",sdf))
}


