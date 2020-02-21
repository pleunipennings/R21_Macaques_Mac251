
SeqDataFiles<-list.files("ProcessedData/",pattern="SeqData.csv$", recursive = TRUE)

for (sdf in SeqDataFiles){
  #sdf=SeqDataFiles[1]
  X<-read.csv(paste0("ProcessedData/",sdf), row.names = 1)
  X$pi<-1-(X$a/X$TotalReads)^2-(X$c/X$TotalReads)^2-(X$g/X$TotalReads)^2-(X$t/X$TotalReads)^2
  write.csv(x = X, file = paste0("ProcessedData/",sdf))
}



