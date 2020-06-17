#Add WTAA MUTAA and TypeOfSite and makesCPG to all files 

for (i in 1:length(SampleSheet$SeqDataFileName)){

  sdf=SampleSheet$SeqDataFileName[i]
  print(sdf)
  X<-read.csv(paste0("ProcessedData/DataCpG/",sdf), row.names = 1)
  X<-getWTAA(X)
  X<-getMUTAA(X)
  X<-synFunction(X)
  X<-CPG_site(X)
  write.csv(x = X, file = paste0("ProcessedData/DataCpG/",sdf))

}