#Add WTAA MUTAA and TypeOfSite and makesCPG to all files 

source("Rscripts/BaseRscript2.R")

r=which(SampleSheet$Monkey=="stockvirus")
sdf=SampleSheet$SeqDataFileName[r]
print(sdf)
SV<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1) #SV = stock virus
dim(SV)

for (sdf in SampleSheet$SeqDataFileName){
  print(sdf)
  X<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1)
  X<-getWTAA(X)
  X<-getMUTAA(X)
  X<-synFunction(X)
  X<-CPG_site(X)
  X$StockMajNt<-SV$MajNt
  X$freqa = X$a/X$TotalReads
  X$freqc = X$c/X$TotalReads
  X$freqg = X$g/X$TotalReads
  X$freqt = X$t/X$TotalReads
  X=getWTAA_Stock(X)
  X=getMUTAA_onemut(X)
  X<-synFunction_stock(X)
  for (i in 1:nrow(X)){
    if (X$StockMajNt[i]=="a") X$Typea[i]="self"
    if (X$StockMajNt[i]=="c") X$Typec[i]="self"
    if (X$StockMajNt[i]=="g") X$Typeg[i]="self"
    if (X$StockMajNt[i]=="t") X$Typet[i]="self"
  }
  write.csv(x = X, file = paste0("ProcessedData/SeqData/",sdf))
}
