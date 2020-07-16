source("Rscripts/BaseRscript2.R")

read.csv("OriginalData/SampleSheetMac251AllSamples.csv", stringsAsFactors = FALSE)->SampleSheet
SampleSheet$Week<-as.numeric(SampleSheet$Week)
SampleSheet<-SampleSheet[order(SampleSheet$Week),]
UniqueMonkeys<-unique(SampleSheet$Monkey)
UniqueMonkeys<-UniqueMonkeys[UniqueMonkeys!="stockvirus"&UniqueMonkeys!="control"]

#pdf(paste0("Output/TRY",Sys.Date(),".pdf"), height=12,width=8)
for (M in UniqueMonkeys){
  #par(mfrow=c(2,2))
  M=UniqueMonkeys[1]
  print(M)
  monkeyrows<-which(SampleSheet$Monkey==M)
  for (mr in monkeyrows){
    #par(mfrow=c(2,2))
    print("")
    sdf=SampleSheet$SeqDataFileName[mr]
    MV<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1)
    #p=which(MV$StockMajNt=="a")
    #plot(MV$pos[p],MV$freqg[p], col=MV$Typeg[p]) 
    print(round(wilcox.test(MV$freqa[MV$Typea=="syn"&MV$StockMajNt=="g"],MV$freqa[MV$Typea=="nonsyn"&MV$StockMajNt=="g"], alternative="gr")$p.value,3))
    print(round(wilcox.test(MV$freqc[MV$Typec=="syn"&MV$StockMajNt=="t"],MV$freqc[MV$Typec=="nonsyn"&MV$StockMajNt=="t"], alternative="gr")$p.value,3))
    print(round(wilcox.test(MV$freqg[MV$Typeg=="syn"&MV$StockMajNt=="a"],MV$freqg[MV$Typeg=="nonsyn"&MV$StockMajNt=="a"], alternative="gr")$p.value,3))
    print(round(wilcox.test(MV$freqt[MV$Typet=="syn"&MV$StockMajNt=="c"],MV$freqt[MV$Typet=="nonsyn"&MV$StockMajNt=="c"], alternative="gr")$p.value,3))
  }
  