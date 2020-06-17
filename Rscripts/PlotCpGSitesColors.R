#I am going to plot frequencies along the gene with 6 colors: SYN / NONSYN / CpG / Non CpG / Nonsense

read.csv("OriginalData/SampleSheetMac251AllSamples.csv", stringsAsFactors = FALSE)->SampleSheet
#SeqDataFiles<-list.files("ProcessedData/",pattern="SeqData.csv$", recursive = TRUE)
SampleSheet$Week<-as.numeric(SampleSheet$Week)
SampleSheet<-SampleSheet[order(SampleSheet$Week),]
UniqueMonkeys<-unique(SampleSheet$Monkey)

#pdf(paste0("Output/CpG_6000totalreads",Sys.Date(),".pdf"))
readscutoff=6000
#for (M in UniqueMonkeys){
  M=UniqueMonkeys[1]
  monkeyrows<-which(SampleSheet$Monkey==M)
  #for (mr in monkeyrows){
    mr = monkeyrows[1]
    sdf=SampleSheet$SeqDataFileName[mr]
    print(sdf)
    X<-read.csv(paste0("ProcessedData/DataCpG/",sdf), row.names = 1)
    X<-X[X$TotalReads>readscutoff,]
    X$makesCpG<-as.factor(X$makesCpG)
    ggplot(X[X$TypeOfSite=="syn",], aes(x=pos, y=trasitionmutFreq, color=makesCpG, shape = MajNt)) +
    geom_point() 
    