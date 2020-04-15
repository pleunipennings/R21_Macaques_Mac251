source("Rscripts/BaseRscript2.R")
#Let's look at total reads

read.csv("OriginalData/SampleSheetMac251AllSamples.csv", stringsAsFactors = FALSE)->SampleSheet
SeqDataFiles<-list.files("ProcessedData/",pattern="SeqData.csv$", recursive = TRUE)
SampleSheet$Week<-as.numeric(SampleSheet$Week)
SampleSheet<-SampleSheet[order(SampleSheet$Week),]

UniqueMonkeys<-unique(SampleSheet$Monkey)

#pdf(paste0("Output/Totalreads",Sys.Date(),".pdf"))
for (M in UniqueMonkeys){
  #M="4016"
  print(M)
  monkeyrows<-which(SampleSheet$Monkey==M)  
  sampleinfo=c()
  for (i in monkeyrows){sampleinfo<-c(sampleinfo, paste0(SampleSheet$Sample[i],"_w_",SampleSheet$Week[i]))}
  counter=1
  for (mr in monkeyrows){
    sdf=SampleSheet$SeqDataFileName[mr]
    print(paste0("    ", SampleSheet$Sample[mr]))
    X<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1)
    counter=counter+1
    print(paste0("       ",length(which(X$TotalReads>10000))))
    print(paste0("       ",mean(X$TotalReads)))
    print("       ")
  }
  }
