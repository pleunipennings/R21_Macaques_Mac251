source("Rscripts/BaseRscript2.R")
#OK, my next goal is to create a DF for each monkey (but starting with one 4016)
#rows are compartments and columns are nucleotide positions

read.csv("OriginalData/SampleSheetMac251AllSamples.csv", stringsAsFactors = FALSE)->SampleSheet
SeqDataFiles<-list.files("ProcessedData/",pattern="SeqData.csv$", recursive = TRUE)
SampleSheet$Week<-as.numeric(SampleSheet$Week)
SampleSheet<-SampleSheet[order(SampleSheet$Week),]

UniqueMonkeys<-unique(SampleSheet$Monkey)
read.csv("ProcessedData/PositionsSameStock.csv")->PositionsToUse

pdf(paste0("Output/SD_vs_Mean",Sys.Date(),".pdf"))
for (M in UniqueMonkeys){
#  M="4016"
  monkeyrows<-which(SampleSheet$Monkey==M)  
  sampleinfo=c()
  for (i in monkeyrows){sampleinfo<-c(sampleinfo, paste0(SampleSheet$Sample[i],"_w_",SampleSheet$Week[i]))}
  
  #read data for one sample
  mr= monkeyrows[1]
  sdf=SampleSheet$SeqDataFileName[mr]
  print(sdf)
  X<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1)
  DF<-data.frame(matrix(ncol=length(X$pos),nrow=length(monkeyrows), dimnames=list(sampleinfo, X$pos)))

  counter=1
  for (mr in monkeyrows){
    sdf=SampleSheet$SeqDataFileName[mr]
    print(sdf)
    X<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1)
    #print(which(X$TotalReads<10000))
    DF[counter,]<-X$pi
    DF[counter, which(X$TotalReads<10000)]<-NA
    counter = counter +1
  }

#calculate rho for each column and ave pi for each colunm
  
  SD <- apply(DF, MARGIN = 2, sd, na.rm = TRUE)
  Mean <- apply(DF, MARGIN = 2, mean, na.rm = TRUE)
  plot(SD/Mean, Mean, main = M)  
}

dev.off()