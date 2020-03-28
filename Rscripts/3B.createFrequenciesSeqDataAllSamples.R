###########################
#
SIVFiles<-list.files("ProcessedData/CSV/",pattern="csv")

read.csv("ProcessedData/PositionsSameStock.csv")->PositionsToUse
no_df<-data.frame("pos"=PositionsToUse$x)
read.csv("ProcessedData/RefPositionsSameStock.csv")->RefPositionsToUse

print(nrow(RefPositionsToUse))

for (i in 1:length(SIVFiles)){
  print(i)
  SeqData<-read.csv(paste("ProcessedData/CSV/",SIVFiles[i],sep=""))
  SeqData<-SeqData[,-c(1,2,9)]
  colnames(SeqData)[1]<-"pos"
  colnames(SeqData)[7:8]<-c("deletion","insertion")
  #remove the reads<100 #PSP 2020 not doin this, I want the same sites in every file
  #determine the majority nucleotide base at each site
  colnames(SeqData)[2:5]<-c("a","c","g","t")
  SeqData$MajNt<-apply(SeqData[,2:5],1,function(x) c("a","c","g","t")[which.max(x)])
  #read the refrence sequence:
  SeqData<-merge(no_df,SeqData,by="pos",all.x=T)
  SeqData$ref<-RefPositionsToUse$x
  head(SeqData,10)
  SeqData$transition.maj<-NA
  SeqData$transition.ref<-NA
  for (j in 1:nrow(SeqData)) SeqData$transition.maj[j]<-transition(SeqData$MajNt[j])
  for (j in 1:nrow(SeqData)) SeqData$transition.ref[j]<-transition(SeqData$ref[j])
  #rearrange the columns
  SeqData<-SeqData[,c("a","c","g","t","deletion","insertion","N","pos","TotalReads","MajNt","ref","transition.maj","transition.ref")]
  SeqDataFilename<-gsub(pattern=".csv", replace="SeqData.csv",x=SIVFiles[i])
  write.csv(SeqData, paste0("ProcessedData/SeqData/",SeqDataFilename))
}

