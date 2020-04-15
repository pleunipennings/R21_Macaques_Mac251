#Here I will try to determine what are synonymous and non-synonymous sites. 

source("Rscripts/BaseRscript2.R")

read.csv("OriginalData/SampleSheetMac251AllSamples.csv", stringsAsFactors = FALSE)->SampleSheet
SeqDataFiles<-list.files("ProcessedData/",pattern="SeqData.csv$", recursive = TRUE)
SampleSheet$Week<-as.numeric(SampleSheet$Week)
SampleSheet<-SampleSheet[order(SampleSheet$Week),]
UniqueMonkeys<-unique(SampleSheet$Monkey)


sdf=SampleSheet$SeqDataFileName[1]
X<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1)

#First 90 nucleotides are in frame :-) 