
read.csv("OriginalData/SampleSheetMac251AllSamples.csv", stringsAsFactors = FALSE)->SampleSheet
SeqDataFiles<-list.files("ProcessedData/",pattern="SeqData.csv$", recursive = TRUE)

#Add SeqDataFile names to the sample sheet
SampleSheet$SeqDataFileName<-list.files("ProcessedData/SeqData/",pattern="SeqData.csv$", recursive = TRUE)
write.csv(x=SampleSheet,file = "OriginalData/SampleSheetMac251AllSamples.csv")
