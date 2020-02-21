#create the bash files to run bbmap and bwa
# read the template command text file

#create the bash files to run bbmap and bwa
# read the template command text file

#choose the fastq files to be prrocessed
read.csv("OriginalData/SampleSheetMac251AllSamples.csv", stringsAsFactors = FALSE)->SampleSheet

stockrows<-which(SampleSheet$Monkey=="stock virus")

#for (i in 1:nrow(SampleSheet)){
for (i in stockrows){
  AnimalNumber= SampleSheet$MiseqSample[i]
  if (AnimalNumber<10) AnimalNumber = paste(c("0",AnimalNumber),collapse="")
  SampleName = paste(c("Run",SampleSheet$SampleSheet[i], AnimalNumber, "Animal", as.character(SampleSheet$Monkey[i])), collapse = "_")
  SampleName<-gsub(pattern=" ", replace="",x=SampleName)
  print(SampleName)
  t1=Sys.time()
  
  FilePath = paste(c(SampleSheet$FilePathPart1[i], SampleSheet$FilePathPart2[i]),collapse="/")
  filelist<-list.files(FilePath, full.names = TRUE)
  FileIn1=filelist[1]
  FileIn2=filelist[2]
  
  #BashscriptFilename = paste0(SampleName,".sh")
  BashLines<-readLines("ProcessedData/BashSampleSheetMac251.sh")
  BashLines = gsub(pattern="FASTQFile1", replace=FileIn1,x=BashLines)
  BashLines = gsub(pattern="FASTQFile2", replace=FileIn2,x=BashLines)
  BashLines = gsub(pattern="SAMPLE", replace=SampleName,x=BashLines)
  bashscriptName<-paste0("ProcessedData/BashScripts/", SampleName,".sh")
  writeLines(BashLines,bashscriptName )
}

bashfiles=list.files("ProcessedData/BashScripts/", pattern = "stockvirus")
#i=1
for (i in 1:length(bashfiles)){
    system(paste("chmod 755 ", "ProcessedData/BashScripts/",bashfiles[i], sep=''))
    system(paste("ProcessedData/BashScripts/",bashfiles[i], sep=''))
}
#To run: ./Data/BashscriptsSampleSheet1/M11.sh