#OK, my next goal is to plot average pi for each sample. 

read.csv("OriginalData/SampleSheetMac251AllSamples.csv", stringsAsFactors = FALSE)->SampleSheet
SeqDataFiles<-list.files("ProcessedData/",pattern="SeqData.csv$", recursive = TRUE)

SampleSheet$Week<-as.numeric(SampleSheet$Week)
SampleSheet<-SampleSheet[order(SampleSheet$Week),]

UniqueMonkeys<-unique(SampleSheet$Monkey)

pdf(paste0("Output/DiversityPiInMonkeys",Sys.Date(),".pdf"))
for (M in UniqueMonkeys){
  #M=UniqueMonkeys[1]
  monkeyrows<-which(SampleSheet$Monkey==M)  
  sampleinfo=c()
  for (i in monkeyrows){sampleinfo<-c(sampleinfo, paste0(SampleSheet$Sample[i],"_w_",SampleSheet$Week[i]))}
  
  par(mar=c(10, 4, 4, 8) + 0.1)
  plot(x=1:length(monkeyrows), y=1:length(monkeyrows), ylim=c(0,0.01),type="n", xaxt="n",ylab="diversity(pi)", xlab = "", 
       main = paste0("Monkey ",M))
       axis(1, at=1:length(monkeyrows),labels=sampleinfo, col.axis="red", las=2)
  
  counter=1
  for (mr in monkeyrows){
    sdf=SampleSheet$SeqDataFileName[mr]
    print(sdf)
    X<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1)
    avePi=mean(X$pi,na.rm=TRUE)
    points(counter,avePi,pch=16)
    counter = counter +1
  }
}
dev.off()

