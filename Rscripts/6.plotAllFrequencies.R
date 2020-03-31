source("Rscripts/BaseRscript2.R")
#OK, my next goal is to plot mutation frequencies for each sample. 

read.csv("OriginalData/SampleSheetMac251AllSamples.csv", stringsAsFactors = FALSE)->SampleSheet
SeqDataFiles<-list.files("ProcessedData/",pattern="SeqData.csv$", recursive = TRUE)

SampleSheet$Week<-as.numeric(SampleSheet$Week)
SampleSheet<-SampleSheet[order(SampleSheet$Week),]

UniqueMonkeys<-unique(SampleSheet$Monkey)
read.csv("ProcessedData/PositionsSameStock.csv")->PositionsToUse

pdf(paste0("Output/MutFreqInMonkeys",Sys.Date(),".pdf"))
for (M in UniqueMonkeys){
  #M=UniqueMonkeys[1]
  monkeyrows<-which(SampleSheet$Monkey==M)  
  sampleinfo=c()
  for (i in monkeyrows){sampleinfo<-c(sampleinfo, paste0(SampleSheet$Sample[i],"_w_",SampleSheet$Week[i]))}
  par(mar=c(4, 4, 4, 4))
  
  plot(x=PositionsToUse$x, y=rep(0,nrow(PositionsToUse)), ylim=c(0.00001,0.5),type="n", xlab="Env Position",ylab="Mutation Frequency",  
       main = paste0("Monkey ",M), log="y")
  #axis(1, at=1:length(monkeyrows),labels=sampleinfo, col.axis="red", las=2)
  
  counter=1
  for (mr in monkeyrows){
    sdf=SampleSheet$SeqDataFileName[mr]
    print(sdf)
    X<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1)
    X<-X[which(X$TotalReads>=10000),]
    points(x=X$pos,y=X$a/X$TotalReads,pch=16, col=alpha(cols[counter],0.5))
    points(x=X$pos,y=X$c/X$TotalReads,pch=16, col=alpha(cols[counter],0.5))
    points(x=X$pos,y=X$g/X$TotalReads,pch=16, col=alpha(cols[counter],0.5))
    points(x=X$pos,y=X$t/X$TotalReads,pch=16, col=alpha(cols[counter],0.5))
    
    if (M!="stockvirus")text(x = 450,y = 10^(-4+counter*0.2), labels = SampleSheet$Sample[mr], col=cols[counter])
    if (M=="stockvirus")text(x = 450,y = 10^(-4+counter*0.2), labels = paste0(SampleSheet$Sample[mr], " ", SampleSheet$Miseq.Date[mr]), col=cols[counter])
    counter = counter +1
  }
}
dev.off()

