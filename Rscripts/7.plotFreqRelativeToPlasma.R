source("Rscripts/BaseRscript2.R")
#OK, my next goal is to plot mutation frequencies for each sample. 

read.csv("OriginalData/SampleSheetMac251AllSamples.csv", stringsAsFactors = FALSE)->SampleSheet
SeqDataFiles<-list.files("ProcessedData/",pattern="SeqData.csv$", recursive = TRUE)

SampleSheet$Week<-as.numeric(SampleSheet$Week)
SampleSheet<-SampleSheet[order(SampleSheet$Week),]

UniqueMonkeys<-unique(SampleSheet$Monkey)
read.csv("ProcessedData/PositionsSameStock.csv")->PositionsToUse

pdf(paste0("Output/MutFreqRelToPlasma",Sys.Date(),".pdf"))
for (M in UniqueMonkeys){
  #M=UniqueMonkeys[1]
  monkeyrows<-which(SampleSheet$Monkey==M)  
  sampleinfo=c()
  for (i in monkeyrows){sampleinfo<-c(sampleinfo, paste0(SampleSheet$Sample[i],"_w_",SampleSheet$Week[i]))}
  par(mar=c(4, 4, 4, 4))
  
  #find which sample is the first plasma sample
  plasmasamples<-grep(pattern ="plasma", x = sampleinfo)
  if (length(plasmasamples)>0){ ###CONTINUE WIT HTHIS THIUNG###
    }
  firstplasmasample<- monkeyrows[plasmasamples[1]] # the reference sample for this monkey
  othersamples<-monkeyrows[-plasmasamples[1]]  # the samples we'll compare with the plasma sample
  
  sdfplasma=SampleSheet$SeqDataFileName[firstplasmasample]
  print(sdfplasma)
  Plasma<-read.csv(paste0("ProcessedData/SeqData/",sdfplasma), row.names = 1)

  plot(x = Plasma$mutFreq, y = Plasma$mutFreq, col = 0, log="xy", xlim=c(0.001,0.5), ylim=c(0.001,0.5),
       main = paste0("Monkey ",M))
  
#  plot(x=PositionsToUse$x, y=rep(0,nrow(PositionsToUse)), ylim=c(0.00001,0.5),type="n", xlab="Env Position",ylab="Mutation Frequency",  
#       main = paste0("Monkey ",M), log="y")

abline(a = 0, b=1)

  counter=1
  for (mr in othersamples){
    sdf=SampleSheet$SeqDataFileName[mr]
    print(sdf)
    Sample<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1)
    points(x=Plasma$mutFreq,y=Sample$mutFreq,pch=16, col=alpha(cols[counter],0.5))
    #text(x = 450,y = 10^(-4+counter*0.2), labels = SampleSheet$Sample[mr], col=cols[counter])
    counter = counter +1
  }
}
dev.off()

