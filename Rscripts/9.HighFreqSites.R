source("Rscripts/BaseRscript2.R")
#OK, my next goal is to create a DF for each monkey (but starting with one 4016)
#rows are compartments and columns are nucleotide positions

read.csv("OriginalData/SampleSheetMac251AllSamples.csv", stringsAsFactors = FALSE)->SampleSheet
SeqDataFiles<-list.files("ProcessedData/",pattern="SeqData.csv$", recursive = TRUE)
SampleSheet$Week<-as.numeric(SampleSheet$Week)
SampleSheet<-SampleSheet[order(SampleSheet$Week),]

UniqueMonkeys<-unique(SampleSheet$Monkey)

#First, let's make a list of high freq sites (or high pi sites, easier)
HighPiSites <- c()
for (M in UniqueMonkeys){
  monkeyrows<-which(SampleSheet$Monkey==M)  
  sampleinfo=c()
  for (i in monkeyrows){sampleinfo<-c(sampleinfo, paste0(SampleSheet$Sample[i],"_w_",SampleSheet$Week[i]))}
  #read data for one sample
  mr= monkeyrows[1]
  sdf=SampleSheet$SeqDataFileName[mr]
  X<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1)
  DF<-data.frame(matrix(ncol=length(X$pos),nrow=length(monkeyrows), dimnames=list(sampleinfo, X$pos)))
  #read for all rows for that monkey
  counter=1
  for (mr in monkeyrows){
    sdf=SampleSheet$SeqDataFileName[mr]
    X<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1)
    DF[counter,]<-X$pi
    DF[counter, which(X$TotalReads<10000)]<-NA
    counter = counter +1
  }
  SD <- apply(DF, MARGIN = 2, sd, na.rm = TRUE)
  Mean <- apply(DF, MARGIN = 2, mean, na.rm = TRUE)
  HighPiSites<-c(HighPiSites, X$pos[which(Mean>0.02)])
}

HighPiSites<-sort(unique(HighPiSites))


pdf(paste0("Output/HighPiSitesFrequencies",Sys.Date(),".pdf"))
plot(x=PositionsToUse$x, y=rep(0,nrow(PositionsToUse)), ylim=c(0.00001,0.5),type="n", xlab="Env Position",ylab="Mutation Frequency",  
     main = paste0("Monkey ",M), log="y")

counter=1
for (M in UniqueMonkeys){
  #M="4016"
  monkeyrows<-which(SampleSheet$Monkey==M)  
  sampleinfo=c()
  for (i in monkeyrows){sampleinfo<-c(sampleinfo, paste0(SampleSheet$Sample[i],"_w_",SampleSheet$Week[i]))}

  for (mr in monkeyrows){
    sdf=SampleSheet$SeqDataFileName[mr]
    print(sdf)
    X<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1)
    X<-X[which(X$TotalReads>=10000),]
    X<-X[which(X$pos %in% HighPiSites), ]
    points(x=X$pos,y=X$a/X$TotalReads,pch=16, col=alpha(cols[counter],0.5))
    points(x=X$pos,y=X$c/X$TotalReads,pch=16, col=alpha(cols[counter],0.5))
    points(x=X$pos,y=X$g/X$TotalReads,pch=16, col=alpha(cols[counter],0.5))
    points(x=X$pos,y=X$t/X$TotalReads,pch=16, col=alpha(cols[counter],0.5))
  }
  if (M!="stockvirus")text(x = 450,y = 10^(-4+counter*0.2), labels = SampleSheet$Monkey[mr], col=cols[counter])
  if (M=="stockvirus")text(x = 450,y = 10^(-4+counter*0.2), labels = paste0(SampleSheet$Sample[mr], " ", SampleSheet$Miseq.Date[mr]), col=cols[counter])
counter = counter + 1
}
dev.off()
