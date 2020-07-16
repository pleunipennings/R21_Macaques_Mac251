#OK, my next goal is to create a DF for each monkey (but starting with one 4016)
#rows are compartments and columns are nucleotide positions

source("Rscripts/BaseRscript2.R")
#SV is stockvirus, Samplesheet, UniqueMonkeys 

head(SV)
plot(sort(SV$mutFreq[SV$TotalReads>readscutoff]), log="y")
HighFreqSites<-SV$pos[which(SV$mutFreq>0.008 & SV$TotalReads>readscutoff)]

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
  HighPiSites<-c(HighPiSites, X$pos[which(Mean>0.015 & SD<0.03)])
}
HighPiSites<-sort(unique(HighPiSites))

#June 2020
#Let's look at the sites that have high freq in stock 
HighPiSites<-HighFreqSites

#Make one plot with all high freq sites in all monkeys
#pdf(paste0("Output/HighPiSitesFrequencies",Sys.Date(),".pdf"))
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
#dev.off()

#Make one plot per monkey? Show sites per week
pdf(paste0("Output/HighFreqSitesFrequencies_Time",Sys.Date(),".pdf"))
for (M in UniqueMonkeys){
  #M="3616"
  print(M)
  monkeyrows<-which(SampleSheet$Monkey==M)  
  sampleinfo=c()
  for (i in monkeyrows){sampleinfo<-c(sampleinfo, paste0(SampleSheet$Sample[i],"_w_",SampleSheet$Week[i]))}
  
  DFmonkey<-data.frame(pos = HighPiSites)
  
  counter=1
  for (mr in monkeyrows){
    sdf=SampleSheet$SeqDataFileName[mr]
    print(sdf)
    X<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1)
    X<-X[which(X$TotalReads>=10000),]
    X<-X[which(X$pos %in% HighPiSites), ]
    if (nrow(X)==length(HighPiSites)){
      colname = paste0(SampleSheet$Sample[mr], "_w_",SampleSheet$Week[mr])
      print(colname)
      DFmonkey[[colname]]<-X$pi
      counter=counter+1
    }}
  
  lengthDF<-ncol(DFmonkey)
  plot(1:(lengthDF-1), DFmonkey[1,2:lengthDF], 
       t="n", pch=20, xlab="", xaxt="n", 
       main=M, ylim=c(0.0001,0.04))
  axis(side = 1, at = 1:length(sampleinfo), labels = FALSE)
  #Not working yet! 
  text(x=1:length(sampleinfo), y=rep(0,length(sampleinfo)),
       labels=names(DFmonkey)[2:(length(sampleinfo)+1)], srt=45, adj=1, xpd=TRUE)
  for (i in 1:nrow(DFmonkey)){
    i=i+1
    print(paste0(i, DFmonkey$pos[i]))
    points(1:(lengthDF-1), DFmonkey[i,2:lengthDF], t="b", pch=2, col=cols[i%%length(cols)])
  }   
}
dev.off()

#Make one plot per monkey? Show sites per week. 
#Show all sites, not just the hifreq sites. 

pdf(paste0("Output/AllSitesFrequencies_Time",Sys.Date(),".pdf"))
for (M in UniqueMonkeys){
  #M="3616"
  print(M)
  monkeyrows<-which(SampleSheet$Monkey==M)  
  sampleinfo=c()
  for (i in monkeyrows){sampleinfo<-c(sampleinfo, paste0(SampleSheet$Sample[i],"_w_",SampleSheet$Week[i]))}
  
  counter=1
  for (mr in monkeyrows){
    sdf=SampleSheet$SeqDataFileName[mr]
    print(sdf)
    X<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1)
    #X<-X[which(X$TotalReads>=10000),]
    if (counter == 1) DFmonkey<-data.frame(pos = X$pos)
    #X<-X[which(X$pos %in% HighPiSites), ]
    colname = paste0(SampleSheet$Sample[mr], "_w_", SampleSheet$Week[mr])
    DFmonkey[[colname]]<-X$pi
    counter=counter+1
    }
  
  lengthDF<-ncol(DFmonkey)
  plot(1:(lengthDF-1), DFmonkey[1,2:lengthDF], 
       t="n", pch=20, xlab="", xaxt="n", 
       main=M, ylim=c(0.0001,0.1))
  axis(side = 1, at = 1:length(sampleinfo), labels = FALSE)
  #Not working yet! 
  text(x=1:length(sampleinfo), y=rep(0,length(sampleinfo)),
       labels=names(DFmonkey)[2:(length(sampleinfo)+1)], srt=45, adj=1, xpd=TRUE)
  for (i in 1:nrow(DFmonkey)){
    #print(i)
    points(1:(lengthDF-1), DFmonkey[i,2:lengthDF], t="b", pch=2, col=cols[i%%length(cols)+1])
  }   
}
dev.off()

