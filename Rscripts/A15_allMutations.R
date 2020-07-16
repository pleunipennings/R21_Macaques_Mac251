
source("Rscripts/BaseRscript2.R")
#SV is stockvirus, Samplesheet, UniqueMonkeys 

pdf(paste0("Output/ItaStyleFigures",Sys.Date(),".pdf"), height=12,width=8)
for (M in UniqueMonkeys){
  par(mfrow=c(2,2))
  #M=UniqueMonkeys[1]
  monkeyrows<-which(SampleSheet$Monkey==M)
  for (mr in monkeyrows){
    #par(mfrow=c(2,2))
    sdf=SampleSheet$SeqDataFileName[mr]
    print(sdf)
    MV<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1)
    levels(MV$Typea)<-c("nonsyn", "syn", "nonsense",  "self")
    levels(MV$Typec)<-c("nonsyn", "syn", "nonsense",  "self")
    levels(MV$Typeg)<-c("nonsyn", "syn", "nonsense",  "self")
    levels(MV$Typet)<-c("nonsyn", "syn", "nonsense",  "self")
    
    colors=c("black","blue","red")
    plot(MV$pos, MV$trasitionmutFreq, xaxt="n", 
         xlab="",
         yaxt="n", ylab="Transition mutation frequency",ylim=c(-0.2,1), 
         bty="n",
         col=0,pch=16,cex=0.7,
         main=paste0(SampleSheet$Monkey[mr], " ",SampleSheet$Sample[mr], " week ",  SampleSheet$Week[mr]))
    V1=c(120,142)
    rect(xleft=(V1[1]-1)*3, ybottom = -0.01, V1[2]*3, 1, col="grey",border="grey")
    text(x = 3*(V1[2]+(V1[1]-1))/2,y = 0.2,labels = "V1")
    V2=c(198,202)
    rect(xleft=(V2[1]-1)*3, ybottom = -0.01, V2[2]*3, 1, col="grey",border="grey")
    text(x = 3*(V2[2]+(V2[1]-1))/2,y = 0.2,labels = "V2")
    
    poslist=c()
    p=which(MV$TotalReads>readscutoff & MV$Typea!="self")
    points(MV$pos[p], MV$freqa[p],col=colors[as.numeric(as.factor(as.character(MV$Typea[p])))],pch=16,cex=0.7)
    poslist<-c(poslist,paste0(MV$pos[which(MV$TotalReads>readscutoff & MV$Typea!="self"&MV$freqa>0.08)],"a"))
    p=which(MV$TotalReads>readscutoff & MV$Typec!="self")
    points(MV$pos[p], MV$freqc[p],col=colors[as.numeric(as.factor(as.character(MV$Typec[p])))],pch=16,cex=0.7)
    poslist<-c(poslist,paste0(MV$pos[which(MV$TotalReads>readscutoff & MV$Typec!="self"&MV$freqc>0.08)],"c"))
    p=which(MV$TotalReads>readscutoff & MV$Typeg!="self")
    points(MV$pos[p], MV$freqg[p],col=colors[as.numeric(as.factor(as.character(MV$Typeg[p])))],pch=16,cex=0.7)
    poslist<-c(poslist,paste0(MV$pos[which(MV$TotalReads>readscutoff & MV$Typeg!="self"&MV$freqg>0.08)],"g"))
    p=which(MV$TotalReads>readscutoff & MV$Typet!="self")
    points(MV$pos[p], MV$freqt[p],col=colors[as.numeric(as.factor(as.character(MV$Typet[p])))],pch=16,cex=0.7)
    poslist<-c(poslist,paste0(MV$pos[which(MV$TotalReads>readscutoff & MV$Typet!="self"&MV$freqt>0.08)],"t"))
    text(300,0.5,paste0(poslist,collapse = " "))
    
    axis(side = 1,at = seq(0,900,30),labels = seq(0,300,10), line=-4.5)
    axis(side = 2,at = seq(0,1,0.1),labels = seq(0,1,0.1), las=2)
    mtext(text = "Env codon position", side = 1,line = -2.)
    legend(x = 440, y=0.8, cex=0.8, legend = c("Non-syn","Syn", "Nonsense"), fill = colors[c(2,3,1)], border = colors[c(2,3,1)], bg = "white")
    
  }
}
dev.off()
