


source("Rscripts/BaseRscript2.R")
read.csv("OriginalData/SampleSheetMac251AllSamples.csv", stringsAsFactors = FALSE)->SampleSheet
SampleSheet$Week<-as.numeric(SampleSheet$Week)
SampleSheet<-SampleSheet[order(SampleSheet$Week),]
UniqueMonkeys<-unique(SampleSheet$Monkey)

#View(SampleSheet)
pdf("Output/DiversityInStockVirus.pdf",width = 6, height = 8)
r=which(SampleSheet$Monkey=="stockvirus")
sdf=SampleSheet$SeqDataFileName[r]
print(sdf)
X<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1)
X$Typea<-factor(X$Typea, levels=c("nonsyn", "syn", "nonsense",  "self"))
X$Typec<-factor(X$Typec, levels=c("nonsyn", "syn", "nonsense",  "self"))
X$Typeg<-factor(X$Typeg, levels=c("nonsyn", "syn", "nonsense",  "self"))
X$Typet<-factor(X$Typet, levels=c("nonsyn", "syn", "nonsense",  "self"))
colors=c("blue","red", "black", "orange")
plot(X$pos, X$trasitionmutFreq, xaxt="n", 
     xlab="",
     yaxt="n", ylab="Transition mutation frequency",ylim=c(-0.1,0.45), 
     bty="n",
     col=0,pch=16,cex=0.7,
     main="SIV_mac_251 stock")
V1=c(120,142)
rect(xleft=(V1[1]-1)*3, ybottom = -0.01, V1[2]*3, 0.45, col="grey",border="grey")
text(x = 3*(V1[2]+(V1[1]-1))/2,y = 0.2,labels = "V1")
V2=c(198,202)
rect(xleft=(V2[1]-1)*3, ybottom = -0.01, V2[2]*3, 0.45, col="grey",border="grey")
text(x = 3*(V2[2]+(V2[1]-1))/2,y = 0.2,labels = "V2")
poslist=c()
p=which(X$TotalReads>readscutoff & X$Typea!="self")
points(X$pos[p], X$freqa[p],col=colors[as.numeric(X$Typea[p])],pch=16,cex=0.5)
poslist<-c(poslist,paste0(X$pos[which(X$TotalReads>readscutoff & X$Typea!="self"&X$freqa>0.08)],"a"))
p=which(X$TotalReads>readscutoff & X$Typec!="self")
points(X$pos[p], X$freqc[p],col=colors[as.numeric(X$Typec[p])],pch=16,cex=0.5)
poslist<-c(poslist,paste0(X$pos[which(X$TotalReads>readscutoff & X$Typec!="self"&X$freqc>0.08)],"c"))
p=which(X$TotalReads>readscutoff & X$Typeg!="self")
points(X$pos[p], X$freqg[p],col=colors[as.numeric(X$Typeg[p])],pch=16,cex=0.5)
poslist<-c(poslist,paste0(X$pos[which(X$TotalReads>readscutoff & X$Typeg!="self"&X$freqg>0.08)],"g"))
p=which(X$TotalReads>readscutoff & X$Typet!="self")
points(X$pos[p], X$freqt[p],col=colors[as.numeric(X$Typet[p])],pch=16,cex=0.5)
poslist<-c(poslist,paste0(X$pos[which(X$TotalReads>readscutoff & X$Typet!="self"&X$freqt>0.08)],"t"))
text(300,0.3,paste0(poslist,collapse = " "))
#points(X$pos[X$TotalReads>readscutoff], X$trasitionmutFreq[X$TotalReads>readscutoff],col=colors[as.factor(X$TypeOfSite[X$TotalReads>readscutoff])],pch=16,cex=0.7)
axis(side = 1,at = seq(0,900,30),labels = seq(0,300,10), line=-4.5)
axis(side = 2,at = seq(0,1,0.1),labels = seq(0,1,0.1), las=2)
mtext(text = "Env codon position", side = 1,line = -2.)
legend(x = 500, y=0.4, cex=0.8, legend = c("Non-syn","Syn", "Nonsense"), fill = colors[c(1,2,3)], border = colors[c(1,2,3)], bg = "white")
#lines(x = c((V1[1]-1),V1[2])*3, y=c(-0.1,-0.1), lwd=10, col=2)#Gp120
dev.off()


#View(SampleSheet)
pdf("Output/DiversityInMonkeyVirus.pdf",width = 6, height = 8)
for (mr in which(SampleSheet$Monkey!="stockvirus")){
sdf=SampleSheet$SeqDataFileName[mr]
print(sdf)
MV<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1)
colors=c("black","blue","red")
plot(X$pos, X$trasitionmutFreq, xaxt="n", 
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
for (p in which(X$MajNt==MV$MajNt & MV$TotalReads>readscutoff)){
  points(MV$pos[p], MV$trasitionmutFreq[p],col=colors[as.factor(MV$TypeOfSite[p])],pch=16,cex=0.7)
}
for (p in which(X$MajNt!=MV$MajNt & MV$TotalReads>readscutoff)){
  print(p)
  points(MV$pos[p], 1-MV$trasitionmutFreq[p],col=colors[as.factor(MV$TypeOfSite[p])],pch=16,cex=0.7)
}
axis(side = 1,at = seq(0,900,30),labels = seq(0,300,10), line=-4.5)
axis(side = 2,at = seq(0,1,0.1),labels = seq(0,1,0.1), las=2)
mtext(text = "Env codon position", side = 1,line = -2.)
legend(x = 500, y=0.8, cex=0.8, legend = c("Non-syn","Syn", "Nonsense"), fill = colors[c(2,3,1)], border = colors[c(2,3,1)], bg = "white")
#lines(x = c((V1[1]-1),V1[2])*3, y=c(-0.1,-0.1), lwd=10, col=2)#Gp120
}
dev.off()

