#OK, my next goal is to determins Fst between lung / LN / plasma in all monkeys

source("Rscripts/BaseRscript2.R")
#SV is stockvirus, Samplesheet, UniqueMonkeys 
nucleotides =c("a", "c", "g", "t")

#I will focus only on the high freq sites (from plasma)
head(SV)
plot(sort(SV$mutFreq[SV$TotalReads>readscutoff]), log="y")
HighFreqSites<-SV$pos[which(SV$mutFreq>0.008 & SV$TotalReads>readscutoff)]

#Plots frequencies
if (FALSE){
for (s in HighFreqSites){
  #for (s in c(347, 345, 351, 348)){  
  #s = HighFreqSites[1] #focus on one site first
  mutnucs = nucleotides[nucleotides!=SV$StockMajNt[SV$pos==s]]
  highfreqnuc = mutnucs[which.max(SV[SV$pos==s,mutnucs])] #which nuc has the highest freq?
  
  #pdf(paste0("Output/HiFreqSite",s,"AllMonkeys",Sys.Date(),".pdf"), height=12,width=8)
  par(mfrow=c(3,2))
  for (M in UniqueMonkeys){
    #M=UniqueMonkeys[1]
    monkeyrows<-which(SampleSheet$Monkey==M)
    
    plot(1:2, 1:2, col=0, ylim=c(0,0.2), xlim=c(0.8,2.8), xaxt="n", xlab = "Stock vs Monkey", ylab= "Freq", 
         main =paste0("monkey ",M, ", site ", s, ", mut ", SV$StockMajNt[SV$pos==s], "->", highfreqnuc, ", ", SV$TypeOfSite[SV$pos==s]))
    axis(side = 1, at = c(1, 1.5, 2, 2.5), labels = c("stock", "plasma", "LN", "lung"))
    points(1, SV[SV$pos==s,highfreqnuc]/SV$TotalReads[SV$pos==s], col=1, pch=16, cex=2)
    abline(h=SV[SV$pos==s,highfreqnuc]/SV$TotalReads[SV$pos==s])
    
    for (mr in monkeyrows){
      #par(mfrow=c(2,2))
      #  mr=monkeyrows[1]
      sdf=SampleSheet$SeqDataFileName[mr]
      print(sdf)
      MV<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1)
      if (MV$TotalReads[MV$pos==s]>readscutoff){ #only plot if sufficient reads
        color=3; xpos=2.5 #for lung
        if (mr %in% grep(pattern ="LN", SampleSheet$Sample)) {color=2; xpos=2} #LN at 2, in red
        if (SampleSheet$Sample[mr]=="plasma") {color=4; xpos=1.5}
        points(xpos, MV[MV$pos==s,highfreqnuc]/MV$TotalReads[MV$pos==s], col=color,pch=1,cex=2) 
      }
    }
    
  }
  #dev.off()
}
}

UniqueMonkeys<-UniqueMonkeys[UniqueMonkeys!="3116" & UniqueMonkeys!="3316"] #because week info missing for this monkey
#Now, let's calculate delta frequency between early plasma and the other points
#for (s in HighFreqSites){

AllDeltaF_plasma<-c()
AllDeltaF_LN<-c()
AllDeltaF_Lung<-c()
DeltaWeek<-c()

pdf(paste0("Output/DeltaFHiFreqSite",Sys.Date(),".pdf"), height=12,width=8)
par(mfrow=c(3,2))
for (count in 1:length(HighFreqSites)){
  s = HighFreqSites[count] #focus on one site first
  mutnucs = nucleotides[nucleotides!=SV$StockMajNt[SV$pos==s]]
  highfreqnuc = mutnucs[which.max(SV[SV$pos==s,mutnucs])] #which nuc has the highest freq?
  
  for (M in UniqueMonkeys){
    print(M)
    #M=UniqueMonkeys[3]
    monkeyrows<-which(SampleSheet$Monkey==M)
    #Determine which of the monkey rows is the early plasma sample
    earlyWeek = min(SampleSheet$Week[monkeyrows])
    print(earlyWeek)
    rowEarlyPlasma<-which(SampleSheet$Sample=="plasma" & SampleSheet$Monkey==M & SampleSheet$Week == earlyWeek) 
    sdf=SampleSheet$SeqDataFileName[rowEarlyPlasma]
    MV<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1)
    #if (MV$TotalReads[MV$pos==s]>readscutoff){ #only plot if sufficient reads
    f_early_plasma = MV[MV$pos==s,highfreqnuc]/MV$TotalReads[MV$pos==s]
    
    lateWeek = max(SampleSheet$Week[SampleSheet$Monkey==M & SampleSheet$Sample=="plasma" ])
    print(lateWeek)
    rowLatePlasma<-which(SampleSheet$Sample=="plasma" & SampleSheet$Monkey==M & SampleSheet$Week == lateWeek) 
    sdf=SampleSheet$SeqDataFileName[rowLatePlasma]
    MV<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1)
    #if (MV$TotalReads[MV$pos==s]>readscutoff){ #only plot if sufficient reads
    f_late_plasma = MV[MV$pos==s,highfreqnuc]/MV$TotalReads[MV$pos==s]
    deltaF_plasma = abs(f_early_plasma-f_late_plasma)
    
    rowsLN<- monkeyrows[which(monkeyrows %in% grep(pattern ="LN", SampleSheet$Sample))]
    f_LN<-c()
    for (r in rowsLN){
      sdf=SampleSheet$SeqDataFileName[r]
      MV<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1)
      if (MV$TotalReads[MV$pos==s]>readscutoff){ #only plot if sufficient reads
        f_LN = c(f_LN, MV[MV$pos==s,highfreqnuc]/MV$TotalReads[MV$pos==s])
      }}
    deltaF_LN = abs(f_LN-f_early_plasma)
    
    rowsLung<- monkeyrows[which(!monkeyrows %in% grep(pattern ="LN", SampleSheet$Sample) & !monkeyrows %in% grep(pattern ="plasma", SampleSheet$Sample))]
    f_Lung<-c()
    for (r in rowsLung){
      sdf=SampleSheet$SeqDataFileName[r]
      MV<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1)
      if (MV$TotalReads[MV$pos==s]>readscutoff){ #only plot if sufficient reads
        f_Lung = c(f_Lung, MV[MV$pos==s,highfreqnuc]/MV$TotalReads[MV$pos==s])
      }}
    deltaF_Lung = abs(f_Lung-f_early_plasma)
    
    if (earlyWeek !=lateWeek){
    AllDeltaF_plasma<-c(AllDeltaF_plasma,deltaF_plasma)
    AllDeltaF_LN<-c(AllDeltaF_LN,deltaF_LN)
    AllDeltaF_Lung<-c(AllDeltaF_Lung,deltaF_Lung)
    DeltaWeek<-c(DeltaWeek, lateWeek-earlyWeek)
    }
    
    plot(1:2, 1:2, col=0, ylim=c(0,1.1*max(deltaF_Lung,deltaF_LN,deltaF_plasma)), xlim=c(1.3,2.8), xaxt="n", xlab = "Stock vs Monkey", ylab= "Freq", 
         main =paste0("monkey ",M, ", site ", s, ", mut ", SV$StockMajNt[SV$pos==s], "->", highfreqnuc, ", ", SV$TypeOfSite[SV$pos==s]))
    axis(side = 1, at = c(1.5, 2, 2.5), labels = c("plasma", "LN", "lung"))
    points(1.5, deltaF_plasma, col=1, pch=count, cex=2)
    points(rep(2,length(deltaF_LN)), deltaF_LN, col=2, pch=count, cex=2)
    points(rep(2.5,length(deltaF_Lung)), deltaF_Lung, col=3, pch=count, cex=2)
  }
}
dev.off()

par(mfrow=c(1,1))

plot(1:2, 1:2, col=0, ylim=c(0.00001,1.1*max(AllDeltaF_plasma,AllDeltaF_LN,AllDeltaF_Lung)), xlim=c(1.3,2.8), 
     xaxt="n", xlab = "Stock vs Monkey", ylab= "Relative difference in frequency", 
     main =paste0("Relative difference in frequency at high frequency sites"),log="y")
axis(side = 1, at = c(1.5, 2, 2.5), labels = c("plasma", "LN", "lung"))
points(rep(1.5,length(AllDeltaF_plasma))+rnorm(n = length(AllDeltaF_plasma),mean = 0,sd = 0.03), AllDeltaF_plasma, col=alpha(1,0.25), pch=16, cex=1)
points(rep(2,length(AllDeltaF_LN))+rnorm(n = length(AllDeltaF_LN),mean = 0,sd = 0.03), AllDeltaF_LN, col=alpha(2,0.25), pch=16, cex=1)
points(rep(2.5,length(AllDeltaF_Lung))+rnorm(n = length(AllDeltaF_Lung),mean = 0,sd = 0.03), AllDeltaF_Lung, col=alpha(4,0.25), pch=16, cex=1)

points(1.5, mean(AllDeltaF_plasma), col=2, pch=21, bg="black",cex=2)
points(2,mean(AllDeltaF_LN), col=4, pch=21, bg="black", cex=2)
points(2.5,mean(AllDeltaF_Lung), col=2, pch=21, bg="black",cex=2)

AllDeltaF_plasma
AllDeltaF_LN
AllDeltaF_Lung
DeltaWeek

wilcox.test(AllDeltaF_plasma,AllDeltaF_LN, alternative = "l")
wilcox.test(AllDeltaF_plasma,AllDeltaF_Lung, alternative = "l")
wilcox.test(AllDeltaF_LN,AllDeltaF_Lung, alternative = "t")

#there is less drift in plasma than between plasma and tissues. 
#no difference between LN and lung. 

