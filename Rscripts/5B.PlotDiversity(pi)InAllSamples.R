#OK, my next goal is to plot average pi for each sample. 
library(dplyr)
source("Rscripts/baseRscript2.R")
#SeqDataFiles<-list.files("ProcessedData/",pattern="SeqData.csv$", recursive = TRUE)
SampleSheet$avePi <-0
SampleSheet$aveTransMutFreq<-0
UniqueMonkeys<-unique(SampleSheet$Monkey)
#UniqueMonkeys<-UniqueMonkeys[UniqueMonkeys!="stockvirus"&UniqueMonkeys!="control"]

#pdf(paste0("Output/DiversityPiInMonkeys",Sys.Date(),".pdf"))
#calculate pi for each sample
for (M in UniqueMonkeys){
  #M=UniqueMonkeys[1]
  monkeyrows<-which(SampleSheet$Monkey==M)  
  sampleinfo=c()
  for (i in monkeyrows){sampleinfo<-c(sampleinfo, paste0(SampleSheet$Sample[i],"_w_",SampleSheet$Week[i]))}
  
  #par(mar=c(10, 4, 4, 8) + 0.1)
  #plot(x=1:length(monkeyrows), y=1:length(monkeyrows), ylim=c(0,0.01),type="n", xaxt="n",ylab="diversity(pi)", xlab = "", 
  #     main = paste0("Monkey ",M))
  #     axis(1, at=1:length(monkeyrows),labels=sampleinfo, col.axis="red", las=2)
  
  counter=1
  for (mr in monkeyrows){
    sdf=SampleSheet$SeqDataFileName[mr]
    print(sdf)
    X<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1)
    X<-X[X$TotalReads>readscutoff,]
    avePi=mean(X$pi,na.rm=TRUE)
    aveTransMutFreq=mean(X$trasitionmutFreq, na.rm=TRUE)
    #Now add average pi for each sample to the DF
    SampleSheet$avePi[mr]<-avePi
    SampleSheet$aveTransMutFreq[mr]<-aveTransMutFreq
    #Add to figure
    #points(counter,avePi,pch=16)
    counter = counter +1
  }
}
#dev.off()

#Next, I want to recreate the graphic that Kaho made, that has the diversity of all samples in one plot. 

#Let me make a dataframe that has data for each monkey. 
#It'll be based on the sample sheet, but with less metadata and more sequence related data

#Order like Kaho did
#but with TB infected monkeys first
UniqueMonkeys <-
c( "16314" , "20615"  , "3516", "30816"  ,  "3116", "3216",  "31316", "3816" , "3616",  "3316" , "4016"  ,   "stockvirus")

c("3816" , "3616",  "3316" , "4016" ) -> SIVonly
c( "16314" , "20615"  , "3516", "30816"  ,  "3116", "3216",  "31316") ->coinfected

if (TRUE){
pdf(paste0("Output/DiversityPiInMonkeys_oneplot",Sys.Date(),".pdf"), width = 16, height = 7)
par(mar=c(10, 4, 4, 4))
plot(1 :nrow(SampleSheet), SampleSheet$avePi, type = "n", 
     xlim = c(3.2,nrow(SampleSheet)-3.15),
     ylim = c(0.0004, max(SampleSheet$avePi, na.rm=TRUE)), 
     xlab="", xaxt="n", ylab = "Genetic diversity", main = "Genetic diversity in all samples")
counter = 1
for (M in UniqueMonkeys){
#M=UniqueMonkeys[1]
  abline(v = counter-0.5)
  monkeyrows<-which(SampleSheet$Monkey==M)  
  if (M!="stockvirus" & M!="control") text(x = counter-0.5 + length(monkeyrows)/2, y=0.002, labels = "Monkey", srt= 0)
  if (M!="stockvirus" & M!="control") text(x = counter-0.5 + length(monkeyrows)/2, y=0.0016, labels = M, srt= 0)
  #if (M!="stockvirus" & M!="control") mtext(text = M, side = 1, at = counter-0.5 + length(monkeyrows)/2, line = 8)
  for (mr in monkeyrows){
    #plot
    col=2
    if (M %in% SIVonly ) col = 4
    #points (counter,SampleSheet$avePi[mr], pch=1, col=col)
    if (SampleSheet$Tissue[mr]=="Plasma") points (counter,SampleSheet$avePi[mr], pch=16, col= col,cex=1.2)
    if (SampleSheet$Tissue[mr]=="LN") points (counter,SampleSheet$avePi[mr], pch=1, col= col)
    if (SampleSheet$Tissue[mr]=="Lung") points (counter,SampleSheet$avePi[mr], pch=2, col= col)
    if (SampleSheet$Sample[mr]=="stockvirus") points (counter,SampleSheet$avePi[mr], pch=4, col= 1)
    if (SampleSheet$Tissue[mr]!="Plasma") mtext(text = SampleSheet$Sample[mr], side = 1, at = counter, las =2, line = 0.5)
    if (SampleSheet$Tissue[mr]=="Plasma") mtext(text = paste0(SampleSheet$Sample[mr],"_w",SampleSheet$Week[mr]), side = 1, at = counter, las =2, line = 0.5)
    counter = counter +1
  }
}
dev.off()
}

unique(SampleSheet$Monkey[grep(pattern = "Latent", x = SampleSheet$Notes)])

#Sign difference between tissue and plasma? TB non TB?
wilcox.test(SampleSheet$avePi[SampleSheet$Tissue=="Plasma"],SampleSheet$avePi[SampleSheet$Tissue=="Lung"|SampleSheet$Tissue=="LN"], alternative ="tw") 
wilcox.test(SampleSheet$avePi[SampleSheet$Tissue=="Plasma"&SampleSheet$Monkey %in% SIVonly], SampleSheet$avePi[SampleSheet$Tissue=="Plasma"& SampleSheet$Monkey %in% coinfected], alternative ="less") 
wilcox.test(SampleSheet$avePi[SampleSheet$Monkey %in% SIVonly], SampleSheet$avePi[SampleSheet$Monkey %in% coinfected], alternative ="less") 
wilcox.test(rep(0,100), rep(1,100), alternative = "less")

mean(SampleSheet$avePi[SampleSheet$Tissue=="Plasma"])*100
mean(SampleSheet$avePi[SampleSheet$Tissue=="Lung"|SampleSheet$Tissue=="LN"])*100
mean(SampleSheet$avePi[SampleSheet$Tissue=="Plasma"&SampleSheet$Monkey %in% SIVonly])*100
mean(SampleSheet$avePi[SampleSheet$Tissue=="Plasma"&SampleSheet$Monkey %in% coinfected])*100

if (FALSE){
  pdf(paste0("Output/TransitionMutFreqInMonkeys_oneplot",Sys.Date(),".pdf"), width = 16, height = 7)
  par(mar=c(10, 4, 4, 4))
  plot(1 :nrow(SampleSheet), SampleSheet$aveTransMutFreq, type = "n", 
       xlim = c(3.2,nrow(SampleSheet)-3.15),
       ylim = c(0.0004, max(SampleSheet$aveTransMutFreq, na.rm=TRUE)), 
       xlab="", xaxt="n", ylab = "average transition mutation frequency",
       main = "Average transition mutation frequency in samples")
  counter = 1
  for (M in UniqueMonkeys){
    #M=UniqueMonkeys[1]
    abline(v = counter-0.5)
    monkeyrows<-which(SampleSheet$Monkey==M)  
    if (M!="stockvirus" & M!="control") text(x = counter-0.5 + length(monkeyrows)/2, y=0.00065, labels = "Monkey", srt= 0)
    if (M!="stockvirus" & M!="control") text(x = counter-0.5 + length(monkeyrows)/2, y=0.0005, labels = M, srt= 0)
    #if (M!="stockvirus" & M!="control") mtext(text = M, side = 1, at = counter-0.5 + length(monkeyrows)/2, line = 8)
    for (mr in monkeyrows){
      #plot
      col=2
      if (M %in% SIVonly ) col = 4
      #points (counter,SampleSheet$aveTransMutFreq[mr], pch=1, col=col)
      if (SampleSheet$Tissue[mr]=="Plasma") points (counter,SampleSheet$aveTransMutFreq[mr], pch=16, col= col,cex=1.2)
      if (SampleSheet$Tissue[mr]=="LN") points (counter,SampleSheet$aveTransMutFreq[mr], pch=1, col= col)
      if (SampleSheet$Tissue[mr]=="Lung") points (counter,SampleSheet$aveTransMutFreq[mr], pch=2, col= col)
      if (SampleSheet$Sample[mr]=="stockvirus") points (counter,SampleSheet$aveTransMutFreq[mr], pch=4, col= 1)
      if (SampleSheet$Tissue[mr]!="Plasma") mtext(text = SampleSheet$Sample[mr], side = 1, at = counter, las =2, line = 0.5)
      if (SampleSheet$Tissue[mr]=="Plasma") mtext(text = paste0(SampleSheet$Sample[mr],"_w",SampleSheet$Week[mr]), side = 1, at = counter, las =2, line = 0.5)
      counter = counter +1
    }
  }
  dev.off()
}

