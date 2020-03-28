#OK, my next goal is to plot average pi for each sample. 
library(dplyr)

read.csv("OriginalData/SampleSheetMac251AllSamples.csv", stringsAsFactors = FALSE)->SampleSheet
SeqDataFiles<-list.files("ProcessedData/",pattern="SeqData.csv$", recursive = TRUE)

SampleSheet$Week<-as.numeric(SampleSheet$Week)
SampleSheet<-SampleSheet[order(SampleSheet$Week),]
SampleSheet$avePi <-0

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
    #Now add average pi for each sample to the DF
    SampleSheet$avePi[mr]<-avePi
    #Add to figure
    points(counter,avePi,pch=16)
    counter = counter +1
  }
}
dev.off()

#Next, I want to recreate the graphic that Kaho made, that has the diversity of all samples in one plot. 

#Let me make a dataframe that has data for each monkey. 
#It'll be based on the sample sheet, but with less metadata and more sequence related data

SampleSheet %>% select(Monkey, Sample, wks.post.SIV, Week, Necropsie, avePi) ->OverviewDF
UniqueMonkeys<-UniqueMonkeys[order(UniqueMonkeys)]

#Order like Kaho did

UniqueMonkeys <-
c( "16314" ,     "20615"  , "3516", "30816"  ,   "3216",  "31316", "3816" , "3616",  "3316" , "4016"  , "3116",   "stockvirus")

c("3816" , "3616",  "3316" , "4016" ) -> SIVonly

pdf(paste0("Output/DiversityPiInMonkeys_oneplot",Sys.Date(),".pdf"), width = 16, height = 7)

par(mar=c(10, 4, 4, 4))

plot(1 :nrow(OverviewDF), OverviewDF$avePi, type = "n", ylim = c(0.004, max(OverviewDF$avePi)), 
     xlab="", xaxt="n")
counter = 1
for (M in UniqueMonkeys){
  #M=UniqueMonkeys[1]
  abline(v = counter-0.5)
  monkeyrows<-which(OverviewDF$Monkey==M)  
  text(x = counter-0.5 + length(monkeyrows)/2, y=0.002, labels = M, srt= 90)
  if (M!="stockvirus")mtext(text = M, side = 1, at = counter-0.5 + length(monkeyrows)/2, line = 8)
  for (mr in monkeyrows){
    #plot
    col=2
    if (M %in% SIVonly ) col = 4
    points (counter,OverviewDF$avePi[mr], pch=1, col=col)
    if (OverviewDF$Sample[mr]=="plasma") points (counter,OverviewDF$avePi[mr], pch=16, col= col)
    mtext(text = OverviewDF$Sample[mr], side = 1, at = counter, las =2)
    counter = counter +1
  }}
dev.off()

