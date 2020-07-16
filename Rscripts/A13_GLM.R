
#Collect one large data frame with all frequencies incl info on plasma / lungs / LN

source("Rscripts/BaseRscript2.R")

read.csv("OriginalData/SampleSheetMac251AllSamples.csv", stringsAsFactors = FALSE)->SampleSheet
SampleSheet$Week<-as.numeric(SampleSheet$Week)
SampleSheet<-SampleSheet[order(SampleSheet$Week),]
UniqueMonkeys<-unique(SampleSheet$Monkey)

sdf=SampleSheet$SeqDataFileName[1]
print(sdf)
X<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1)
Y<-subset(x=X, TotalReads>readscutoff & !is.na(TypeOfSite), select = c(TotalReads, MajNt, trasitionmutFreq, TypeOfSite,makesCpG))
LongDF<-Y[1,]

for (sdf in SampleSheet$SeqDataFileName[SampleSheet$Monkey!="stockvirus"&SampleSheet$Monkey!="control"]){
  print(sdf)
  X<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1)
  Y<-subset(x=X, TotalReads>readscutoff & !is.na(TypeOfSite), select = c(TotalReads, MajNt, trasitionmutFreq, TypeOfSite,makesCpG))
  LongDF<-rbind(LongDF,Y)
  print(dim(LongDF))
}

pdf("Output/LongDF.pdf")
par(mfrow=c(2,2))
for (nuc in c("a","c","g","t")){

LongDF_x = LongDF[LongDF$MajNt==nuc,]

plot(c(1,2,3), c(0,0,0), ylim=c(0,0.02), col=0, xlim=c(0.5, 3.5), 
     main = paste0("all samples together, nuc = ", nuc),
     xaxt="n", xlab="type of site", ylab = "mutation frequency")
mtext(text=c("syn", "non-syn", "nonsense"), side=1, at = 1:3)

types=c("syn","nonsyn","nonsense")
for (t in 1:3){
  typeofsite=types[t]
  rows= which(LongDF_x$TypeOfSite==typeofsite)
  points(rep(t,length(rows))+rnorm(length(rows),mean=0,sd=0.1),LongDF_x$trasitionmutFreq[rows], col=alpha(cols[t],0.1),pch=16)
  #add mean
  points(t, median(LongDF_x$trasitionmutFreq[rows]), col=1, pch="--", cex=2)
}
syn_array<-LongDF_x$trasitionmutFreq[LongDF_x$TypeOfSite=="syn"]
nonsyn_array<-X$trasitionmutFreq[LongDF_x$TypeOfSite!="syn"]
if (length(syn_array)>0&length(nonsyn_array)>0){
  
  text(x=2, y=0.018, label=paste0("pvalue = ",
                                  round(wilcox.test(syn_array, nonsyn_array, alternative = "gr")$p.value,3)))
  #SampleSheet$pvalueSynNonSyn[mr]<-wilcox.test(syn_array, nonsyn_array, alternative = "gr")$p.value
}
}
dev.off()



LongDF_x = LongDF[LongDF$MajNt=="a",]
plot(c(1,2,3,4), c(0,0,0,0), ylim=c(0,0.01), col=0, xlim=c(0.5, 4.5), 
     main = paste0("CpG"), 
     xaxt="n", xlab="type of site")
mtext(text=c("syn non CpG", "syn CpG", "nonsyn non CpG", "nonsyn nonCpG"), side=1, at = 1:4)

syn_nonCpG_array<-LongDF_x$trasitionmutFreq[LongDF_x$TypeOfSite=="syn" & LongDF_x$makesCpG==0 ]
syn_CpG_array<-LongDF_x$trasitionmutFreq[LongDF_x$TypeOfSite=="syn" & LongDF_x$makesCpG==1 ]
nonsyn_nonCpG_array<-LongDF_x$trasitionmutFreq[LongDF_x$TypeOfSite!="syn" & LongDF_x$makesCpG==0 ]
nonsyn_CpG_array<-LongDF_x$trasitionmutFreq[LongDF_x$TypeOfSite!="syn" & LongDF_x$makesCpG==1 ]

t=1; arr= syn_nonCpG_array
points(rep(t,length(arr))+rnorm(length(arr),mean=0,sd=0.1),arr, col=alpha(cols[t],0.5))
points(t, median(arr), col=1, pch="--", cex=2)
t=2; arr= syn_CpG_array
points(rep(t,length(arr))+rnorm(length(arr),mean=0,sd=0.1),arr, col=alpha(cols[t],0.5))
points(t, median(arr), col=1, pch="--", cex=2)
t=3; arr= nonsyn_nonCpG_array
points(rep(t,length(arr))+rnorm(length(arr),mean=0,sd=0.1),arr, col=alpha(cols[t],0.5))
points(t, median(arr), col=1, pch="--", cex=2)
t=4; arr= nonsyn_CpG_array
points(rep(t,length(arr))+rnorm(length(arr),mean=0,sd=0.1),arr, col=alpha(cols[t],0.5))
points(t, median(arr), col=1, pch="--", cex=2)

if (length(syn_nonCpG_array)>0&length(syn_CpG_array)>0){
  text(x=1.5, y=0.008, label=paste0("1-sided pvalue = ",
                                    round(wilcox.test(syn_nonCpG_array, syn_CpG_array, alternative = "gr")$p.value,3)))
  
}
if (length(nonsyn_nonCpG_array)>0&length(nonsyn_CpG_array)>0){
  text(x=3.5, y=0.008, label=paste0("1-sided pvalue = ",
                                    round(wilcox.test(nonsyn_nonCpG_array, nonsyn_CpG_array, alternative = "gr")$p.value,3)))
  
}

#library (betareg)
LongDF<-LongDF[LongDF$MajNt=="a",]
model1<-glm(formula = LongDF_x$trasitionmutFreq ~ LongDF_x$TypeOfSite*LongDF_x$makesCpG, family = "quasi")
summary(model1)


