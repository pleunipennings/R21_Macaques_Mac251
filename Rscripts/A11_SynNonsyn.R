#Here I will try to determine what are synonymous and non-synonymous sites. 

source("Rscripts/BaseRscript2.R")

#read.csv("OriginalData/SampleSheetMac251AllSamples.csv", stringsAsFactors = FALSE)->SampleSheet
#SampleSheet$Week<-as.numeric(SampleSheet$Week)
#SampleSheet<-SampleSheet[order(SampleSheet$Week),]
#UniqueMonkeys<-unique(SampleSheet$Monkey)
SampleSheet$pvalueSynNonSyn<-NA

pdf(paste0("Output/SynNonSyn_500totalreads",Sys.Date(),".pdf"))
#plot frequencies
for (M in UniqueMonkeys){
  #M=UniqueMonkeys[1]
  monkeyrows<-which(SampleSheet$Monkey==M)
  #monkeyrows<-monkeyrows[monkeyrows%in%grep(pattern ="plasma", x = SampleSheet$Sample)]
  #sampleinfo=c()
  #for (i in monkeyrows){sampleinfo<-c(sampleinfo, paste0(SampleSheet$Sample[i],"_w_",SampleSheet$Week[i]))}
  for (mr in monkeyrows){
    sdf=SampleSheet$SeqDataFileName[mr]
    print(sdf)
    X<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1)
    #X<-X[X$MajNt=="t",]
    
    #print(paste0(X$WTAA[seq(3,nrow(X)-2,3)],collapse = ""))
    
    plot(c(1,2,3), c(0,0,0), ylim=c(0,0.01), col=0, xlim=c(0.5, 3.5), 
         main = paste0("Monkey ",M," week ", SampleSheet$Week[mr], " comp ", SampleSheet$Sample[mr]),
          xaxt="n", xlab="type of site", ylab = "mutation frequency")
    mtext(text=c("syn", "non-syn", "nonsense"), side=1, at = 1:3)

    types=c("syn","nonsyn","nonsense")
    for (t in 1:3){
      typeofsite=types[t]
      rows= which(X$TypeOfSite==typeofsite & !is.na(X$a) & X$TotalReads>readscutoff)
      points(rep(t,length(rows))+rnorm(length(rows),mean=0,sd=0.1),X$trasitionmutFreq[rows], col=alpha(cols[t],0.5))
      #add mean
      points(t, median(X$trasitionmutFreq[rows]), col=cols[t], pch=15, cex=2)
    }
    syn_array<-X$trasitionmutFreq[X$TypeOfSite=="syn" & !is.na(X$a) & X$TotalReads>readscutoff ]
    nonsyn_array<-X$trasitionmutFreq[X$TypeOfSite!="syn" & !is.na(X$a) & X$TotalReads>readscutoff ]
    if (length(syn_array)>0&length(nonsyn_array)>0){
    
    text(x=2, y=0.008, label=paste0("pvalue = ",
                                    round(wilcox.test(syn_array, nonsyn_array, alternative = "gr")$p.value,3)))
  
    SampleSheet$pvalueSynNonSyn[mr]<-wilcox.test(syn_array, nonsyn_array, alternative = "gr")$p.value
    }
  }
}
dev.off()

#pdf(paste0("Output/SynNonSynPvaluePlasmaVsRest_6000Totalreads",Sys.Date(),".pdf"))

plot(c(1,2), c(0,0), ylim=c(0,0.6), col=0, xlim=c(0.5, 2.5), 
     main = "pvalues comparison syn vs nonsyn sites")
for (i in 1:nrow(SampleSheet)){
  print(paste0(c(SampleSheet$Sample[i],SampleSheet$pvalueSynNonSyn[i])))
  x=2
  if (SampleSheet$Sample[i]=="plasma")x=1
  points(x, SampleSheet$pvalueSynNonSyn[i])
}
points(1, mean(SampleSheet$pvalueSynNonSyn[SampleSheet$Sample=="plasma"]), col=2,pch=16, cex=2)
points(2, mean(SampleSheet$pvalueSynNonSyn[SampleSheet$Sample!="plasma"], na.rm=TRUE), col=2,pch=16, cex=2)

wilcox.test(SampleSheet$pvalueSynNonSyn[SampleSheet$Sample=="plasma"], SampleSheet$pvalueSynNonSyn[SampleSheet$Sample!="plasma"])
#dev.off()

write.csv(x=SampleSheet,file = "OriginalData/SampleSheetMac251AllSamples.csv", row.names = FALSE)

#how many plasma significant
print(paste0("plasma: ", 
length(which(SampleSheet$pvalueSynNonSyn[SampleSheet$Sample=="plasma"]<0.05)), 
" out of ", 
length(which(SampleSheet$Sample=="plasma"& !is.na(SampleSheet$pvalueSynNonSyn))), " significant"))

print(length(which(SampleSheet$pvalueSynNonSyn[SampleSheet$Sample=="plasma"]<0.05))/length(which(SampleSheet$Sample=="plasma"& !is.na(SampleSheet$pvalueSynNonSyn))))

print(paste0("non-plasma: ", 
             length(which(SampleSheet$pvalueSynNonSyn[SampleSheet$Sample!="plasma" & SampleSheet$Sample!="stockvirus"]<0.05)), 
             " out of ", 
             length(which(SampleSheet$Sample!="plasma" & SampleSheet$Sample!="stockvirus" & !is.na(SampleSheet$pvalueSynNonSyn))), " significant"))

print(length(which(SampleSheet$pvalueSynNonSyn[SampleSheet$Sample!="plasma"]<0.05))/length(which(SampleSheet$Sample!="plasma"& !is.na(SampleSheet$pvalueSynNonSyn))))

SampleSheetLN<-SampleSheet[grep(pattern ="LN", SampleSheet$Sample),]
print(length(which(SampleSheetLN$pvalueSynNonSyn<0.05))/length(which(!is.na(SampleSheetLN$pvalueSynNonSyn))))

SampleSheetLung<-SampleSheet[-grep(pattern ="LN", SampleSheet$Sample),]
SampleSheetLung<-SampleSheetLung[which(SampleSheetLung$Sample!="plasma"&SampleSheetLung$Sample!="stockvirus"&SampleSheetLung$Sample!="control"),]
print(length(which(SampleSheetLung$pvalueSynNonSyn<0.05))/length(which(!is.na(SampleSheetLung$pvalueSynNonSyn))))

