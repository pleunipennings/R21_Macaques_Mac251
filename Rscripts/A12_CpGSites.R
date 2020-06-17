#Here I will try to determine whether CpG creating sites have lower frequencies. 

source("Rscripts/BaseRscript2.R")

read.csv("OriginalData/SampleSheetMac251AllSamples.csv", stringsAsFactors = FALSE)->SampleSheet
#SeqDataFiles<-list.files("ProcessedData/",pattern="SeqData.csv$", recursive = TRUE)
SampleSheet$Week<-as.numeric(SampleSheet$Week)
SampleSheet<-SampleSheet[order(SampleSheet$Week),]
UniqueMonkeys<-unique(SampleSheet$Monkey)
SampleSheet$pvalueNonsyn_CpG<-NA
SampleSheet$pvalueSyn_CpG<-NA

#Calculate P values for CpG effect and plot 
if (TRUE){
pdf(paste0("Output/CpG_6000totalreads_DataCpG_OnlyAT",Sys.Date(),".pdf"))
readscutoff=6000
for (M in UniqueMonkeys){
#  M=UniqueMonkeys[1]
  monkeyrows<-which(SampleSheet$Monkey==M)
  for (mr in monkeyrows){
    sdf=SampleSheet$SeqDataFileName[mr]
    print(sdf)
    X<-read.csv(paste0("ProcessedData/DataCpG/",sdf), row.names = 1)
    X<-X[X$MajNt%in%c("a","t"),] #Only focus on A and T sites , becyase they can make CpG sites
    if (!is.na(SampleSheet$pvalueSynNonSyn[mr])&SampleSheet$pvalueSynNonSyn[mr]<0.05){ #Make plot
    plot(c(1,2,3,4), c(0,0,0,0), ylim=c(0,0.01), col=0, xlim=c(0.5, 4.5), 
         main = paste0("Monkey ",M," week ", SampleSheet$Week[mr], " comp ", SampleSheet$Sample[mr]), 
         xaxt="n", xlab="type of site")
    mtext(text=c("syn non CpG", "syn CpG", "nonsyn non CpG", "nonsyn nonCpG"), side=1, at = 1:4)
    
    text(x=2.5, y=0.009, label=paste0("pvalue Syn/NonSyn = ",
                                      round(SampleSheet$pvalueSynNonSyn[mr], 3)))
    
    syn_nonCpG_array<-X$trasitionmutFreq[X$TypeOfSite=="syn" & X$makesCpG==0 & !is.na(X$a) & X$TotalReads>readscutoff ]
    syn_CpG_array<-X$trasitionmutFreq[X$TypeOfSite=="syn" & X$makesCpG==1 & !is.na(X$a) & X$TotalReads>readscutoff]
    nonsyn_nonCpG_array<-X$trasitionmutFreq[X$TypeOfSite!="syn" & X$makesCpG==0 & !is.na(X$a) & X$TotalReads>readscutoff]
    nonsyn_CpG_array<-X$trasitionmutFreq[X$TypeOfSite!="syn" & X$makesCpG==1 & !is.na(X$a) & X$TotalReads>readscutoff]

    t=1; arr= syn_nonCpG_array
    points(rep(t,length(arr))+rnorm(length(arr),mean=0,sd=0.1),arr, col=alpha(cols[t],0.5))
    points(t, median(arr), col=cols[t], pch=15, cex=2)
    points(t, mean(arr), col=cols[t], pch=16, cex=2)
    t=2; arr= syn_CpG_array
    points(rep(t,length(arr))+rnorm(length(arr),mean=0,sd=0.1),arr, col=alpha(cols[t],0.5))
    points(t, median(arr), col=cols[t], pch=15, cex=2)
    points(t, mean(arr), col=cols[t], pch=16, cex=2)
    t=3; arr= nonsyn_nonCpG_array
    points(rep(t,length(arr))+rnorm(length(arr),mean=0,sd=0.1),arr, col=alpha(cols[t],0.5))
    points(t, median(arr), col=cols[t], pch=15, cex=2)
    t=4; arr= nonsyn_CpG_array
    points(rep(t,length(arr))+rnorm(length(arr),mean=0,sd=0.1),arr, col=alpha(cols[t],0.5))
    points(t, median(arr), col=cols[t], pch=15, cex=2)
    #make plot
    if (length(syn_nonCpG_array)>0&length(syn_CpG_array)>0){
      text(x=1.5, y=0.008, label=paste0("2-sided pvalue = ",
                                      round(wilcox.test(syn_nonCpG_array, syn_CpG_array, alternative = "two")$p.value,3)))
      
      SampleSheet$pvalueSyn_CpG[mr]<-wilcox.test(syn_nonCpG_array, syn_CpG_array, alternative = "two")$p.value
    }
    if (length(nonsyn_nonCpG_array)>0&length(nonsyn_CpG_array)>0){
      text(x=3.5, y=0.008, label=paste0("2-sided pvalue = ",
                                        round(wilcox.test(nonsyn_nonCpG_array, nonsyn_CpG_array, alternative = "two")$p.value,3)))
      
      SampleSheet$pvalueNonsyn_CpG[mr]<-wilcox.test(nonsyn_nonCpG_array, nonsyn_CpG_array, alternative = "two")$p.value
    }
    }
    write.csv(x = X, file = paste0("ProcessedData/DataCpG/",sdf))
  }
}
dev.off()
}


#how many plasma significant
if (TRUE){
print(paste0("plasma: ", 
             length(which(SampleSheet$pvalueSyn_CpG[SampleSheet$Sample=="plasma"]<0.05)), 
             " out of ", 
             length(which(SampleSheet$Sample=="plasma"& !is.na(SampleSheet$pvalueSyn_CpG))), " significant"))
print(paste0("plasma, nonSyn: ", 
             length(which(SampleSheet$pvalueNonsyn_CpG[SampleSheet$Sample=="plasma"]<0.05)), 
             " out of ", 
             length(which(SampleSheet$Sample=="plasma"& !is.na(SampleSheet$pvalueNonsyn_CpG))), " significant"))

print(paste0("non-plasma Syn: ", 
             length(which(SampleSheet$pvalueSyn_CpG[SampleSheet$Sample!="plasma"]<0.05)), 
             " out of ", 
             length(which(SampleSheet$Sample!="plasma" & !is.na(SampleSheet$pvalueSyn_CpG))), " significant"))

print(paste0("non-plasma NonSyn: ", 
             length(which(SampleSheet$pvalueNonsyn_CpG[SampleSheet$Sample!="plasma"]<0.05)), 
             " out of ", 
             length(which(SampleSheet$Sample!="plasma" & !is.na(SampleSheet$pvalueNonsyn_CpG))), " significant"))
}

#Write Sample sheet with CpG Pvalues
if (TRUE){
write.csv(x=SampleSheet,file = "OriginalData/SampleSheetMac251AllSamples.csv", row.names = FALSE)
}
