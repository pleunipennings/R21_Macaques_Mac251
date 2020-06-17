
read.csv("Run_4_01_Kaho.csv")->K
read.csv("ProcessedData/DataCpG/Run_4_01_Animal_16314SeqData.csv")->P

#Only look at A and T sites (only those can make CpG sites)
P<-P[P$MajNt%in%c("a","t"),]

head(P)

library(dplyr)

P %>%                                        # Specify data frame
  group_by(makesCpG) %>%                         # Specify group indicator
  summarize(mean_size = mean(trasitionmutFreq, na.rm = TRUE))
P %>%                                        # Specify data frame
  group_by(makesCpG) %>%                         # Specify group indicator
  summarize(median = median(trasitionmutFreq, na.rm = TRUE))

P %>%                                        # Specify data frame
  group_by(TypeOfSite, makesCpG) %>%                         # Specify group indicator
  summarize(mean_size = mean(trasitionmutFreq, na.rm = TRUE))

P %>%                                        # Specify data frame
  group_by(TypeOfSite, makesCpG) %>%                         # Specify group indicator
  summarize(median = median(trasitionmutFreq, na.rm = TRUE))

P[!is.na(P$a)&P$TotalReads>readscutoff,]%>%
group_by(TypeOfSite, makesCpG) %>%
  summarize(mean_size = mean(trasitionmutFreq, na.rm = TRUE))

pdf("Run4_01_June12_B.pdf")
par(mfrow=c(2,2))
X<-P[P$trasitionmutFreq<0.5,]
X<-P
readscutoff=6000
if (TRUE){
plot(c(1,2,3,4), c(0,0,0,0), ylim=c(0,.01), col=0, xlim=c(0.5, 4.5), 
     #main = paste0("Monkey ",M," week ", SampleSheet$Week[mr], " comp ", SampleSheet$Sample[mr]), 
     main = "Run_4_01_Animal_16314SeqData.csv _ freq<0.5",
     xaxt="n", xlab="type of site")
mtext(text=c("syn non CpG", "syn CpG", "nonsyn non CpG", "nonsyn nonCpG"), side=1, at = 1:4)

#text(x=2.5, y=0.009, label=paste0("pvalue Syn/NonSyn = ",
 #                                 round(SampleSheet$pvalueSynNonSyn[mr], 3)))

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
  text(x=1.5, y=0.08, label=paste0("1-sided pvalue = ",
                                    round(wilcox.test(syn_nonCpG_array, syn_CpG_array, alternative = "l")$p.value,3)))
  
#  SampleSheet$pvalueSyn_CpG[mr]<-wilcox.test(syn_nonCpG_array, syn_CpG_array, alternative = "two")$p.value
}
if (length(nonsyn_nonCpG_array)>0&length(nonsyn_CpG_array)>0){
  text(x=3.5, y=0.08, label=paste0("1-sided pvalue = ",
                                    round(wilcox.test(nonsyn_nonCpG_array, nonsyn_CpG_array, alternative = "l")$p.value,3)))
  
#  SampleSheet$pvalueNonsyn_CpG[mr]<-wilcox.test(nonsyn_nonCpG_array, nonsyn_CpG_array, alternative = "two")$p.value
}
}
dev.off()
