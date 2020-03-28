
SeqData_stock5<-read.csv("ProcessedData/SeqData/Run_5_01_Animal_stockvirusSeqData.csv", row.names = 1)
SeqData_stock6<-read.csv("ProcessedData/SeqData/Run_6_01_Animal_stockvirusSeqData.csv", row.names = 1)



SeqData_stock5$MaxMutFreq<-0
columns=which(names(SeqData_stock5)%in%c("a","c","g","t"))
nts<-c("a","c","g","t")

for (i in 1:nrow(SeqData_stock5)){
#i=1
#which column is the most common mut?
mutcolumns<-columns[which(names(SeqData_stock5)[columns]!=SeqData_stock5$MajNt[i])]
maxmutcol<-which(nts==names(which.max(SeqData_stock5[i,mutcolumns])))
SeqData_stock5$MaxMut[i]<-names(SeqData_stock5)[maxmutcol]
SeqData_stock5$MaxMutFreq[i]<-SeqData_stock5[i,maxmutcol]/SeqData_stock5$TotalReads[i]
}


SeqData_stock6$MutFreqMaxMutStock5 <-0
for (i in 1:nrow(SeqData_stock6)){
maxmutcol<-which(nts==SeqData_stock5$MaxMut[i])
SeqData_stock6$MutFreqMaxMutStock5[i]<-SeqData_stock6[i,maxmutcol]/SeqData_stock6$TotalReads[i]
}

pdf (paste0("Output/StockComparison", Sys.Date(), ".pdf"))
plot(SeqData_stock5$MaxMutFreq,SeqData_stock6$MutFreqMaxMutStock, log="xy", col=alpha(cols,0.5), pch=16,
          main = "Comparison between two stock samples, R21", xlab = "Stock 5 Mutation frequencies",ylab = "Stock 6 Mutation frequencies"
     , xlim = c(0.00001,0.05), ylim = c(0.00001,0.05) )
abline(a = 0, b = 1)
cor.test(SeqData_stock5$MaxMutFreq,SeqData_stock6$MutFreqMaxMutStock, method="spearman")
corcoef= round(cor(SeqData_stock5$MaxMutFreq,SeqData_stock6$MutFreqMaxMutStock, method="spearman"),3)
pvalue=cor.test(SeqData_stock5$MaxMutFreq,SeqData_stock6$MutFreqMaxMutStock, method="spearman")$p.value
text(x = 10^-4,y = 10^-2,labels = paste0("Spearman Corr = ", corcoef))
text(x = 10^-4,y = 0.5*10^-2,labels = paste0("P-value = ", round(pvalue,80)))
dev.off()

pdf (paste0("Output/StockComparison_linearscale", Sys.Date(), ".pdf"))
plot(SeqData_stock5$MaxMutFreq,SeqData_stock6$MutFreqMaxMutStock, col=alpha(cols,0.7), pch=16,
     main = "Comparison between two stock samples, R21", xlab = "Stock 5 Mutation frequencies",ylab = "Stock 6 Mutation frequencies"
     ,xlim = c(0.00001,0.05), ylim = c(0.00001,0.05) )
abline(a = 0, b = 1)
cor.test(SeqData_stock5$MaxMutFreq,SeqData_stock6$MutFreqMaxMutStock, method="spearman")
corcoef= round(cor(SeqData_stock5$MaxMutFreq,SeqData_stock6$MutFreqMaxMutStock, method="spearman"),3)
pvalue=cor.test(SeqData_stock5$MaxMutFreq,SeqData_stock6$MutFreqMaxMutStock, method="spearman")$p.value
text(x = 10^-4,y = 10^-2,labels = paste0("Spearman Corr = ", corcoef))
text(x = 10^-4,y = 0.5*10^-2,labels = paste0("P-value = ", round(pvalue,80)))
dev.off()


