library(tidyverse)
source("Rscripts/BaseRscript2.R")

###########################
#
SIVFiles<-list.files("ProcessedData/CSV/",pattern="csv")
stockrows=grep(pattern = "stock", x = SIVFiles)

coding.start<-95
coding.end<-1000
no_df<-data.frame("pos"=c(coding.start:coding.end))
for (i in stockrows){
        print(i)
        SeqData<-read.csv(paste("ProcessedData/CSV/",SIVFiles[i],sep=""))
        SeqData<-SeqData[,-c(1,2,9)]
        colnames(SeqData)[1]<-"pos"
        colnames(SeqData)[7:8]<-c("deletion","insertion")
        #remove the reads<1000 
        SeqData<-SeqData[-which(SeqData$TotalReads<10000),] #PSP March 23 2020 now removing <1000 instead of <100 #PSP I don't think this is needed FEB 2020 now I do!
        #determine the majority nucleotide base at each site
        colnames(SeqData)[2:5]<-c("a","c","g","t")
        SeqData$MajNt<-apply(SeqData[,2:5],1,function(x) c("a","c","g","t")[which.max(x)])
        #read the refrence sequence:
        SeqData<-merge(no_df,SeqData,by="pos",all.x=T)
        reference<-read.dna("OriginalData/SIV_ENV.fasta", format = "fasta",as.character=TRUE)
        ref.code<-reference[coding.start:coding.end]
        SeqData$ref<-ref.code[1:length(SeqData[,1])]
        head(SeqData,100)
        SeqData$transition.maj<-NA
        SeqData$transition.ref<-NA
        for (j in 1:nrow(SeqData)) SeqData$transition.maj[j]<-transition(SeqData$MajNt[j])
        for (j in 1:nrow(SeqData)) SeqData$transition.ref[j]<-transition(SeqData$ref[j])
        #rearrange the columns
        SeqData<-SeqData[,c("a","c","g","t","deletion","insertion","N","pos","TotalReads","MajNt","ref","transition.maj","transition.ref")]
        SeqDataFilename<-gsub(pattern=".csv", replace="SeqData.csv",x=SIVFiles[i])
        write.csv(SeqData, paste0("ProcessedData/SeqData/",SeqDataFilename))
}
     
#Now that we have the stock viruses, we can see what positions are the same and then get only those positions for all the other samples

SeqData_stock5<-read.csv("ProcessedData/SeqData/Run_5_01_Animal_stockvirusSeqData.csv")
SeqData_stock6<-read.csv("ProcessedData/SeqData/Run_6_01_Animal_stockvirusSeqData.csv")

#My first question is: Is the majority the same for these two stocks? 

length(which(SeqData_stock5$MajNt != SeqData_stock6$MajNt)) # 5 positions are different

length(which(SeqData_stock5$MajNt == SeqData_stock6$MajNt)) #475 are the same #March 25 only 449 left because of more stringent filtering

SeqData_stock5[which(SeqData_stock5$MajNt != SeqData_stock6$MajNt),]
SeqData_stock6[which(SeqData_stock5$MajNt != SeqData_stock6$MajNt),]

#I think I should remove these sites from analysis. 
#Focus analysis on positions that are the same in the two Stock virus samples. 
#Writing the positions to a csv file. 

write.csv(x = SeqData_stock5$pos[which(SeqData_stock5$MajNt == SeqData_stock6$MajNt)], 
          file = "ProcessedData/PositionsSameStock.csv", row.names = FALSE)

write.csv(x = SeqData_stock5$ref[which(SeqData_stock5$MajNt == SeqData_stock6$MajNt)], 
          file = "ProcessedData/RefPositionsSameStock.csv", row.names = FALSE)

#Pleuni May 1st 2020
#Alternative to not filter out the 5 sites that are different between the stock samples: 
#write.csv(x = SeqData_stock5$pos, file = "ProcessedData/PositionsSameStock.csv", row.names = FALSE)
#write.csv(x = SeqData_stock5$ref, file = "ProcessedData/RefPositionsSameStock.csv", row.names = FALSE)

