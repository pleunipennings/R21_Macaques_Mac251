library(ape)
library(seqinr)
library(pegas)
library(sfsmisc)
library(ggplot2)
library(scales)
library(plotrix)
library(RColorBrewer)

source("Rscripts/CPG_Function.R")
source("Rscripts/WTAA_consensus.R")
source("Rscripts/MUTAA.R")
source("Rscripts/SynNonSyn.R")
source("Rscripts/CPG_Function.R")

readscutoff=500
nucleotides =c("a", "c", "g", "t")

read.csv("OriginalData/SampleSheetMac251AllSamples.csv", stringsAsFactors = FALSE)->SampleSheet
SampleSheet$Week<-as.numeric(SampleSheet$Week)
SampleSheet$Tissue<-"Lung"
SampleSheet$Tissue[SampleSheet$Sample=="stockvirus"|SampleSheet$Sample=="control"]<-""
SampleSheet$Tissue[grep(pattern ="LN", SampleSheet$Sample)]<-"LN"
SampleSheet$Tissue[SampleSheet$Sample=="plasma"]<-"Plasma"
SampleSheet<-SampleSheet[order(SampleSheet$Week),]
SampleSheet<-SampleSheet[order(SampleSheet$Tissue),]

UniqueMonkeys<-unique(SampleSheet$Monkey)
UniqueMonkeys<-UniqueMonkeys[UniqueMonkeys!="stockvirus"&UniqueMonkeys!="control"]

r=which(SampleSheet$Monkey=="stockvirus")
sdf=SampleSheet$SeqDataFileName[r]
print(sdf)
SV<-read.csv(paste0("ProcessedData/SeqData/",sdf), row.names = 1)
dim(SV)

#* Transition function*
transition<-function(nuc){
    if (is.na(nuc) ==T ){ 
                return (NA)
                next}
    else if (nuc=="a") return("g")
    else if (nuc=="g") return("a")
    else if (nuc=="c") return("t")
    else if (nuc=="t") return("c")
    else if (nuc=='n') return('n')
}

transv1<-function(nuc){
        if (is.na(nuc) ==T ){ 
                return (NA)
                next}
        else if (nuc=="a") return("c")
        else if (nuc=="g") return("c")
        else if (nuc=="c") return("a")
        else if (nuc=="t") return("a")
        else if (nuc=='n') return('n')
}

transv2<-function(nuc){
        if (is.na(nuc) ==T ){ 
                return ("NA")
                next}
        else if (nuc=="a") return("t")
        else if (nuc=="g") return("t")
        else if (nuc=="c") return("g")
        else if (nuc=="t") return("g")
        else if (nuc=='n') return('n')
}






typeofsitefunction<-function(WTcodon, mutantcodon){
        WTAA<-seqinr::translate(WTcodon)
        MUTAA<-seqinr::translate(mutantcodon)
        if (WTAA == "X"){ 
                return ("NA")
                next}
        else if (WTAA == MUTAA) return ("syn")
        else if (MUTAA == "*") return ("stop")
        else return ("nonsyn")
}


EstimatedS <- function(mu, meanfreq){
        if (is.na(mu)|is.na(meanfreq) ==T){ 
                return ("NA")
                next}
         else if (meanfreq == 0) return (1)
         else return (min(c(mu/meanfreq,1)))
}

#Amino acid changes
pos <- "R|H|K"
neg <- "D|E"
unc <- "S|T|N|Q"
spe <- "C|U|G|P"
hyd <- "A|I|L|F|M|W|Y|V"
amCat <- function(AA){
    if(regexpr(pos, AA) > 0){ return(0) }
    if(regexpr(neg, AA) > 0){ return(1) }
    if(regexpr(unc, AA) > 0){ return(2) }
    if(regexpr(spe, AA) > 0){ return(3) }
    if(regexpr(hyd, AA) > 0){ return(4) }
    return(5)
}

#Colors
cols <- c("#66CCEE","#228833","#CCBB44","#EE6677","#AA3377","#4477AA","#BBBBBB")
cols<-c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695')
cols<-c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99')
# from https://personal.sron.nl/~pault/

