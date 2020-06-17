###Get Wildtype amino acid###

getWTAA<-function(df){
  #Assign consensus to a variable
  df$MajNt<-as.character(df$MajNt)
  df$MajNt[is.na(df$MajNt)]<-"n"
  cons =  as.character(df$MajNt)

  #Create "WTAA" column if not already
  if (length(which(names(df)=="WTAA"))==0){
    df$WTAA =0}
  
  #Loop for translating consensus
  frame = 3
  for(x in seq(frame, length(cons) - 2, 3)){ #PSP June 2020 start translating at pos 3
    codon <- c(cons[x], cons[x+1], cons[x+2])
    #print(paste0(x, codon))
    new_AA <- seqinr::translate(codon)
    df$WTAA [x] <- new_AA
    df$WTAA [x+1] <- new_AA
    df$WTAA [x+2] <- new_AA
  }

  #Insert value into column
  #df$WTAA_consensus <-WTAA_consensus 

  #Return the data frame
  return(df)
}
