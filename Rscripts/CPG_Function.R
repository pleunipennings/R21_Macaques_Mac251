######## Investigating makesCpG sites: 
library(stringr)
## Note : this function assumes that wtnt column is present in the 
##        dataframe and your dataframe that is set as df

# Runs makesCpG_site function which takes in the arguement/input dataframe
CPG_site<-function(df){
  
  # Collapses wtnt_consensus vector sequence into a character string and sets to variable STRING
  paste(df$MajNt, collapse = '') -> STRING

  # create new column name makesCpG with values zero
  df$makesCpG<- 0
  #df$ca<- 0
  #df$tg<- 0

  #Looks for pattern tg in the data STRING, return location of TG sites within the STRING and store in variable TG
  gregexpr(pattern ='tg',STRING ) -> TG
  #df$tg<-sum(str_count(STRING, "tg"))
  
  # Since TG is a list, insert 'list' TG into a data frame BELL
  BELL <- data.frame(matrix(unlist(TG)))
 
  # Inserting value 1 in the column "makesCpG" and using the values found dataframe BELL as the row #
  df[BELL[,1],"makesCpG"] <- 1

  # For CA sites: 
  
  #Looks for pattern ca in the data STRING, return location of CA sites within the STRING and store in variable CA
  gregexpr(pattern ='ca',STRING ) -> CA
  #df$ca<-sum(str_count(STRING, "ca"))
  # Since TG is a list, insert 'list' TG into a data frame BELL
  CASITES <- data.frame(matrix(unlist(CA)))
  
  # Inserting value 1 in the column "makesCpG" and using the values found dataframe BELL as the row #
  df[CASITES[,1]+1,"makesCpG"] <- 1
  #returns new dataframe
  return(df)
}

# Calls the makesCpG_site function with the dataframe as the argument Old dataframe becomes the new dataframe


