#function for syn/non/nonsence
#This function intakes a dataframe with the columns MUTAA, WTAA , and TypeOfSite created. Then evaluates the values in the WTAA and MUTAA columns to determine the value for the TypeOfSite column
synFunction <- function(df) {
  for (h in 1:nrow(df)){ #looks at each row in the dataframe df
    if(df[h,"MUTAA"]== df[h,"WTAA"]){ #if the value in the MUTAA column for the row of interest is equal to the value in the WTAA of the same row
      df[h,"TypeOfSite"] = "syn" #then insert "syn" for the value in the TypeOfSite column for this row
    }
    if(df[h,"MUTAA"] != df[h,"WTAA"]){ #if the value in the MUTAA column for the row of interest is NOT equal to the value in the WTAA of the same row
      if(df[h,"MUTAA"]=="*"){ #and it the MUTAA value equal "*"
        df[h,"TypeOfSite"] = "nonsense"} #then insert "nonsense" for the value in the TypeOfSite column for this row
      else {
        df[h,"TypeOfSite"] = "nonsyn" #else insert "nonsyn" for the value in the TypeOfSite column for this row
      }
    }
    if (df[h,"WTAA"]==0 | df[h,"WTAA"]=="X") df[h,"TypeOfSite"] = NA
  }
return(df)
}


synFunction_stock <- function(df) {
  df$Typea=0; df$Typec=0;df$Typeg=0; df$Typet=0
  for (h in 1:nrow(df)){ #looks at each row in the dataframe df
    if(df[h,"MUTAAa"]== df[h,"WTAA_stock"]){ #if the value in the MUTAA column for the row of interest is equal to the value in the WTAA of the same row
      df[h,"Typea"] = "syn" #then insert "syn" for the value in the TypeOfSite column for this row
    }
    if(df[h,"MUTAAa"] != df[h,"WTAA_stock"]){ #if the value in the MUTAA column for the row of interest is NOT equal to the value in the WTAA of the same row
      if(df[h,"MUTAAa"]=="*"){ #and it the MUTAA value equal "*"
        df[h,"Typea"] = "nonsense"} #then insert "nonsense" for the value in the TypeOfSite column for this row
      else {
        df[h,"Typea"] = "nonsyn" #else insert "nonsyn" for the value in the TypeOfSite column for this row
      }
    }
    if (df[h,"WTAA_stock"]==0 | df[h,"WTAA_stock"]=="X") df[h,"Typea"] = NA
  } #a
  for (h in 1:nrow(df)){ #looks at each row in the dataframe df
    if(df[h,"MUTAAc"]== df[h,"WTAA_stock"]){ #if the value in the MUTAA column for the row of interest is equal to the value in the WTAA of the same row
      df[h,"Typec"] = "syn" #then insert "syn" for the value in the TypeOfSite column for this row
    }
    if(df[h,"MUTAAc"] != df[h,"WTAA_stock"]){ #if the value in the MUTAA column for the row of interest is NOT equal to the value in the WTAA of the same row
      if(df[h,"MUTAAc"]=="*"){ #and it the MUTAA value equal "*"
        df[h,"Typec"] = "nonsense"} #then insert "nonsense" for the value in the TypeOfSite column for this row
      else {
        df[h,"Typec"] = "nonsyn" #else insert "nonsyn" for the value in the TypeOfSite column for this row
      }
    }
    if (df[h,"WTAA_stock"]==0 | df[h,"WTAA_stock"]=="X") df[h,"Typec"] = NA
  } #c
  for (h in 1:nrow(df)){ #looks at each row in the dataframe df
    if(df[h,"MUTAAg"]== df[h,"WTAA_stock"]){ #if the value in the MUTAA column for the row of interest is equal to the value in the WTAA of the same row
      df[h,"Typeg"] = "syn" #then insert "syn" for the value in the TypeOfSite column for this row
    }
    if(df[h,"MUTAAg"] != df[h,"WTAA_stock"]){ #if the value in the MUTAA column for the row of interest is NOT equal to the value in the WTAA of the same row
      if(df[h,"MUTAAg"]=="*"){ #and it the MUTAA value equal "*"
        df[h,"Typeg"] = "nonsense"} #then insert "nonsense" for the value in the TypeOfSite column for this row
      else {
        df[h,"Typeg"] = "nonsyn" #else insert "nonsyn" for the value in the TypeOfSite column for this row
      }
    }
    if (df[h,"WTAA_stock"]==0 | df[h,"WTAA_stock"]=="X") df[h,"Typeg"] = NA
  } #g
  for (h in 1:nrow(df)){ #looks at each row in the dataframe df
    if(df[h,"MUTAAt"]== df[h,"WTAA_stock"]){ #if the value in the MUTAA column for the row of interest is equal to the value in the WTAA of the same row
      df[h,"Typet"] = "syn" #then insert "syn" for the value in the TypeOfSite column for this row
    }
    if(df[h,"MUTAAt"] != df[h,"WTAA_stock"]){ #if the value in the MUTAA column for the row of interest is NOT equal to the value in the WTAA of the same row
      if(df[h,"MUTAAt"]=="*"){ #and it the MUTAA value equal "*"
        df[h,"Typet"] = "nonsense"} #then insert "nonsense" for the value in the TypeOfSite column for this row
      else {
        df[h,"Typet"] = "nonsyn" #else insert "nonsyn" for the value in the TypeOfSite column for this row
      }
    }
    if (df[h,"WTAA_stock"]==0 | df[h,"WTAA_stock"]=="X") df[h,"Typet"] = NA
  }#t
    return(df)
}

