
###############################################################################################
###############################################################################################
####	Get the variables needed for each of the formulas.
####	
####	Returns a charater vector of the variable names
####
####	Programmed by Scott in the first half of 2022
####
###############################################################################################
###############################################################################################

getVarNames <- function( fo, bfo, foList){

  #object for the variable names
  varNames <- NULL
  #for the distribution formula
  if( !is.null( fo)){
    #just add the new unique ones
    varNames <- setdiff( as.character( attr( terms( fo), "variables")), "list")
  }
  #for the bias formula (PO data)
  if( !is.null( bfo)){
    #just add the new ones
    varNames <- c( varNames, setdiff( as.character( attr( terms( bfo), "variables")), "list"))
  }
  #for the artefact formulas
  for( ii in names( foList))
    if( !is.null( foList[[ii]])){
      #just add the new ones
      varNames <- c( varNames, setdiff( as.character( attr( terms( foList[[ii]]), "variables")), "list"))
    }

  #make sure that they are unique between formulas too
  varNames <- unique( varNames)
  return( varNames)
}

