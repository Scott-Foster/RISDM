
###############################################################################################
###############################################################################################
####	Get the set of variables used in the model, from any artefact sub-component.
####	
####	Returns a list with unique variable names in each sub-component.
####
####	Programmed by Scott in the first half of 2022 
####
###############################################################################################
###############################################################################################

uniqueVarNames <- function( obsList, arteForm, DCsurvID){

  #container for the altered artefact formuals
  newForm <- arteForm
  #container for the altered observation lists (names to match newForm)
  newObs <- obsList
  #cycle through each of the formulas (artefact)
  for( ii in names( newForm)){
    #terms used in the formula
    tt <- terms( newForm[[ii]])
    #the variable labels
    tmp.labels <- attr( tt, "term.labels")
    if( length( tmp.labels) > 0){
      #append the data type to label
      tmp <- paste( tmp.labels, ii, sep="_")
      dataname <- paste0(ii,"dat")
      #change data variable names
      colnames( newObs[[dataname]])[colnames( newObs[[dataname]]) %in% attr( tt, "term.labels")] <- tmp
      #change it back to a formula
      newForm[[ii]] <- stats::reformulate( tmp)
      #fix up environments
      environment( newForm[[ii]]) <- environment( arteForm[[ii]])
      environment( newObs[[dataname]]) <- environment( obsList[[dataname]])

      #a bit extra for DC data -- just to make it extra unique :-)
      if( ii == "DC")
        if( DCsurvID %in% attr( tt, "term.labels"))
          DCsurvID <- paste( DCsurvID,"DC",sep="_")
    }
  }
  
  #the return object
  res <- list( obsList=newObs, arteForm=newForm, DCsurvID=DCsurvID)
}

