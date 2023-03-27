
###############################################################################################
###############################################################################################
####	
####	Combine formulas together to get a formula for each data type. Note that not all
####	formulas need to be supplied.
####
####	Returns a formula for passing to INLA
####
####	Programmed by Scott in March 2023
####
###############################################################################################
###############################################################################################

makePartialFormulas <- function( dform, bform, aforms, addRE=TRUE) {
  #not going to use a combined intercept.  Rather let each data type have its own
  
  #try/do add the bias specific terms to each of the different data types.  Do so by including a nested effect within intercept types
  formmy <- list()
  for( ii in c("PO","PA","AA","DC")){
    tmpForm <- dform
    if( ii == "PO"){
      if( !is.null( bform)){
	#belt and braces to remove the PO intercept. PO intercept will get added later (under a differnt name)
        tmpForm <- update( bform, "~.-0-1")  
      }
    }
    #add the terms from each component formula
    if( ii != "PO"){
      addTerms <- as.character( tmpForm)
      addTerms <- addTerms[length(addTerms)]
      #the intercept under a component label.
      tmpForm <- update( tmpForm, paste0("~.+Intercept.",ii,"/(",addTerms,")"))
    }
    formmy[[ii]] <- tmpForm
  }

  #add the random effects, if asked to.  Note that the label is already specified internally.
  if( addRE){
    #index is always assumed to be "isdm.spat.XXX" and spde model is always "my.spde"
    formmy <- lapply( formmy, function(xx) update( xx, "~ . + f(isdm.spat.XXX,model=my.spde)")) 
  }
  #make a common outcome name
  formmy <- lapply( formmy, function(xx) update( xx, "resp~."))
  
  return( formmy)
}

