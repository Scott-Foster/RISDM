
###############################################################################################
###############################################################################################
####	
####	Combine formulas together to get a single formula for use in INLA. Note that not all
####	formulas need to be supplied.
####
####	Returns a formula for passing to INLA
####
####	Programmed by Scott in the first half of 2022
####
###############################################################################################
###############################################################################################

makeFormula <- function( dform, bform, aforms, addRE=TRUE) {
  #not going to use a combined intercept.  Rather let each data type have its own
  
  #try/do add the bias specific terms to each of the different data types.  Do so by including a nested effect within intercept types
  formmy  <- dform
  for( ii in c("PO","PA","AA","DC")){
    tmpForm <- NULL
    if( ii == "PO"){
      if( !is.null( bform)){
	#belt and braces to remove the PO intercept. PO intercept will get added later (under a differnt name)
        tmpForm <- update( bform, "~.-0-1")  
      }
    }
    else{
      if( !is.null( aforms[[ii]])){
	#belt and braces to remove intercept.  Will get added later (under a different name)
        tmpForm <- update( aforms[[ii]], "~.-0-1")  
      }
    }
    #add the terms from each component formula
    if( !is.null( tmpForm)){
      addTerms <- as.character( tmpForm)
      addTerms <- addTerms[length(addTerms)]
      #the intercept under a component label.
      formmy <- update( formmy, paste0("~.+Intercept.",ii,"/(",addTerms,")"))
    }
  }

  #add the random effects, if asked to.  Note that the label is already specified internally.
  if( addRE){
    #index is always assumed to be "isdm.spat.XXX" and spde model is always "my.spde"
    formmy <- update( formmy, "~ . + f(isdm.spat.XXX,model=my.spde)") 
  }
  #just in case the formula still has a (user-specified) outcome.
  if( length( formmy)==3)
    warning( "Ignoring distributionFormula's outcome and replacing with 'resp' (data generated internally)" )
  formmy <- update( formmy, "resp~.")
  
  return( formmy)
}

