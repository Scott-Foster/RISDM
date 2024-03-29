
###############################################################################################
###############################################################################################
####	
####	Combine formulas together to get a single formula for use in INLA. Note that not all
####	formulas need to be supplied.
####
####	Returns a formula for passing to INLA
####
####	Based on old makeFormula, but utlises new design matrix structure etc...
####	This version programmed in March 2023
####
###############################################################################################
###############################################################################################

combineFormulae <- function( forms, addRE=FALSE) {

  fullForm <- ~-1+removeMe
  for( jj in 1:length( forms)) {
    if( inherits( forms[[jj]], "list")) {
      for( ii in 1:length( forms[[jj]])) {
	addTerms <- as.character( forms[[jj]][[ii]])
	addTerms <- addTerms[length( addTerms)]
	fullForm <- stats::update.formula( fullForm, paste0("~.+",addTerms))
      }
    }
    else {
      addTerms <- as.character( forms[[jj]])
      addTerms <- addTerms[length( addTerms)]
      if( length( addTerms) > 0 )
	fullForm <- stats::update.formula( fullForm, paste0("~.+",addTerms))
    }
  }
  fullForm <- stats::update.formula( fullForm, "~.-removeMe")
    
  #add the random effects, if asked to.  Note that the label is already specified internally.
  if( addRE){
    #index is always assumed to be "isdm.spat.XXX" and spde model is always "my.spde"
    fullForm <- stats::update.formula( fullForm, "~ . + f(isdm.spat.XXX,model=my.spde)") 
  }
  #just in case the formula still has a (user-specified) outcome.
  fullForm <- stats::update.formula( fullForm, "resp~.")
  
  return( fullForm)
}

