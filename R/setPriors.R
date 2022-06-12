
###############################################################################################
###############################################################################################
####	
####	Set the priors on the combined stack based on what is in contarol.
####
####	Returns a list of priors
####
####	Programmed by Scott in the first half of 2022
####
###############################################################################################
###############################################################################################

setPriors <- function( control, stck){

  #container for the result
  tmp <- list( mean=control$prior.mean)
  tmp$prec <- list( default=control$other.prec)
  #for the intercept terms, in turn of data type.
  if( "Intercept.PO" %in% colnames( stck$effects$data))
    tmp$prec$Intercept.PO <- control$int.prec
  if( "Intercept.PA" %in% colnames( stck$effects$data))
    tmp$prec$Intercept.PA <- control$int.prec
  if( "Intercept.AA" %in% colnames( stck$effects$data))
    tmp$prec$Intercept.AA <- control$int.prec
  if( "Intercept.DC" %in% colnames( stck$effects$data))
    tmp$prec$Intercept.DC <- control$int.prec
  return( tmp)
}

