
###############################################################################################
###############################################################################################
####	
####	Set the priors on the combined stack based on what is in control.
####
####	Returns a list of priors
####
####	Programmed by Scott in the first half of 2022
####
###############################################################################################
###############################################################################################

setPriors <- function( control, stck){

  if( !is.null( control$prior.list))
    return( control$prior.list)

  #container for the result
  tmp <- list( mean=control$prior.mean)
  tmp$prec <- list( default=control$other.prec)
  #for the intercept terms, in turn of data type.
  if( "PO_Intercept" %in% colnames( stck$effects$data))
    tmp$prec$Intercept.PO <- control$int.prec
  if( "PA_Intercept" %in% colnames( stck$effects$data))
    tmp$prec$Intercept.PA <- control$int.prec
  if( "AA_Intercept" %in% colnames( stck$effects$data))
    tmp$prec$Intercept.AA <- control$int.prec
  if( "DC_Intercept" %in% colnames( stck$effects$data))
    tmp$prec$Intercept.DC <- control$int.prec
  return( tmp)
}

