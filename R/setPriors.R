
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

  if( !is.null( control$prior.list)){
    if( ! names( control$prior.list) %in% c("prec","mean"))
      stop( "Priors for fixed effects contain elements that are not 'mean' nor 'prec'. Please check (see ?isdm too)")
    if( !is.null(control$prior.list$mean) & !names( control$prior.list$mean) %in% stck$effects$names)
      message( "Priors set for the mean of at least one fixed effect that is not in the formula.  Please check.")
    if( !is.null(control$prior.list$prec) & !names( control$prior.list$prec) %in% stck$effects$names)
      message( "Priors set for the precision of at least one fixed effect that is not in the formula. Please check.")
    return( control$prior.list)
  }

  #container for the result
  tmp <- list( mean=control$prior.mean)
  tmp$prec <- list( default=control$other.prec)
  #for the intercept terms, in turn of data type.
  if( "PO_Intercept" %in% colnames( stck$effects$data))
    tmp$prec$PO_Intercept <- control$int.prec
  if( "PA_Intercept" %in% colnames( stck$effects$data))
    tmp$prec$PA_Intercept <- control$int.prec
  if( "AA_Intercept" %in% colnames( stck$effects$data))
    tmp$prec$AA_Intercept <- control$int.prec
  if( "DC_Intercept" %in% colnames( stck$effects$data))
    tmp$prec$DC_Intercept <- control$int.prec
  return( tmp)
}

