
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
    if( ! any( names( control$prior.list) %in% c("prec","mean")))
      stop( "Priors for fixed effects contain elements that are not 'mean' nor 'prec'. Please check (see ?isdm too)")
    tmpNames <- names( control$prior.list$mean)
    tmpNames <- setdiff( tmpNames, "default")
    if( !is.null(control$prior.list$mean) & !all( tmpNames %in% stck$effects$names))
      message( "Priors set for the mean of at least one fixed effect that is not in the formula.  Please check.")
    tmpNames <- names( control$prior.list$prec)
    tmpNames <- setdiff( tmpNames, "default")
    if( !is.null(control$prior.list$prec) & !all( tmpNames %in% stck$effects$names))
      message( "Priors set for the precision of at least one fixed effect that is not in the formula. Please check.")
#    return( control$prior.list)
    tmp <- control$prior.list
    if( !"default" %in% names( tmp$mean))
      tmp$mean$default <- control$prior.mean
    if( !"default" %in% names( tmp$prec))
      tmp$prec$default <- control$other.prec #intercepts will be updated separately and presently
  }
  else
    tmp <- list( mean=list( default=control$prior.mean), prec=list( default=control$other.prec)) #intercepts will be updated separately and presently

  if( ("PO_Intercept" %in% stck$effects$names) & (!"PO_Intercept" %in% names( tmp$prec)))
    tmp$prec$PO_Intercept <- control$int.prec
  if( ("PA_Intercept" %in% stck$effects$names) & (!"PA_Intercept" %in% names( tmp$prec)))
    tmp$prec$PA_Intercept <- control$int.prec
  if( ("AA_Intercept" %in% stck$effects$names) & (!"AA_Intercept" %in% names( tmp$prec)))
    tmp$prec$AA_Intercept <- control$int.prec
  if( ("DC_Intercept" %in% stck$effects$names) & (!"DC_Intercept" %in% names( tmp$prec)))
    tmp$prec$DC_Intercept <- control$int.prec
    
  return( tmp)
}

