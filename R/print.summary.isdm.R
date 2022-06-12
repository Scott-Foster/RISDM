
###############################################################################################
###############################################################################################
####	
####	print from a summary.isdm object
####
####	Returns the object invisibly
####
####	Programmed by Scott in the first half of 2022
####
###############################################################################################
###############################################################################################

print.summary.isdm <- function( x, digits = max(3L, getOption("digits") - 3L), ...){
  
  #print the distribution terms
  if( !is.null( x$DISTRIBUTION)){
    cat("\nCoefficients for distribution model (all data)\n")
    print( round( as.matrix( x$DISTRIBUTION), digits))
  }
  #print the PO bias
  if( !is.null( x$PO_BIAS)){
    cat("\nCoefficients for sampling bias model (PO data)\n")
    print( round( as.matrix( x$PO_BIAS), digits))
  }
  #print the DC artefact object
  if( !is.null( x$DC_ARTEFACT)){
    cat("\nCoefficients for artefect model for DC data\n")
    print( round( as.matrix( x$DC_ARTEFACT), digits))
  }
  #print( the AA artefact
  if( !is.null( x$AA_ARTEFACT)){
    cat("\nCoefficients for artefact model for AA data\n")
    print( round( as.matrix( x$AA_ARTEFACT), digits))
  }
  #print the PA artecfact
  if( !is.null( x$PA_ARTEFACT)){
    cat("\nCoefficients for artefact model for PA data\n")
    print( round( as.matrix( x$PA_ARTEFACT), digits))
  }
  #print the spatial terms, if present
  if( !is.null( x$SPATIAL)){
    cat("\nHyperparameters for Spatial Effects\n")
    print( round( as.matrix( x$SPATIAL), digits))
  }
  #print the logl.
  cat("\nMarginal log-likelihood\n")
  cat( round( x$marg.lik, digits))
  
  cat("\n\n")
  invisible(x)
}
