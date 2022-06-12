
summary.isdm <- function( object, ...){
  
  allFixedNames <- object$mod$names.fixed
  
  #labels for the PO-only model
  biasTerms <- allFixedNames[grep( "Intercept.PO", allFixedNames)]
  allFixedNames <- setdiff( allFixedNames, biasTerms)

  #labels for the DC-only model
  DCTerms <- allFixedNames[grep( "Intercept.DC", allFixedNames)]
  allFixedNames <- setdiff( allFixedNames, DCTerms)

  #labels for the AA-only model
  AATerms <- allFixedNames[grep( "Intercept.AA", allFixedNames)]
  allFixedNames <- setdiff( allFixedNames, AATerms)
  
  #labels for the PA-only model
  PATerms <- allFixedNames[grep( "Intercept.PA", allFixedNames)]
  allFixedNames <- setdiff( allFixedNames, PATerms)

  distTerms <- allFixedNames #what is left *must*? be distribution...
  rm( allFixedNames)

  res <- list()
  res$DISTRIBUTION <- if( length( distTerms)!=0) object$mod$summary.fixed[distTerms,1:5] else NULL
  res$PO_BIAS <- if( length( biasTerms)!=0) object$mod$summary.fixed[biasTerms,1:5] else NULL
  if( length( DCTerms)!=0) {
    alpha_id <- grep( "alpha", DCTerms)
    ord <- order( DCTerms[alpha_id])
    tmp <- c(DCTerms[-alpha_id], DCTerms[alpha_id[ord]])
    res$DC_ARTEFACT <- object$mod$summary.fixed[DCTerms,1:5]      
  }
  else
    res$DC_ARTEFACT <- NULL
  res$AA_ARTEFACT <- if( length( AATerms)!=0) object$mod$summary.fixed[AATerms,1:5]  else NULL
  res$PA_ARTEFACT <- if( length( PATerms)!=0) object$mod$summary.fixed[PATerms,1:5] else NULL

  res$SPATIAL <- if( length( object$mod$summary.hyperpar) != 0) object$mod$summary.hyperpar[,1:5] else NULL
  
  res$marg.lik <- object$mod$mlik[2]

  class( res) <- "summary.isdm"
  return( res)
          
}


print.summary.isdm <- function( x, digits = max(3L, getOption("digits") - 3L), ...){
  if( !is.null( x$DISTRIBUTION)){
    cat("\nCoefficients for distribution model (all data)\n")
    print( round( as.matrix( x$DISTRIBUTION), digits))
  }
  if( !is.null( x$PO_BIAS)){
    cat("\nCoefficients for sampling bias model (PO data)\n")
    print( round( as.matrix( x$PO_BIAS), digits))
  }
  if( !is.null( x$DC_ARTEFACT)){
    cat("\nCoefficients for artefect model for DC data\n")
    print( round( as.matrix( x$DC_ARTEFACT), digits))
  }
  if( !is.null( x$AA_ARTEFACT)){
    cat("\nCoefficients for artefact model for AA data\n")
    print( round( as.matrix( x$AA_ARTEFACT), digits))
  }
  if( !is.null( x$PA_ARTEFACT)){
    cat("\nCoefficients for artefact model for PA data\n")
    print( round( as.matrix( x$PA_ARTEFACT), digits))
  }
  
  if( !is.null( x$SPATIAL)){
    cat("\nHyperparameters for Spatial Effects\n")
    print( round( as.matrix( x$SPATIAL), digits))
  }

  cat("\nMarginal log-likelihood\n")
  cat( round( x$marg.lik, digits))
  
  cat("\n\n")
  invisible(x)
}
