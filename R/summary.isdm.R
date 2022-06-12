
###############################################################################################
###############################################################################################
####	
####	Producing (hopefully) easier to read printed summary of the isdm
####
####	Returns a summary.isdm object
####
####	Programmed by Scott in the first half of 2022
####
###############################################################################################
###############################################################################################

summary.isdm <- function( object, ...){
  
  #find names
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

  #find labels for the distribution terms
  distTerms <- allFixedNames #what is left *must*? be distribution...
  rm( allFixedNames)

  #bundling together for return
  res <- list()
  #distribution terms
  res$DISTRIBUTION <- if( length( distTerms)!=0) object$mod$summary.fixed[distTerms,1:5] else NULL
  #PO bias terms
  res$PO_BIAS <- if( length( biasTerms)!=0) object$mod$summary.fixed[biasTerms,1:5] else NULL
  #DC terms
  if( length( DCTerms)!=0) {
    alpha_id <- grep( "alpha", DCTerms)
    ord <- order( DCTerms[alpha_id])
    tmp <- c(DCTerms[-alpha_id], DCTerms[alpha_id[ord]])
    res$DC_ARTEFACT <- object$mod$summary.fixed[DCTerms,1:5]      
  }
  else
    res$DC_ARTEFACT <- NULL
  #AA artefact terms
  res$AA_ARTEFACT <- if( length( AATerms)!=0) object$mod$summary.fixed[AATerms,1:5]  else NULL
  #PA artefact terms
  res$PA_ARTEFACT <- if( length( PATerms)!=0) object$mod$summary.fixed[PATerms,1:5] else NULL
  #spatial terms
  res$SPATIAL <- if( length( object$mod$summary.hyperpar) != 0) object$mod$summary.hyperpar[,1:5] else NULL
  #marginal logl
  res$marg.lik <- object$mod$mlik[2]

  class( res) <- "summary.isdm"
  return( res)
          
}


