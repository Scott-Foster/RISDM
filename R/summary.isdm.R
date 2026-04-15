
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
  biasTerms <- allFixedNames[grep( "PO_", allFixedNames)]
  allFixedNames <- setdiff( allFixedNames, biasTerms)

  #labels for the DC-only model
  DCTerms <- allFixedNames[grep( "DC_", allFixedNames)]
  allFixedNames <- setdiff( allFixedNames, DCTerms)

  #labels for the AA-only model
  AATerms <- allFixedNames[grep( "AA_", allFixedNames)]
  allFixedNames <- setdiff( allFixedNames, AATerms)
  
  #labels for the PA-only model
  PATerms <- allFixedNames[grep( "PA_", allFixedNames)]
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
#    alpha_id <- grep( "alpha", DCTerms)
#    ord <- order( DCTerms[alpha_id])
#    tmp <- c(DCTerms[-alpha_id], DCTerms[alpha_id[ord]])
    res$DC_ARTEFACT <- object$mod$summary.fixed[DCTerms,1:5]   
    logPisNames <- rownames( res$DC_ARTEFACT)[grep( "logDetectPi",rownames( res$DC_ARTEFACT))]
    if( length( logPisNames) > 0){
      piSum <- matrix( NA, nrow=length( logPisNames), ncol=5, dimnames= list(logPisNames, colnames(object$mod$summary.fixed)[1:5]))
      for( ii in logPisNames){
        #mean and variance of pi (not log pi)
        m1 <- INLA::inla.emarginal( function(xx) exp( xx), object$mod$marginals.fixed[[ii]])
        piSum[ii,"mean"] <- m1
        m2 <- INLA::inla.emarginal( function(xx) exp( 2*xx), object$mod$marginals.fixed[[ii]])
        piSum[ii,"sd"] <- sqrt( m2-m1^2)
        #quantiles, remember that quantiles are preserved under monotonic transformations.
        piSum[ii,c("0.025quant","0.5quant","0.975quant")] <- exp( INLA::inla.qmarginal( c(0.025, 0.5, 0.975), object$mod$marginals.fixed[[ii]]))
        res$DC_ARTEFACT[rownames( res$DC_ARTEFACT)==ii,] <- piSum[ii,]
        rownames( res$DC_ARTEFACT)[rownames( res$DC_ARTEFACT)==ii] <- sub( "logDetect", "Detect", ii)
      }
    }
    else{
      tmp <- matrix( NA, nrow=length( attr( object, "MoM_pi")), ncol=5)
      tmp[,1] <- attr( object, "MoM_pi")
      rownames( tmp) <- paste("MoM_pi",1:nrow( tmp))
      colnames( tmp) <- colnames( res$DC_ARTEFACT)
      res$DC_ARTEFACT <- rbind( res$DC_ARTEFACT, tmp)
    }
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
  #if a barrier model has been used, then transform spatial parameters
  if( !is.null( object$mesh$barrier$triBarrier.xy)){
    spatSum <- matrix( NA, nrow=2, ncol=5, dimnames= list(NULL, colnames(object$mod$summary.fixed)[1:5]))
    ronam <- names( object$mod$marginals.hyperpar)
    rownames( spatSum) <- gsub( "Theta2", "Stdev", gsub( "Theta1", "Range", ronam))
    for( ii in 1:2){
      #mean and variance in terms of raw moments
      m1 <- INLA::inla.emarginal( function(xx) exp( xx), object$mod$marginals.hyperpar[[ii]])
      spatSum[ii,"mean"] <- m1
      m2 <- INLA::inla.emarginal( function(xx) exp( 2*xx), object$mod$marginals.hyperpar[[ii]])
      spatSum[ii,"sd"] <- sqrt( m2-m1^2)
      #quantiles, remember that quantiles are preserved under monotonic transformations.
      spatSum[ii,c("0.025quant","0.5quant","0.975quant")] <- exp( INLA::inla.qmarginal( c(0.025, 0.5, 0.975), object$mod$marginals.hyperpar[[ii]]))
    }
    res$SPATIAL <- spatSum
  }

  class( res) <- "summary.isdm"
  return( res)
          
}


