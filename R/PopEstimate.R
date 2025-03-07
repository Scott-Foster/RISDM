
###############################################################################################
###############################################################################################
####	
####	Predict a population estimate based on an isdm
####
####	Returns a list of predictions
####
####	Programmed by Scott in the first half of 2022
####
###############################################################################################
###############################################################################################

PopEstimate <- function( preds, intercept.terms=NULL, control=NULL){

  control <- makePopEstiControl( control)

  #checking terms
  if( !all( intercept.terms %in% rownames( preds$fixedSamples)))
    stop( "intercept.terms specified are not part of the model. Please check call and model.")

  #assuming, if needed, or bundling if we have to
  if( is.null( intercept.terms)){
    message( "Assuming that Intercept terms have already been included in predictions.  Is this what you want? See ?predict.isdm for how to include them.")
  }
  else{
    int.contr <- preds$fixedSamples[rownames( preds$fixedSamples) %in% intercept.terms,,drop=FALSE]
    int.contr <- colSums( int.contr)
    int.contr <- exp( int.contr)
    tmp <- sweep( x=preds$cell.samples, MARGIN=2, STATS=int.contr, FUN="*")
    preds$cell.samples <- tmp
  }
    
  #Winsorisation, if requested (default)
  if( control$winsor){
    my.quants <- apply( preds$cell.samples, 1, stats::quantile, probs=c(control$percent, 1-control$percent))
    if( control$tail %in% c( "both", "upper")){
      tmp <- sweep( preds$cell.samples, 1, my.quants[2,], "-")
      tmp1 <- tmp > 0
      tmp2 <- sweep( tmp1, 1, my.quants[2,], "*")
      preds$cell.samples[tmp1==1] <- tmp2[tmp1==1]
    }
    if( control$tail %in% c( "both","lower")){
      tmp <- sweep( preds$cell.samples, 1, my.quants[1,], "-")  #preds$cell.samples should all be non-negative
      tmp1 <- tmp < 0
      tmp2 <- sweep( tmp1, 1, my.quants[1,], "*")
      preds$cell.samples[tmp1==1] <- tmp2[tmp1==1]
    }
  }
 
  #summary over cells.
  samplePopN <- colSums( preds$cell.samples)
  quants <- stats::quantile( samplePopN, probs=control$probs, na.rm=TRUE)
  
  #adding some sampling noise (this is prediction after all).
  samplePopN.pred <- stats::rpois( n=length( samplePopN), lambda=samplePopN)
  quants.preds <- stats::quantile( samplePopN.pred, probs=control$probs, na.rm=TRUE)
  
  res <- list( mean=mean( samplePopN, na.rm=TRUE), median=stats::median( samplePopN, na.rm=TRUE), interval=quants,
               mean.pred=mean( samplePopN.pred, na.rm=TRUE), median.pred=stats::median( samplePopN.pred, na.rm=TRUE), interval.preds=quants.preds)
  
  return( res)
}


