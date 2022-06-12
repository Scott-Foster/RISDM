
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

PopEstimate <- function( preds, probs=c(0.025,0.975), intercept.terms=NULL){

  #checking terms
  if( !all( intercept.terms %in% preds$fixed.names))
    stop( "intercept.terms specified are not part of the model. Please check call and model.")

  #assuming, if needed, or bundling if we have to
  if( is.null( intercept.terms)){
    message( "Assuming that Intercept terms have already been included in predictions.  Is this what you want? See ?predict.isdm for how to include them.")
  }
  else{
    int.contr <- preds$fixedSamples[preds$fixed.names %in% intercept.terms,,drop=FALSE]
    int.contr <- colSums( int.contr)
    int.contr <- exp( int.contr)
    tmp <- sweep( x=preds$cell.samples, MARGIN=2, STATS=int.contr, FUN="*")
    preds$cell.samples <- tmp
  }

  #summary over cells.
  samplePopN <- colSums( preds$cell.samples)
  quants <- stats::quantile( samplePopN, probs=probs)
  
  #adding some sampling noise (this is prediction after all).
  samplePopN.pred <- stats::rpois( n=length( samplePopN), lambda=samplePopN)
  quants.preds <- stats::quantile( samplePopN.pred, probs=probs)
  
  res <- list( mean=mean( samplePopN), median=stats::median( samplePopN), interval=quants,
               mean.pred=mean( samplePopN.pred), median.pred=stats::median( samplePopN.pred), interval.preds=quants.preds)
  
  return( res)
}


