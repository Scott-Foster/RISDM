
plot.isdm <- function( object, covarRaster, ...){
  message( paste0("Generating ",n.predSamps," samples to form prediction (with distribution, random and bias effects)."))
  if( !hasArg( S))
    S <- 250
  
  preds <- predict( object, covarRaster, intercept.terms=NULL, type='intensity', S=S, includeFixed=TRUE, includeRandom=TRUE, includeBias=TRUE)
  
  
  
}
  
