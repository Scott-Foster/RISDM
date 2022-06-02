
plot.isdm <- function( object, covarRaster, ...){
  message( paste0("Generating ",n.predSamps," samples to form prediction"))
  preds <- predict( object, covarRaster, intercept.terms=NULL, type='intensity', S=S)
}
  