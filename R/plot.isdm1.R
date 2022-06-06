
plot.isdm <- function( object, covarRaster, ...){
  message( paste0("Generating ",n.predSamps," samples to form prediction (with distribution, random and bias effects)."))
  if( !hasArg( S))
    S <- 250
  
  preds <- predict( object, covarRaster, intercept.terms=NULL, type='intensity', S=S, includeFixed=TRUE, includeRandom=TRUE, includeBias=TRUE)
  
  POspP <- sp::spatialPoints( object$observationList$PO[,attr( object, "coord.names")], proj4string=crs( covarRaster))
  
  raster::rasterize( POspP, covarRaster, fun='count')
  
  
  
#  inla.stack.index( fm$VIC$stack, "AA")$data  #from Krainski et al.
  
}
  
