
plot.isdm <- function( object, covarRaster, ...){
 
  if( !hasArg( S))
    S <- 250
  message( paste0("Generating ",S," samples to form prediction (with distribution, random and bias effects)."))
   
  preds <- predict( object, covarRaster, intercept.terms=NULL, type='intensity', S=S, includeFixed=TRUE, includeRandom=TRUE, includeBias=TRUE)
  
  POspP <- sp::SpatialPoints( object$observationList$PO[,attr( object, "coord.names")], proj4string=crs( covarRaster))
  
  rasCount <- raster::rasterize( POspP, covarRaster, fun='count')
  
  
  
  
#  inla.stack.index( fm$VIC$stack, "AA")$data  #from Krainski et al.
  
}
  
