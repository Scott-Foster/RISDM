
plot.isdm <- function( object, covarRaster, ...){
 
  if( !hasArg( S))
    S <- 250
  message( paste0("Generating ",S," samples to form prediction (with distribution, random and bias effects)."))
   
  preds <- predict( object, covarRaster, intercept.terms=NULL, type='intensity', S=S, includeFixed=TRUE, includeRandom=TRUE, includeBias=TRUE)
  
  POspP <- sp::SpatialPoints( object$observationList$PO[,attr( object, "coord.names")], proj4string=crs( covarRaster))
  
  rasCount <- raster::rasterize( POspP, covarRaster, fun='count', background=0)
  rasCount <- raster::mask( rasCount, covarRaster[[1]])
  
  resids <- stats::ppois( raster::values( rasCount), raster::values( preds$mean.field$mu.mean))
  tmp <- stats::ppois( pmax( raster::values( rasCount)-1, 0), raster::values( preds$mean.field$mu.mean))
  tmp[values( rasCount)==0] <- 0
  tmp1 <- stats::runif( n=length( resids), min=resids, max=tmp)
  
  
  
  
  
  resids <- (values( rasCount) - preds$mean.field$mu.mean) / sqrt(preds$mean.field$mu.mean)
  resids <- raster::rasterFromXYZ( cbind( resids, raster::coordinates( rasCount)))
  
  
  
#  inla.stack.index( fm$VIC$stack, "AA")$data  #from Krainski et al.
  
}
  
