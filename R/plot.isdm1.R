
plot.isdm <- function( object, covarRaster, ...){
 
  if( !hasArg( S))
    S <- 250
  message( paste0("Generating ",S," samples to form prediction (with distribution, random and bias effects)."))
   
  preds <- predict( object, covarRaster, intercept.terms=NULL, type='intensity', S=S, includeFixed=TRUE, includeRandom=TRUE, includeBias=TRUE)
  
  POspP <- sp::SpatialPoints( object$observationList$PO[,attr( object, "coord.names")], proj4string=crs( covarRaster))
  
  rasCount <- raster::rasterize( POspP, covarRaster, fun='count', background=0)
  rasCount <- raster::mask( rasCount, covarRaster[[1]])
  
  tmp1 <- stats::ppois( raster::values( rasCount), raster::values( preds$mean.field$mu.mean))
  tmp2 <- stats::ppois( pmax( raster::values( rasCount)-1, 0), raster::values( preds$mean.field$mu.mean))
  tmp2[values( rasCount)==0] <- 0
  suppressWarnings( tmp3 <- stats::runif( n=length( tmp1), max=tmp1, min=tmp2))
  POresids <- pnorm( tmp3)
  
  POresids <- raster::rasterFromXYZ( cbind( raster::coordinates( rasCount), POresids))

  par( mfrow=c(1,2))
  raster::plot( POresids)
  stats::qqnorm( values( POresids))
  stats::qqline( values( POresids), col='red')



  
  
#  inla.stack.index( fm$VIC$stack, "AA")$data  #from Krainski et al.
  
}
  
