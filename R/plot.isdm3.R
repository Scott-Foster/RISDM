
plot.isdm <- function( object, covarRaster, ...){
 
  #number of data types
  numTypes <- 0
  #number of columns in the plot
  ncolly <- 2
 
  DCresids <- AAresids <- PAresids <- POresids <- NULL
  if( "Intercept.DC" %in% object$mod$names.fixed){
    numTypes <- numTypes+1
    preds <- object$mod$summary.fitted.values[INLA::inla.stack.index( object$stack, "DC")$data,"mean"]
    outcomes <- object$observationList$DCdat$DCcountDC
    tmp1 <- stats::ppois( outcomes, lambda=preds)
    tmp2 <- stats::ppois( pmax( outcomes-1, 0), lambda=preds)
    tmp2[outcomes==0] <- 0
    tmp3 <- stats::runif( n=length( tmp1), min=tmp2, max=tmp1)
    DCresids <- data.frame( fitted=preds, observed=outcomes, residual=qnorm( tmp3))
  }
  if( "Intercept.AA" %in% object$mod$names.fixed){
    numTypes <- numTypes+1
    preds <- object$mod$summary.fitted.values[INLA::inla.stack.index( object$stack, "AA")$data,"mean"]
    outcomes <- object$observationList$AAdat[,object$responseNames["AA"]]
    tmp1 <- stats::ppois( outcomes, lambda=preds)
    tmp2 <- stats::ppois( pmax( outcomes-1, 0), lambda=preds)
    tmp2[outcomes==0] <- 0
    tmp3 <- stats::runif( n=length( tmp1), min=tmp2, max=tmp1)
    AAresids <- data.frame( fitted=preds, observed=outcomes, residual=qnorm( tmp3))
  }
  if( "Intercept.PA" %in% object$mod$names.fixed){
    numTypes <- numTypes+1
    preds <- object$mod$summary.fitted.values[INLA::inla.stack.index( object$stack, "PA")$data,"mean"]
    outcomes <- object$observationList$AAdat[,object$responseNames["PA"]]
    tmp1 <- stats::pbinom( outcomes, size=1, prob=preds)
    tmp2 <- stats::pbinom( pmax( outcomes-1, 0), size=1, prob=preds)
    tmp2[outcomes==0] <- 0
    tmp3 <- stats::runif( n=length( tmp1), min=tmp2, max=tmp1)
    PAresids <- data.frame( fitted=preds, observed=outcomes, residual=qnorm( tmp3))
  }
  
  if( "Intercept.PO" %in% object$mod$names.fixed){
    numTypes <- numTypes+1
    ncolly <- 3
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
    tmp3 <- stats::runif( n=length( tmp1), max=tmp1, min=tmp2)
    POresids <- list()
    ressy <- qnorm( tmp3)
    POresids$ras <- raster::rasterFromXYZ( cbind( raster::coordinates( rasCount), ressy))
    POresids$POresids <- data.frame( fitted=raster::values( preds$mean.field$mu.mean), observed=values( rasCount), residual=ressy)
  }
  
  graphics::par( mfrow=c(numTypes,ncolly))
  if( "Intercept.DC" %in% object$mod$names.fixed){
    if( ncolly==3)
      graphics::plot.new()
    plot( DCresids$fitted, DCresids$residual, pch=20, ylab="DC residuals", xlab="DC fitted", main="Double Count")
    graphics::abline( h=0, col='green')
    stats::qqnorm( DCresids$residual, pch=20, ylab="DC quantile", main="Double Count")
    stats::qqline( DCresids$residual, col='green')
  }
  if( "Intercept.AA" %in% object$mod$names.fixed){
    if( ncolly==3)
      graphics::plot.new()
    plot( AAresids$fitted, AAresids$residual, pch=20, ylab="AA residuals", xlab="AA fitted", main="Abundance-Absence")
    graphics::abline( h=0, col='green')
    stats::qqnorm( AAresids$residual, pch=20, ylab="AA quantile", main="Abundance-Absence")
    stats::qqline( AAresids$residual, col='green')
  }
  if( "Intercept.PA" %in% object$mod$names.fixed){
    if( ncolly==3)
      graphics::plot.new()
    plot( PAresids$fitted, PAresids$residual, pch=20, ylab="PA residuals", xlab="PA fitted", main="Presence-Absence")
    graphics::abline( h=0, col='green')
    stats::qqnorm( values( PAresids$residual), pch=20, ylab="PA quantile", main="Presence-Absence")
    stats::qqline( values( PAresids$residual), col='green')
  }
  if( "Intercept.PO" %in% object$mod$names.fixed){
    raster::plot( POresids$ras, main="Presence-Only Residuals")
#    plot( POresids$POresids$fitted, POresids$POresids$residual, col=c("black","blue")[(POresids$POresids$observed>0)+1], pch=20)
    plot( POresids$POresids$fitted, POresids$POresids$residual, pch=20, ylab="PO residuals", xlab="PO fitted", main="Prsence-Only")
    graphics::abline( h=0, col='green')
#    stats::qqnorm( POresids$POresids$residual, col=c("black","blue")[(POresids$POresids$observed>0)+1], pch=20)
    stats::qqnorm( POresids$POresids$residual, pch=20, ylab="PO quantile", main="Presence Only")
    stats::qqline( POresids$POresids$residual, col='green')
#    legend( "bottomright", legend=c("Zero","Non-zero"), pch=20, col=c("black","blue"))
  }

#  res <- list( DC=DCresids, AA=AAresids, PA=PAresids, PO=POresids)
#  invisible( res)
  invisible()
  
}
  
