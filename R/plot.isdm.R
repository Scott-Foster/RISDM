
###############################################################################################
###############################################################################################
####	
####	create diagnostic plots for isdm objects.
####
####	Returns NULL but does so invisibily
####
####	Programmed by Scott in the first half of 2022
####
###############################################################################################
###############################################################################################

plot.isdm <- function( x, covarRaster, nFigRow=1, ask=TRUE, ...){
 
  #number of data types
  numTypes <- 0
  #number of columns in the plot
  ncolly <- 2
 
  #containers for the residuals
  DCresids <- AAresids <- PAresids <- POresids <- NULL
  #DC residuals
  if( "Intercept.DC" %in% x$mod$names.fixed){
    #increment the number of data types.
    numTypes <- numTypes+1
    #get the fitted values as predictions
    preds <- x$mod$summary.fitted.values[INLA::inla.stack.index( x$stack, "DC")$data,"mean"]
    #get the outcomes and arrange them appropriately
    outcomes <- x$observationList$DCdat$DCcountDC
    #calculate the two probs for RQR (lower and upper)
    tmp1 <- stats::ppois( outcomes, lambda=preds)
    tmp2 <- stats::ppois( pmax( outcomes-1, 0), lambda=preds)
    #make sure that lower isn't too low...
    tmp2[outcomes==0] <- 0
    #do the random part of RQR
    tmp3 <- stats::runif( n=length( tmp1), min=tmp2, max=tmp1)
    #bundle the residuals
    DCresids <- data.frame( fitted=preds, observed=outcomes, residual=stats::qnorm( tmp3))
  }
  if( "Intercept.AA" %in% x$mod$names.fixed){
    #increment the number of data types.
    numTypes <- numTypes+1
    #get the fitted values as predictions
    preds <- x$mod$summary.fitted.values[INLA::inla.stack.index( x$stack, "AA")$data,"mean"]
    #get the outcomes and arrange them appropriately
    outcomes <- x$observationList$AAdat[,x$responseNames["AA"]]
    #calculate the two probs for RQR (lower and upper)
    tmp1 <- stats::ppois( outcomes, lambda=preds)
    tmp2 <- stats::ppois( pmax( outcomes-1, 0), lambda=preds)
    #make sure that lower isn't too low...
    tmp2[outcomes==0] <- 0
    #do the random part of RQR
    tmp3 <- stats::runif( n=length( tmp1), min=tmp2, max=tmp1)
    #bundle the residuals
    AAresids <- data.frame( fitted=preds, observed=outcomes, residual=stats::qnorm( tmp3))
  }
  if( "Intercept.PA" %in% x$mod$names.fixed){
    #increment the number of data types.
    numTypes <- numTypes+1
    #get the fitted values as predictions
    preds <- x$mod$summary.fitted.values[INLA::inla.stack.index( x$stack, "PA")$data,"mean"]
    #provide minor adjustments to preds that are identically 0 or 1.
    preds[preds==0] <- min( min( preds[preds!=0])/2, 1e-4)
    preds[preds==1] <- max( max( (1-preds[preds!=1])/2), 1-1e-4)
    #get the outcomes and arrange them appropriately
    outcomes <- x$observationList$PAdat[,x$responseNames["PA"]]
    #calculate the two probs for RQR (lower and upper)
    tmp1 <- stats::pbinom( outcomes, size=1, prob=preds)
    tmp2 <- stats::pbinom( pmax( outcomes-1, 0), size=1, prob=preds)
    #make sure that lower isn't too low...
    tmp2[outcomes==0] <- 0
    #do the random part of RQR
    tmp3 <- stats::runif( n=length( tmp1), min=tmp2, max=tmp1)
    #bundle the residuals
    PAresids <- data.frame( fitted=preds, observed=outcomes, residual=stats::qnorm( tmp3))
  }
  
  if( "Intercept.PO" %in% x$mod$names.fixed){
    #increment the number of data types.
    numTypes <- numTypes+1
    #increase the plotting columns to 3 (default is 2)
    ncolly <- 3

    
    POspP <- sp::SpatialPoints( x$observationList$PO[,attr( x, "coord.names")], proj4string=crs( covarRaster))
    #get the outcomes and arrange them appropriately
    rasCount <- raster::rasterize( POspP, covarRaster, fun='count', background=0)
    rasCount <- raster::mask( rasCount, covarRaster[[1]])

#####  General solution, but use the quick one for now...
    #number of draws to use in the simulation
#    if( !methods::hasArg( S))
#      S <- 250
#    message( paste0("Generating ",S," samples to form prediction (with distribution, random and bias effects)."))
    #get the fitted values as predictions
#    preds1 <- predict( x, covarRaster, intercept.terms=NULL, type='intensity', S=S, includeFixed=TRUE, includeRandom=TRUE, includeBias=TRUE)

#####	Quick running solution
    preds <- list()
    preds$mean.field <- list()
    preds$mean.field$mu.mean <- raster::rasterFromXYZ( cbind( coordinates( covarRaster), x$mod$summary.fitted.values[inla.stack.index(x$stack,"PO")$data,"mean"]))
    preds$mean.field$mu.mean <- raster::mask( preds$mean.field$mu.mean, rasCount)
    
    #calculate the two probs for RQR (lower and upper)
    tmp1 <- stats::ppois( raster::values( rasCount), raster::values( preds$mean.field$mu.mean))
    tmp2 <- stats::ppois( pmax( raster::values( rasCount)-1, 0), raster::values( preds$mean.field$mu.mean))
    #make sure that lower isn't too low...
    tmp2[raster::values( rasCount)==0] <- 0
    #trim the NAs
    na.id <- is.na( tmp1) | is.na( tmp2)
    tmp1 <- tmp1[!na.id]
    tmp2 <- tmp2[!na.id]
    #do the random part of RQR
    tmp3 <- stats::runif( n=length( tmp1), max=tmp1, min=tmp2)
    ressy <- rep( NA, length( raster::values( rasCount)))
    ressy[!na.id] <- stats::qnorm( tmp3)
    #bundle the residuals
    POresids <- list()
    POresids$ras <- raster::rasterFromXYZ( cbind( raster::coordinates( rasCount), ressy))
#    #re-extend the residual field to match the original field
#    POresids$ras <- raster::extend( POresids$ras, preds$mean.field$mu.mean)
    POresids$POresids <- data.frame( fitted=raster::values( preds$mean.field$mu.mean), observed=values( rasCount), residual=ressy)
  }
  
  #plotting up the residuals. RQR versus fitted and for PO data the raster.
  if( ask){
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit( grDevices::devAskNewPage(FALSE))
  }
  graphics::par( mfrow=c(nFigRow,3))
  
#  oask <- grDevices::devAskNewPage(FALSE) #don't ask for the first page of plots
  if( "Intercept.DC" %in% x$mod$names.fixed){
    plot.new()
    plot( DCresids$fitted, DCresids$residual, pch=20, ylab="DC residuals", xlab="DC fitted", main="Double Count")
    graphics::abline( h=0, col='green')
    stats::qqnorm( DCresids$residual, pch=20, ylab="DC quantile", main="Double Count")
    stats::qqline( DCresids$residual, col='green')
  }
  if( "Intercept.AA" %in% x$mod$names.fixed){
    plot.new()
    plot( AAresids$fitted, AAresids$residual, pch=20, ylab="AA residuals", xlab="AA fitted", main="Abundance-Absence")
    graphics::abline( h=0, col='green')
    stats::qqnorm( AAresids$residual, pch=20, ylab="AA quantile", main="Abundance-Absence")
    stats::qqline( AAresids$residual, col='green')
  }
  if( "Intercept.PA" %in% x$mod$names.fixed){
    plot.new()
    plot( PAresids$fitted, PAresids$residual, pch=20, ylab="PA residuals", xlab="PA fitted", main="Presence-Absence")
    graphics::abline( h=0, col='green')
    stats::qqnorm( PAresids$residual, pch=20, ylab="PA quantile", main="Presence-Absence")
    stats::qqline( PAresids$residual, col='green')
  }
  if( "Intercept.PO" %in% x$mod$names.fixed){
    raster::plot( POresids$ras, main="Presence-Only")
    plot( POresids$POresids$fitted, POresids$POresids$residual, pch=20, ylab="PO residuals", xlab="PO fitted", main="Presence-Only")
    graphics::abline( h=0, col='green')
    stats::qqnorm( POresids$POresids$residual, pch=20, ylab="PO quantile", main="Presence Only")
    stats::qqline( POresids$POresids$residual, col='green')
  }

  invisible( NULL)
  
}
  
