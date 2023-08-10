
###############################################################################################
###############################################################################################
####	
####	create residuals for isdm objects.
####
####	Returns a list with the different sets of residuals
####
####	Programmed by Scott in the May half of 2023
####
###############################################################################################
###############################################################################################

residuals.isdm <- function( object, ...){

  #number of data types
  numTypes <- 0
  #number of columns in the plot
  ncolly <- 2
 
  #containers for the residuals
  DCresids <- AAresids <- PAresids <- POresids <- NULL
  #DC residuals
  if( "DC_Intercept" %in% object$mod$names.fixed){
    #increment the number of data types.
    numTypes <- numTypes+1
    #get the fitted values as predictions
    preds <- object$mod$summary.fitted.values[INLA::inla.stack.index( object$stack, "DC")$data,"mean"]
    #get the outcomes and arrange them appropriately
    outcomes <- object$observationList$DCdat$DCcountDC
    #calculate the two probs for RQR (lower and upper)
    tmp1 <- stats::ppois( outcomes, lambda=preds)
    tmp2 <- stats::ppois( pmax( outcomes-1, 0), lambda=preds)
    #make sure that lower isn't too low...
    tmp2[outcomes==0] <- 0
    #guard against numerical erros
    tmp2[tmp2>tmp1] <- tmp1[tmp2>tmp1]
    #do the random part of RQR
    tmp3 <- stats::runif( n=length( tmp1), min=tmp2, max=tmp1)
    #make sure that there aren't any zeros or ones...
    tmp3 <- (1-1e-8) * (tmp3-0.5) + 0.5
    #bundle the residuals
    DCresids <- data.frame( fitted=preds, observed=outcomes, residual=stats::qnorm( tmp3))
  }
  if( "AA_Intercept" %in% object$mod$names.fixed){
    #increment the number of data types.
    numTypes <- numTypes+1
    #get the fitted values as predictions
    preds <- object$mod$summary.fitted.values[INLA::inla.stack.index( object$stack, "AA")$data,"mean"]
    #get the outcomes and arrange them appropriately
    outcomes <- object$observationList$AAdat[,object$responseNames["AA"]]
    #calculate the two probs for RQR (lower and upper)
    tmp1 <- stats::ppois( outcomes, lambda=preds)
    tmp2 <- stats::ppois( pmax( outcomes-1, 0), lambda=preds)
    #make sure that lower isn't too low...
    tmp2[outcomes==0] <- 0
    #guard against numerical erros
    tmp2[tmp2>tmp1] <- tmp1[tmp2>tmp1]
    #do the random part of RQR
    tmp3 <- stats::runif( n=length( tmp1), min=tmp2, max=tmp1)
    #make sure that there aren't any zeros or ones...
    tmp3 <- (1-1e-8) * (tmp3-0.5) + 0.5
    #bundle the residuals
    AAresids <- data.frame( fitted=preds, observed=outcomes, residual=stats::qnorm( tmp3))
  }
  if( "PA_Intercept" %in% object$mod$names.fixed){
    #increment the number of data types.
    numTypes <- numTypes+1
    #get the fitted values as predictions
    preds <- object$mod$summary.fitted.values[INLA::inla.stack.index( object$stack, "PA")$data,"mean"]
    #provide minor adjustments to preds that are identically 0 or 1.
    preds[preds==0] <- min( min( preds[preds!=0])/2, 1e-4)
    preds[preds==1] <- max( max( (1-preds[preds!=1])/2), 1-1e-4)
    #get the outcomes and arrange them appropriately
    outcomes <- object$observationList$PAdat[,object$responseNames["PA"]]
    #calculate the two probs for RQR (lower and upper)
    tmp1 <- stats::pbinom( outcomes, size=1, prob=preds)
    tmp2 <- stats::pbinom( pmax( outcomes-1, 0), size=1, prob=preds)
    #make sure that lower isn't too low...
    tmp2[outcomes==0] <- 0
    #guard against numerical erros
    tmp2[tmp2>tmp1] <- tmp1[tmp2>tmp1]
    #do the random part of RQR
    tmp3 <- stats::runif( n=length( tmp1), min=tmp2, max=tmp1)
    #make sure that there aren't any zeros or ones...
    tmp3 <- (1-1e-8) * (tmp3-0.5) + 0.5
    #bundle the residuals
    PAresids <- data.frame( fitted=preds, observed=outcomes, residual=stats::qnorm( tmp3))
  }
  
  if( "PO_Intercept" %in% object$mod$names.fixed){
    #increment the number of data types.
    numTypes <- numTypes+1
    #increase the plotting columns to 3 (default is 2)
    ncolly <- 3
    
    POspP <- sp::SpatialPoints( object$observationList$PO[,attr( object, "coord.names")], proj4string=crs( object$data$covars))
    #get the outcomes and arrange them appropriately
    rasCount <- raster::rasterize( POspP, object$data$covars, fun='count', background=0)
    rasCount <- raster::mask( rasCount, object$data$covars[[1]])
    ###Something to go in here...

#####  General solution, but use the quick one for now...
    #number of draws to use in the simulation
#    if( !methods::hasArg( S))
#      S <- 250
#    message( paste0("Generating ",S," samples to form prediction (with distribution, random and bias effects)."))
    #get the fitted values as predictions
#    preds1 <- predict( object, covars, intercept.terms=NULL, type='intensity', S=S, includeFixed=TRUE, includeRandom=TRUE, includeBias=TRUE)

#####	Quick running solution
    preds <- list()
    preds$mean.field <- list()
    preds$mean.field$mu.mean <- raster::rasterFromXYZ( cbind( raster::coordinates( rasCount)[!is.na( raster::values( rasCount)),], object$mod$summary.fitted.values[inla.stack.index(object$stack,"PO")$data,"mean"]))
    preds$mean.field$mu.mean <- raster::extend( preds$mean.field$mu.mean, rasCount)
    preds$mean.field$mu.mean <- raster::crop( preds$mean.field$mu.mean, rasCount)
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
    #guard against numerical erros
    tmp2[tmp2>tmp1] <- tmp1[tmp2>tmp1]
    #do the random part of RQR
    tmp3 <- stats::runif( n=length( tmp1), max=tmp1, min=tmp2)
    #make sure that there aren't any zeros or ones...
    tmp3 <- (1-1e-8) * (tmp3-0.5) + 0.5
    #return object
    ressy <- rep( NA, length( raster::values( rasCount)))
    ressy[!na.id] <- stats::qnorm( tmp3)
    #bundle the residuals
    POresids <- list()
    POresids$ras <- raster::rasterFromXYZ( cbind( raster::coordinates( rasCount), ressy))
#    #re-extend the residual field to match the original field
#    POresids$ras <- raster::extend( POresids$ras, preds$mean.field$mu.mean)
    POresids$POresids <- data.frame( fitted=raster::values( preds$mean.field$mu.mean), observed=values( rasCount), residual=ressy)
  }
  
  res <- list( DC=DCresids, AA=AAresids, PA=PAresids, PO=POresids)
  attr( res, "numTypes") <- numTypes
  attr( res, "ncolly") <- ncolly
  
  return( res)
  
}
