
#Function to get prediction from a fitted INLA model.
predict.isdm <- function( object, covarRaster, S=500, intercept.terms=NULL, n.threads=NULL, includeRandom=TRUE, includeFixed=TRUE, includeBias=FALSE, type="intensity", ...){
  
  if( !includeRandom & !includeFixed & !includeBias)
    stop( "Neither fixed, random, nor bias included in model predictions. Please choose one or more, but probably fixed and random.")
  
  ####determine the number of threads to use.  Default is to use the same as the fit
  if( is.null( n.threads))
    n.threads <- attr( object, 'n.threads')
    
  if( !all( intercept.terms %in% object$mod$names.fixed))
    stop( "One or more of the intercept.terms supplied is not in the model.  Please check.")
  
  ####Get (grid of) predictions
  #build predictions
  samples <- draw.posterior.samps(object$mod, B=S, what="effects", field="isdm.spat.XXX")
  allFixedEffectSamples <- samples$fixedEffects

  #a data.frame containing prediciton points (no NAs).
  #add cell areas first
  covarRaster <- raster::addLayer( covarRaster, raster::raster( terra::cellSize( terra::rast( covarRaster))))
  names( covarRaster)[raster::nlayers( covarRaster)] <- "myCellAreas"
  #get the coordinates of the prediction points
  predcoords <- raster::coordinates( covarRaster)
  #extract the covariates
  rasterVarNames <- getVarNames( object$distributionFormula, object$biasFormula, NULL)
  covarData <- as.data.frame( raster::extract( covarRaster, predcoords[,1:2])[,rasterVarNames, drop=FALSE])
  myCellAreas <- as.matrix( raster::extract( covarRaster, predcoords[,1:2])[,"myCellAreas", drop=FALSE])
  #cut down to just those areas without NAs.
  noNAid <- apply( covarData, 1, function(x) !any( is.na( x)))
  predcoords <- predcoords[noNAid,]
  covarData <- covarData[noNAid,, drop=FALSE]
  myCellAreas <- myCellAreas[noNAid,,drop=FALSE]
  
  if( !is.null( intercept.terms)){
    intercept.terms.legal <- gsub( pattern=":", replacement="XCOLONX", x=intercept.terms)
    for( ii in 1:length( intercept.terms.legal)){
      covarData$tmptmptmptmp <- 1
      colnames( covarData)[ ncol( covarData)] <- intercept.terms.legal[ii]
    }
  }
  
  #the model matrix
  myForm <- object$distributionFormula
  if( !is.null( intercept.terms))
    myForm <- update( myForm, paste0("~.-1+",paste( intercept.terms.legal, collapse="+")))
  
  X <- stats::model.matrix( myForm, data=covarData)
  #undoing the hack from before.
  colnames( X) <- gsub( "XCOLONX", ":", colnames( X))
  
  #sorting the design matrix and the effects so that they match
  fix.names <- object$mod$names.fixed
  #those (factor levels) that are made when producing (constrained/unconstrained) design matrix
  #update 30/5/2022 -- this shouldn't do anything?
  missedLevels <- setdiff( colnames( X), fix.names)
  newSampsFixedDist <- samples$fixedEffects
  if( length( missedLevels)>0){
    fix.names <- c( fix.names, missedLevels)
    newSampsFixedDist <- rbind( samples$fixedEffects, matrix( 0, nrow=length( missedLevels), ncol=ncol( samples$fixedEffects)))
  }
  
  fix.subset <- which( fix.names %in% colnames( X))
  fix.names.ord <- order( fix.names[fix.subset])
  newSampsFixedDist <- newSampsFixedDist[fix.subset[fix.names.ord],]
  X <- X[,order( colnames( X))]
  
  #predictions start with cell area
  #always included no matter what components are included.
  eta <- log( as.numeric( myCellAreas))
  
  #predictions due to only fixed effects
  if( includeFixed==TRUE)
    eta <- eta + X %*% newSampsFixedDist
  
  #adding in the random effects, if present and wanted
  if( length( samples$fieldAtNodes[[1]])!=0 & includeRandom==TRUE){
    #projector matrix( linking prediction points to mesh)
    A.prd <- INLA::inla.spde.make.A( object$mesh, loc=predcoords)
    eta <- eta + A.prd %*% samples$fieldAtNodes
  }
  
  #adding in the bias terms, if wanted.
  if( includeBias==TRUE){
    bform <- object$biasFormula
    bX <- stats::model.matrix( bform, data=covarData)
    colnames( bX)[grep( "Intercept", colnames( bX))] <- "Intercept.PO"

    
    eta <- eta + 0
  }
  
  mu.all <- NULL
  if( type=='intensity')
    mu.all <- as.matrix( exp( eta))
  if( type=='probability')
    mu.all <- as.matrix( 1-exp( -exp( eta)))
  #if not 'intensity' or 'probability' then must(?!) be link
  if( is.null( mu.all) & type != "link")
    stop( "unknown type.  Must be 'intensity', 'probability' or 'link'. Please check function call.")
  mu.mean <- rowMeans( mu.all)
#  mu.median <- apply( mu.all, 1, stats::median)
  mu.sd <- apply( mu.all, 1, stats::sd)
  mu.lower <- apply( mu.all, 1, stats::quantile, probs=0.025)
  mu.upper <- apply( mu.all, 1, stats::quantile, probs=0.975)
      
  muRaster <- raster::rasterFromXYZ( cbind( predcoords, mu.mean), crs=raster::crs( covarRaster))
  muRaster <- raster::addLayer(muRaster, raster::rasterFromXYZ( cbind( predcoords,mu.sd), crs=raster::crs( covarRaster)))
#  muRaster <- raster::addLayer(muRaster, raster::rasterFromXYZ( cbind( predcoords,mu.median), crs=raster::crs( covarRaster)))
  muRaster <- raster::addLayer(muRaster, raster::rasterFromXYZ( cbind( predcoords,mu.lower), crs=raster::crs( covarRaster)))
  muRaster <- raster::addLayer(muRaster, raster::rasterFromXYZ( cbind( predcoords,mu.upper), crs=raster::crs( covarRaster)))
  
  res <- list( mean.field=muRaster, cell.samples=mu.all, fixedSamples=allFixedEffectSamples, fixed.names=object$mod$names.fixed, predLocats=predcoords)
  
  return( res)
}

PopEstimate <- function( preds, probs=c(0.025,0.975), intercept.terms=NULL){

  if( !all( intercept.terms %in% preds$fixed.names))
    stop( "intercept.terms specified are not part of the model. Please check call and model.")
  
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

  samplePopN <- colSums( preds$cell.samples)
  quants <- stats::quantile( samplePopN, probs=probs)
  
  samplePopN.pred <- stats::rpois( n=length( samplePopN), lambda=samplePopN)
  quants.preds <- stats::quantile( samplePopN.pred, probs=probs)
  
  res <- list( mean=mean( samplePopN), median=stats::median( samplePopN), interval=quants,
               mean.pred=mean( samplePopN.pred), median.pred=stats::median( samplePopN.pred), interval.preds=quants.preds)
  
  return( res)
}

