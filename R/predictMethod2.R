
#Function to get prediction from a fitted INLA model.
predict.isdm <- function( fit, covarRaster, S=500, intercept.terms=NULL, n.threads=NULL){
  
  ####determine the number of threads to use.  Default is to use the same as the fit
  if( is.null( n.threads))
    n.threads <- attr( fit, 'n.threads')
  
  ####Get (grid of) predictions
  #build predictions
  message( "INLA::draw.posterior.samps will sometimes produce warnings. These seem to be internally correct -- please ignore (for now).")
  samples <- draw.posterior.samps(fit$mod, B=S, what="effects", field="isdm.spat.XXX")
  message( "Any warnings from now on should be taken more seriously.")
  
  #a data.frame containing prediciton points (no NAs).
  #add cell areas first
  covarRaster <- raster::addLayer( covarRaster, raster::raster( terra::cellSize( terra::rast( covarRaster))))
  names( covarRaster)[raster::nlayers( covarRaster)] <- "myCellAreas"
  #get the coordinates of the prediction points
  predcoords <- raster::coordinates( covarRaster)
  #extract the covariates
  rasterVarNames <- getVarNames( fit$distributionFormula, NULL, NULL)
  covarData <- as.data.frame( raster::extract( covarRaster, predcoords[,1:2])[,rasterVarNames, drop=FALSE])
  myCellAreas <- as.matrix( raster::extract( covarRaster, predcoords[,1:2])[,"myCellAreas", drop=FALSE])
  #cut down to just those areas without NAs.
  noNAid <- apply( covarData, 1, function(x) !any( is.na( x)))
  predcoords <- predcoords[noNAid,]
  covarData <- covarData[noNAid,, drop=FALSE]
  myCellAreas <- myCellAreas[noNAid,,drop=FALSE]
  
  if( !is.null( intercept.terms))
    for( ii in 1:length( intercept.terms)){
      covarData$tmptmptmptmp <- 1
      colnames( covarData)[ ncol( covarData)] <- intercept.terms[ii]
    }
  else{
    if( "Intercept.DC" %in% fit$mod$names.fixed) 
      covarData$Intercept.DC <- 1
    else{
      if( "Intercept.AA" %in% fit$mod$names.fixed)
	covarData$Intercept.AA <- 1
      else {
	if( "Intercept.PA" %in% fit$mod$names.fixed)
	  covarData$Intercept.PA <- 1
	else {
	  if( "Intercept.PO" %in% fit$mod$names.fixed)
	    covarData$Intercept.PO <- 1
	}
      }
    }
  }
  
  #the model matrix
  myForm <- fit$distributionFormula
#  if( "Intercept.DC" %in% colnames( covarData))
#    myForm <- update( myForm, paste0("~.-1+Intercept.DC/",attr(fit,"DCobserverInfo")$SurveyID))
#  else
  myForm <- update( myForm, paste0("~.-1+",paste( intercept.terms, collapse="+")))  #colnames( covarData)[grep( "Intercept.", colnames( covarData))]))
  
  X <- model.matrix( myForm, data=covarData)
  
  #sorting the design matrix and the effects so that they match
  fix.names <- fit$mod$names.fixed
  #those (factor levels) that are made when producing (constrained/unconstrained) design matrix
  #update 30/5/2022 -- this shouldn't do anything?
  missedLevels <- setdiff( colnames( X), fix.names)
  if( length( missedLevels)>0){
    fix.names <- c( fix.names, missedLevels)
    samples$fixedEffects <- rbind( samples$fixedEffects, matrix( 0, nrow=length( missedLevels), ncol=ncol( samples$fixedEffects)))
  }
  
  fix.subset <- which( fix.names %in% colnames( X))
  fix.names.ord <- order( fix.names[fix.subset])
  samples$fixedEffects <- samples$fixedEffects[fix.subset[fix.names.ord],]
  X <- X[,order( colnames( X))]
  
  #the size of a grid cell
#  areaGrid <- prod( raster::res( covarRaster))  #assumed to be the same for all cells (not lat long)

  #predictions due to only fixed effects
  eta <- X %*% samples$fixedEffects + log( as.numeric( myCellAreas))
  
  #adding in the random effects, if present
  if( length( samples$fieldAtNodes[[1]])!=0){
    #projector matrix( linking prediction points to mesh)
    A.prd <- INLA::inla.spde.make.A( fit$mesh, loc=predcoordsExpand)
    eta <- eta + A.prd %*% samples$fieldAtNodes
  }
  mu.all <- as.matrix( exp( eta))
#  if( nDClevs>0){
#    cellIDs <- rep( 1:(nrow( covarData) %/% 3), each=nDClevs)
#    mu.cell.mean <- apply( mu.all, 2, function(xx) tapply( xx, cellIDs, mean))
#  }
#  else
#    mu.cell.mean <- mu.all
  mu.mean <- rowMeans( mu.all)  #mu.cell.mean)
  mu.sd <- apply( mu.all, 1, sd)  #mu.cell.mean, 1, sd)
    
  muRaster <- raster::rasterFromXYZ( cbind( predcoords, mu.mean), crs=raster::crs( covarRaster))
  muRaster <- raster::addLayer(muRaster, raster::rasterFromXYZ( cbind( predcoords,mu.sd), crs=raster::crs( covarRaster)))
  
  res <- list( mean.field=muRaster, samples.all=mu.all, predLocats=predcoords)
  
  return( res)
}

PopEstimate <- function( preds, probs=c(0.025,0.975)){
  samplePopN <- colSums( preds$sample.cell.mean)
  quants <- quantile( samplePopN, probs=probs)
  
#  sample.preds <- matrix( rpois(n=prod( dim( preds$sample.cell.mean)), lambda=preds$sample.cell.mean), 
#                          nrow=nrow( preds$sample.cell.mean), ncol=ncol( preds$sample.cell.mean))
  samplePopN.pred <- rpois( n=length( samplePopN), lambda=samplePopN)
  quants.preds <- quantile( samplePopN.pred, probs=probs)
  
  res <- list( mean=mean( samplePopN), median=median( samplePopN), interval=quants,
               mean.pred=mean( samplePopN.pred), median.pred=median( samplePopN.pred), interval.preds=quants.preds)
  
  return( res)
}

