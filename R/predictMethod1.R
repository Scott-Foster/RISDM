
#Function to get prediction from a fitted INLA model.
predict.isdm <- function( fit, covarRaster, S=500, DCaverage=TRUE){
  
  ####Get (grid of) predictions
  #build predictions
  samples <- draw.posterior.samps(fit$mod, B=S, what="effects", field="isdm.spat.XXX")
  
  #a data.frame containing prediciton points (no NAs).
  #add cell areas first
  covarRaster <- addLayer( covarRaster, raster::raster( terra::cellSize( terra::rast( covarRaster))))
  names( covarRaster)[nlayers( covarRaster)] <- "myCellAreas"
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
  if( any( grepl( paste0("Intercept.DC:",attr(fit,"DCobserverInfo")$SurveyID), fit$mod$names.fixed))){
    covarData$Intercept.DC <- 1
    nDClevs <- length( attr( fit, "DCSurveyIDLevels"))
    tmp <- covarData[rep( 1:nrow( covarData), each=nDClevs),, drop=FALSE]
    tmp$tmptmptmp <- as.factor( rep( attr( fit, "DCSurveyIDLevels"), each=nrow( covarData)))
    colnames(tmp)[colnames( tmp)=="tmptmptmp"] <- attr(fit,"DCobserverInfo")$SurveyID
    covarData <- tmp
    rm(tmp)
    predcoordsExpand <- predcoords[rep( 1:nrow( predcoords), each=nDClevs),]
  }
  else {
    nDClevs <- 0
    predcoordsExpand <- predcoords
    if( "Intercept.AA" %in% fit$mod$names.fixed)
      covarData$Intercept.AA <- 1
    else {
      if( "Intercept.PA" %in% fit$mod$names.fixed)
        covarData$Intercept.PA <- 1
      else {
        if( "Intercept.PO" %in% fit$mod$names.fixed)
          covarData$Intercept.PO <- 1
      }}}
  #the model matrix
  myForm <- fit$distributionFormula
  if( "Intercept.DC" %in% colnames( covarData))
    myForm <- update( myForm, paste0("~.-1+Intercept.DC/",attr(fit,"DCobserverInfo")$SurveyID))
  else
    myForm <- update( myForm, paste0("~.-1+",colnames( covarData)[grep( "Intercept.", colnames( covarData))]))
  
  X <- model.matrix( myForm, data=covarData)
  
  #sorting the design matrix and the effects so that they match
  fix.names <- fit$mod$names.fixed
  fix.subset <- which( fix.names %in% colnames( X))
  fix.names.ord <- order( fix.names[fix.subset])
  samples$fixedEffects <- samples$fixedEffects[fix.subset[fix.names.ord],]
  X <- X[,order( colnames( X))]
  
  #the size of a grid cell
  areaGrid <- prod( raster::res( covarRaster))  #assumed to be the same for all cells (not lat long)

  #predictions due to only fixed effects
  eta <- X %*% samples$fixedEffects + log( as.numeric( myCellAreas))
  
  #adding in the random effects, if present
  if( length( samples$fieldAtNodes[[1]])!=0){
    #projector matrix( linking prediction points to mesh)
    A.prd <- INLA::inla.spde.make.A( fit$mesh, loc=predcoordsExpand)
    eta <- eta + A.prd %*% samples$fieldAtNodes
  }
  mu.all <- as.matrix( exp( eta))
  if( nDClevs>0){
    cellIDs <- rep( 1:(nrow( covarData) %/% 3), each=nDClevs)
    mu.cell.mean <- apply( mu.all, 2, function(xx) tapply( xx, cellIDs, mean))
  }
  else
    mu.cell.mean <- mu.all
  mu.mean <- rowMeans( mu.cell.mean)
  mu.sd <- apply( mu.cell.mean, 1, sd)
  
  
  muRaster <- raster::rasterFromXYZ( cbind( predcoords, mu.mean), crs=raster::crs( covarRaster))
  muRaster <- raster::addLayer(muRaster, raster::rasterFromXYZ( cbind( predcoords,mu.sd), crs=raster::crs( covarRaster)))
  
  res <- list( mean.field=muRaster, sample.cell.mean=mu.cell.mean, samples.all=mu.all, predLocats=predcoords)
  
  return( res)
}

PopEstimate <- function( preds, probs=c(0.025,0.975)){
  samplePopN <- colSums( preds$sample.cell.mean)
  quants <- quantile( samplePopN, probs=probs)
  res <- list( mean=mean( samplePopN), median=median( samplePopN), interval=quants)
  
  return( res)
}

