
###############################################################################################
###############################################################################################
####	
####	Predict distribution and/or bias and/or artefact from an isdm
####
####	Returns a list of predictions
####
####	Programmed by Scott in the first half of 2022
####
###############################################################################################
###############################################################################################

#Function to get prediction from a fitted INLA model.
predict.isdm <- function( object, covarRaster, S=500, intercept.terms=NULL, n.threads=NULL, includeRandom=TRUE, includeFixed=TRUE, includeBias=FALSE, type="intensity", ...){
  
  #check if there's anything to do.
  if( !includeRandom & !includeFixed & !includeBias)
    stop( "Neither fixed, random, nor bias included in model predictions. Please choose one or more, but probably fixed and random.")
  
  #determine the number of threads to use.  Default is to use the same as the fit
  if( is.null( n.threads))
    n.threads <- attr( object, 'n.threads')
  
  #check to see if the intercept supplied is useable.  
  if( !all( intercept.terms %in% object$mod$names.fixed))
    stop( "One or more of the intercept.terms supplied is not in the model.  Please check.")
  
  ####Get (grid of) predictions
  #build predictions
  samples <- draw.posterior.samps(object$mod, B=S, what="effects", field="isdm.spat.XXX")

  #a data.frame containing prediciton points (no NAs).
  #add cell areas first
  covarRaster <- raster::addLayer( covarRaster, raster::area( covarRaster))#raster::raster( terra::cellSize( terra::rast( covarRaster))))
  names( covarRaster)[raster::nlayers( covarRaster)] <- "myCellAreas"
  #get the coordinates of the prediction points
  predcoords <- raster::coordinates( covarRaster)
  #extract the covariates
  rasterVarNames <- getVarNames( object$distributionFormula, object$biasFormula, NULL)
  #edit from Andrew to handle as.factor()...  Edit made 4/11/22
  tf2 <- grepl("as.factor",rasterVarNames)
  if(any(tf2)){rasterVarNames[tf2] <- gsub("as.factor[(]","",rasterVarNames[tf2])
			   rasterVarNames[tf2] <- gsub("[)]","",rasterVarNames[tf2])}


  covarData <- as.data.frame( raster::extract( covarRaster, predcoords[,1:2])[,rasterVarNames, drop=FALSE])
  myCellAreas <- as.matrix( raster::extract( covarRaster, predcoords[,1:2])[,"myCellAreas", drop=FALSE])
  #cut down to just those areas without NAs.
  noNAid <- apply( covarData, 1, function(x) !any( is.na( x)))
  predcoords <- predcoords[noNAid,]
  covarData <- covarData[noNAid,, drop=FALSE]
  myCellAreas <- myCellAreas[noNAid,,drop=FALSE]
  
  #predictions start with cell area
  #always included no matter what components are included.
  eta <- log( as.numeric( myCellAreas))
  
  #predictions due to only fixed effects
  if( includeFixed==TRUE){
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
      newSampsFixedDist <- rbind( newSampsFixedDist, matrix( 0, nrow=length( missedLevels), ncol=ncol( newSampsFixedDist)))
    }
    #the fixed effects involved in this term
    fix.subset <- which( fix.names %in% colnames( X))
    fix.names.ord <- order( fix.names[fix.subset])
    #ordering
    newSampsFixedDist <- newSampsFixedDist[fix.subset[fix.names.ord],]
    X <- X[,order( colnames( X)),drop=FALSE]
    #the addition to the linear predictor
    eta <- eta + X %*% newSampsFixedDist
  }
  
  #adding in the random effects, if present and wanted
  if( length( samples$fieldAtNodes[[1]])!=0 & includeRandom==TRUE){
    #projector matrix( linking prediction points to mesh)
    A.prd <- INLA::inla.spde.make.A( object$mesh, loc=predcoords)
    #the addition to the linear predictor
    eta <- eta + A.prd %*% samples$fieldAtNodes
  }
  
  #adding in the bias terms, if wanted.
  if( includeBias==TRUE){
    #sampling bias formula
    bform <- object$biasFormula
    #sampling bias model.matrix
    bX <- stats::model.matrix( bform, data=covarData)
    #fixing up intercept name
    tmpID <- grepl( "(Intercept)", colnames( bX))
    colnames( bX)[tmpID] <- "Intercept.PO"
    colnames( bX)[!tmpID] <- paste0( "Intercept.PO:",colnames( bX)[!tmpID])
    bX <- bX[,order( colnames( bX))]
    #sorting samples and design matrix
    tmpID1 <- grep( "Intercept.PO", object$mod$names.fixed)
    bfixedNames <- object$mod$names.fixed[tmpID1][order( tmpID1)]
    newSampsFixedBias <- samples$fixedEffects[tmpID1,]
    newSampsFixedBias <- newSampsFixedBias[order( tmpID1),]
    #the addition to the linear predictor
    eta <- eta + bX %*% newSampsFixedBias
  }
  
  #putting together on prediction scale
  mu.all <- NULL
  if( type=='intensity')
    mu.all <- as.matrix( exp( eta))
  if( type=='probability')
    mu.all <- as.matrix( 1-exp( -exp( eta)))
  if( type=='link')
    mu.all <- as.matrix( eta)
  if( is.null( mu.all) & type != "link")
    stop( "unknown type.  Must be 'intensity', 'probability' or 'link'. Please check function call.")
  #summaries
  mu.median <- apply( mu.all, 1, stats::quantile, probs=0.5, na.rm=TRUE)
  mu.lower <- apply( mu.all, 1, stats::quantile, probs=0.025, na.rm=TRUE)
  mu.upper <- apply( mu.all, 1, stats::quantile, probs=0.975, na.rm=TRUE)
  mu.mean <- rowMeans( mu.all)
  mu.sd <- apply( mu.all, 1, stats::sd)
  
  #raster format
  muRaster <- raster::rasterFromXYZ( cbind( predcoords, mu.median), crs=raster::crs( covarRaster))
  muRaster <- raster::addLayer(muRaster, raster::rasterFromXYZ( cbind( predcoords,mu.lower), crs=raster::crs( covarRaster)))
  muRaster <- raster::addLayer(muRaster, raster::rasterFromXYZ( cbind( predcoords,mu.upper), crs=raster::crs( covarRaster)))
  muRaster <- raster::addLayer( muRaster, raster::rasterFromXYZ(cbind( predcoords, mu.mean), crs=raster::crs( covarRaster)))
  muRaster <- raster::addLayer(muRaster, raster::rasterFromXYZ( cbind( predcoords,mu.sd), crs=raster::crs( covarRaster)))

  #sort out extent in case...
  muRaster <- raster::extend( muRaster, covarRaster)  #just in case it is needed -- could be dropped throughout the creation of the raster.
  
  res <- list( mean.field=muRaster, cell.samples=mu.all, fixedSamples=samples$fixedEffects, fixed.names=object$mod$names.fixed, predLocats=predcoords)
  
  return( res)
}

