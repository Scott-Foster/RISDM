
###############################################################################################
###############################################################################################
####	
####	Predict distribution and/or bias and/or artefact from an isdm
####  	When the prediction space is large or the number of samples is large.
####
####	Returns a list of predictions
####
####	Programmed by Scott in Nov/Dec 2022
####
###############################################################################################
###############################################################################################

#Function to get prediction from a fitted INLA model.
predict.isdm <- function( object, covars, habitatArea=NULL, S=500, intercept.terms=NULL, n.threads=NULL, n.batches=1, includeRandom=TRUE, includeFixed=TRUE, includeBias=FALSE, type="intensity", ...){
  
  #check if there's anything to do.
  if( !is.logical(includeFixed) & !all( includeFixed %in% names( covars)))
    stop( "Fixed effect terms specified are not in covariate layers")
  if( !includeRandom & any(includeFixed == FALSE) & !includeBias)
    stop( "Neither fixed, random, nor bias included in model predictions. Please choose one or more, but probably fixed and random.")
  
  #determine the number of threads to use.  Default is to use the same as the fit
  if( is.null( n.threads))
    n.threads <- attr( object, 'n.threads')
  
  #check to see if the intercept supplied is useable.  
  if( !is.null( intercept.terms) & !all( intercept.terms %in% object$mod$names.fixed))
    stop( "One or more of the intercept.terms supplied is not in the model.  Please check.")

#  batchStartEnd <- seq( from=0, to=S, by=ceiling( S/n.batches))
  batchStartEnd <- round( seq( from=0, to=S, by=S/n.batches),0)
#  if( tail( batchStartEnd,1) < S)
#    batchStartEnd <- c( batchStartEnd, S) #[length( batchStartEnd)] <- S  #last batch is smaller
  batchSize <- diff( batchStartEnd)

  #### function for drawing posterior sampls in parallel.
  funny <- function( ii){
    tmp.samples <- draw.posterior.samps(object$mod, B=batchSize[ii], what="effects", field="isdm.spat.XXX", n.threads=n.threads)
    return( tmp.samples)
  }

  samp.element <- funny(1)

  samples <- list()
  if( object$control$addRandom){
    samples$fieldAtNodes <- matrix( NA, nrow=nrow( samp.element$fieldAtNodes), ncol=S)
    samples$fieldAtNodes[,(batchStartEnd[1]+1):batchStartEnd[2]] <- samp.element$fieldAtNodes    
  }
  samples$fixedEffects <- matrix( NA, nrow=nrow( samp.element$fixedEffects), ncol=S)

  samples$fixedEffects[,(batchStartEnd[1]+1):batchStartEnd[2]] <- samp.element$fixedEffects
  
  if( n.batches > 1){
    for( ii in 2:n.batches){
      samp.element <- funny(ii)
      samples$fieldAtNodes[,(batchStartEnd[ii]+1):batchStartEnd[ii+1]] <- samp.element$fieldAtNodes
      samples$fixedEffects[,(batchStartEnd[ii]+1):batchStartEnd[ii+1]] <- samp.element$fixedEffects
      rm( samp.element)
      gc()
    }
  }
  
  #a data.frame containing prediciton points (no NAs).
  #add cell areas first
  if( !is.null( habitatArea)){
    tmp <- covars[[habitatArea]]
    names( tmp) <- "blah"  #will be changed in a bit anyway
    covars <- c( covars, tmp)#raster::addLayer( covars, tmp)  #this wastes memory a bit, temporarily (only really doing to rename things easily)
    covars <- covars[[names( covars) != habitatArea]]
  }
  else{
#    if( terra::is.lonlat( covars))
      covars <- c( covars, terra::cellSize( covars))
#    else{
#      tmp <- covars[[1]]
#      names( tmp) <- "tmpName"
#      terra::values( tmp) <- prod( terra::res( tmp))
#      covars <- c( covars, tmp)
#    }
  }
  names( covars)[terra::nlyr( covars)] <- "myCellAreas"
  covars[["myCellAreas"]] <- terra::mask( covars[["myCellAreas"]], covars[[1]])

  #get the coordinates of the prediction points
  predcoords <- terra::crds( covars, na.rm=FALSE)
  #extract the covariates
  
  #Get expanded data (model matrix) and corresponding formulae
  newInfo <- uniqueVarNames( obsList=list(), covarBrick=covars, distForm=object$distributionFormula, biasForm=object$biasFormula, arteForm=list(), habitatArea="myCellAreas", DCsurvID=attr( object, "DCobserverInfo"), coord.names=attr( res, "coord.names"), responseNames=object$responseNames, sampleAreaNames=NULL, stdCovs=object$control$standardiseCovariates, na.action=object$control$na.action)
  #putting it into a data frame
  covarData <- as.data.frame( terra::extract( newInfo$covarBrick, predcoords[,1:2]))#[,names( newInfo$covarBrick), drop=FALSE])
  myCellAreas <- as.matrix( terra::extract( covars, predcoords[,1:2])[,"myCellAreas", drop=FALSE])
    
  #cut down to just those areas without NAs.
  noNAid <- apply( covarData, 1, function(x) !any( is.na( x)))
  predcoords <- predcoords[noNAid,]
  covarData <- covarData[noNAid,, drop=FALSE]
  myCellAreas <- myCellAreas[noNAid,,drop=FALSE]

  #predictions start with cell area
  #always included no matter what components are included.
  eta <- matrix( NA, nrow=length( myCellAreas), ncol=S)
  eta[,] <- rep( log( as.numeric( myCellAreas)), times=S)
  
  #container for names of fixed effects
  fix.names <- object$mod$names.fixed

  #adding the intercepts, if any
  if( !is.null( intercept.terms))
    for( jj in intercept.terms)
      eta <- sweep( eta, MARGIN=2, STATS=samples$fixedEffects[fix.names == jj,], FUN="+")
  
  #predictions due to only fixed effects
  if( any( includeFixed!=FALSE)){
    #the model matrix
    myForm <- newInfo$distForm#object$distributionFormula
  
    X <- stats::model.matrix( myForm, data=covarData)
  
    #sorting the design matrix and the effects so that they match

    #the fixed effects involved in this term
    fix.subset <- which( fix.names %in% colnames( X))
    fix.names.ord <- order( fix.names[fix.subset])
    #ordering
    fixedSamps <- samples$fixedEffects[fix.subset[fix.names.ord],]
    X <- X[,order( colnames( X)),drop=FALSE]
    #zero-ing out effects other than requested.
    if( !is.logical( includeFixed)){
      colID <- unlist( lapply( includeFixed, function(xx) grep( xx, colnames( X))))
      X[,-colID] <- 0
    }
    #the addition to the linear predictor
    eta <- eta + X %*% fixedSamps
  }
  
  #adding in the random effects, if present and wanted
  if( object$control$addRandom & includeRandom==TRUE){
    #projector matrix( linking prediction points to mesh)
    A.prd <- as.matrix( INLA::inla.spde.make.A( object$mesh, loc=predcoords))
    #the addition to the linear predictor
#    Allow batch matrix multiplication for case where memory limits reached -- Dave's code
#    eta <- eta + A.prd %*% samples$fieldAtNodes
    if ((n.batches>1)&(!is.null( dim(eta))))
    {
      eta[,(batchStartEnd[1] + 1):batchStartEnd[2]] <- eta[,(batchStartEnd[1] + 1):batchStartEnd[2]] + as.matrix(A.prd %*% samples$fieldAtNodes[,(batchStartEnd[1] + 1):batchStartEnd[2]])
      for (ii in 2:n.batches) {
        eta[,(batchStartEnd[ii] + 1):batchStartEnd[ii + 1]] <- eta[,(batchStartEnd[ii] + 1):batchStartEnd[ii + 1]] + as.matrix(A.prd %*% samples$fieldAtNodes[,(batchStartEnd[ii] + 1):batchStartEnd[ii + 1]])         
      }  
    }
    else
    {
      eta <- eta + A.prd %*% samples$fieldAtNodes
    }
  }
  
  #adding in the bias terms, if wanted.
  if( includeBias==TRUE){
    #sampling bias formula
    bform <- newInfo$biasForm
    bform <- update( bform, "~.-1")  #belt and braces
    #sampling bias model.matrix
    bX <- stats::model.matrix( bform, data=covarData)
#    #fixing up intercept name
#    tmpID <- grepl( "(Intercept)", colnames( bX))
#    colnames( bX)[tmpID] <- "Intercept.PO"
#    colnames( bX)[!tmpID] <- paste0( "Intercept.PO:",colnames( bX)[!tmpID])
    bX <- bX[,order( colnames( bX))]
    #sorting samples and design matrix
    tmpID1 <- grep( "PO_", object$mod$names.fixed)
    bfixedNames <- object$mod$names.fixed[tmpID1]
    newSampsFixedBias <- samples$fixedEffects[tmpID1,][order( bfixedNames),]
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
  muRaster <- terra::rast( cbind( predcoords, mu.median), crs=terra::crs( covars), type='xyz')
  muRaster <- c(muRaster, terra::rast( cbind( predcoords, mu.lower), crs=terra::crs( covars), type='xyz'))
  muRaster <- c(muRaster, terra::rast( cbind( predcoords, mu.upper), crs=terra::crs( covars), type='xyz'))
  muRaster <- c(muRaster, terra::rast( cbind( predcoords, mu.mean), crs=terra::crs( covars), type='xyz'))
  muRaster <- c(muRaster, terra::rast( cbind( predcoords, mu.sd), crs=terra::crs( covars), type='xyz'))

  #sort out extent in case...
  muRaster <- terra::extend( muRaster, terra::ext( covars))  #just in case it is needed -- could be dropped throughout the creation of the raster.

  res <- list( field=muRaster, cell.samples=mu.all, fixedSamples=samples$fixedEffects, fixed.names=object$mod$names.fixed, predLocats=predcoords)
  
  return( res)
}

