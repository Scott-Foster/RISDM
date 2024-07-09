
###############################################################################################
###############################################################################################
####	
####	Predict distribution and/or bias and/or artefact from an isdm
####  	When the prediction space is large or the number of samples is large.
####
####	Returns a list of predictions
####
####	Programmed by Scott in Nov/Dec 2022
####	Refactored by Scott July 2024
####
###############################################################################################
###############################################################################################

#Function to get prediction from a fitted INLA model.
predict.isdm <- function( object, covars, habitatArea=NULL, S=500, intercept.terms=NULL, n.threads=NULL, n.batches=1, includeRandom=TRUE, includeFixed=TRUE, includeBias=FALSE, type="intensity", confidence.level=0.95, ...){
  
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

  #set up the batch details
  batchStartEnd <- round( seq( from=0, to=S, by=S/n.batches),0)
  batchSize <- diff( batchStartEnd)

  #### function for drawing posterior sampls in parallel.
  funny <- function( ii){
    tmp.samples <- draw.posterior.samps(object$mod, B=batchSize[ii], what="effects", field="isdm.spat.XXX", n.threads=n.threads)
    return( tmp.samples)
  }

  #initiate (get dimensions etc)
  samp.element <- funny(1)

  #container for the posterior samples
  samples <- list()
  # first batch: the random field
  if( object$control$addRandom){
    samples$fieldAtNodes <- matrix( NA, nrow=nrow( samp.element$fieldAtNodes), ncol=S)
    samples$fieldAtNodes[,(batchStartEnd[1]+1):batchStartEnd[2]] <- samp.element$fieldAtNodes    
  }
  # first batch: the posterior fixed effects
  samples$fixedEffects <- matrix( NA, nrow=nrow( samp.element$fixedEffects), ncol=S)
  samples$fixedEffects[,(batchStartEnd[1]+1):batchStartEnd[2]] <- samp.element$fixedEffects
  
  #if there are multiple batches  (don't want too many batches given for loop)
  if( n.batches > 1){
    for( ii in 2:n.batches){
      samp.element <- funny(ii)
      if( object$control$addRandom)
	samples$fieldAtNodes[,(batchStartEnd[ii]+1):batchStartEnd[ii+1]] <- samp.element$fieldAtNodes
      samples$fixedEffects[,(batchStartEnd[ii]+1):batchStartEnd[ii+1]] <- samp.element$fixedEffects
      rm( samp.element)
      gc()
    }
  }

  # move from parameter samples to predictions
  #
  if( is.null( habitatArea)){
    covars <- c( covars, terra::cellSize( covars))  #could possibly use units, somehow...
    habitatArea <- "habArea"
    names( covars)[terra::nlyr( covars)] <- habitatArea
    covars[[habitatArea]] <- terra::mask( covars[[habitatArea]], covars[[1]])
  }

  #get the coordinates of the prediction points
  predcoords <- terra::crds( covars, na.rm=FALSE)
  #extract the covariates
  
  #Get expanded data (model matrix) and corresponding formulae
  newInfo <- uniqueVarNames( obsList=list(), covarBrick=covars, distForm=object$distributionFormula, biasForm=object$biasFormula, arteForm=list(), habitatArea=habitatArea, DCsurvID=attr( object, "DCobserverInfo"), coord.names=attr( res, "coord.names"), responseNames=object$responseNames, sampleAreaNames=NULL, stdCovs=object$control$standardiseCovariates, na.action=object$control$na.action)
  #putting it into a data frame
  covarData <- as.data.frame( terra::extract( newInfo$covarBrick, predcoords[,1:2]))
#  myCellAreas <- as.matrix( terra::extract( covars[[habitatArea]], predcoords[,1:2]))
    
  #cut down to just those areas without NAs.
  noNAid <- apply( covarData, 1, function(x) !any( is.na( x)))
  predcoords <- predcoords[noNAid,,drop=FALSE]
  covarData <- covarData[noNAid,,drop=FALSE]
#  myCellAreas <- myCellAreas[noNAid,,drop=FALSE]

  #predictions start with cell area
  #always included no matter what components are included.
  eta <- matrix( NA, nrow=length( covarData[,habitatArea]), ncol=S)
  eta[,] <- rep( log( as.numeric( covarData[,habitatArea])), times=S)
  
  #container for names of fixed effects
  fix.names <- object$mod$names.fixed

  #adding the intercepts, if any
  if( !is.null( intercept.terms))
    for( jj in intercept.terms)
      eta <- sweep( eta, MARGIN=2, STATS=samples$fixedEffects[fix.names == jj,,drop=FALSE], FUN="+")
  
  #predictions due to only fixed effects
  if( any( includeFixed!=FALSE)){
    #the model matrix
    myForm <- newInfo$distForm
  
    X <- covarData[,attr( terms( myForm), "term.labels"),drop=FALSE]#stats::model.matrix( myForm, data=covarData)
  
    #sorting the design matrix and the effects so that they match

    #the fixed effects involved in this term
    fix.subset <- which( fix.names %in% colnames( X))
    fix.names.ord <- order( fix.names[fix.subset])
    #ordering
    fixedSamps <- samples$fixedEffects[fix.subset[fix.names.ord],,drop=FALSE]
    X <- X[,order( colnames( X)),drop=FALSE]
    #zero-ing out effects other than requested. Generally all will be requested
    if( !is.logical( includeFixed)){
      colID <- unlist( lapply( includeFixed, function(xx) grep( xx, colnames( X))))
      X[,-colID] <- 0
    }
    #the addition to the linear predictor
    eta <- eta + as.matrix( X) %*% fixedSamps
  }
  
  #adding in the random effects, if present and wanted
  if( object$control$addRandom & includeRandom==TRUE){
    #projector matrix( linking prediction points to mesh)
    A.prd <- INLA::inla.spde.make.A( object$mesh, loc=predcoords)
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
      eta <- eta + as.matrix( A.prd %*% samples$fieldAtNodes)
    }
  }
  
  #adding in the bias terms, if wanted.
  if( includeBias==TRUE){
    #sampling bias formula
    bform <- newInfo$biasForm
    bform <- update( bform, "~.-1")  #belt and braces
    #sampling bias model.matrix
    #bX <- stats::model.matrix( bform, data=covarData)
    bX <- covarData[,attr( terms( bform), "term.labels"), drop=FALSE]
    #ordering shouldn't be needed, but it won't hurt!?
    bX <- bX[,order( colnames( bX)),drop=FALSE]
    #sorting samples and design matrix
    tmpID1 <- grep( "PO_", object$mod$names.fixed)
    bfixedNames <- object$mod$names.fixed[tmpID1]
    #ordering shouldn't be needed, but it won't hurt!?
    newSampsFixedBias <- samples$fixedEffects[tmpID1,][order( bfixedNames),,drop=FALSE]
    #the addition to the linear predictor
    eta <- eta + as.matrix( bX) %*% newSampsFixedBias
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
  limitty <- c( (1-confidence.level)/2, 1-(1-confidence.level)/2)
  lambda.median <- apply( mu.all, 1, stats::quantile, probs=0.5, na.rm=TRUE)
  lambda.lower <- apply( mu.all, 1, stats::quantile, probs=limitty[1], na.rm=TRUE)
  lambda.upper <- apply( mu.all, 1, stats::quantile, probs=limitty[2], na.rm=TRUE)
  lambda.mean <- rowMeans( mu.all)
  lambda.sd <- apply( mu.all, 1, stats::sd)
  
  #raster format
  lambdaRaster <- terra::rast( cbind( predcoords, lambda.median), crs=terra::crs( covars), type='xyz')
  lambdaRaster <- c(lambdaRaster, terra::rast( cbind( predcoords, lambda.lower), crs=terra::crs( covars), type='xyz'))
  lambdaRaster <- c(lambdaRaster, terra::rast( cbind( predcoords, lambda.upper), crs=terra::crs( covars), type='xyz'))
  lambdaRaster <- c(lambdaRaster, terra::rast( cbind( predcoords, lambda.mean), crs=terra::crs( covars), type='xyz'))
  lambdaRaster <- c(lambdaRaster, terra::rast( cbind( predcoords, lambda.sd), crs=terra::crs( covars), type='xyz'))
  names( lambdaRaster) <- c("Median","Lower","Upper","Mean","SD")

  #sort out extent in case...
  lambdaRaster <- terra::extend( lambdaRaster, terra::ext( covars))  #just in case it is needed -- could be dropped throughout the creation of the raster.

  res <- list( field=lambdaRaster, cell.samples=mu.all, fixedSamples=samples$fixedEffects, fixed.names=object$mod$names.fixed, predLocats=predcoords, confidence.limits=limitty)
  
  return( res)
}

