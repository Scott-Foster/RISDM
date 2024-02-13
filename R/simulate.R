
###############################################################################################
###############################################################################################
####	
####	Simulate data that is useful for checking isdm approach
####
####	Returns some data, within a list
####
####	Programmed by Scott in the first half of 2022 (rejigged again in late 2023)
####
###############################################################################################
###############################################################################################

simulateData.isdm <- function( expected.pop.size=10000, expected.n.PO=300, n.PA=150, n.AA=50, n.DC=50,
                           coefs=list(dist=c(NA,0.5,-0.75), bias=c(-2,-0.75)), 
#                           DC.pis=matrix( c(0.8,0.76, 0.9,0.85, 0.82,0.87), nrow=3, ncol=2, byrow=TRUE),
                           DC.pis=c(0.8, 0.9, 0.75),
                           transect.size = 0.125, #a proportion of cell size.
                           rasterBoundary=NULL,
			   rasterCovars=NULL,
			   rasterBiasCovar=NULL,
                           control=list()){
			     
  if( sum( is.null( rasterCovars), is.null( rasterBiasCovar)) == 1)  #one but not both layers specified
    stop( "You must specify either both rasterCovars and rasterBiasCovar, or neither.")

#  if( is.null( rasterCovars)){
#    thiscall <- match.call( expand.dots=TRUE)
#    thiscall[[1]] <- as.name("simulateData.isdm")
#    tmp <- eval.parent( thiscall)
#    return( tmp)
#  }

  #note that coefs$dist[1] will be ignored to match expected.n.PO  
  #set the control for the simulation
  control <- set.sim.control(control)
  #set the seed, if requested
  if( control$set.random.seed)
    set.seed( control$random.seed)

  #mask and standardise the covariates and bias layers if supplied
  #assumes that if rasterCovars is present then rasterBiasCovar will also be present
  if( !is.null( rasterCovars)){
    effortLayerName <- names( rasterBiasCovar)[1]
    tmp <- apply( cbind( terra::values( rasterCovars, na.rm=FALSE), terra::values( rasterBiasCovar)), 1, function(xx) any( is.na( xx)))
    tmp1 <- apply( cbind( terra::values( rasterCovars, na.rm=FALSE), terra::values( rasterBiasCovar)), 1, function(xx) all( is.na( xx)))
    if( !( all( tmp-tmp1)==0))
      warning( "Some, but not all, NAs in one raster layer are absent in the other. Altering data to maintain NA pattern.")
    rasterBoundary <- rasterCovars[[1]]
    terra::values( rasterBoundary) <- as.numeric( !is.na( terra::values( rasterCovars[[1]])))
    terra::values( rasterBoundary)[tmp,] <- NA

    rasterCovars <- terra::scale( terra::mask( rasterCovars, rasterBoundary))
    rasterBiasCovar <- terra::scale( terra::mask( rasterBiasCovar, rasterBoundary))
    names( rasterBiasCovar) <- effortLayerName#"effortLayer"
    
    #getting raster parameters for simulation
    X <- as.data.frame( terra::crds( rasterBoundary, na.rm=FALSE))
    xSeq <- sort( unique( X[,1]))
    ySeq <- sort( unique( X[,2]))
    if( is.na( control$range))  #arbitrary!
      control$range <- 5*max( terra::res(rasterBoundary))
#    rm( X)
  }
  else{
    #  Define survey area based on square
    if( is.null( rasterBoundary)){
      xSeq <- seq( from=0, to=10, length=control$raster.dim[1])
      ySeq <- seq( from=0, to=10, length=control$raster.dim[2])
      X <- expand.grid( xSeq, ySeq)
      
      my.scale <- sqrt( 100+100) / 15 #1/15 of max observed
      if( is.na( control$range))
	control$range <- my.scale / 3 #to get quite variable random effects
    }
    #or survey area based on supplied boundary raster
    else{
      X <- as.data.frame( terra::crds( rasterBoundary, na.rm=FALSE))
      X <- X[order( X[,1],X[,2]),]
      xSeq <- sort( unique( X[,1]))
      ySeq <- sort( unique( X[,2]))
      if( is.na( control$range))  #arbitrary!
	control$range <- 5*max( terra::res(rasterBoundary))
#      rm( X)
    }
    simmy1 <- fftGPsim2( x=xSeq, y=ySeq, sig2 = 1, rho = control$range, nu = 1/2, nugget = 0.01)
    simmy2 <- fftGPsim2( x=xSeq, y=ySeq, sig2 = 1, rho = control$range, nu = 1/2, nugget = 0.01)
   
    X <- cbind( X, as.numeric( t( simmy1)), as.numeric( t( simmy2)))#, as.numeric( t( simmy3)))
    X[,-(1:2)] <- apply( X[,-(1:2)], 2, scale, center=TRUE, scale=TRUE)
    colnames( X) <- c("x","y","Altitude","Temperature")#,"bias")
    rasterCovars <- c( terra::rast( X[,c(1,2,3)], type='xyz', crs="epsg:3857"), terra::rast( X[,c(1,2,4)], type='xyz', crs="epsg:3857"))
    rm( X)   
    
    #  Now deal with the biassing layer
    #  Assumed bias layer is NULL as covariates are null
    #distance to city, or a proxy (something quite like it)...
    rasterBiasCovar <- rasterCovars[[1]]
    myCRDS <- terra::crds( rasterBiasCovar, na.rm=FALSE)
    midPT <- c( mean( myCRDS[,1], na.rm=TRUE), mean( myCRDS[,2], na.rm=TRUE))
    ref <- terra::cellFromXY( rasterBiasCovar, matrix( midPT, nrow=1))  
    myCRDS[,1] <- (myCRDS[,1] - myCRDS[ref,1])^2
    myCRDS[,2] <- (myCRDS[,2] - myCRDS[ref,2])^2
    myCRDS <- sqrt( rowSums( myCRDS))
    terra::values( rasterBiasCovar) <- myCRDS
    rasterBiasCovar <- terra::mask( rasterBiasCovar, rasterCovars[[1]])
    terra::values( rasterBiasCovar) <- scale( terra::values( rasterBiasCovar, na.rm=FALSE))
    effortLayerName <- names( rasterBiasCovar) <- "effortLayer"
  }

  dataBrick <- c( rasterCovars[[1]], rasterCovars[[2]], rasterBiasCovar)
  
  #spatial random effect
  if( is.na( control$range))
     control$range <- 5*max( terra::res( dataBrick))
  
  #random effect for the log-gauss process
  if( control$addRandom){
   REff <- fftGPsim2( x=xSeq, y=ySeq, sig2=control$sd^2, rho=control$range, nu=1/2)
   REff <- as.numeric( REff)
  }
  else
    REff <- rep(0, terra::ncell( rasterCovars))
  
  #data brick for the covariates and random effect
  tmp <- dataBrick[[1]]
  names( tmp) <- "REff"
  terra::values( tmp) <- REff
  dataBrick <- c( dataBrick, tmp)#terra::rast( cbind( terra::crds(rasterCovars), REff), type='xyz', crs=terra::crs( dataBrick)))
  if( !is.null( rasterBoundary)){
    terra::crs( dataBrick) <- terra::crs( rasterBoundary)
    dataBrick <- terra::mask( dataBrick, rasterBoundary)
  }
  #a data frame model.matrix( ok a df, but still)
  X <- terra::as.data.frame( dataBrick, xy=TRUE)
  #cellsizes
#  if( terra::crs( dataBrick) != "")
    myCellSize <- terra::cellSize( dataBrick)
#  else
#    myCellSize <- rep( prod( terra::res( dataBrick)), terra::ncell( dataBrick))
  #linear predictor
  LinPred <- #coefs$dist[1] + 
              dataBrick[[1]]*coefs$dist[2] + 
              dataBrick[[2]]*coefs$dist[3] + 
              dataBrick$REff +
	      log(myCellSize)
  
  #Intensity for log-gauss process
  Intensity <- exp( LinPred)
  #the total expectation and observed number from unbiassed (actual distribution)
  tmp.Lambda <- sum( terra::values( Intensity), na.rm=TRUE)
  Intensity <- Intensity / tmp.Lambda  #Intensity should sum to 1
  Intensity <- Intensity * expected.pop.size  #per unit area, not per cell nor transect nor...
  
  #scale lin pred
  LinPred <- log( Intensity)
  
  #add to databrick
  dataBrick <- c(dataBrick,LinPred,Intensity)
  
  names( dataBrick)[-(1:4)] <- c("logIntensity","Intensity") 

  #get the boundary again, if needed
  if( is.null( rasterBoundary))
    tmpBoundary <- !is.na( dataBrick[[1]])
  else
    tmpBoundary <- rasterBoundary

  ####  Create the PA data
  #first create the design -- an equal probability spat-bal design
  if( control$useMBHdesign)#nchar( system.file( package="MBHdesign")) > 0)
    PAdata <- MBHdesign::quasiSamp.raster( n=n.PA, tmpBoundary)
  else{
    tmp <- terra::as.data.frame( tmpBoundary, xy=TRUE, cells=TRUE)
    tmp <- tmp[sample.int(n=nrow(tmp), size=n.PA, replace=FALSE),]
    PAdata <- data.frame( x=tmp$x, y=tmp$y, inclusion.probabilities=NA, ID=tmp$cell)
  }
  #add the transect area
#  if( terra::crs( dataBrick)!="")
      mySizes <- terra::values( terra::cellSize( dataBrick))[PAdata$ID]
#  else
#    mySizes <- rep( prod( terra::res( dataBrick)), length( PAdata$ID))
  PAdata$transectArea <- transect.size * mySizes

  tmpIntensity <- ( terra::values( dataBrick$Intensity)[PAdata$ID,] / terra::values( terra::cellSize( dataBrick))[PAdata$ID,]) * PAdata$transectArea
  #now simulate the PA response at those sites (based on point process)
  PAdata$PA <- stats::rbinom( n=n.PA, size=1, prob=pmax( 0, pmin( 1-exp( -tmpIntensity),1)))
  
  ####  Create the count data (AA)
  if( control$useMBHdesign)#nchar( system.file( package="MBHdesign")) > 0)
    AAData <- MBHdesign::quasiSamp.raster( n=n.AA, tmpBoundary)
  else{
    tmp <- terra::as.data.frame( tmpBoundary, xy=TRUE, cells=TRUE)
    tmp <- tmp[sample.int(n=nrow(tmp), size=n.AA, replace=FALSE),]
    AAData <- data.frame( x=tmp$x, y=tmp$y, inclusion.probabilities=NA, ID=tmp$cell)
  }
  #add the transect area
  AAData$transectArea <- transect.size * terra::values( terra::cellSize( dataBrick))[AAData$ID]#prod( terra::res( dataBrick))
  #intensity over transect area
  tmpIntensity <- ( terra::values( dataBrick$Intensity)[AAData$ID,] / terra::values( terra::cellSize( dataBrick))[AAData$ID,]) * AAData$transectArea
  
  #now simulate the AA response at those sites (based on point process)
  AAData$AA <- stats::rpois( n=n.AA, lambda=tmpIntensity)
  
  ####  Create the double count data (DC)
  if( control$useMBHdesign)#nchar( system.file( package="MBHdesign")) > 0)
    DCData <- MBHdesign::quasiSamp.raster( n=n.DC, tmpBoundary)
  else{
    tmp <- terra::as.data.frame( tmpBoundary, xy=TRUE, cells=TRUE)
    tmp <- tmp[sample.int(n=nrow(tmp), size=n.DC, replace=FALSE),]
    DCData <- data.frame( x=tmp$x, y=tmp$y, inclusion.probabilities=NA, ID=tmp$cell)
  }
  #add the transect area
  DCData$transectArea <- transect.size * terra::values( terra::cellSize( dataBrick))[DCData$ID]#prod( terra::res( dataBrick))
  #intensity over transect area
  tmpIntensity <- ( terra::values( dataBrick$Intensity)[DCData$ID,] / terra::values( terra::cellSize( dataBrick))[DCData$ID,]) * DCData$transectArea

  #now simulate the DC response at those sites (based on point process)
  #total number of animals available along a transect
  tmpTotCount <- stats::rpois( n=n.DC, lambda=tmpIntensity)
  #container for the DC data
  tmptmptmp <- as.data.frame( matrix( NA, ncol=3, nrow=n.DC))
  colnames( tmptmptmp) <- c( "Obs1","Obs2","Both")
  #the survey ID
#  tmpID <- rep( 1:nrow( DC.pis), each=n.DC %/% nrow( DC.pis))
  tmpID <- rep( 1:length( DC.pis), each=n.DC %/% length( DC.pis))
#  tmpID <- c( tmpID, rep( nrow( DC.pis), n.DC %% nrow( DC.pis))) #make the last survey biggest, if needed
  tmpID <- c( tmpID, rep( length( DC.pis), n.DC %% length( DC.pis))) #make the last survey biggest, if needed
  
  #thin into different categories for the observers
#  cat.pis <- cbind( Obs1=DC.pis[,1]*(1-DC.pis[,2]), Obs2=(1-DC.pis[,1])*DC.pis[,2], Both=DC.pis[,1]*DC.pis[,2], Neither=(1-DC.pis[,1])*(1-DC.pis[,2]))
  cat.pis <- cbind( Obs1=DC.pis*(1-DC.pis), Obs2=(1-DC.pis)*DC.pis, Both=DC.pis*DC.pis, Neither=(1-DC.pis)*(1-DC.pis))
#  observs <- matrix( NA, nrow=n.DC, ncol=nrow( DC.pis))
  observs <- matrix( NA, nrow=n.DC, ncol=length( DC.pis))
  colnames( observs) <- colnames( tmptmptmp)[1:3]
  for( ii in 1:n.DC)
    observs[ii,] <- stats::rmultinom( n=1, size=tmpTotCount[ii], prob=cat.pis[tmpID[ii],])[1:3]  #don't want the last one.
  
  DCData <- cbind( DCData, observs)
  
  #survey labels
  DCData$Survey <- paste0("DCsurvey_",tmpID)
    
  ####	PO data (as per the model fitted)
  #linear predictor
  LinPred <- #coefs$dist[1] + 
              dataBrick[[1]]*coefs$dist[2] + 
              dataBrick[[2]]*coefs$dist[3] + 
              dataBrick$REff +
	      dataBrick[[effortLayerName]]*coefs$bias[2] + 
	      log(terra::cellSize( dataBrick))
  
  #Intensity for log-gauss process
  Intensity <- exp( LinPred)
  #the total expectation and observed number from unbiassed (actual distribution)
  tmp.Lambda <- sum( terra::values( Intensity), na.rm=TRUE)
  Intensity <- Intensity / tmp.Lambda  #Intensity should sum to 1
  if( !control$exact.n.PO)
    N <- stats::rpois( n=1, lambda=expected.n.PO)
  else
    N <- expected.n.PO
  Intensity <- Intensity * N  #per unit area, not per cell nor transect nor...
  
  #scale lin pred
  LinPred <- log( Intensity)
  
  #add to databrick
  dataBrick <- c(dataBrick,LinPred,Intensity)
  
  names( dataBrick)[-(1:6)] <- c("biasLogIntensity","biasIntensity") 

  if( !is.null( rasterBoundary))
    dataBrick <- terra::mask( dataBrick, rasterBoundary)
  
  ####  Take the sample of po data
  
#  tmp.Lambda <- sum( terra::values( dataBrick$biasIntensity), na.rm=TRUE)
#  probs <- dataBrick$biasIntensity * expected.n.PO / tmp.Lambda  #probably not needed as sample() internally rescales...
  
#  NAid <- is.na( terra::values( probs))
  tmpProbs <- terra::values( dataBrick$biasIntensity) / N  #so that the probs sum to 1 (will be done in sample() anyway, but making it clear here)
  tmpProbs[is.na( tmpProbs)] <- 0
  
  sampleID <- sample.int( n=length( tmpProbs), size=N, prob=tmpProbs, replace=TRUE)
  presences <- as.matrix( terra::crds( dataBrick, na.rm=FALSE)[sampleID,])

  if( control$doPlot){
    oldPar <- graphics::par( mfrow=c(4,2))
    #PA
    terra::plot( dataBrick$Intensity, main="Intensity")
    terra::plot( !is.na( dataBrick[[1]]), col=grey(c(1,0.95)), legend=FALSE, main="PA survey")
#    terra::plot( PAdata[,1:2], type='n', main="PA survey", xlim=range( xSeq), ylim=range( ySeq), asp=1, add=TRUE)
    points( jitter( PAdata[PAdata$PA==0,1]), jitter( PAdata[PAdata$PA==0,2]), cex=0.5, pch=20, col=grDevices::grey(0.7))
    points( jitter( PAdata[PAdata$PA==1,1]), jitter( PAdata[PAdata$PA==1,2]), col='red', pch=20, cex=0.5)
    legend( x="bottomleft", pch=c(20, 20), col=c(grDevices::grey(0.7),'red'), legend=c("Absence", "Presence"))
    #AA
    terra::plot( dataBrick$Intensity, main="Intensity")
    terra::plot( !is.na( dataBrick[[1]]), col=grey(c(1,0.95)), legend=FALSE, main="AA survey")
#    terra::plot( AAData[,1:2], type='n', main="AA survey", xlim=range( xSeq), ylim=range( ySeq), asp=1)
    for( ii in 0:max( AAData$AA))
      points( AAData[AAData$AA==ii,1:2], col=ii+1, pch=20, cex=1.5)
    legend( x="bottomleft", pch=20, col=1:(max( AAData$AA)+1), legend=0:max( AAData$AA))
    #DC
    terra::plot( dataBrick$Intensity, main="Intensity")
    terra::plot( !is.na( dataBrick[[1]]), col=grey(c(1,0.95)), legend=FALSE, main="DC survey")
#    terra::plot( DCData[,1:2], type='n', main="DC survey", xlim=range( xSeq), ylim=range( ySeq), asp=1)
    tmpCount <- rowSums( DCData[,c("Obs1","Obs2","Both")])
    for( ii in 0:max( tmpCount))
      points( DCData[tmpCount==ii,1:2], col=ii+1, pch=20, cex=1.5)
    legend( x="bottomleft", pch=20, col=1:(max( tmpCount)+1), legend=0:max( tmpCount))
    
    #PO  
    terra::plot( dataBrick$biasIntensity, main="Biassed Intensity")
    terra::plot( !is.na( dataBrick[[1]]), col=grey(c(1,0.95)), legend=FALSE, main="PO observations")
    points( jitter( as.matrix( presences)), pch=20, cex=0.5)
    
    #reset graphics back to user's
    graphics::par( oldPar)
  }
  
#  res <- list( PO=sf::st_as_sf( as.data.frame( cbind( 1:nrow( presences), presences)), coords=c("x","y")),#sp::SpatialPoints( presences), 
#               PA=sf::st_as_sf( PAdata, coords=c("x","y")),#sp::SpatialPointsDataFrame(coords = PAdata[,1:2], data = PAdata[,-(1:2)]), 
#               AA=sf::st_as_sf( AAData, coords=c("x","y")),#sp::SpatialPointsDataFrame( coords=AAData[,1:2], data=AAData[,-(1:2)]), 
#               DC=sf::st_as_sf( DCData, coords=c("x","y")),#sp::SpatialPointsDataFrame( coords=DCData[,1:2], data=DCData[,-(1:2)]),
#               covarBrick=dataBrick,
#               expected.pop.size=expected.pop.size)
	       
  res <- list( PO=presences,
               PA=PAdata,
               AA=AAData,
               DC=DCData,
               covarBrick=dataBrick,
               expected.pop.size=expected.pop.size)
  
  class( res) <- "simISDMdata"
  
  return(res)
}
