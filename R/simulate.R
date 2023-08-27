
###############################################################################################
###############################################################################################
####	
####	Simulate data that is useful for checking isdm approach
####
####	Returns some data, within a list
####
####	Programmed by Scott in the first half of 2022
####
###############################################################################################
###############################################################################################

simulateData.isdm <- function( expected.pop.size=10000, expected.n.PO=300, n.PA=150, n.AA=50, n.DC=50,
                           coefs=list(dist=c(-1.5,0.5,-0.75), bias=c(-2,-0.75)), 
                           DC.pis=matrix( c(0.8,0.76, 0.9,0.85, 0.82,0.87), nrow=3, ncol=2, byrow=TRUE),
                           transect.size = 0.125, #a proportion of cell size.
                           rasterBoundary=NULL,
			   rasterCovars=NULL,
			   rasterBiasCovar=NULL,
                           control=list()){
			     
  #note that coefs$dist[1] will be largely ignored to match expected.n.PO  
  #set the control for the simulation
  control <- set.sim.control(control)
  #set the seed, if requested
  if( control$set.random.seed)
    set.seed( control$random.seed)
  
  if( is.null( rasterCovars)){
    #  Define survey area based on supplied raster
    if( is.null( rasterBoundary)){
      xSeq <- seq( from=0, to=10, length=control$raster.dim[1])
      ySeq <- seq( from=0, to=10, length=control$raster.dim[2])
      X <- expand.grid( x=xSeq, y=ySeq)
      my.scale <- sqrt( 100+100) / 15 #1/10 of observed
      if( is.na( control$range))
	control$range <- my.scale / 3 #to get quite variable random effects
    }
    #or survey area based on unit square
    else{
      X <- as.data.frame( terra::crds( rasterBoundary, na.rm=FALSE))
      xSeq <- sort( unique( X[,1]))
      ySeq <- sort( unique( X[,2]))
      my.scale <- sqrt( (utils::tail( xSeq,1) - utils::head( xSeq,1))^2 + (utils::tail( ySeq,1) - utils::head( ySeq,1))^2) / 15  #arbitrary
      if( is.na( control$range))
	control$range <- my.scale / 3
    }
    
#    simmy1 <- fftGPsim( x=xSeq, y=ySeq, sig2 = 1, rho = my.scale, nu = 5/2, nugget = 0.01)  #5/2 as the a big value -- most gaussian...
#    simmy2 <- fftGPsim( x=xSeq, y=ySeq, sig2 = 1, rho = my.scale, nu = 5/2, nugget = 0.01)  #5/2 as the a big value -- most gaussian...
       
    commonSim <- GPMaternSim(x=xSeq, y=ySeq, sig2 = 1, rho = my.scale, nu = 5/2, nugget = 0.01, n=2)  #5/2 as the a big value -- most gaussian...
    simmy1 <- commonSim[,3]
    simmy2 <- commonSim[,4]
    
    X <- cbind( X, as.numeric( t( simmy1)), as.numeric( t( simmy2)))
    X[,-(1:2)] <- apply( X[,-(1:2)], 2, scale)
    colnames( X) <- c("x","y","Altitude","Temperature")
    
    dataBrick <- c( terra::rast( X[,c(1,2,3)], type='xyz'), terra::rast( X[,c(1,2,4)], type='xyz'))
  }
  else{
#    simmy1 <- rasterCovars[[1]]
#    simmy2 <- rasterCovars[[2]]
    xSeq <- terra::xFromCol( rasterCovars, 1:terra::ncol(rasterCovars))
    ySeq <- terra::yFromRow( rasterCovars, 1:terra::nrow(rasterCovars))
    X <- as.data.frame( terra::crds( rasterCovars, na.rm=FALSE))
    rasterBoundary <- !is.na( rasterCovars[[1]])
    terra::values( rasterBoundary) <- ifelse( terra::values( rasterBoundary)==1, 1, NA)
    dataBrick <- c( rasterCovars[[1]], rasterCovars[[2]])
    if( is.na( control$range))
      control$range <- sqrt( (utils::tail( xSeq,1) - utils::head( xSeq,1))^2 + (utils::tail( ySeq,1) - utils::head( ySeq,1))^2) / (3*15)  #arbitrary
  }
  #random effect for the log-gauss process
  if( control$addRandom){
#   REff <- fftGPsim( x=xSeq, y=ySeq, sig2=control$sd^2, rho=control$range / 2, nu=3/2)  #3/2 is what is used to estimate (in RISDM) 14/8/23. Previously nu=1...?
#    REff <- fftGPsim( x=xSeq, y=ySeq, sig2=control$sd^2, rho=sqrt( 1.5 * 8) * control$range, nu=3/2)  
   REff <- GPMaternSim( x=xSeq, y=ySeq, sig2=control$sd^2, rho=control$range / 2, nu=3/2)[,3]  #3/2 is what is used to estimate (in RISDM) 14/8/23. Previously nu=1...?
   REff <- as.numeric( t( REff))
  }
  else
    REff <- rep(0,nrow( X))
  
  #scale covariates
  dataBrick <- terra::scale( dataBrick, center=TRUE, scale=FALSE)

  #data brick for the covariates and random effect
  dataBrick <- c( dataBrick, terra::rast( cbind( X[,1:2], REff), type='xyz', crs=terra::crs( dataBrick)))
  if( !is.null( rasterBoundary)){
    terra::crs( dataBrick) <- terra::crs( rasterBoundary)
    dataBrick <- terra::mask( dataBrick, rasterBoundary)
  }
  #a data frame model.matrix( ok a df, but still)
  X <- terra::as.data.frame( dataBrick, xy=TRUE)
  #linear predictor
  LinPred <- coefs$dist[1] + 
              dataBrick[[1]]*coefs$dist[2] + 
              dataBrick[[2]]*coefs$dist[3] + 
              dataBrick$REff
  
  #Intensity for log-gauss process
  Intensity <- exp( LinPred)
  #the total expectation and observed number from unbiassed (actual distribution)
  tmp.Lambda <- sum( terra::values( Intensity), na.rm=TRUE)
  Intensity <- Intensity / tmp.Lambda  #Intensity should sum to 1
  Intensity <- Intensity * expected.pop.size  #now sum to expected.pop.size
  
  #add to databrick
  dataBrick <- c(dataBrick,LinPred,Intensity)
  
  names( dataBrick)[-(1:3)] <- c("logIntensity","Intensity") 

  #get the boundary again, if needed
  if( is.null( rasterBoundary))
    tmpBoundary <- !is.na( dataBrick[[1]])
  else
    tmpBoundary <- rasterBoundary

  ####  Create the PA data
  #first create the design -- an equal probability spat-bal design
  PAdata <- MBHdesign::quasiSamp.raster( n=n.PA, tmpBoundary)
  
  #now simulate the PA response at those sites (based on point process)
  PAdata$PA <- stats::rbinom( n=n.PA, size=1, prob=pmax( 0, pmin( 1-exp( -terra::values( dataBrick$Intensity)[PAdata$ID]*transect.size),1)))
  
  #add the transect area
  PAdata$transectArea <- transect.size * prod( terra::res( dataBrick))
  
  ####  Create the count data (AA)
  AAData <- MBHdesign::quasiSamp.raster( n=n.AA, tmpBoundary)
  
  #now simulate the AA response at those sites (based on point process)
  AAData$AA <- stats::rpois( n=n.AA, lambda=terra::values( dataBrick$Intensity)[AAData$ID]*transect.size)
  
  #add the transect area
  AAData$transectArea <- transect.size * prod( terra::res( dataBrick))

  ####  Create the double count data (DC)
  DCData <- MBHdesign::quasiSamp.raster( n=n.DC, tmpBoundary)
  
  #now simulate the DC response at those sites (based on point process)
  #total number of animals available along a transect
  tmpTotCount <- stats::rpois( n=n.DC, lambda=terra::values( dataBrick$Intensity)[DCData$ID]*transect.size)
  #container for the DC data
  tmptmptmp <- as.data.frame( matrix( NA, ncol=3, nrow=n.DC))
  colnames( tmptmptmp) <- c( "Obs1","Obs2","Both")
  #the survey ID
  tmpID <- rep( 1:nrow( DC.pis), each=n.DC %/% nrow( DC.pis))
  tmpID <- c( tmpID, rep( nrow( DC.pis), n.DC %% nrow( DC.pis))) #make the last survey biggest, if needed
  
  #thin into different categories for the observers
  cat.pis <- cbind( Obs1=DC.pis[,1]*(1-DC.pis[,2]), Obs2=(1-DC.pis[,1])*DC.pis[,2], Both=DC.pis[,1]*DC.pis[,2], Neither=(1-DC.pis[,1])*(1-DC.pis[,2]))
  observs <- matrix( NA, nrow=n.DC, ncol=nrow( DC.pis))
  colnames( observs) <- colnames( tmptmptmp)[1:3]
  for( ii in 1:n.DC)
    observs[ii,] <- stats::rmultinom( n=1, size=tmpTotCount[ii], prob=cat.pis[tmpID[ii],])[1:3]  #don't want the last one.
  
  DCData <- cbind( DCData, observs)
  
  #survey labels
  DCData$Survey <- paste0("DCsurvey_",tmpID)
  
  #add the transect area
  DCData$transectArea <- transect.size * prod( terra::res( dataBrick))
  
  ####  Now deal with the biassed observations
  if( is.null( rasterBiasCovar)){
    #distance to city, or a proxy (something quite like it)...
    tmptmptmp <- dataBrick[[1]]
    myCRDS <- terra::crds( tmptmptmp, na.rm=FALSE)
    midPT <- c( mean( myCRDS[,1], na.rm=TRUE), mean( myCRDS[,2], na.rm=TRUE))
    ref <- terra::cellFromXY( dataBrick, matrix( midPT, nrow=1))  
    myCRDS[,1] <- (myCRDS[,1] - myCRDS[ref,1])^2
    myCRDS[,2] <- (myCRDS[,2] - myCRDS[ref,2])^2
    myCRDS <- sqrt( rowSums( myCRDS))
    terra::values( tmptmptmp) <- myCRDS
  #  terra::values( tmptmptmp) <- max( myCRDS, na.rm=TRUE) - myCRDS
    tmptmptmp <- terra::mask( tmptmptmp, dataBrick[[1]])
    terra::values( tmptmptmp) <- scale( terra::values( tmptmptmp, na.rm=FALSE))
    dataBrick <- c( dataBrick, tmptmptmp)
    names( dataBrick)[length( names( dataBrick))] <- "dist2City"
  }
  else
    dataBrick <- c( dataBrick, terra::scale( rasterBiasCovar))
  
  #probability of being observed (thinning process) -- cloglog link for prob
  dataBrick <- c( dataBrick, 1-exp(- exp( coefs$bias[1]+coefs$bias[2]*dataBrick[[terra::nlyr( dataBrick)]])))
  names( dataBrick)[length( names( dataBrick))] <- "obsProb"
  #The biassed intensity after accounting for observation bias.
  dataBrick <- c( dataBrick, biasIntensity=dataBrick$Intensity * dataBrick$obsProb)
  names( dataBrick)[length( names( dataBrick))] <- "biasIntensity"
  
  if( !is.null( rasterBoundary))
    dataBrick <- terra::mask( dataBrick, rasterBoundary)
  
  ####  Take the sample of po data
  N <- stats::rpois( n=1, lambda=expected.n.PO)
  tmp.Lambda <- sum( terra::values( dataBrick$biasIntensity), na.rm=TRUE)
  probs <- dataBrick$biasIntensity * expected.n.PO / tmp.Lambda  #probably not needed as sample() internally rescales...
  
  NAid <- is.na( terra::values( probs))
  tmpProbs <- terra::values( probs)
  tmpProbs[NAid] <- 0
  
  sampleID <- sample.int( n=length( tmpProbs), size=N, prob=tmpProbs, replace=TRUE)
  presences <- terra::crds( dataBrick, na.rm=FALSE)[sampleID,]

  if( control$doPlot){
    graphics::par( mfrow=c(4,2))
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
    points( jitter( presences), pch=20, cex=0.5)
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
