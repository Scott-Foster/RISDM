
###############################################################################################
###############################################################################################
####	
####	Simulate data that is useful for checking isdm approach
####
####	Returns some data, within a list
####
####	Programmed by Scott in the first half of 2022 (rejigged again in late 2023)
####	Refactored in mid 2024
####
###############################################################################################
###############################################################################################

simulateData.isdm <- function( pop.size=10000,
				      distForm=~-1+var1, biasForm=~1,
				      Intercept=NULL, distCoefs=NULL, biasCoefs=NULL, DC.pi=NULL,
				      n.PO=300, n.PA=150, n.AA=50, n.DC=50,
				      rasterBoundary=NULL, covarBrick=NULL,
				      transect.size=0.125,
				      Intensity=NULL,
				      control=list()){

  #set the control for the simulation
  control <- set.sim.control(control)
  #set the seed, if requested
  if( control$set.random.seed)
    set.seed( control$random.seed)
    
  if( !is.null( Intensity)){
    if( is.null( covarBrick))
      covarBrick <- Intensity
    else
      covarBrick <- c( covarBrick, Intensity)
    distForm <- stats::reformulate( termlabels=names( Intensity))
    Intercept <- NULL
  }
				      
  #variable names
  my.allVars <- unique( c( all.vars( distForm), all.vars( biasForm)))
#  if( !is.null( habitatArea))
#    my.allVars <- unique( c( my.allVars, habitatArea))

  #check raster contents
  if( !is.null( covarBrick)){
    if( !all( all.vars( distForm) %in% names( covarBrick)))
      stop( "Distribution formula terms are not present in covarBrick")
    if( !all( all.vars( biasForm) %in% names( covarBrick)))
      stop( "Bias distribution formula terms are not present in covarBrick")
    #generate design matrices, just as isdm() would do.
    newInfo <- uniqueVarNames( obsList=NULL, covarBrick=covarBrick, distForm=distForm, biasForm=biasForm, arteForm=list(), habitatArea=NULL, DCsurvID=list(), coord.names=NULL, responseNames=NULL, sampleAreaNames=NULL, stdCovs=TRUE, na.action=na.omit)
    
    xSeq <- terra::xFromCol( newInfo$covarBrick)#sort( unique( terra::crds( covarBrick)[,1]))
    ySeq <- terra::yFromRow( newInfo$covarBrick)#sort( unique( terra::crds( covarBrick)[,2]))
    
    rasterBoundary <- newInfo$covarBrick[[1]]
    tmpID <- is.na( values( rasterBoundary))
    values( rasterBoundary)[tmpID,] <- NA
    values( rasterBoundary)[!tmpID,] <- 1
  }
  else{
    message( "No raster provided. Variables will be simulated.")
    if( is.null( rasterBoundary)){
      message( "Simulation over square")
      xSeq <- seq( from=0, to=10, length=control$raster.dim[1])
      ySeq <- seq( from=0, to=10, length=control$raster.dim[2])
      my.scale <- sqrt( 100+100) / 15 #1/15 of max observed
      if( is.na( control$range))
	control$range <- my.scale / 3 #to get quite variable random effects
    }
    #or survey area based on supplied boundary raster
    else{
      message( "Simulation over supplied raster boundary")
      xSeq <- terra::xFromCol( rasterBoundary)#sort( unique( terra::crds( covarBrick)[,1]))
      ySeq <- terra::yFromRow( rasterBoundary)#sort( unique( terra::crds( covarBrick)[,2]))
      my.scale <- sqrt( (diff( range( xSeq)))^2 + ( diff( range( ySeq)))^2) / 15
      if( is.na( control$range))  #arbitrary! #1/15 of max observed 
	control$range <- my.scale / 3 #to get quite variable random effects
    }
    newVars <- matrix( NA, nrow=length(xSeq)*length(ySeq), ncol=length( my.allVars))
    for( ii in 1:length( my.allVars))
      newVars[,ii] <- as.numeric( t( fftGPsim2( x=xSeq, y=ySeq, sig2 = 1, rho = my.scale, nu = 1/2, nugget = 0.01)))
    colnames( newVars) <- my.allVars
    if( !is.null( rasterBoundary)){
      covarBrick <- terra::rast( rasterBoundary, nlyrs=ncol( newVars), names=colnames( newVars), vals=newVars)
#      names( covarBrick) <- colnames( newVars)
    }
    else{
      covarBrick <- terra::rast( nrows=length( xSeq), ncols=length( ySeq), nlyrs=ncol( newVars), xmin=min( xSeq), xmax=max( xSeq), ymin=min( ySeq), ymax=max( ySeq), names=colnames( newVars))#, type='xyz', crs="epsg:3857")
      terra::values( covarBrick) <- newVars
      rasterBoundary <- !is.na( covarBrick[[1]])
    }
    newInfo <- uniqueVarNames( obsList=NULL, covarBrick=covarBrick, distForm=distForm, biasForm=biasForm, arteForm=list(), habitatArea=NULL, DCsurvID=list(), coord.names=NULL, responseNames=NULL, sampleAreaNames=NULL, stdCovs=TRUE, na.action=na.omit)
  }
  
  #random effect for the log-gauss process
  if( is.null( Intensity)){
    if( control$addRandom){
      if( is.na( control$range))
	control$range <- 5*max( terra::res( covarBrick))
    REff <- fftGPsim2( x=xSeq, y=ySeq, sig2=control$sd^2, rho=control$range, nu=1/2)
    REff <- as.numeric( REff)
    }
    else
      REff <- rep(0, terra::ncell( covarBrick))
   
    #data brick for the covariates and random effect
    tmp <- terra::rast( newInfo$covarBrick, nlyrs=1, names="REff")
    values( tmp) <- REff
    
    tmp <- terra::mask( tmp, rasterBoundary)
    names( tmp) <- "REff"
    newInfo$covarBrick <- c( newInfo$covarBrick, tmp)
  }
  tmp <- terra::cellSize( newInfo$covarBrick)
  tmp <- terra::mask( tmp, rasterBoundary)
  names( tmp) <- "myCellSize"
  newInfo$covarBrick <- c( newInfo$covarBrick, tmp)
  
  newInfo$Xall <- as.data.frame( newInfo$covarBrick)
#  if( is.null( habitatArea))
#    newInfo$Xall[[habitatArea]] <- terra::values( terra::cellSize( newInfo$covarBrick))
#  else
#    newInfo$Xall[[habitatArea]] <- terra::values( newInfor$covarBrick[[habitatArea]])
  if( is.null( Intensity)){
    Xdist <- stats::model.matrix( newInfo$distForm, newInfo$Xall)
  
    #set up coefficients (if not specified)
    if( is.null( distCoefs)){
      distCoefs <- stats::rnorm( ncol( Xdist), mean=0, sd=1)
      names( distCoefs) <- colnames( Xdist)
    }
    if( ncol( Xdist) != length( distCoefs))
      stop( "There are a different number of distribution coefficients specified for the formula provided. Please check.")
  }
  newInfo$biasForm <- stats::update.formula( newInfo$biasForm, ~.-1)
  Xbias <- stats::model.matrix( newInfo$biasForm, newInfo$Xall)
  if( is.null( biasCoefs)){
    biasCoefs <- stats::rnorm( ncol( Xbias), mean=0, sd=1)
    names( biasCoefs) <- colnames( Xbias)
  }
  if( ncol( Xbias) != length( biasCoefs))
    stop( "There are a different number of bias coefficients specified for the formula provided. Please check.")

  if( is.null( Intensity)){
    #linear predictor
    LinPred <- Xdist %*% distCoefs +
		newInfo$Xall$REff +
		log( newInfo$Xall$myCellSize)
  
    if( !is.null( Intercept))
      LinPred <- LinPred + Intercept
      
    #Intensity for log-gauss process
    Intensity <- exp( LinPred)
  }
  else{  #just to make format consistent
    Intensity <- terra::mask( Intensity, rasterBoundary)
    Intensity <- terra::values( Intensity, na.rm=TRUE)
  }
  
  if( is.null( Intercept)){
    #the total expectation and observed number from unbiassed (actual distribution)
    tmp.Lambda <- sum( Intensity, na.rm=TRUE)
    Intensity <- Intensity / tmp.Lambda  #Intensity should sum to 1
    Intensity <- Intensity * pop.size  #per unit area, not per cell nor transect nor...
    Intercept <- log( pop.size) - log( tmp.Lambda)
  }
  else{
    pop.size <- sum( Intensity, na.rm=TRUE)
  }
  
  #scale lin pred
  LinPred <- log( Intensity)
  
  #add to databrick
  tmpR <- terra::rast( newInfo$covarBrick, nlyrs=2, names=c("LinPred","Intensity"), val=NA)
  terra::values( tmpR)[as.numeric( rownames( newInfo$Xall)),] <- cbind( LinPred, Intensity)
  newInfo$covarBrick <- c( newInfo$covarBrick, tmpR)
  
  ####################################
  ####	Create observations
  ####################################

  ####	PA data
  if( control$useMBHdesign)#nchar( system.file( package="MBHdesign")) > 0)
    PAdata <- MBHdesign::quasiSamp.raster( n=n.PA, rasterBoundary)
  else{
    tmp <- terra::as.data.frame( rasterBoundary, xy=TRUE, cells=TRUE)
    tmp <- tmp[sample.int(n=nrow(tmp), size=n.PA, replace=FALSE),]
    PAdata <- data.frame( x=tmp$x, y=tmp$y, inclusion.probabilities=NA, ID=tmp$cell)
  }
  #add the transect area
  mySizes <- terra::values( newInfo$covarBrick$myCellSize)[PAdata$ID]
  PAdata$transectArea <- transect.size * mySizes

  #intensity over transect area
  tmpIntensity <- terra::values( newInfo$covarBrick$Intensity)[PAdata$ID,] / mySizes * PAdata$transectArea
  #now simulate the PA response at those sites (based on point process)
  PAdata$PA <- stats::rbinom( n=n.PA, size=1, prob=pmax( 0, pmin( 1-exp( -tmpIntensity),1)))
  
  ####  Create the count data (AA)
  if( control$useMBHdesign)#nchar( system.file( package="MBHdesign")) > 0)
    AAdata <- MBHdesign::quasiSamp.raster( n=n.AA, rasterBoundary)
  else{
    tmp <- terra::as.data.frame( rasterBoundary, xy=TRUE, cells=TRUE)
    tmp <- tmp[sample.int(n=nrow(tmp), size=n.AA, replace=FALSE),]
    AAdata <- data.frame( x=tmp$x, y=tmp$y, inclusion.probabilities=NA, ID=tmp$cell)
  }
  #add the transect area
  mySizes <- terra::values( newInfo$covarBrick$myCellSize)[AAdata$ID]
  AAdata$transectArea <- transect.size * mySizes

  #intensity over transect area
  tmpIntensity <- ( terra::values( newInfo$covarBrick$Intensity)[AAdata$ID,] / mySizes) * AAdata$transectArea
  #now simulate the AA response at those sites (based on point process)
  AAdata$AA <- stats::rpois( n=n.AA, lambda=tmpIntensity)
  
  ####  Create the double count data (DC)
  
  #number of DC surveys is 1
  if( is.null( DC.pi))
    DC.pi <- stats::runif( 1, min=0.6, max=0.99)
    
  if( control$useMBHdesign)#nchar( system.file( package="MBHdesign")) > 0)
    DCdata <- MBHdesign::quasiSamp.raster( n=n.DC, rasterBoundary)
  else{
    tmp <- terra::as.data.frame( rasterBoundary, xy=TRUE, cells=TRUE)
    tmp <- tmp[sample.int(n=nrow(tmp), size=n.DC, replace=FALSE),]
    DCdata <- data.frame( x=tmp$x, y=tmp$y, inclusion.probabilities=NA, ID=tmp$cell)
  }
  #add the transect area
  mySizes <- terra::values( newInfo$covarBrick$myCellSize)[DCdata$ID]
  DCdata$transectArea <- transect.size * mySizes
  #intensity over transect area
  tmpIntensity <- ( terra::values( newInfo$covarBrick$Intensity)[DCdata$ID,] / mySizes) * DCdata$transectArea

  #now simulate the DC response at those sites (based on point process)
  #total number of animals available along a transect
  tmpTotCount <- stats::rpois( n=n.DC, lambda=tmpIntensity)
  #container for the DC data
  tmptmptmp <- as.data.frame( matrix( NA, ncol=3, nrow=n.DC))
  colnames( tmptmptmp) <- c( "Obs1","Obs2","Both")
  
  #thin into different categories for the observers
  cat.pis <- c( Obs1=DC.pi*(1-DC.pi), Obs2=(1-DC.pi)*DC.pi, Both=DC.pi*DC.pi, Neither=(1-DC.pi)*(1-DC.pi))
  observs <- matrix( NA, nrow=n.DC, ncol=3)
  colnames( observs) <- colnames( tmptmptmp)[1:3]
  for( ii in 1:n.DC)
    observs[ii,] <- stats::rmultinom( n=1, size=tmpTotCount[ii], prob=cat.pis)[1:3]  #don't want the last one.
  
  DCdata <- cbind( DCdata, observs)

  ####	PO data (as per the model fitted)
  #linear predictor
  LinPred_PO <- LinPred + Xbias %*% biasCoefs
  
  #Intensity for log-gauss process
  Intensity_PO <- exp( LinPred_PO)
  #the total expectation and observed number from unbiassed (actual distribution)
  tmp.Lambda <- sum( Intensity_PO, na.rm=TRUE)
  Intensity_PO <- Intensity_PO / tmp.Lambda  #Intensity should sum to 1
  Intensity_PO <- Intensity_PO * n.PO  #per unit area, not per cell nor transect nor...
  
  #scale lin pred
  LinPred_PO <- log( Intensity_PO)
  
  #add to databrick
  tmpR <- terra::rast( newInfo$covarBrick, nlyrs=2, names=c("biasLinPred","biasIntensity"), val=NA)
  terra::values( tmpR)[as.numeric( rownames( newInfo$Xall)),] <- cbind( LinPred_PO, Intensity_PO)
  newInfo$covarBrick <- c( newInfo$covarBrick, tmpR)
  
  if( !is.null( rasterBoundary))
    newInfo$covarBrick <- terra::mask( newInfo$covarBrick, rasterBoundary)
  
  ####  Take the sample of po data
  
#  tmp.Lambda <- sum( terra::values( dataBrick$biasIntensity), na.rm=TRUE)
#  probs <- dataBrick$biasIntensity * expected.n.PO / tmp.Lambda  #probably not needed as sample() internally rescales...
  
#  NAid <- is.na( terra::values( probs))
  tmpProbs <- terra::values( newInfo$covarBrick$biasIntensity) / n.PO  #so that the probs sum to 1 (will be done in sample() anyway, but making it clear here)
  tmpProbs[is.na( tmpProbs)] <- 0
  
  sampleID <- sample.int( n=length( tmpProbs), size=n.PO, prob=tmpProbs, replace=TRUE)
  presences <- as.matrix( terra::crds( newInfo$covarBrick, na.rm=FALSE)[sampleID,])

  if( control$doPlot){
    oldPar <- graphics::par( mfrow=c(4,2))
    #PA
    terra::plot( newInfo$covarBrick$Intensity, main="Intensity")
    terra::plot( !is.na( newInfo$covarBrick[[1]]), col=grey(c(1,0.95)), legend=FALSE, main="PA survey")
    points( jitter( PAdata[PAdata$PA==0,1]), jitter( PAdata[PAdata$PA==0,2]), cex=0.5, pch=20, col=grDevices::grey(0.7))
    points( jitter( PAdata[PAdata$PA==1,1]), jitter( PAdata[PAdata$PA==1,2]), col='red', pch=20, cex=0.5)
    legend( x="bottomleft", pch=c(20, 20), col=c(grDevices::grey(0.7),'red'), legend=c("Absence", "Presence"))
    #AA
    terra::plot( newInfo$covarBrick$Intensity, main="Intensity")
    terra::plot( !is.na( newInfo$covarBrick[[1]]), col=grey(c(1,0.95)), legend=FALSE, main="AA survey")
    for( ii in 0:max( AAdata$AA))
      points( AAdata[AAdata$AA==ii,1:2], col=ii+1, pch=20, cex=1.5)
    legend( x="bottomleft", pch=20, col=1:(max( AAdata$AA)+1), legend=0:max( AAdata$AA))
    #DC
    terra::plot( newInfo$covarBrick$Intensity, main="Intensity")
    terra::plot( !is.na( newInfo$covarBrick[[1]]), col=grey(c(1,0.95)), legend=FALSE, main="DC survey")
    tmpCount <- rowSums( DCdata[,c("Obs1","Obs2","Both")])
    for( ii in 0:max( tmpCount))
      points( DCdata[tmpCount==ii,1:2], col=ii+1, pch=20, cex=1.5)
    legend( x="bottomleft", pch=20, col=1:(max( tmpCount)+1), legend=0:max( tmpCount))
    
    #PO  
    terra::plot( newInfo$covarBrick$biasIntensity, main="Biassed Intensity")
    terra::plot( !is.na( newInfo$covarBrick[[1]]), col=grey(c(1,0.95)), legend=FALSE, main="PO observations")
    points( jitter( as.matrix( presences)), pch=20, cex=0.5)
    
    #reset graphics back to user's
    graphics::par( oldPar)
  }
	       
  res <- list( PO=presences,
               PA=PAdata,
               AA=AAdata,
               DC=DCdata,
               covarBrick=newInfo$covarBrick,
	       distCoefs=distCoefs, biasCoefs=biasCoefs, Intercept=Intercept, DC.pi=DC.pi,
	       pop.size=pop.size, n.PO=n.PO, n.PA=n.PA, n.AA=n.AA, n.DC=n.DC)
  
  class( res) <- "simISDMdata"
  
  return(res)
}
