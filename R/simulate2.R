
set.sim.control <- function( contr){

  if( ! "set.random.seed" %in% names( contr))
    contr$set.random.seed <- FALSE
  if( ! "random.seed" %in% names( contr))  
    contr$random.seed <- 787
  if( ! "raster.dim" %in% names( contr))
    contr$raster.dim <- rep(100,2)
  if( ! "do.plot" %in% names( contr))
    contr$do.plot <- TRUE
  if( ! "addRandom" %in% names( contr))
    contr$addRandom <- TRUE
  if( ! "sd" %in% names( contr))
    contr$sd <- 0.5
  if( ! "range" %in% names( contr))
    contr$range <- 4
  
  #add to as we go along
  if( !all( names( contr) %in% c("set.random.seed","random.seed","raster.dim","do.plot","addRandom", "sd", "range")))
    warning( "There are control parameters specified that are not used.")
  return( contr)
  
}


simulate.isdm <- function( expected.pop.size=400000, expected.n.PO=300, n.PA=150, n.AA=50, n.DC=50,
                           coefs=list(dist=c(-1,-0.25,0.75), bias=c(-2,-1.5)), 
                           DC.pis=matrix( c(0.8,0.76, 0.7,0.73, 0.82,0.67), nrow=3, ncol=2, byrow=TRUE),
                           transect.size = 0.125, #a proportion of cell size.
                           rasterBoundary=NULL,
                           control=list()){
  #note that coefs$dist[1] will be largely ignored to match expected.n.PO
  
  control <- set.sim.control(control)
  if( control$set.random.seed)
    set.seed( control$random.seed)
  
  ####  Define survey area
  if( is.null( rasterBoundary)){
    xSeq <- seq( from=0, to=10, length=control$raster.dim[1])
    ySeq <- seq( from=0, to=10, length=control$raster.dim[2])
    X <- expand.grid( x=xSeq, y=ySeq)
    my.scale <- sqrt( 100+100) / 10 #1/10 of observed
  }
  else{
    X <- as.data.frame( coordinates( rasterBoundary))
    xSeq <- sort( unique( X[,1]))
    ySeq <- sort( unique( X[,2]))
    my.extent <- as.matrix( extent( rasterBoundary))
    my.scale <- sqrt( (tail( xSeq,1) - head( xSeq,1))^2 + (tail( ySeq,1) - head( ySeq,1))^2) / 10
  }
  
  #define two covariates within the study area
  #not sure why these next two lines throw warnings...
  #  somewhere in the assignment, and the printing of the summary.  Odd.
  Mod1 <- RandomFields::RMgauss( var=1, scale=my.scale) + RandomFields::RMnugget( var=0.01)
  Mod2 <- RandomFields::RMgauss( var=1, scale=my.scale) + RandomFields::RMnugget( var=0.01)
  simmy1 <- RandomFields::RFsimulate( Mod1, x=xSeq, y=ySeq)
  simmy2 <- RandomFields::RFsimulate( Mod2, x=xSeq, y=ySeq)
  X <- cbind( X, as.numeric( as.matrix( simmy1)), as.numeric( as.matrix( simmy2)))
  X[,-(1:2)] <- apply( X[,-(1:2)], 2, scale)
  colnames( X) <- c("x","y","Altitude","Temperature")
  
  #random effect for the log-gauss process
  if( control$addRandom){
    Mod3 <- RandomFields::RMmatern( nu=1, var=control$sd^2, scale= control$range / (2))
    REff <- RandomFields::RFsimulate( Mod3, x=xSeq, y=ySeq)
    REff <- as.numeric( as.matrix( REff))
  }
  else
    REff <- rep(0,nrow( X))
  
  dataBrick <- raster::brick( raster::rasterFromXYZ( X[,c(1,2,3)]), raster::rasterFromXYZ( X[,c(1,2,4)]), raster::rasterFromXYZ( cbind( X[,1:2], REff)))
  if( !is.null( rasterBoundary))
    dataBrick <- mask( dataBrick, rasterBoundary)
  dataBrick <- raster::scale( dataBrick, center=TRUE, scale=FALSE)
  X <- raster::as.data.frame( dataBrick, xy=TRUE)
  #linear predictor
  LinPred <- coefs$dist[1] + 
              dataBrick$Altitude*coefs$dist[2] + 
              dataBrick$Temperature*coefs$dist[3] + 
              dataBrick$REff#with( X[,-5], model.matrix( ~1+Altitude+Temperature+REff)) %*% c( coefs$dist,1)
  
  #Intensity for log-gauss process
  Intensity <- exp( LinPred)
  #the total expectation and observed number from unbiassed (actual distribution)
  tmp.Lambda <- sum( raster::values( Intensity), na.rm=TRUE)
  Intensity <- Intensity / tmp.Lambda  #Intensity should sum to 1
  Intensity <- Intensity * expected.pop.size  #now sum to expected.pop.size
  
  dataBrick <- addLayer( dataBrick, LinPred)
  dataBrick <- addLayer( dataBrick, Intensity)
  
  names( dataBrick) <- c("Altitude","Temperature","REff","logIntensity","Intensity") 
  
  if( is.null( rasterBoundary))
    tmpBoundary <- !is.na( dataBrick$Altitude)
  else
    tmpBoundary <- rasterBoundary

  ####  Create the PA data
  #first create the design -- an equal probability spat-bal design
  PAdata <- MBHdesign::quasiSamp.raster( n=n.PA, tmpBoundary)
  
  #now simulate the PA response at those sites (based on point process)
  PAdata$PA <- stats::rbinom( n=n.PA, size=1, prob=1-exp( -dataBrick$Intensity[PAdata$ID]*transect.size))
  
  #add the transect area
  PAdata$transectArea <- transect.size * prod( raster::res( dataBrick))
  
  ####  Create the count data (AA)
  AAData <- MBHdesign::quasiSamp.raster( n=n.AA, tmpBoundary)
  
  #now simulate the AA response at those sites (based on point process)
  AAData$AA <- stats::rpois( n=n.AA, lambda=dataBrick$Intensity[AAData$ID]*transect.size)
  
  #add the transect area
  AAData$transectArea <- transect.size * prod( raster::res( dataBrick))

  ####  Create the double count data (DC)
  DCData <- MBHdesign::quasiSamp.raster( n=n.DC, tmpBoundary)
  
  #now simulate the DC response at those sites (based on point process)
  #total number of animals available along a transect
  tmpTotCount <- stats::rpois( n=n.DC, lambda=dataBrick$Intensity[DCData$ID]*transect.size)
  #container for the DC data
  tmptmptmp <- as.data.frame( matrix( NA, ncol=3, nrow=n.DC))
  colnames( tmptmptmp) <- c( "Obs1","Obs2","Both")
  #the survey ID
  tmpID <- rep( 1:nrow( DC.pis), each=n.DC %/% nrow( DC.pis))
  tmpID <- c( tmpID, rep( nrow( DC.pis), n.DC %% nrow( DC.pis))) #make the last survey biggest, if needed
  
  #thin into different categories  
  cat.pis <- cbind( Obs1=DC.pis[,1]*(1-DC.pis[,2]), Obs2=(1-DC.pis[,1])*DC.pis[,2], Both=DC.pis[,1]*DC.pis[,2], Neither=(1-DC.pis[,1])*(1-DC.pis[,2]))
  observs <- matrix( NA, nrow=n.DC, ncol=nrow( DC.pis))
  colnames( observs) <- colnames( tmptmptmp)[1:3]
  for( ii in 1:n.DC)
    observs[ii,] <- rmultinom( n=1, size=tmpTotCount[ii], prob=cat.pis[tmpID[ii],])[1:3]  #don't want the last one.
  
  DCData <- cbind( DCData, observs)
  
  #survey labels
  DCData$Survey <- paste0("DCsurvey_",tmpID)
  
  #add the transect area
  DCData$transectArea <- transect.size * prod( raster::res( dataBrick))
  
  ####  Now deal with the biassed observations
  
  #distance to city, or a proxy (something quite like it)...
  maxTmp <- extent( dataBrick)[2] + extent( dataBrick)[4]#max( xSeq) + max( ySeq)  #sqrt( max( xSeq)^2+max( ySeq)^2)
  dataBrick <- raster::addLayer( dataBrick, 
        raster::rasterFromXYZ( cbind( raster::coordinates( dataBrick)[,1:2], dist2City=(maxTmp-rowSums( raster::coordinates(dataBrick))) / maxTmp)))
  names( dataBrick)[length( names( dataBrick))] <- "dist2City"
  
  #probability of being observed (thinning process) -- cloglog link for prob
  dataBrick <- raster::addLayer( dataBrick, 1-exp(- exp( coefs$bias[1]+coefs$bias[2]*dataBrick$dist2City)))
  names( dataBrick)[length( names( dataBrick))] <- "obsProb"
  #The biassed intensity after accounting for observation bias.
  dataBrick <- raster::addLayer( dataBrick, biasIntensity=dataBrick$Intensity * dataBrick$obsProb)
  names( dataBrick)[length( names( dataBrick))] <- "biasIntensity"
  
  if( !is.null( rasterBoundary))
    dataBrick <- raster::mask( dataBrick, rasterBoundary)
  dataBrick$dist2City <- raster::scale( dataBrick$dist2City, center=TRUE, scale=FALSE)
  
  ####  Take the sample of po data
  N <- stats::rpois( n=1, lambda=expected.n.PO)
  tmp.Lambda <- sum( raster::values( raster::raster( dataBrick, layer="biasIntensity")), na.rm=TRUE)
  probs <- dataBrick$biasIntensity * expected.n.PO / tmp.Lambda  #probably not needed as sample() internally rescales...
  
  NAid <- is.na( raster::values( probs))
  tmpProbs <- raster::values( probs)
  tmpProbs[NAid] <- 0
  
  sampleID <- sample.int( n=length( tmpProbs), size=N, prob=tmpProbs, replace=TRUE)
  presences <- raster::coordinates( dataBrick)[sampleID,]

  if( control$do.plot){
    par( mfrow=c(4,2))
    #PA
    raster::plot( dataBrick$Intensity, main="Intensity")
    raster::plot( PAdata[,1:2], type='n', main="PA survey", xlim=range( xSeq), ylim=range( ySeq), asp=1)
    points( jitter( PAdata[PAdata$PA==0,1]), jitter( PAdata[PAdata$PA==0,2]), cex=0.5, pch=20, col=grey(0.7))
    points( jitter( PAdata[PAdata$PA==1,1]), jitter( PAdata[PAdata$PA==1,2]), col='red', pch=20, cex=0.5)
    legend( x="bottomleft", pch=c(20, 20), col=c(grey(0.7),'red'), legend=c("Absence", "Presence"))
    #AA
    raster::plot( dataBrick$Intensity, main="Intensity")
    raster::plot( AAData[,1:2], type='n', main="AA survey", xlim=range( xSeq), ylim=range( ySeq), asp=1)
    for( ii in 0:max( AAData$AA))
      points( AAData[AAData$AA==ii,1:2], col=ii+1, pch=20, cex=1.5)
    legend( x="bottomleft", pch=20, col=1:(max( AAData$AA)+1), legend=0:max( AAData$AA))
    #DC
    raster::plot( dataBrick$Intensity, main="Intensity")
    raster::plot( DCData[,1:2], type='n', main="DC survey", xlim=range( xSeq), ylim=range( ySeq), asp=1)
    tmpCount <- rowSums( DCData[,c("Obs1","Obs2","Both")])
    for( ii in 0:max( tmpCount))
      points( DCData[tmpCount==ii,1:2], col=ii+1, pch=20, cex=1.5)
    legend( x="bottomleft", pch=20, col=1:(max( tmpCount)+1), legend=0:max( tmpCount))
    
    #PO  
    raster::plot( dataBrick$biasIntensity, main="Biassed Intensity")
    plot( jitter( presences), pch=20, cex=0.5, main="Presences", asp=1)
  }
  
  res <- list( PO=sp::SpatialPoints( presences), 
               PA=sp::SpatialPointsDataFrame(coords = PAdata[,1:2], data = PAdata[,-(1:2)]), 
               AA=sp::SpatialPointsDataFrame( coords=AAData[,1:2], data=AAData[,-(1:2)]), 
               DC=sp::SpatialPointsDataFrame( coords=DCData[,1:2], data=DCData[,-(1:2)]),
               covarBrick=dataBrick,
               expected.pop.size=expected.pop.size)
  
  return(res)
}
