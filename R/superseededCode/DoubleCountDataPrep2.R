
DC_loglikFun <- function( par, dat){
  par <- exp( par)
  lambda <- matrix( NA, nrow=nrow( dat), ncol=3)
  lambda[,1] <- par[1]*(1-par[2])
  lambda[,2] <- (1-par[1])*par[2]
  lambda[,3] <- par[1]*par[2]
 
  logl_ij <- dat * log( lambda) #multinomial 
  # logl_ij <- dat*log( lambda) - lambda
  logl_i <- rowSums( logl_ij)
  logl <- sum( logl_i)
  
  return( -logl)
}

estimatePisDoubleCount <- function( dat){
  #dat is an nx3 matrix with colnames c("obs1","obs2","both")
#data should already be ordered by the way that this function is called...
  #dat <- as.matrix( dat[,c("Obs1","Obs2","Both")])  #make sure that data is in 'right' order and format
  start.vals <- cbind( (dat[,1]+dat[,3]), dat[,2]+dat[,3]) / rowSums( dat, na.rm=TRUE)
  start.vals <- colMeans( start.vals, na.rm=TRUE)
  start.vals <- log( start.vals + 0.01)
  
  tmp <- stats::optim( par=start.vals, fn=DC_loglikFun, dat=dat)
  if( tmp$convergence != 0)
    warning("Nelder-Mead optimisation for expansion points/probabilities has not worked properly (convergence code != 0)")
  expans.pts <- exp( tmp$par)  
  
  return( expans.pts)
}

prepareSingleSurvey <- function( singleDataSet, datasetID, DCcols, survIDname, sampAreaDC, DCmethod){
  
  n <- nrow( singleDataSet)
  
  #the bits and pieces for approximating log( 1-pi)
  expansPis <- estimatePisDoubleCount( singleDataSet[,DCcols])
  singleDataSet$Obs1_pi <- expansPis[1]
  singleDataSet$Obs2_pi <- expansPis[2]
  
  #reformatting and filling in
  county <- unlist( singleDataSet[,DCcols])
  index <- rep( 1:nrow( singleDataSet), times=3)  #once for each double count type
  singleDataSet <- singleDataSet[index,]
  singleDataSet$DCcountDC <- county #this variable name shouldn't get confused with others...?

  #identifiers
  singleDataSet$observer <- rep( DCcols, each=n)
  
  #removing unwanted cols
  singleDataSet <- singleDataSet[,-(1:3)]
  
  #storing the original offset before it is altered.
  singleDataSet$originalSampleArea <- singleDataSet[,sampAreaDC]
  
  #the offset for plugin observer prob.
  if( DCmethod=="plugin"){
    offy <- rep( NA, 3*n)
    offy[1:n] <- expansPis[1]*(1-expansPis[2])
    offy[n+1:n] <- (1-expansPis[1]) * expansPis[2]
    offy[2*n + 1:n] <- expansPis[1] * expansPis[2]
    singleDataSet$observerOffy <- offy
    singleDataSet[,sampAreaDC] <- singleDataSet[,sampAreaDC] * offy
    
    return( singleDataSet)
  }

  #the offset for TaylorLinApprox of log( 1-pi) around log( pi)  
  singleDataSet$expansOffset <- singleDataSet$alpha1Coef <- singleDataSet$alpha2Coef <- NA
  #Obs1 (not Obs2)
  singleDataSet$expansOffset[1:n] <- log( 1-expansPis[2]) - (expansPis[2] / (expansPis[2]-1)) * log( expansPis[2])
  singleDataSet$alpha1Coef[1:n] <- 1
  singleDataSet$alpha2Coef[1:n] <- expansPis[2] / (expansPis[2]-1)
  #Obs2 (not Obs1)
  singleDataSet$expansOffset[n+1:n] <- log( 1-expansPis[1]) - (expansPis[1] / (expansPis[1]-1)) * log( expansPis[1])
  singleDataSet$alpha1Coef[n+1:n] <- expansPis[1] / (expansPis[1]-1)
  singleDataSet$alpha2Coef[n+1:n] <- 1
  #Both
  singleDataSet$expansOffset[2*n+1:n] <- 0
  singleDataSet$alpha1Coef[2*n+1:n] <- 1
  singleDataSet$alpha2Coef[2*n+1:n] <- 1

  singleDataSet[,sampAreaDC] <- exp( log( singleDataSet[,sampAreaDC]) + singleDataSet$expansOffset)
    
  return( singleDataSet)
  
}

prepareDCdata <- function( DCdat, DCobserverInfo, sampAreaDC, DCmethod){
  nsurvey <- length( unique( DCdat[,DCobserverInfo$SurveyID]))
  DCcols <- unlist( DCobserverInfo[c("Obs1","Obs2","Both")])
  tmp <- list()  #this way of compiling information is going to be slow and memory hungry for high numbers of datasets and possibly for large data.  Should be OK for almost all real applications though...
  kount <- 1
  for( ss in unique( DCdat[,DCobserverInfo$SurveyID])){
    ssID <- which( DCdat[,DCobserverInfo$SurveyID]==ss)
    tmp[[kount]] <- prepareSingleSurvey( singleDataSet=DCdat[ssID,], datasetID=ss, DCcols=DCcols, survIDname=DCobserverInfo$SurveyID, sampAreaDC=sampAreaDC, DCmethod=DCmethod)
    tmp[[kount]] <- cbind( ss, tmp[[kount]])
    colnames( tmp[[kount]])[1] <- DCobserverInfo$SurveyID
    kount <- kount + 1
  }
  DCdatExpand <- do.call( "rbind", tmp)
  DCdatExpand <- DCdatExpand[,!colnames( DCdatExpand) %in% DCcols]
  
  return( DCdatExpand)
}


