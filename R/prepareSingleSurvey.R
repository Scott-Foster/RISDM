
###############################################################################################
###############################################################################################
####	
####	For a single DC survey, arrange expansion stuff, pre-estimate
####
####	Returns a matrix for each survey full of good bits 'n' pieces for DC data estimation 
####
####	Programmed by Scott in the first half of 2022
####
###############################################################################################
###############################################################################################

prepareSingleSurvey <- function( singleDataSet, datasetID, DCcols, survIDname, sampAreaDC, DCmethod){
  
  n <- nrow( singleDataSet)
  
  #check if estimation can be done (without penalties and boundary effects etc)
  tmpDat <- singleDataSet[,DCcols]
  if( ! all( colSums( tmpDat) > 0))
    stop( "At least one of the DC data sets has an observer that hasn't seen a single individual. This is currently outside of RISDM's scope. Sorry. Consider adjusting model.")
  
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
  
#  #removing unwanted cols
#  singleDataSet <- singleDataSet[,-(1:3)]
  
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

