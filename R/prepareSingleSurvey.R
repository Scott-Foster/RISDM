
###############################################################################################
###############################################################################################
####	
####	For a single DC survey, arrange expansion stuff, pre-estimate
####
####	Returns a matrix for each survey full of good bits 'n' pieces for DC data estimation 
####
####	Programmed by Scott in the first half of 2022
####	Re-factored by Scott Feb 2024
#### 		to have only 1 prob, shared for both observers (sometimes the labels are
####		swapped and observers change and...)
####
###############################################################################################
###############################################################################################

prepareSingleSurvey <- function( singleDataSet, datasetID, DCcols, survIDname, sampAreaDC, DCmethod){
  
  n <- nrow( singleDataSet)
  
  #check if estimation can be done (without penalties and boundary effects etc)
  tmpDat <- singleDataSet[,DCcols]
  tmpSums <- colSums( tmpDat)
  if( ( tmpSums[1]==0 & tmpSums[3]==0) | ( tmpSums[2]==0 & tmpSums[3]==0))
    stop( "One (or more) of the DC data sets has not seen a single individual. This is currently outside of RISDM's scope. Sorry.")
#  if( ! all( colSums( tmpDat) > 0))
#    stop( "At least one of the DC data sets has an observer that hasn't seen a single individual. This is currently outside of RISDM's scope. Sorry. Consider adjusting model.")
  
  #the bits and pieces for approximating log( 1-pi)
  expansPis <- estimatePisDoubleCount( singleDataSet[,DCcols])
  singleDataSet$pi <- expansPis
#  singleDataSet$Obs1_pi <- expansPis[1]
#  singleDataSet$Obs2_pi <- expansPis[2]
  
  #reformatting and filling in
  county <- unlist( singleDataSet[,DCcols])
  index <- rep( 1:nrow( singleDataSet), times=3)  #once for each double count type
  singleDataSet <- singleDataSet[index,]
  singleDataSet$DCcountDC <- county #this variable name shouldn't get confused with others...?

  #identifiers
  singleDataSet$observer <- rep( DCcols, each=n)
  
  #storing the original offset before it is altered.
  singleDataSet$originalSampleArea <- singleDataSet[,sampAreaDC]
  
  #the offset for plugin observer prob.
  if( DCmethod=="plugin"){
    offy <- rep( NA, 3*n)
#    offy[1:n] <- expansPis[1]*(1-expansPis[2])
#    offy[n+1:n] <- (1-expansPis[1]) * expansPis[2]
#    offy[2*n + 1:n] <- expansPis[1] * expansPis[2]
    offy[1:(2*n)] <- expansPis*(1-expansPis)
    offy[2*n+1:n] <- expansPis^2
    singleDataSet$observerOffy <- offy
    singleDataSet[,sampAreaDC] <- singleDataSet[,sampAreaDC] * offy
    
    return( singleDataSet)
  }

  #the offset for TaylorLinApprox of log( 1-pi) around log( pi)  
#  singleDataSet$expansOffset <- singleDataSet$alpha1Coef <- singleDataSet$alpha2Coef <- NA
  singleDataSet$expansOffset <- singleDataSet$logDetectPi <- NA
  #Obs1 (not Obs2)
#  singleDataSet$expansOffset[1:n] <- log( 1-expansPis[2]) - (expansPis[2] / (expansPis[2]-1)) * log( expansPis[2])
#  singleDataSet$alpha1Coef[1:n] <- 1
#  singleDataSet$alpha2Coef[1:n] <- expansPis[2] / (expansPis[2]-1)
  singleDataSet$expansOffset[1:n] <- log( 1-expansPis) + expansPis * log( expansPis) / (1-expansPis)
  singleDataSet$logDetectPi[1:n] <- -expansPis / (1-expansPis)
  #Obs2 (not Obs1)
#  singleDataSet$expansOffset[n+1:n] <- log( 1-expansPis[1]) - (expansPis[1] / (expansPis[1]-1)) * log( expansPis[1])
#  singleDataSet$alpha1Coef[n+1:n] <- expansPis[1] / (expansPis[1]-1)
#  singleDataSet$alpha2Coef[n+1:n] <- 1
  singleDataSet$expansOffset[n+1:n] <- log( 1-expansPis) + expansPis * log( expansPis) / (1-expansPis)
  singleDataSet$logDetectPi[n+1:n] <- -expansPis / (1-expansPis)
  #Both
#  singleDataSet$expansOffset[2*n+1:n] <- 0
#  singleDataSet$alpha1Coef[2*n+1:n] <- 1
#  singleDataSet$alpha2Coef[2*n+1:n] <- 1
  singleDataSet$expansOffset[2*n+1:n] <- 0
  singleDataSet$logDetectPi[2*n+1:n] <- 2


  singleDataSet[,sampAreaDC] <- exp( log( singleDataSet[,sampAreaDC]) + singleDataSet$expansOffset)
    
  return( singleDataSet)
  
}

