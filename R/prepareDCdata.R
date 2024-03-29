
###############################################################################################
###############################################################################################
####	
####	Get the DC data ready for inclusion into the isdm framework.
####
####	Returns a matrix full of good bits 'n' pieces for DC data estimation 
####
####	Programmed by Scott in the first half of 2022
####	Re-factored by Scott Feb 2024
#### 		to have only 1 prob, shared for both observers (sometimes the labels are
####		swapped and observers change and...)
####
###############################################################################################
###############################################################################################

prepareDCdata <- function( DCdat, DCobserverInfo, sampAreaDC, DCmethod, coord.names){

  #the number of DC surveys to include (multiplied by 2 to get the number of observer probs).
  nsurvey <- length( unique( DCdat[,DCobserverInfo$SurveyID]))
  #identify the right columns.
  DCcols <- unlist( DCobserverInfo[c("Obs1","Obs2","Both")])
  #this way of compiling information might be slow and memory hungry for high numbers of datasets and possibly for large data.  
  #Should be OK for almost all real applications though...
  tmp <- list()  
  kount <- 1
  for( ss in unique( DCdat[,DCobserverInfo$SurveyID])){
    #subset the survey
    ssID <- which( DCdat[,DCobserverInfo$SurveyID]==ss)
    #prepare each survey's data
    tmp[[kount]] <- prepareSingleSurvey( singleDataSet=DCdat[ssID,], datasetID=ss, DCcols=DCcols, survIDname=DCobserverInfo$SurveyID, sampAreaDC=sampAreaDC, DCmethod=DCmethod)
    #append it.
    tmp[[kount]] <- cbind( ss, tmp[[kount]])
    #corretly name
    colnames( tmp[[kount]])[1] <- DCobserverInfo$SurveyID
    kount <- kount + 1
  }
  #combine information
  DCdatExpand <- do.call( "rbind", tmp)
  #remove excess
  DCdatExpand <- DCdatExpand[,!colnames( DCdatExpand) %in% DCcols]
#locations should already be part of the data, in the right order too!
##  #add locations
##  ntimes <- nrow( DCdatExpand) / nrow( DCdat)
##  DCdatExpand <- cbind( kronecker( matrix(1,nrow=ntimes,ncol=1), as.matrix( DCdat[,coord.names])), DCdatExpand)
##  colnames( DCdatExpand)[1:2] <- coord.names

  attr( DCdatExpand, "nsurvey") <- nsurvey
  
  return( DCdatExpand)
}

