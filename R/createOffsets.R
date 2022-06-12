
###############################################################################################
###############################################################################################
####	Each non-PO data set will need to have an offset.  This function extracts them and 
####	uses the (dumb) value of 1 (or log(1)) as default.
####
####	Returns a list that contains the objects for each daa type.
####
####	Programmed by Scott in the first half of 2022
####
###############################################################################################
###############################################################################################

createOffsets <- function( sampleAreaNames, observationList){
  
  #Are areanames already specified by the user?  If not, then create them as dumb values.
  if( is.null( sampleAreaNames)){
    sampleAreaNames <- rep("area1",3)
    names( sampleAreaNames) <- c("PA","AA","DC")
    if( is.null( observationList$PAdat))
      sampleAreaNames <- sampleAreaNames[names( sampleAreaNames) != "PA"]
    else
      observationList$PAdat$area1 <- 1
    
    if( is.null( observationList$AAdat))
      sampleAreaNames <- sampleAreaNames[names( sampleAreaNames) != "AA"]
    else
      observationList$AAdat$area1 <- 1
    
    if( is.null( observationList$DCdat))
      sampleAreaNames <- sampleAreaNames[names( sampleAreaNames) != "DC"]
    else
      observationList$DCdat$area1 <- 1
  }
  
  res <- list( names=sampleAreaNames, dat=observationList)
  
}
