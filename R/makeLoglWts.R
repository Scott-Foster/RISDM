
###############################################################################################
###############################################################################################
####	
####	Make weights for the approximation of the ppp logl.
####
####	Returns a vector double.
####
####	Programmed by Scott in the first half of 2022
####
####	Currently unused. It was a hack to start with...
####
###############################################################################################
###############################################################################################

makeLoglWts <- function( stck){
  wts <- rep( 1, nrow( stck$data$data))
  
  wts[which( stck$data$data$resp==0 & stck$data$data$e==0)] <- 0
  
  return( wts)
}

