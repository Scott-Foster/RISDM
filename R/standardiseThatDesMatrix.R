
###############################################################################################
###############################################################################################
####	Take a matrix or data.frame (or similar) and standardises all numerics
####	So that they have mean 0 and sd 1.
####	
####	Returns the same format with changed data.
####
####	Programmed by Scott in March/April 2023
####
###############################################################################################
###############################################################################################

standardiseThatDesMatrix <- function( matty){
  if( is.null( dim( matty)))
    matty <- matrix( matty, ncol=1)
    
    for( jj in 1:ncol( matty)){ #loop through cols
      if( inherits( matty[,jj], "numeric"))  #standardise only if a numeric (not a factor or char)
	if( ! all( matty[,jj] %in% c(0,1,NA))) #a design matrix factor (another interpretation of a factor)
#	  if( ! any( grepl( "alpha1Coef", colnames( matty)[jj]), grepl( "alpha2Coef", colnames( matty[jj]))))  #a special DC variable
	  if( ! grepl( "logDetectPi", colnames( matty)[jj]))  #a special DC variable
	    matty[,jj] <- scale( matty[,jj])
    }
  
  return( matty)
}
