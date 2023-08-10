
###############################################################################################
###############################################################################################
####	
####	Get the covariates at mesh locations.
####
####	Returns a data.frame of covariates at nice locations.
####
####	Programmed by Scott in the first half of 2022
####
###############################################################################################
###############################################################################################

ExtractCovarsAtNodes <- function( mesh=NULL, covars=NULL){
  
  #mesh is the INLA mesh for the region
  #covars is a raster/brick containing all the covariate values.  Mush contain spatial area for all the mesh nodes
  
  locs <- as.matrix( mesh$loc)
  #blinear to follow Krainski et al and Flagg + Koegh.  Well almost consistent (different packages and methods)
  #FWIW, I'd use 'simple' to avoid Berkson errors.
  covarAtLocs <- terra::extract( covars, locs[,1:2], cells=TRUE)#, small=TRUE, buffer=min( res( covars)), fun='median', methd='simple')

  return( covarAtLocs)
}

