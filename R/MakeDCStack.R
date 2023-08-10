
###############################################################################################
###############################################################################################
####	Set up the INLA stacks for the DC component of an isdm.
####	Assumes that there has already been a fair bit of processing to get things into the
####	the right, consistent formats.
####
####	Returns an INLA stack object with correctly set-up multiple responses etc...
####
####	Programmed by Scott in the first half of 2022
####	Re-jigged in March 2023 to make explicit for INLA the particular pattern
####		of constraints for the factor level effects
####
###############################################################################################
###############################################################################################


MakeDCstack=function( observs, mesh = NULL, DCname = "DCabund", tag = "DC", sampleAreaName, ind) {

  # Function to create stack for DC abundance absence points
  
  #casting to numeric
  y.pp <- as.integer( terra::as.data.frame( observs[,DCname])[,1])
  
  #the area sampled
  offy <- terra::as.data.frame( observs[,sampleAreaName])[,1]
  
#  #expanding the artefact formulas, just to make sure that they get expanded properly in the INLA call (funny parsing for aliasing?)
  covars <- terra::as.data.frame( observs[,-(1:2),drop=FALSE])[,1:(ncol( observs)-3), drop=FALSE]
  
  #building the response matrix.
  resp <- matrix( NA, nrow=nrow( covars), ncol=sum( ind))
  colnames( resp) <- names( ind[ind!=0])  #name the variables
  resp[,"DC"] <- y.pp

  # Projector matrix from mesh to data.
  projmat <- INLA::inla.spde.make.A(mesh, as.matrix( sf::st_coordinates( observs))) # from mesh to point observations
  
  #make the stack
  stk.DCabund <- INLA::inla.stack(data=list(resp=resp, e=rep( 1, nrow( covars)), offy=log( offy)),
                                A=list(1,projmat), tag=tag,
                                effects=list( covars, list(isdm.spat.XXX=1:mesh$n)))
				
#  attr( stk.DCabund, "newArtForm") <- newArtForm
  
  return( stk.DCabund)
}



