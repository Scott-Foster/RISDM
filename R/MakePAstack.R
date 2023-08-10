
###############################################################################################
###############################################################################################
####	Set up the INLA stacks for the PA component of an isdm.
####	Assumes that there has already been a fair bit of processing to get things into the
####	the right, consistent formats.
####
####	Returns an INLA stack object with correctly set-up multiple responses etc...
####
####	Programmed by Scott in the first half of 2022
####	Genesis of code started with Isaac's et al. (2020) R functions
####	####	Re-jigged in March 2023 to make explicit for INLA the particular pattern
####		of constraints for the factor level effects
####
###############################################################################################
###############################################################################################

#MakePAstack=function( observs, distForm=NULL, mesh = NULL, presname = "GroupSize", tag = "PA", sampleAreaName, ind) {#, varNames
MakePAstack=function( observs, mesh = NULL, presname = "GroupSize", tag = "PA", sampleAreaName, ind) {#, varNames

  # Function to create stack for presence absence points  NOT binomial points
  
  #casting to numeric
  y.pp <- as.integer( terra::as.data.frame( observs[,presname])[,1])
  
  #the number of binomial trials (Bernoulli)
  ntrials <- rep(1, nrow(observs))

  #the area sampled
  offy <- terra::as.data.frame( observs[,sampleAreaName])[,1]

  #expanding the artefact formulas, just to make sure that they get expanded properly in the INLA call (funny parsing for aliasing?)
  covars <- terra::as.data.frame( observs[,-(1:2),drop=FALSE])[,1:(ncol( observs)-3), drop=FALSE]
  
  #building the response matrix.
  resp <- matrix( NA, nrow=nrow( covars), ncol=sum( ind))
  colnames( resp) <- names( ind[ind!=0])  #name the variables
  resp[,"PA"] <- y.pp

  # Projector matrix from mesh to data.
  projmat <- INLA::inla.spde.make.A( mesh, as.matrix( sf::st_coordinates( observs))) # from mesh to point observations

  #make the stack
  stk.binom <- INLA::inla.stack( data=list( resp=resp, Ntrials=ntrials, offy=log( offy)),
                          A=list(1,projmat), tag=tag,
                          effects=list( covars, list(isdm.spat.XXX=1:mesh$n)))

#  attr( stk.binom, "newArtForm") <- newArtForm
  
  return( stk.binom)
}

