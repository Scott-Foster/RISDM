
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
####
###############################################################################################
###############################################################################################

MakePAstack=function( observs, mesh = NULL, presname = "GroupSize", tag = "PA", varNames, sampleAreaName, ind) {

  # Function to create stack for presence absence points  NOT binomial points
  
  # observs SpatialPointsDataFrame of covariates and observations. Must contain a column corresponding to presename, which must be numeric (0 or 1) or boolean
  # mesh INLA mesh for the SPDE model
  # presname Name of presences column in observs.
  # offsetname Name of the offset variable for PA data.
  # tag Name for tag for the stack.
  
  # An INLA stack with binomial data: include Ntrials, which is the number of trials
  
  #casting to numeric
  y.pp <- as.integer( observs@data[,presname])

  #the number of binomial trials (Bernoulli)
  ntrials <- rep(1, nrow(observs@data))

  #the area sampled
  offy <- observs@data[,sampleAreaName]

  #covariates 'n' stuff
  tmp <- varNames %in% colnames( observs@data)
#  if( ! all( tmp))
#    warning( "Not all bias covariates in PA data. Missing: ", paste( varNames[!tmp], " "), "Creating variable and padding with NAs.")
  tmptmp <- matrix( NA, nrow=nrow( observs@data), ncol=sum( !tmp), dimnames=list( NULL, varNames[!tmp]))
  observs@data <- cbind( observs@data, as.data.frame( tmptmp))
  
  covars <- as.data.frame( observs@data[,varNames])
  colnames( covars) <- varNames

  #intercept added for each data type
  covars$Intercept.PA <- 1

  resp <- matrix( NA, nrow=nrow( covars), ncol=sum( ind))
  colnames( resp) <- names( ind[ind!=0])  #name the variables
  resp[,"PA"] <- y.pp

  # Projector matrix from mesh to data.
  projmat <- INLA::inla.spde.make.A( mesh, as.matrix( raster::coordinates( observs))) # from mesh to point observations

  #make the stack
  stk.binom <- INLA::inla.stack( data=list( resp=resp, Ntrials=ntrials, offy=log( offy)),
                          A=list(1,projmat), tag=tag,
                          effects=list( covars, list(isdm.spat.XXX=1:mesh$n)))

  return( stk.binom)
}

