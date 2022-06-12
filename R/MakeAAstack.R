
###############################################################################################
###############################################################################################
####	Set up the INLA stacks for the AA component of an isdm.
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

MakeAAStack=function( observs, mesh = NULL, abundname = "abund", tag = "AA", varNames, sampleAreaName, ind) {

  # Function to create stack for ABUNDANCE absence points
    
  #casting to numeric
  y.pp <- as.integer( observs@data[,abundname])

  #the area sampled
  offy <- observs@data[,sampleAreaName]

  #covariates 'n' stuff
  tmp <- varNames %in% colnames( observs@data)
  if( ! all( tmp))
    warning( "Not all bias covariates in AA data. Missing: ", paste( varNames[!tmp], " "), "Creating variable and padding with NAs.")
  tmptmp <- matrix( NA, nrow=nrow( observs@data), ncol=sum( !tmp), dimnames=list( NULL, varNames[!tmp]))
  observs@data <- cbind( observs@data, as.data.frame( tmptmp))
    
  covars <- as.data.frame( observs@data[,varNames])
  colnames( covars) <- varNames
  #  covars$Intercept <- 1

  #still need to think about what is going to be the reference level for the overall prevalence
  #  if( sum( ind) > 1)
  covars$Intercept.AA <- 1

  resp <- matrix( NA, nrow=nrow( covars), ncol=sum( ind))
  colnames( resp) <- names( ind[ind!=0])  #name the variables
  resp[,"AA"] <- y.pp

  # Projector matrix from mesh to data.
  projmat <- INLA::inla.spde.make.A( mesh, as.matrix( raster::coordinates( observs))) # from mesh to point observations

  #make the stack
  stk.abund <- INLA::inla.stack(data=list(resp=resp, e=rep( 1, nrow( covars)), offy=log( offy)),
                          A=list(1,projmat), tag=tag,
                          effects=list( covars, list(isdm.spat.XXX=1:mesh$n)))

  return( stk.abund)
}

