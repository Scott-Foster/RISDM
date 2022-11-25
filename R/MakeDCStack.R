
###############################################################################################
###############################################################################################
####	Set up the INLA stacks for the DC component of an isdm.
####	Assumes that there has already been a fair bit of processing to get things into the
####	the right, consistent formats.
####
####	Returns an INLA stack object with correctly set-up multiple responses etc...
####
####	Programmed by Scott in the first half of 2022
####
###############################################################################################
###############################################################################################


MakeDCStack=function( observs, mesh = NULL, DCname = "DCabund", tag = "DC", varNames, sampleAreaName, ind) {
  
  #casting to numeric
  y.pp <- as.integer( observs@data[,DCname])
  
  #the area sampled
  offy <- observs@data[,sampleAreaName]
  
  #covariates 'n' stuff
  tmp <- varNames %in% colnames( observs@data)
#  if( ! all( tmp))
#    warning( "Not all bias covariates in DC data. Missing: ", paste(varNames[!tmp], sep=" "), "Creating variable and padding with NAs.")
  tmptmp <- matrix( NA, nrow=nrow( observs@data), ncol=sum( !tmp), dimnames=list( NULL, varNames[!tmp]))
  observs@data <- cbind( observs@data, as.data.frame( tmptmp))
  
  covars <- as.data.frame( observs@data[,varNames])
  colnames( covars) <- varNames
  
  covars$Intercept.DC <- 1
  
  resp <- matrix( NA, nrow=nrow( covars), ncol=sum( ind))
  colnames( resp) <- names( ind[ind!=0])  #name the variables
  resp[,"DC"] <- y.pp
  
  # Projector matrix from mesh to data.
  projmat <- INLA::inla.spde.make.A(mesh, as.matrix( raster::coordinates( observs))) # from mesh to point observations
  
  #make the stack
  stk.DCabund <- INLA::inla.stack(data=list(resp=resp, e=rep( 1, nrow( covars)), offy=log( offy)),
                                A=list(1,projmat), tag=tag,
                                effects=list( covars, list(isdm.spat.XXX=1:mesh$n)))
  
  return( stk.DCabund)
}
