
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
  y.pp <- as.integer( observs@data[,DCname])
  
  #the area sampled
  offy <- observs@data[,sampleAreaName]
  
#  #expanding the artefact formulas, just to make sure that they get expanded properly in the INLA call (funny parsing for aliasing?)
  covars <- observs@data[,-(1:2),drop=FALSE]#stats::model.matrix( artForm, as.data.frame( observs))
#  tmptmp <- stats::model.matrix( artForm, as.data.frame( observs))
#  colnames( tmptmp) <- paste0( "Intercept.DC_", colnames( tmptmp))
#  colnames( tmptmp) <- gsub( "_(Intercept)", "", colnames( tmptmp), fixed=TRUE)
#  #without this substituion the ordering of the compound variable names can sometimes be reversed
#  #best to remove things that are parsed by formulas
#  colnames( tmptmp) <- gsub( ":", "_", colnames( tmptmp), fixed=TRUE)
#  observs@data <- cbind( observs@data, as.data.frame( tmptmp))
  
#  #a new formula for the new data labels
#  newArtForm <- reformulate( colnames( tmptmp))
  
#  #a new formula for the expanded data
#  newFullForm.DC <- reformulate( c( strsplit( deparse1( distForm[[2]]), " + ", fixed=TRUE)[[1]], strsplit( deparse1(newArtForm[[2]]), " + ", fixed=TRUE)[[1]]))
  
#  #just those covariates needed (with sensible labels)
#  covars <- as.data.frame( observs@data[,attr( terms( newFullForm.DC), "term.labels")])
#  colnames( covars) <- attr( terms( newFullForm.DC), "term.labels")
  
  #building the response matrix.
  resp <- matrix( NA, nrow=nrow( covars), ncol=sum( ind))
  colnames( resp) <- names( ind[ind!=0])  #name the variables
  resp[,"DC"] <- y.pp

  # Projector matrix from mesh to data.
  projmat <- INLA::inla.spde.make.A(mesh, as.matrix( raster::coordinates( observs))) # from mesh to point observations
  
  #make the stack
  stk.DCabund <- INLA::inla.stack(data=list(resp=resp, e=rep( 1, nrow( covars)), offy=log( offy)),
                                A=list(1,projmat), tag=tag,
                                effects=list( covars, list(isdm.spat.XXX=1:mesh$n)))
				
#  attr( stk.DCabund, "newArtForm") <- newArtForm
  
  return( stk.DCabund)
}



