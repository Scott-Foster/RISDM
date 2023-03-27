
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
  y.pp <- as.integer( observs@data[,presname])

  #the number of binomial trials (Bernoulli)
  ntrials <- rep(1, nrow(observs@data))

  #the area sampled
  offy <- observs@data[,sampleAreaName]

  #expanding the artefact formulas, just to make sure that they get expanded properly in the INLA call (funny parsing for aliasing?)
  covars <- observs@data[,-(1:2),drop=FALSE]#stats::model.matrix( artForm, as.data.frame( observs))
  
#  colnames( tmptmp) <- paste0( "Intercept.PA_", colnames( tmptmp))
#  colnames( tmptmp) <- gsub( "_(Intercept)", "", colnames( tmptmp), fixed=TRUE)
#  #without this substituion the ordering of the compound variable names can sometimes be reversed
#  #best to remove things that are parsed by formulas
#  colnames( tmptmp) <- gsub( ":", "_", colnames( tmptmp), fixed=TRUE)
#  observs@data <- cbind( observs@data, as.data.frame( tmptmp))
  
#  #a new formula for the new data labels
#  newArtForm <- reformulate( colnames( tmptmp))
  
#  #a new formula for the expanded data
#  newFullForm.PA <- reformulate( c( strsplit( deparse1( distForm[[2]]), " + ", fixed=TRUE)[[1]], strsplit( deparse1(newArtForm[[2]]), " + ", fixed=TRUE)[[1]]))
##  newFullForm.PA <- reformulate( attr( terms( reformulate( c( attr( terms( distForm), "term.labels"), attr( terms( newArtForm), "term.labels")))),"term.labels"))
  
#  #just those covariates needed (with sensible labels)
#  covars <- as.data.frame( observs@data[,attr( terms( newFullForm.PA), "term.labels")])
#  colnames( covars) <- attr( terms( newFullForm.PA), "term.labels")

  #building the response matrix.
  resp <- matrix( NA, nrow=nrow( covars), ncol=sum( ind))
  colnames( resp) <- names( ind[ind!=0])  #name the variables
  resp[,"PA"] <- y.pp

  # Projector matrix from mesh to data.
  projmat <- INLA::inla.spde.make.A( mesh, as.matrix( raster::coordinates( observs))) # from mesh to point observations

  #make the stack
  stk.binom <- INLA::inla.stack( data=list( resp=resp, Ntrials=ntrials, offy=log( offy)),
                          A=list(1,projmat), tag=tag,
                          effects=list( covars, list(isdm.spat.XXX=1:mesh$n)))

#  attr( stk.binom, "newArtForm") <- newArtForm
  
  return( stk.binom)
}

