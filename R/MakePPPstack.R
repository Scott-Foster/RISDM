
###############################################################################################
###############################################################################################
####	Set up the INLA stacks for PO data (point process)
####	Assumes that there has already been a fair bit of processing to get things into the
####	the right, consistent formats.
####	Uses the on-grid approximation.
####
####	Returns an INLA stack object for PPP data.with correctly set-up multiple responses etc...
####
####	Programmed by Scott in the first half of 2022
####	Genesis of code started with Isaac's et al. (2020) R function but significiant modifactions
####
###############################################################################################
###############################################################################################

MakePPPstack <- function( observs, covarBrick, Mesh, presname="presence", tag="PO", varNames, ind) {

  ##update 17-May-22  Going on the grid -- increase robustness (hopefully) with little increase in computation (for coarse enough grid)

  #create a raster to perform approximation and analysis on.  Assumed to be the same as the covariate brick.
  tmp <- raster::rasterize( x=coordinates( observs), y=covarBrick, background=0, fun='count')
  tmp <- raster::mask( tmp, covarBrick[[1]])
  #the cell areas will act as an offset for the PPP
  tmp <- raster::addLayer( tmp, raster::area(tmp))#raster( terra::cellSize( terra::rast( tmp))))
  names( tmp) <- c( "count", "cellArea")

  #change to spatial data for some ease.
  observs <- SpatialPointsDataFrame( coords=coordinates( tmp), data=cbind( raster::as.data.frame( tmp), raster::as.data.frame( covarBrick)))
  abundname <- 'count'
  sampleAreaName <- "cellArea"
  
  ### This is copied from makeAAStack.  Changes here or there will need to be transfered
  y.pp <- as.integer( observs@data[,abundname])
  
  #the area sampled
  offy <- observs@data[,sampleAreaName]
  
  #covariates 'n' stuff
  tmp <- varNames %in% colnames( observs@data)
#  if( ! all( tmp))
#    warning( "Not all bias covariates in PO data. Missing: ", paste( varNames[!tmp], " "), "Creating variable and padding with NAs.")
  tmptmp <- matrix( NA, nrow=nrow( observs@data), ncol=sum( !tmp), dimnames=list( NULL, varNames[!tmp]))
  observs@data <- cbind( observs@data, as.data.frame( tmptmp))
  
  #Convert back to non-spatial
  covars <- as.data.frame( observs@data[,varNames])
  colnames( covars) <- varNames

  #add intercept, whose name is not prespecified.
  covars$Intercept.PO <- 1

  #Getting the response ready
  resp <- matrix( NA, nrow=nrow( covars), ncol=sum( ind))
  colnames( resp) <- names( ind[ind!=0])  #name the variables
  resp[,"PO"] <- y.pp
  
  # Projector matrix from mesh to data.
  projmat <- INLA::inla.spde.make.A( Mesh$mesh, as.matrix( raster::coordinates( observs))) # from mesh to point observations
  
  #make the stack
  stk.po <- INLA::inla.stack(data=list(resp=resp, e=rep( 1, nrow( covars)), offy=log( offy)),
                                A=list(1,projmat), tag=tag,
                                effects=list( covars, list(isdm.spat.XXX=1:Mesh$mesh$n)))

  return( stk.po)

}
