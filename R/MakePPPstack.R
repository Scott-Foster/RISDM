
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

MakePPPstack <- function( observs, covarBrick, habitatArea, Mesh, presname="presence", tag="PO", ind){#varNames, ind) {

  ##update 17-May-22  Going on the grid -- increase robustness (hopefully) with little increase in computation (for coarse enough grid)

  #create a raster to perform approximation and analysis on.  Assumed to be the same as the covariate brick.
  tmp <- terra::rasterize( x=sf::st_coordinates( observs)[,1:2], y=covarBrick, background=0, fun='count')
  tmp <- terra::mask( tmp, covarBrick[[1]])
  #the cell areas will act as an offset for the PPP
  if( !is.null( habitatArea)){
    tmp <- c(tmp, covarBrick[[habitatArea]])#raster::addLayer( tmp, covarBrick[[habitatArea]])
    if( any( unique( terra::values( tmp[[1]]>0 & !(tmp[[2]] > 0))), na.rm=TRUE)){
      xxz <- sum( unique( terra::values( tmp[[1]]>0 & !(tmp[[2]] > 0))), na.rm=TRUE)
      warning( paste0( "Presences occur in cells with zero (or NA) habitat area. Please check POdat and habitatArea arguments.\n For now, these observations will be removed. There are ",xxz," of them."))
      terra::values( tmp[[1]])[terra::values( tmp[[1]]>0 & !(tmp[[2]] > 0))] <- 0
    }
  }
  else{
    if( terra::is.lonlat( tmp))
      tmp <- c( tmp, terra::cellSize(tmp))
    else{
      tmp1 <- prod( terra::res( tmp))
      tmp2 <- tmp[[1]]
      terra::values( tmp2) <- tmp1
      tmp <- c( tmp, tmp2)
      tmp <- terra::mask( tmp, covarBrick[[1]])
    }
  }  
  names( tmp) <- c( "count", "cellArea")

  #change to spatial data for some ease.
#  observs <- SpatialPointsDataFrame( coords=terra::crds( tmp), data=cbind( ter::as.data.frame( tmp), raster::as.data.frame( covarBrick)))
  #or non-spatial df
  observs <- cbind( terra::as.data.frame( tmp, na.rm=FALSE), terra::as.data.frame( covarBrick, na.rm=FALSE), terra::crds( tmp, na.rm=FALSE))

  abundname <- 'count'
  sampleAreaName <- "cellArea"
  #remove the cells that are not within the study area (masked out)
  observs <- observs[!is.na( observs$count),]  
  #remove the cells that have zero habitat
  observs <- observs[!is.na( observs[,sampleAreaName]) & observs[,sampleAreaName]>0,]
  
  ### This is copied from makeAAStack.  Changes here or there will need to be transfered
  y.pp <- as.integer( observs[,abundname])
  
  #the area sampled
  offy <- observs[,sampleAreaName]
  
#  #covariates 'n' stuff
#  tmp <- varNames %in% colnames( observs@data)
##  if( ! all( tmp))
##    warning( "Not all bias covariates in PO data. Missing: ", paste( varNames[!tmp], " "), "Creating variable and padding with NAs.")
#  tmptmp <- matrix( NA, nrow=nrow( observs@data), ncol=sum( !tmp), dimnames=list( NULL, varNames[!tmp]))
#  observs@data <- cbind( observs@data, as.data.frame( tmptmp))
  
  #Convert back to non-spatial
#  as.data.frame( observs@data[,varNames])
#  colnames( covars) <- varNames

#  #add intercept, whose name is not prespecified.
#  observs$PO_Intercept <- 1
#  covars$Intercept.PO <- 1

  #Getting the response ready
  resp <- matrix( NA, nrow=nrow( observs), ncol=sum( ind))
  colnames( resp) <- names( ind[ind!=0])  #name the variables
  resp[,"PO"] <- y.pp
  
  # Projector matrix from mesh to data.
  projmat <- INLA::inla.spde.make.A( Mesh$mesh, as.matrix( observs[,ncol( observs)- (1:0)])) # from mesh to point observations
  
  #make the stack
  stk.po <- INLA::inla.stack(data=list(resp=resp, e=rep( 1, nrow( observs)), offy=log( offy)),
                                A=list(1,projmat), tag=tag,
                                effects=list( observs[,-c( 1:2, ncol(observs)-(1:0))], list(isdm.spat.XXX=1:Mesh$mesh$n)))

  return( stk.po)

}
