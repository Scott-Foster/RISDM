
###############################################################################################
###############################################################################################
####	Set up the INLA stacks for an isdm.
####	Assumes that there has already been a fair bit of processing to get things into the
####	the right, consistent formats.
####
####	Returns an INLA stack object with correctly set-up multiple responses etc...
####
####	Programmed by Scott in the first half of 2022
####	Genesis of code started with the corresponding R functions of
####		Dambly et al. (2019) https://doi.org/10.5281/zenodo.3363936 
####	Re-jigged in March 2023 to make explicit for INLA the particular pattern
####		of constraints for the factor level effects
####
###############################################################################################
###############################################################################################

makeCombinedStack <- function( obsList, covarBrick, habitatArea, distForm, artForms, Mesh, contr, responseNames, ind, sampleAreaNames) {#varNames, 
  
#  #indicator for presence of data sets
#  ind <- c( PO=0, PA=0, AA=0, DC=0) #for PO, PA, AA  Others to come later...
  
  #convert each of the observationList elements into a spatial points data.frame.  Can be easier to deal with.
  if( !is.null( obsList$POdat)){
#    ind["PO"] <- 1
    responseNames <- c( responseNames, PO="anybloodything")  #temporary
#    tmp <- as.data.frame( cbind( obsList$POdat, raster::extract( x=covarBrick, obsList$POdat[,contr$coord.names])))
#    obsList$POdat <- sp::SpatialPointsDataFrame( coords=tmp[,contr$coord.names], data=tmp)
  }
#  if( !is.null( obsList$PAdat)){
#    ind["PA"] <- 1
#    tmp <- as.data.frame( cbind( obsList$PAdat, raster::extract( x=covarBrick, obsList$PAdat[,contr$coord.names])))
#    obsList$PAdat <- sp::SpatialPointsDataFrame( coords=tmp[,contr$coord.names], data=tmp)
#  }
#  if( !is.null( obsList$AAdat)){
#    ind["AA"] <- 1
#    tmp <- as.data.frame( cbind( obsList$AAdat, raster::extract( x=covarBrick, obsList$AAdat[,contr$coord.names])))
#    obsList$AAdat <- sp::SpatialPointsDataFrame( coords=tmp[,contr$coord.names], data=tmp)
#  }
#  if( !is.null( obsList$DCdat)){
#    ind["DC"] <- 1
#    tmp <- as.data.frame( cbind( obsList$DCdat, raster::extract( x=covarBrick, obsList$DCdat[,contr$coord.names])))
#    obsList$DCdat <- sp::SpatialPointsDataFrame( coords=tmp[,contr$coord.names], data=tmp)
#  }

  #elementrary check to see if there are necessary components
  if( ! all( names( ind)[ind==1] %in% names( responseNames)))
    stop( "Data supplied with no response variable.  Please check responseNames argument and make sure it matches the data input.")
  if( ! all( names( responseNames) %in% names( ind)[ind==1]))
    stop( "Response variables named, but no data supplied.  Please check responseNames argument and make sure it matches data input.")
  
  #more elementary checks
  nOutcomes <- sum( ind)
  if( nOutcomes==0)
    stop("Error: No species data supplied.")
    
  #check if habitat area is in the covarBrick
  if( !is.null( habitatArea))
    if( !habitatArea %in% names( covarBrick))
      stop( "habitatArea name supplied but not present in covarBrick. Please check.")
  
  stck <- NULL
  CombinedArtForms <- list()
  #create each of the component stacks and then add them to a combined stack, sequentially.
  for( ss in which( ind==1)){
    tmp <- switch( names( ind)[ss],  #obtuse way to get it, but it preserves the names
                   PO = MakePPPstack( observs=obsList$POdat, covarBrick=covarBrick, habitatArea=habitatArea, Mesh=Mesh, presname=contr$presence, ind=ind),
                   PA = MakePAstack( observs=obsList$PAdat, mesh=Mesh$mesh, presname=responseNames["PA"], tag = "PA", sampleAreaName=sampleAreaNames["PA"], ind=ind),
		   AA = MakeAAstack( observs=obsList$AAdat, mesh=Mesh$mesh, abundname=responseNames["AA"], tag="AA", sampleAreaName=sampleAreaNames["AA"], ind=ind),
                   DC = MakeDCstack( observs=obsList$DCdat, mesh=Mesh$mesh, DCname=responseNames["DC"], tag="DC", sampleAreaName=sampleAreaNames["DC"], ind=ind))
    if( is.null( stck))
      stck <- tmp
    else
      stck <- INLA::inla.stack( stck, tmp)
    
#    #building pieces for the formula
#    if( names( ind)[ss] != "PO")
#      CombinedArtForms[[names( ind)[ss]]] <- attr( tmp, "newArtForm")
  }
  
  #add another thing indicating which data is present.
  attr( stck, "ind") <- ind
  
  #add another thing giving the number of each observation
  tmp <- apply( stck$data$data[,1:nOutcomes,drop=FALSE], 2, function(xx) sum( !is.na( xx)))
  attr( stck, "nObs") <- tmp
  
#  #add the updated artefact formulae
#  attr( stck, "newArtForms") <- CombinedArtForms
  
  return( stck)
}
