
###############################################################################################
###############################################################################################
####	Get the set of variables used in the model, from any artefact sub-component.
####	Also generates design matrices (or associated rasters)
####	
####	Returns a list with unique variable names in each sub-component,
####			 along with dataframes or rasters
####
####	Programmed by Scott in the first half of 2022 
####	Touched up by Scott in March/April 2023 (significantly!)
####	And again July 2024
####
###############################################################################################
###############################################################################################
uniqueVarNames <- function( obsList, covarBrick, distForm, biasForm, arteForm, habitatArea, DCsurvID, coord.names, responseNames, sampleAreaNames, stdCovs, na.action){

  ####	Distribution formula for entire region -- individual (artefact) datasets will be taken from this
  #Design matrix/raster for distribution
  tmpXX <- as.data.frame( terra::values( covarBrick))
  #get rid of intercept in distribution formula (and make sure of it)
  if( !is.null( distForm))
    distForm <- stats::update.formula( distForm, ~.-1+0)
  #make sure that all variables in dist form are available
  if( !all( all.vars( distForm) %in% colnames( tmpXX)))
    stop( "There is a variable in the distribution formula that is not present in the covariate raster brick.")
  #prune variables not needed in formulas. Both done to ensure NA pattern is consistent (hence basis expansion too)
  #make sure that biasForm is available and has the variables
  if( !all( all.vars( biasForm) %in% colnames( tmpXX))){
    warning( "Bias variables missing from covariates. Continuting without these variables. For calls to isdm(), this will(?!) cause a later error. For calls to predict.isdm(), this may cause an error but most likely it will mean that the basis for the distribution formula will be slightly different than the basis for the estimation -- be wary.")
    myVars <- all.vars( distForm)
  }
  else
    myVars <- unique( c( all.vars( distForm), all.vars( biasForm)))
  if( !is.null( habitatArea))
    myVars <- unique( c( myVars, habitatArea))
  tmpXX <- tmpXX[,myVars,drop=FALSE]
  #issue summary/warning if some covariates have a different NA pattern than others. Code stolen from isdm.model.matrix.
  allNAs.id <- apply( tmpXX, 1, function(zz) all( is.na(zz)))
  anyNAs.id <- apply( tmpXX, 1, function(zz) any( is.na(zz)))
  if( !all( allNAs.id == anyNAs.id))
    warning( paste0("Covariate rasters: Missing data (NA) present in a partial number of covariates for an observation (or cell). There are ",sum(anyNAs.id)-sum(allNAs.id)," such observations (cells). They are excluded from analysis."))
  #actually unifying the NAs
  tmpXX[anyNAs.id,] <- NA
  zeroEffort <- !is.na( tmpXX[,habitatArea]) & tmpXX[,habitatArea] == 0
  if( any( zeroEffort == TRUE))
    warning( paste0("Habitat area raster: There are ", sum( zeroEffort)," non-NA observations with zero (0) effort. They are excluded from analysis."))
  tmpXX[zeroEffort,] <- NA

  #the model frame for the distribution data
  if( !is.null( distForm)){
    XX <- isdm.model.matrix( formmy=distForm, obsy=tmpXX, namy=NULL, includeNA=TRUE) #includeNA to get raster pattern
    #make new formula based on expanded names
    newDistForm <- stats::reformulate( colnames( XX))
    newDistForm <- stats::update.formula( newDistForm, ~.-1)
    #put it in the 'correct' environment
    environment( newDistForm) <- environment( distForm)
    #assign the design matrix spatially
    newCovarBrick <- terra::rast( covarBrick[[1]], nlyrs=ncol( XX))
    terra::values( newCovarBrick) <- XX
    names( newCovarBrick) <- colnames( XX)
    #put it in the 'correct' environent
    environment( newCovarBrick) <- environment( covarBrick)
    if( stdCovs)
      terra::values( newCovarBrick) <- standardiseThatDesMatrix( terra::values( newCovarBrick))
  }
  else {
    newDistForm <- NULL
    newCovarBrick <- terra::rast( covarBrick[[1]], nlyrs=0)
    environment( newCovarBrick) <- environment( covarBrick)
  }
  
  ####	Bias formula, if present
  if( !is.null( biasForm)){
    #Design matrix/raster for bias
    XX <- isdm.model.matrix( formmy=biasForm, obsy=tmpXX, namy="PO", includeNA=TRUE)
    #make new formula
    newBiasForm <- stats::reformulate( colnames( XX))
    #put it in the 'correct' environment
    environment( newBiasForm) <- environment( biasForm)
    #assign the design matrix
    newBiasCovarBrick <- terra::rast( covarBrick[[1]], nlyrs=ncol(XX))
    terra::values( newBiasCovarBrick) <- XX
    names( newBiasCovarBrick) <- colnames( XX)  
    if( stdCovs)
      terra::values( newBiasCovarBrick) <- standardiseThatDesMatrix( terra::values( newBiasCovarBrick))
  
    ####  Combine distribution and bias	(there may be some duplicated data but it should have undergone the same standardisation in dist and in PO...
    newCovarBrick <- c( newCovarBrick, newBiasCovarBrick)
  }
  else
    newBiasForm <- NULL

  ####	HabitatArea variable too
  if( !is.null( habitatArea)){
    if( !is.null( newDistForm)){#terra::nlyr( newCovarBrick) > 0){
      newCovarBrick <- c( newCovarBrick, covarBrick[[habitatArea]]) #addLayer( newCovarBrick, covarBrick[[habitatArea]])
      names( newCovarBrick)[terra::nlyr( newCovarBrick)] <- habitatArea
    }
    else
      newCovarBrick <- covarBrick[[habitatArea]]
    tmp <- terra::values( newCovarBrick[[habitatArea]])
    tmp[anyNAs.id] <- NA
    tmp[tmp==0] <- NA
    terra::values( newCovarBrick)[,habitatArea] <- tmp
#    terra::values( newCovarBrick[[habitatArea]])[anyNAs.id] <- NA  #to match other variables.
#    terra::values( newCovarBrick[[habitatArea]])[terra::values( newCovarBrick[[habitatArea]])==0] <- NA  #to match other variables
  }
  #put it in the 'correct' environent
  environment( newCovarBrick) <- environment( covarBrick)

  #### Artefact data (not standardised)
  #container for the altered observation lists (names to match newForm)
  #note using the unexpanded covarBrick rather than newCovarBrick
  newObs <- lapply( obsList, function(xx) as.data.frame( cbind( xx, terra::extract( x=covarBrick, as.matrix( xx[,coord.names]), method='simple'))))
  #indicator for if there is that type of data
  ind <- c( PO=0, PA=0, AA=0, DC=0) #for PO, PA, AA, DC

  ####	Artefact formulas 'n' data
  #container for the altered artefact formuals
  newForm <- arteForm
  #cycle through each of the formulas (artefact)
  for( ii in names( arteForm)){
    ind[ii] <- 1
    dataname <- paste0( ii,'dat')
    #a design matrix for the survey data.  No scaling. No alteration of names
    XX <- isdm.model.matrix( formmy=newForm[[ii]], obsy=as.data.frame( newObs[[dataname]]), namy=ii, includeNA=TRUE)
    #make new formula
    newForm[[ii]] <- stats::reformulate( colnames( XX))
    #put it in the 'correct' environment
    environment( newForm[[ii]]) <- environment( arteForm[[ii]])
    #add sample areas
    if( !is.null( sampleAreaNames[ii])){
      XX <- cbind( newObs[[dataname]][,sampleAreaNames[ii]], XX)
      colnames( XX)[1] <- sampleAreaNames[ii]
    }
    #add outcome, if appropriate (not DC)
    if( length( responseNames[ii]) != 0){
      XX <- cbind( newObs[[dataname]][,responseNames[ii]], XX)
      colnames( XX)[1] <- responseNames[ii]
    }
    #convert to sf
    XX <- cbind( XX, as.matrix( newObs[[dataname]][,coord.names]))
    XX <- sf::st_as_sf( as.data.frame( XX), coords=coord.names)
    #rename variables to something that will parse
    names( XX) <- removeParsingChars( names( XX))
    #if there are variable in both the survey data *and* the distribution+bias variables
    for( jj in setdiff( colnames( XX), "geometry")){
      if( jj %in% paste0(ii,"_",names( covarBrick))){
	tmpID <- which( jj == paste0(ii,"_",names( covarBrick)))
	tmpCovarjj <- terra::extract( x=covarBrick[[tmpID]], y=as.matrix( newObs[[dataname]][,coord.names]))[,1]
      }
    }
    #add in the distribution variables -- data
    XXdist <- terra::extract( newCovarBrick, XX, ID=FALSE)
    #remove the bias terms from the data -- they will have "PO_" as a prefix.
    XXdist <- XXdist[,!grepl( "PO_", colnames( XXdist)),drop=FALSE]
    #combine
    XX <- cbind( XX, XXdist)

    #assign the design matrix
    newObs[[dataname]] <- XX
    #put it in the 'correct' environent
    environment( newObs[[dataname]]) <- environment( obsList[[dataname]])
    
    #a bit extra for DC data -- just to make it extra unique :-)  IS THIS REALLY NEEDED?
    if( ii == "DC")
      if( DCsurvID %in% attr( terms( newForm[[ii]]), "term.labels"))
        DCsurvID <- paste( DCsurvID,"DC",sep="_")
  }

  ####	Casting PO data too, if it is there
  if( "POdat" %in% names( newObs)){
    ind['PO'] <- 1
    if( !inherits( newObs$POdat, "sf"))
      newObs$POdat <- sf::st_multipoint( as.matrix( newObs$POdat[,coord.names]))
  }
   
  ###	the return object
  res <- list( obsList=newObs, covarBrick=newCovarBrick, arteForm=newForm, distForm=newDistForm, biasForm=newBiasForm, DCsurvID=DCsurvID, ind=ind)
}

