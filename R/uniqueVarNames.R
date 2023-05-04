
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
####
###############################################################################################
###############################################################################################

uniqueVarNames <- function( obsList, covarBrick, distForm, biasForm, arteForm, habitatArea, DCsurvID, coord.names, responseNames, sampleAreaNames, stdCovs){

  ####	Distribution formula for entire region -- individual (artefact) datasets will be taken from this
  #Design matrix/raster for distribution
  tmpXX <- as.data.frame( raster::values( covarBrick))
  #the model frame for the data
  XX <- isdm.model.matrix( formmy=distForm, obsy=tmpXX, namy=NULL)
  #make new formula
  newDistForm <- stats::reformulate( colnames( XX))
  newDistForm <- stats::update.formula( newDistForm, ~.-1)
  #put it in the 'correct' environment
  environment( newDistForm) <- environment( distForm)
  #assign the design matrix spatially
  newCovarBrick <- raster::brick( covarBrick[[rep(1,ncol( XX))]])
  raster::values( newCovarBrick) <- XX
  #put it in the 'correct' environent
  environment( newCovarBrick) <- environment( covarBrick)
  if( stdCovs)
    raster::values( newCovarBrick) <- standardiseThatDesMatrix( raster::values( newCovarBrick))

  ####	Bias formula, if present
  if( !is.null( biasForm)){
  #Design matrix/raster for distribution
    XX <- isdm.model.matrix( formmy=biasForm, obsy=as.data.frame( raster::values( covarBrick)), namy="PO")
    #make new formula
    newBiasForm <- stats::reformulate( colnames( XX))
    #put it in the 'correct' environment
    environment( newBiasForm) <- environment( biasForm)
    #assign the design matrix
    newBiasCovarBrick <- covarBrick[[rep(1,ncol( XX))]]
    raster::values( newBiasCovarBrick) <- XX
    names( newBiasCovarBrick) <- colnames( XX)  
    if( stdCovs)
      raster::values( newBiasCovarBrick) <- standardiseThatDesMatrix( raster::values( newBiasCovarBrick))
  
    ####	Combine distribution and bias
    tmpNames <- c( names( newCovarBrick), names( newBiasCovarBrick))
    raster::values( newCovarBrick) <- cbind( raster::values( newCovarBrick), raster::values( newBiasCovarBrick))
    names( newCovarBrick) <- tmpNames
  }
  else
    newBiasForm <- NULL
  
  ####	HabitatArea variable too
  if( !is.null( habitatArea))
    newCovarBrick <- addLayer( newCovarBrick, covarBrick[[habitatArea]])
  #put it in the 'correct' environent
  environment( newCovarBrick) <- environment( covarBrick)

  #### Artefact data (not standardised)
  #container for the altered observation lists (names to match newForm)
  #note using the unexpanded covarBrick rather than newCovarBrick
  newObs <- lapply( obsList, function(xx) as.data.frame( cbind( xx, raster::extract( x=covarBrick, xx[,coord.names]))))
  #indicator for if there is that type of data
  ind <- c( PO=0, PA=0, AA=0, DC=0) #for PO, PA, AA, DC

  ####	Artefact formulas 'n' data
  #container for the altered artefact formuals
  newForm <- arteForm
  #cycle through each of the formulas (artefact)
  for( ii in names( newForm)){
    ind[ii] <- 1
    dataname <- paste0( ii,'dat')
    #a design matrix for the survey data.  No scaling. No alteration of names
    XX <- model.matrix( newForm[[ii]], newObs[[dataname]])
    if( "(Intercept)" %in% colnames( XX))
      colnames( XX)["(Intercept)" == colnames( XX)] <- paste0( ii,"_Intercept")
    #make sure other variable names are also unique
    colnames( XX)[ !grepl( "_Intercept", colnames( XX))] <- paste0(ii,"_",colnames( XX)[!grepl( "_Intercept",colnames( XX))])
    #remove special characters for parsing in inla()
    colnames( XX) <- removeParsingChars( colnames( XX))
    #make new formula
    newForm[[ii]] <- stats::reformulate( colnames( XX))
    #put it in the 'correct' environment
    environment( newForm[[ii]]) <- environment( arteForm[[ii]])
    #add sample areas
    XX <- cbind( newObs[[dataname]][,sampleAreaNames[ii]], XX)
    colnames( XX)[1] <- sampleAreaNames[ii]
    #add outcome, if appropriate (not DC)
    if( length( responseNames[ii]) != 0){
      XX <- cbind( newObs[[dataname]][,responseNames[ii]], XX)
      colnames( XX)[1] <- responseNames[ii]
    }
    #convert to SpatialPoints*
    XX <- sp::SpatialPointsDataFrame( coords=newObs[[dataname]][,coord.names], data=as.data.frame( XX))
    #rename variables to something that will parse
    names( XX) <- removeParsingChars( names( XX))
    #if there are variable in both the survey data *and* the distribution+bias variables
    for( jj in colnames( XX)){
      if( jj %in% paste0(ii,"_",names( covarBrick))){
	tmpID <- which( jj %in% paste0(ii,"_",names( covarBrick)))
	XX[,jj] <- raster::extract( x=newCovarBrick[[tmpID]], coords=newObs[[dataname]][,coord.names])
      }
    }
    #add in the distribution variables -- data
    XXdist <- raster::extract( newCovarBrick, XX)
    #remove the bias terms from the data -- they will have "PO_" as a prefix.
    XXdist <- XXdist[,!grepl( "PO_", colnames( XXdist))]
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
    newObs$POdat <- sp::SpatialPoints( coords=newObs$POdat[,coord.names])
  }
   
  ###	the return object
  res <- list( obsList=newObs, covarBrick=newCovarBrick, arteForm=newForm, distForm=newDistForm, biasForm=newBiasForm, DCsurvID=DCsurvID, ind=ind)
}

