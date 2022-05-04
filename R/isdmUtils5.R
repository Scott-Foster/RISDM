
checkInput <- function( respNames, biasForm, arteForm, DCobsInfo, obsList){
  
  #checking the input for obvious inconsistencies
  
  if( all( is.null( respNames)))
    stop("No response specified.  Please do so through responseNames argument.")
  if( length( respNames)<1 | length( respNames)>3)
    stop("Please specify responseNames as a named vector (no response required for PO or DC, see DCobserverInfo).  See help for guidance.")
  
  if( length( arteForm)<1 | length( arteForm)>3)
    stop("Please specify artefactFormulas as a named list of length 3 or less.  See help for guidance.")
  
  if( length( DCobsInfo) != 4)
    stop("Please specify DCobserverInfo as a named list of exactly length 4.  See help for guidance.")
  
  if( all( is.null( obsList)))
    stop("No data specified.  Please do so through observationList argument.")
  
  #are the pieces available for PO?
  if( any( c( is.null( obsList$POdat), is.null( biasForm))))
    if ( !all( c( is.null( obsList$POdat), is.null( biasForm))))
      stop("To include PO data, you *must* supply observationList$POdat (argument) AND a biasFormula")
  #are the pieces available for PA?
  if( any( c( is.null( obsList$PAdat), is.na( respNames["PA"]), is.null( arteForm[["PA"]]))))
    if( !all( c( is.null( obsList$PAdat), is.na( respNames["PA"]), is.null( arteForm[["PA"]]))))
      stop("To include PA data, you *must* supply observationList$PAdat (argument) AND an entry in responseNames AND an entry in artefactFormulas")
  #are the pieces available for AA?
  if( any( c( is.null( obsList$AAdat), is.na( respNames["AA"]), is.null( arteForm[["AA"]]))))
    if( !all( c( is.null( obsList$AAdat), is.na( respNames["AA"]), is.null( arteForm[["AA"]]))))
      stop("To include AA data, you *must* supply observationList$AAdat (argument) AND an entry in responseNames AND an entry in artefactFormulas")
  #are the pieces available for DC?
  if( any( c( is.null( obsList$DCdat), is.null( arteForm[["DC"]])))){#, is.null( DCobsInfo)))){
    if( !all( c( is.null( obsList$DCdat), is.null( arteForm[["DC"]]))))#, is.null( DCobsInfo))))
      stop("To include DC data, you *must* supply observationList$DCdat (argument) AND an entry in artefactFormulas AND the information about observers (DCobserverInfo")
    else
      rm( DCobsInfo)
  }
  
  return( TRUE)  #if execution makes it this far
}

makeFormula <- function( dform, bform, aforms, addRE=TRUE) {
  #not going to use a combined intercept.  Rather let each data type have its own
  
  #try/do add the bias specific terms to each of the different data types.  Do so by including a nested effect within intercept types
  formmy  <- dform
  for( ii in c("PO","PA","AA","DC")){
    tmpForm <- NULL
    if( ii == "PO"){
      if( !is.null( bform))
        tmpForm <- update( bform, "~.-0-1")
    }
    else{
      if( !is.null( aforms[[ii]]))
        tmpForm <- update( aforms[[ii]], "~.-0-1")
    }
    if( !is.null( tmpForm)){
      addTerms <- as.character( tmpForm)
      addTerms <- addTerms[length(addTerms)]
      #formmy <- update( formmy, paste0("~.+Intercept.",ii,"/",addTerms))
#      if( ii != "DC")
        formmy <- update( formmy, paste0("~.+Intercept.",ii,"/(",addTerms,")"))
#      else
#        formmy <- update( formmy, paste0("~.+Intercept.",ii,":(",addTerms,")"))
    }
  }
  
  #   
  #   
  #   if( !is.null( bforms[[ii]])){
  #     tmpForm <- update( bforms[[ii]], "~.-0-1")
  #     addTerms <- as.character( tmpForm)
  #     addTerms <- addTerms[length(addTerms)]
  #     #formmy <- update( formmy, paste0("~.+Intercept.",ii,"/",addTerms))
  #     formmy <- update( formmy, paste0("~.+Intercept.",ii,"/(",addTerms,")"))
  #   }
  # }
  # #DC doesn't want/need an intercept at all.  Treat it as a special case
  # if( !is.null( bforms[["DC"]])){
  #   tmpForm <- update( bforms[["DC"]], "~.-0-1")
  #   addTerms <- as.character( tmpForm)
  #   addTerms <- addTerms[length(addTerms)]
  #   formmy <- update( formmy, paste0("~.+Intercept.DC:(",addTerms,")"))
  # }
  # 
  
  if( addRE)
    formmy <- update( formmy, "~ . + f(isdm.spat.XXX,model=my.spde)") #index is always assumed to be "isdm.spat.XXX" and spde model is always "my.spde"
  if( length( formmy)==3)
    warning( "Ignoring distributionFormula's outcome and replacing with 'resp' (data generated internally)" )
  formmy <- update( formmy, "resp~.")
  
  return( formmy)
}

getVarNames <- function( fo, bfo, foList){

  varNames <- NULL
  if( !is.null( fo))
    varNames <- attr( terms( fo), "term.labels")
  if( !is.null( bfo))
    varNames <- c( varNames, setdiff( as.character( attr( terms( bfo), "variables")), "list"))
  for( ii in names( foList))
    if( !is.null( foList[[ii]]))
      varNames <- c( varNames, setdiff( as.character( attr( terms( foList[[ii]]), "variables")), "list"))

  varNames <- unique( varNames)
  return( varNames)
}

createOffsets <- function( sampleAreaNames, observationList){
  if( is.null( sampleAreaNames)){
    sampleAreaNames <- rep("area1",3)
    names( sampleAreaNames) <- c("PA","AA","DC")
    if( is.null( observationList$PAdat))
      sampleAreaNames <- sampleAreaNames[names( sampleAreaNames) != "PA"]
    else
      observationList$PAdat$area1 <- 1
    
    if( is.null( observationList$AAdat))
      sampleAreaNames <- sampleAreaNames[names( sampleAreaNames) != "AA"]
    else
      observationList$AAdat$area1 <- 1
    
    if( is.null( observationList$DCdat))
      sampleAreaNames <- sampleAreaNames[names( sampleAreaNames) != "DC"]
    else
      observationList$DCdat$area1 <- 1
  }
  
  res <- list( names=sampleAreaNames, dat=observationList)
  
}

makeCombinedStack <- function( obsList, covarBrick, Mesh, contr, varNames, rasterVarNames, sampleAreaNames, responseNames) {
  
  #indicator for presence of data sets
  ind <- c( PO=0, PA=0, AA=0, DC=0) #for PO, PA, AA  Others to come later...
  
  if( !is.null( obsList$POdat)){
    ind["PO"] <- 1
    responseNames <- c( responseNames, PO="anybloodything")  #temporary
    tmp <- as.data.frame( cbind( obsList$POdat, raster::extract( x=covarBrick, obsList$POdat[,contr$coord.names])))
    obsList$POdat <- sp::SpatialPointsDataFrame( coords=tmp[,contr$coord.names], data=tmp)
  }
  if( !is.null( obsList$PAdat)){
    ind["PA"] <- 1
    tmp <- as.data.frame( cbind( obsList$PAdat, raster::extract( x=covarBrick, obsList$PAdat[,contr$coord.names])))
    obsList$PAdat <- sp::SpatialPointsDataFrame( coords=tmp[,contr$coord.names], data=tmp)
  }
  if( !is.null( obsList$AAdat)){
    ind["AA"] <- 1
    tmp <- as.data.frame( cbind( obsList$AAdat, raster::extract( x=covarBrick, obsList$AAdat[,contr$coord.names])))
    obsList$AAdat <- sp::SpatialPointsDataFrame( coords=tmp[,contr$coord.names], data=tmp)
  }
  if( !is.null( obsList$DCdat)){
    ind["DC"] <- 1
    tmp <- as.data.frame( cbind( obsList$DCdat, raster::extract( x=covarBrick, obsList$DCdat[,contr$coord.names])))
    obsList$DCdat <- sp::SpatialPointsDataFrame( coords=tmp[,contr$coord.names], data=tmp)
  }
  
  if( ! all( names( ind)[ind==1] %in% names( responseNames)))
    stop( "Data supplied with no response variable.  Please check responseNames argument and make sure it matches the data input.")
  if( ! all( names( responseNames) %in% names( ind)[ind==1]))
    stop( "Response variables named, but no data supplied.  Please check responseNames argument and make sure it matches data input.")
  
  nOutcomes <- sum( ind)
  if( nOutcomes==0)
    stop("Error: No species data supplied.")
  
  stck <- NULL
  for( ss in which( ind==1)){
    tmp <- switch( names( ind)[ss],  #obtuse way to get it, but it preserves the names
                   PO = MakePPPstack( observs=obsList$POdat, Mesh=Mesh, presname=contr$presence, ind=ind, 
                                      varNames=rasterVarNames),
                   PA = MakePAstack( observs=obsList$PAdat, mesh=Mesh$mesh, presname=responseNames["PA"], tag = "PA", 
                                     varNames=varNames, sampleAreaName=sampleAreaNames["PA"], ind=ind),
                   AA = MakeAAStack( observs=obsList$AAdat, mesh=Mesh$mesh, abundname=responseNames["AA"], tag="AA", 
                                     varNames=varNames, sampleAreaName=sampleAreaNames["AA"], ind=ind),
                   DC = MakeDCStack( observs=obsList$DCdat, mesh=Mesh$mesh, DCname=responseNames["DC"], tag="DC", 
                                     varNames=varNames, sampleAreaName=sampleAreaNames["DC"], ind=ind))
    if( is.null( stck))
      stck <- tmp
    else
      stck <- INLA::inla.stack( stck, tmp)
  }
  #  pred.stck <- MakePredictionStack( mesh=Mesh$mesh, covar_raster_data=covarBrick, tag='pred', varNames=rasterVarNames, ind=ind)
  #  stck <- INLA::inla.stack( stck, pred.stck$stk)
  
  attr( stck, "ind") <- ind
  #  attr( stck, "predcoords") <- pred.stck$predcoords
  
  return( stck)
}

makeFamLink <- function( ind) {

  #the names (datatype) contained within the data
  nammy <- names( ind)

  #vectors for the families and links
  fammy <- rep( NA, length( nammy))
  names( fammy) <- nammy
  linkky <- list()

  if( ind["PO"]==1){
    fammy["PO"] <- "poisson"
    linkky[["PO"]] <- list(link="log")
  }
  if( ind["PA"]==1){
    fammy["PA"] <- "binomial"
    linkky[["PA"]] <- list(link="cloglog")
  }
  if( ind["AA"]==1){
    fammy["AA"] <- "poisson"
    linkky[["AA"]] <- list(link="log")
  }
  if( ind["DC"]==1){
    fammy["DC"] <- "poisson"
    linkky[["DC"]] <- list(link="log")
  }
  #remove the NA entries
  fammy <- fammy[!is.na( fammy)]

  linkID <- NA  #the index of a log link, if PO AA data present.  Otherwise its gotta be cloglog.
  if( is.na( linkID) & ( "PO" %in% names( fammy)))
    linkID <- which( names( fammy) == "PO")
  if( is.na( linkID) & ( "AA" %in% names( fammy)))
    linkID <- which( names( fammy) == "AA")
  if( is.na( linkID) & ( "PA" %in% names( fammy)))
    linkID <- which( names( fammy) == "PA")
  if( is.na( linkID) & ( "DC" %in% names( fammy)))
    linkID <- which( names( fammy) == "DC")
  #  if( length( fammy)==1)
  #    linkky <- linkky[[1]]  #remove a level of the listing for INLA call.

  names( linkky) <- NULL  #to satisfy one of inla's internal checks.  Sigh...

  res <- list( fam=fammy, link=linkky, linkID=linkID)

  return( res)
}

setPriors <- function( control, stck){
  tmp <- list( mean=control$prior.mean)
  tmp$prec <- list( default=control$other.prec)
  if( "Intercept.PO" %in% colnames( stck$effects$data))
    tmp$prec$Intercept.PO <- control$int.prec
  if( "Intercept.PA" %in% colnames( stck$effects$data))
    tmp$prec$Intercept.PA <- control$int.prec
  if( "Intercept.AA" %in% colnames( stck$effects$data))
    tmp$prec$Intercept.AA <- control$int.prec
  if( "Intercept.DC" %in% colnames( stck$effects$data))
    tmp$prec$Intercept.DC <- control$int.prec
  return( tmp)
}

makeControl <- function( contr) {

  if( ! "n.threads" %in% names( contr))  #number of cores to use.  Default is greedy but not super-super greedy
    contr$n.threads <- parallel::detectCores()-1
  if( ! "coord.names" %in% names( contr))
    contr$coord.names <- c("Easting","Northing")
  if( ! "prior.mean" %in% names( contr))
    contr$prior.mean <- 0
  if( ! "int.prec" %in% names( contr))
    contr$int.prec <- 0.0001
  if( ! "other.prec" %in% names( contr))
    contr$other.prec <- 0.0001
  if( ! "calcICs" %in% names( contr))
    contr$calcICs <- FALSE
  if( !"prior.range" %in% names( contr))
    contr$prior.range <- c(0.2, 0.1)  #suitable for a unit square geometry!
  if( !"prior.space.sigma" %in% names( contr))
    contr$prior.space.sigma <- c(1, 0.01)
  if( !"verbose" %in% names( contr))
    contr$verbose <- FALSE
  if( !"addRandom" %in% names( contr))
    contr$addRandom <- TRUE
  if( !"returnStack" %in% names( contr))
    contr$returnStack <- FALSE
  if( !"DCmethod" %in% names( contr))
    contr$DCmethod <- "plugin"  #the other option is "TaylorsLinApprox"
  
  #add to as we go along
  if( !all( names( contr) %in% c("n.threads","tag.pred","spat.index", "coord.names", "verbose",
                                 "prior.mean","int.prec","other.prec", "calcICs", "prior.range", 
                                 "prior.space.sigma", "addRandom", "returnStack", "DCmethod")))
    warning( "There are control parameters specified that are not used.")
  return( contr)
}

MakeDCStack=function( observs, mesh = NULL, DCname = "DCabund", tag = "DC", varNames, sampleAreaName, ind) {
  
  #casting to numeric
  y.pp <- as.integer( observs@data[,DCname])
  
  #the area sampled
  offy <- observs@data[,sampleAreaName]
  
  #covariates 'n' stuff
  tmp <- varNames %in% colnames( observs@data)
  if( ! all( tmp))
    warning( "Not all bias covariates in DC data. Missing: ", paste(varNames[!tmp], sep=" "), "Creating variable and padding with NAs.")
  tmptmp <- matrix( NA, nrow=nrow( observs@data), ncol=sum( !tmp), dimnames=list( NULL, varNames[!tmp]))
  observs@data <- cbind( observs@data, as.data.frame( tmptmp))
  
  covars <- as.data.frame( observs@data[,varNames])
  
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

MakeAAStack=function( observs, mesh = NULL, abundname = "abund", tag = "AA", varNames, sampleAreaName, ind) {

  # Function to create stack for ABUNDANCE absence points
  
  # observs SpatialPointsDataFrame of covariates and observations. Must contain a column corresponding to abundName, which must be non-negative integer (0, 1, 2, ...)
  # mesh INLA mesh for the SPDE model
  # presname Name of abundance column in observs.
  # offsetname Name of the offset variable for AA data.
  # tag Name for tag for the stack.
  
  # An INLA stack with binomial data: include Ntrials, which is the number of trials
  
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
  if( ! all( tmp))
    warning( "Not all bias covariates in PA data. Missing: ", paste( varNames[!tmp], " "), "Creating variable and padding with NAs.")
  tmptmp <- matrix( NA, nrow=nrow( observs@data), ncol=sum( !tmp), dimnames=list( NULL, varNames[!tmp]))
  observs@data <- cbind( observs@data, as.data.frame( tmptmp))
  
  covars <- as.data.frame( observs@data[,varNames])

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

MakePPPstack <- function( observs, Mesh, presname="presence", tag="PO", ind=1, varNames) {

  #observs is a SpatialPointsDataFrame containing a row for each presence and quadrature point. Data must contain all covariates needed
  #for the model.
  #mesh is the mesh for spde model containing mesh object, the spde object (not used here though), the area of the in-region Voronoi polygons, and the covariate data
  #                 note that the names of the covariate data MUST match the names within observs.
  #tag is the tag for the INLA stack.
  #ind is a 0,1 vector indicating which outcomes/likelihoods will be needed.  Default will produce a single outcome stack.
  
  #following Simpson et al 2016 (off the grid), as described in https://becarioprecario.bitbucket.io/spde-gitbook/ch-lcox.html
  #originally following Isaac's et al code, but they don't follow the same recipe quite...  This seems much more in line with Simpson et al,
  #as it incorporates mesh nodes (background points) into the observation model too (doing the point process integration).
  
  #Note that the tesselation areas (e.pp from Mesh$w are calcualted using a different method (deldir))
  #This shouldn't matter, and is faster, and based on Voronoi polygons rather than GIS polygons 
  #   (see above website Fig 4.3 for some probably effective but odd shapes...)

  #the working variate
  y.pp <- c( rep( 0, Mesh$mesh$n), rep( 1, nrow( observs)))

  #its weights
  e.pp <- c( Mesh$w, rep( 0, nrow( observs)))

  #projection matrices
  #mesh nodes to mesh nodes
  imat <- Matrix::Diagonal( Mesh$mesh$n, rep( 1, Mesh$mesh$n))
  #observations to mesh nodes
  lmat <- INLA::inla.spde.make.A( Mesh$mesh, as.matrix( raster::coordinates( observs)))
  #entire projection
  A.pp <- rbind( imat, lmat)

  #covariates 'n' stuff
  covars <- as.data.frame( rbind( Mesh$covars, as.data.frame( observs)[,varNames]))

  #each data type gets its own intercept
  covars$Intercept.PO <- 1

  resp <- matrix( NA, nrow=nrow( covars), ncol=sum( ind))
  colnames( resp) <- names( ind[ind!=0])  #name the variables
  resp[,"PO"] <- y.pp

  #make the stack
  stk.pp <- INLA::inla.stack( data = list( resp = resp, e = e.pp, offy=0),
                        A = list(1, A.pp),
                        effects = list(covars, list(isdm.spat.XXX = 1:Mesh$mesh$n)),
                        tag = tag)
  return( stk.pp)

}

MakeSpatialRegion <- function (data = NULL, coords = c("X", "Y"), mesh, bdry, proj = sp::CRS("+proj=utm"), dataBrick=NULL, varNames=NULL) {

  #largely taken from Isaacs et al.  But large chunks removed.
  #Note that the approach to calculate the weights, and do the quadrature is quite different
  #Here, we separate the SPDE mesh and the quadrature points (unlike Isaacs et al). We calculate the weights using ppmData.
  
#  #make the boundary into something that INLA will understand.
#  region.bdry <- INLA::inla.sp2segment( raster::union( bdry))
#  #make the mesh, may take a while if lots of points and/or mesh configuration that is not well suited to bdry
#  mesh <- INLA::inla.mesh.2d(boundary = region.bdry, cutoff = meshpars$cutoff, max.edge = meshpars$max.edge, offset = meshpars$offset, max.n = meshpars$max.n)
  #set up the parameters for the spatial random field.

  #finding the areas of the Voronoi polygons around each of the mesh points.
  dd <- deldir::deldir(mesh$loc[, 1], mesh$loc[, 2], eps=1e-6)
  tiles <- deldir::tile.list(dd)  #the polygons

  #convert the bdry to a single dissolved polygon
  GPC.polies <- list()
  lenny <- length( bdry@polygons[[1]]@Polygons)
  for( ii in 1:lenny)
    GPC.polies[[ii]] <- as(bdry@polygons[[1]]@Polygons[[ii]]@coords, "gpc.poly")
  if( lenny > 1){
    warning( "currently on first to polygon boundary segments are used. FIX IN MakeSpatialRegion")
    poly.gpc <- rgeos::union( GPC.polies[[1]], GPC.polies[[2]]) #warning
  }
  else
    poly.gpc <- GPC.polies[[1]]

  #Figure out how much of each tile is within bdry.
  #This could be parallelised.
  suppressWarnings( w <- sapply(tiles, function(p) rgeos::area.poly(rgeos::intersect(as(cbind(p$x, p$y), "gpc.poly"), poly.gpc))))

  #find the covariate values under the mesh
  covars <- NULL
  if( !is.null( dataBrick)){
    if( ! all( varNames %in% names( dataBrick)))
      stop( "Covariates in formula not supplied in raster data")
    tmp <- ExtractCovarsAtNodes( mesh=mesh, covars=dataBrick)
    covars <- tmp[,varNames]
  }

  return(list(mesh = mesh, w=w, covars=covars))
}

ExtractCovarsAtNodes <- function( mesh=NULL, covars=NULL){
  
  #mesh is the INLA mesh for the region
  #covars is a raster/brick containing all the covariate values.  Mush contain spatial area for all the mesh nodes
  
  locs <- as.matrix( mesh$loc)
  covarAtLocs <- raster::extract( covars, locs[,1:2])

  return( covarAtLocs)
}

