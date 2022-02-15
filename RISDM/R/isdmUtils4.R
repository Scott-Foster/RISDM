getVarNames <- function( fo, foList){

  varNames <- NULL
  if( !is.null( fo))
    varNames <- attr( terms( fo), "term.labels")
  for( ii in names( foList))
    if( !is.null( foList[[ii]]))
      varNames <- c( varNames, attr( terms( foList[[ii]]), "term.labels"))

    varNames <- unique( varNames)

    return( varNames)
}

makeFormula <- function( dform, bforms) {
  #try/do add the bias specific terms to each of the different data types.  Do so by including a nested effect within intercept types
  formmy  <- dform
  for( ii in c("PO","PA","AA")){
    if( !is.null( bforms[[ii]])){
      tmpForm <- update( bforms[[ii]], "~.-0-1")
      addTerms <- as.character( tmpForm)
      addTerms <- addTerms[length(addTerms)]
      #formmy <- update( formmy, paste0("~.+Intercept.",ii,"/",addTerms))
      formmy <- update( formmy, paste0("~.+Intercept.",ii,"/(",addTerms,")"))
    }
  }
  #not going to use a combined intercept.  Rather let each data type have its own
  #formmy <- update( formmy, "~ Intercept + .")
  formmy <- update( formmy, "~ . + f(i,model=my.spde)") #index is always assumed to be "i" and spde model is always "my.spde"
  if( length( formmy)==3)
    warning( "Ignoring distributionFormula's outcome and replacing with 'resp' (data generated internally)" )
  formmy <- update( formmy, "resp~.")

  return( formmy)
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
  #remove the NA entries
  fammy <- fammy[!is.na( fammy)]

  linkID <- NA  #the index of a log link, if PO AA data present.  Otherwise its gotta be cloglog.
  if( is.na( linkID) & ( "PO" %in% names( fammy)))
    linkID <- which( names( fammy) == "PO")
  if( is.na( linkID) & ( "AA" %in% names( fammy)))
    linkID <- which( names( fammy) == "AA")
  if( is.na( linkID) & ( "PA" %in% names( fammy)))
    linkID <- which( names( fammy) == "PA")

  #  if( length( fammy)==1)
  #    linkky <- linkky[[1]]  #remove a level of the listing for INLA call.

  names( linkky) <- NULL  #to satisfy one of inla's internal checks.  Sigh...

  res <- list( fam=fammy, link=linkky, linkID=linkID)

  return( res)
}

makeCombinedStack <- function( POdat, PAdat, AAdat, covarBrick, Mesh, contr, varNames, rasterVarNames, sampleAreaNames, responseNames) {

  #indicator for presence of data sets
  ind <- c( PO=0, PA=0, AA=0) #for PO, PA, AA  Others to come later...

  if( !is.null( POdat)){
    ind["PO"] <- 1
    tmp <- as.data.frame( cbind( POdat, raster::extract( x=covarBrick, POdat[,contr$coord.names])))
#    nonRasterVars <- setdiff( varNames, rasterVarNames)
#    if( length( nonRasterVars) > 0){
#      tmpDF <- as.data.frame( matrix( NA, ncol=length( nonRasterVars), nrow=nrow( tmp)))
#      colnames( tmpDF) <- nonRasterVars
#      tmp <- cbind( tmp, tmpDF)
#    }
    POdat <- sp::SpatialPointsDataFrame( coords=tmp[,contr$coord.names], data=tmp)
  }
  if( !is.null( PAdat)){
    ind["PA"] <- 1
    tmp <- as.data.frame( cbind( PAdat, raster::extract( x=covarBrick, PAdat[,contr$coord.names])))
    PAdat <- sp::SpatialPointsDataFrame( coords=tmp[,contr$coord.names], data=tmp)
  }
  if( !is.null( AAdat)){
    ind["AA"] <- 1
    tmp <- as.data.frame( cbind( AAdat, raster::extract( x=covarBrick, AAdat[,contr$coord.names])))
    AAdat <- sp::SpatialPointsDataFrame( coords=tmp[,contr$coord.names], data=tmp)
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
                   PO = MakePPPstack( observs=POdat, Mesh=Mesh, presname=contr$presence, ind=ind, varNames=rasterVarNames),
                   PA = MakePAstack( observs=PAdat, mesh=Mesh$mesh, presname=responseNames["PA"], tag = "PA", varNames=varNames, sampleAreaName=sampleAreaNames["PA"], ind=ind),
                   AA = MakeAAStack( observs=AAdat, mesh=Mesh$mesh, abundname=responseNames["AA"], tag="AA", varNames=varNames, sampleAreaName=sampleAreaNames["AA"], ind=ind))
    if( is.null( stck))
      stck <- tmp
    else
      stck <- INLA::inla.stack( stck, tmp)
  }
  pred.stck <- MakePredictionStack( mesh=Mesh$mesh, covar_raster_data=covarBrick, tag='pred', varNames=rasterVarNames, ind=ind)
  stck <- INLA::inla.stack( stck, pred.stck$stk)

  attr( stck, "ind") <- ind
  attr( stck, "predcoords") <- pred.stck$predcoords

  return( stck)
}

makeControl <- function( contr) {

  if( ! "n.threads" %in% names( contr))  #number of cores to use.  Default is greedy but not super-super greedy
    contr$n.threads <- parallel::detectCores()-1
  if( ! "tag.pred" %in% names( contr))  #how to reference the predictions in the output
    contr$tag.pred <- "pred"
  if( ! "spat.index" %in% names( contr))  #name for the spatial term
    contr$spat.index <- "i"
  if( ! "coord.names" %in% names( contr))
    contr$coord.names <- c("Easting","Northing")
#  if( ! "Meshpars" %in% names( contr))   #constraints on mesh creation
#    contr$Meshpars <- list( max.edge = 0.1, offset = 0, cutoff = 0.001, max.n=1000)
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

  #add to as we go along
  if( !all( names( contr) %in% c("n.threads","tag.pred","spat.index", "coord.names", "verbose",
                                 #"Meshpars",
                                 "prior.mean","int.prec","other.prec", "calcICs", "prior.range", "prior.space.sigma")))
    warning( "There are control parameters specified that are not used.")
  return( contr)
}

# Function to create stack for predictions
# mesh INLA mesh
# covar_raster_data Raster (brick) of covariates used in the model.  Must be on the same measurement scale as the model's covariates.
# tag Name for tag for the stack
# multiLike a boolean to indicate if there are two outcomes or just 1.  It doesn't matter her whether it is PO or PA.

# An INLA stack onto which new data can be projected

MakePredictionStack <- function( mesh, covar_raster_data, tag='pred', varNames, ind) {

  #get the coordinates of the prediction points
  predcoords <- raster::coordinates( covar_raster_data)
  #extract the covariates
  covarData <- raster::extract( covar_raster_data, predcoords[,1:2])[,varNames]

  #find the na positions and remove them.
  na.id <- apply( covarData, 1, function(x) any( is.na( x)))
  predcoords <- predcoords[!na.id,]
  covarData <- as.data.frame( covarData[!na.id,])

  #  covarData$Intercept <- 1
  covarData <- cbind( covarData, 1)
  for( ii in c("PO","AA","PA")){  #ordered in decreasing preference for prediction levels
    if( ind[ii]==1)
      colnames( covarData)[ncol(covarData)] <- paste0("Intercept.",ii)
  }

  resp <- matrix( NA, nrow=nrow( covarData), ncol=sum( ind))

  projgrid <- INLA::inla.mesh.projector( mesh, predcoords)

  # stack the predicted data
  stk <- INLA::inla.stack(list(resp=resp, e=rep(0, nrow(covarData))),
                    A=list(1,projgrid$proj$A), tag=tag, effects=list(covarData, list(i=1:mesh$n)))
  pred=list(stk=stk, predcoords=predcoords)
  return( pred)
}

# Function to create stack for ABUNDANCE absence points

# observs SpatialPointsDataFrame of covariates and observations. Must contain a column corresponding to abundName, which must be non-negative integer (0, 1, 2, ...)
# mesh INLA mesh for the SPDE model
# presname Name of abundance column in observs.
# offsetname Name of the offset variable for AA data.
# tag Name for tag for the stack.

# An INLA stack with binomial data: include Ntrials, which is the number of trials


MakeAAStack=function( observs, mesh = NULL, abundname = "abund", tag = "AA", varNames, sampleAreaName, ind) {

  #casting to numeric
  y.pp <- as.integer( observs@data[,abundname])

  #the area sampled
  offy <- observs@data[,sampleAreaName]

  #covariates 'n' stuff
  covars <- as.data.frame( observs@data[,varNames])
  #  covars$Intercept <- 1

  #still need to think about what is going to be the reference level for the overall prevalence
  #  if( sum( ind) > 1)
  covars$Intercept.AA <- 1

  resp <- matrix( NA, nrow=nrow( covars), ncol=sum( ind))
  colnames( resp) <- names( ind[ind!=0])  #name the variables
  resp[,"AA"] <- y.pp

  # Projector matrix from mesh to data.
  projmat <- INLA::inla.spde.make.A(mesh, as.matrix( coordinates( observs))) # from mesh to point observations

  #make the stack
  stk.abund <- INLA::inla.stack(data=list(resp=resp, e=rep( 1, nrow( covars)), offy=log( offy)),
                          A=list(1,projmat), tag=tag,
                          effects=list( covars, list(i=1:mesh$n)))

  return( stk.abund)
}

# Function to create stack for presence absence points  NOT binomial points

# observs SpatialPointsDataFrame of covariates and observations. Must contain a column corresponding to presename, which must be numeric (0 or 1) or boolean
# mesh INLA mesh for the SPDE model
# presname Name of presences column in observs.
# offsetname Name of the offset variable for PA data.
# tag Name for tag for the stack.

# An INLA stack with binomial data: include Ntrials, which is the number of trials

MakePAstack=function( observs, mesh = NULL, presname = "GroupSize", tag = "PA", varNames, sampleAreaName, ind) {

  #casting to numeric
  y.pp <- as.integer( observs@data[,presname])

  #the number of binomial trials (Bernoulli)
  ntrials <- rep(1, nrow(observs@data))

  #the area sampled
  offy <- observs@data[,sampleAreaName]

  #covariates 'n' stuff
  covars <- as.data.frame( observs@data[,varNames])
  #  covars$Intercept <- 1  #added as datatype specific terms

  #still need to think about what is going to be the reference level for the overall prevalence
  #  if( sum( ind) > 1)
  covars$Intercept.PA <- 1

  resp <- matrix( NA, nrow=nrow( covars), ncol=sum( ind))
  colnames( resp) <- names( ind[ind!=0])  #name the variables
  resp[,"PA"] <- y.pp

  # Projector matrix from mesh to data.
  projmat <- INLA::inla.spde.make.A(mesh, as.matrix( coordinates( observs))) # from mesh to point observations

  #make the stack
  stk.binom <- INLA::inla.stack(data=list(resp=resp, Ntrials=ntrials, offy=log( offy)),
                          A=list(1,projmat), tag=tag,
                          effects=list( covars, list(i=1:mesh$n)))

  return( stk.binom)
}

#observs is a SpatialPointsDataFrame containing a row for each presence and quadrature point. Data must contain all covariates needed
#for the model.
#mesh is the mesh for spde model containing mesh object, the spde object (not used here though), the area of the in-region Voronoi polygons, and the covariate data
#                 note that the names of the covariate data MUST match the names within observs.
#tag is the tag for the INLA stack.
#ind is a 0,1 vector indicating which outcomes/likelihoods will be needed.  Default will produce a single outcome stack.

MakePPPstack <- function( observs, Mesh, presname="presence", tag="PO", ind=1, varNames) {

  #following Simpson et al 2016 (off the grid), as described in https://becarioprecario.bitbucket.io/spde-gitbook/ch-lcox.html and Isaac's et al code

  #the working variate
  y.pp <- c( rep( 0, Mesh$mesh$n), rep( 1, nrow( observs)))

  #its weights
  e.pp <- c( Mesh$w, rep( 0, nrow( observs)))

  #projection matrices
  #mesh nodes to mesh nodes
  imat <- Matrix::Diagonal( Mesh$mesh$n, rep( 1, Mesh$mesh$n))
  #observations to mesh nodes
  lmat <- INLA::inla.spde.make.A( Mesh$mesh, as.matrix( coordinates( observs)))
  #entire projection
  A.pp <- rbind( imat, lmat)

  #covariates 'n' stuff
  covars <- as.data.frame( rbind( Mesh$covars, as.data.frame( observs)[,varNames]))

  #put in an intercept to use in model -- but not here (added as data-type specific terms)
  #  covars$Intercept <- 1
  #  if( sum( ind) > 1)  #put in a constant term for PO data (difference in intensity from PA data)
  covars$Intercept.PO <- 1

  resp <- matrix( NA, nrow=nrow( covars), ncol=sum( ind))
  colnames( resp) <- names( ind[ind!=0])  #name the variables
  resp[,"PO"] <- y.pp

  #make the stack
  stk.pp <- INLA::inla.stack( data = list( resp = resp, e = e.pp, offy=0),
                        A = list(1, A.pp),
                        effects = list(covars, list(i = 1:Mesh$mesh$n)),
                        tag = tag)
  return( stk.pp)

}

#largely taken from Isaacs et al.  But large chunks removed.
#Note that the approach to calculate the weights, and do the quadrature is quite different
#Here, we separate the SPDE mesh and the quadrature points (unlike Isaacs et al). We calculate the weights using ppmData.

MakeSpatialRegion <- function (data = NULL, coords = c("X", "Y"), mesh, bdry, proj = sp::CRS("+proj=utm"), dataBrick=NULL, varNames=NULL) {

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
  if( !is.null( dataBrick))
    covars <- ExtractCovarsAtNodes( mesh=mesh, covars=dataBrick)[,varNames]

  return(list(mesh = mesh, w=w, covars=covars))
}

#mesh is the INLA mesh for the region
#covars is a raster/brick containing all the covariate values.  Mush contain spatial area for all the mesh nodes

ExtractCovarsAtNodes <- function( mesh=NULL, covars=NULL){
  locs <- as.matrix( mesh$loc)
  covarAtLocs <- raster::extract( covars, locs[,1:2])

  return( covarAtLocs)
}
