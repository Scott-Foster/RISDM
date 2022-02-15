
#isdm( POdat=observationData$PO, PAdat=observationData$PA, AAdat=observationData$AA, predLocs=coordinates( tmp_predPoints), responseNames=c(PO="loc",PA="PA",AA="abund"), distributionFormula= ~ 0 + Altitude + Temperature, biasFormula=list(PO=~0+dist2City, PA=~1, AA=~1), covarBrick=dataBrick, boundary=simBorder, control=list(coord.names=c("x","y")))

#reponseNames: a character vector of max length 3.  Elements are names "PO","PA", and/or "AA".
#               The contents of vector gives the name of the response variable in each of the data sets.
#               Note that the
#               No default.  If it is not specified an error will be thrown.
isdm <- function( POdat=NULL, PAdat=NULL, AAdat=NULL,
                                covarBrick=NULL,
                                mesh=NULL,
                                boundary=NULL,
                                responseNames=NULL,
                                sampleAreaNames=NULL,  #a character vector of max length 2. Elements are named "PA" and/or "AA". Refers to one of the cols in PAdat, AAdat respectively
                                distributionFormula=NULL,
                                biasFormulas=list( PO=~1, PO=~1, AA=~1),  #an intercept for each data type.
                                control=list()) {

  #preliminary checks.  Very loose ATM
  if( is.null( responseNames))
    stop("No response specified.  Please do so through responseNames argument.")
  if( length( responseNames)<1 | length( responseNames)>3)
    stop("Please specify responseNames as a named vector of length 3 or less.  See help for guidance.")
  if( length( biasFormulas)<1 | length( biasFormulas)>3)
    stop("Please specify biasFormulas as a named list of length 3 or less.  See help for guidance.")

  #are the pieces available for PO?
  if( any( c( is.null( POdat), is.na( responseNames["PO"]), is.null( biasFormulas[["PO"]]))))
    if (!all( c( is.null( POdat), is.na( responseNames["PO"]), is.null( biasFormulas[["PO"]]))))
    stop("To include PO data, you *must* supply POdat (argument) AND an entry in responseNames AND an entry in biasFormulas")
  #are the pieces available for PA?
  if( any( c( is.null( PAdat), is.na( responseNames["PA"]), is.null( biasFormulas[["PA"]]))))
    if( !all( c( is.null( PAdat), is.na( responseNames["PA"]), is.null( biasFormulas[["PA"]]))))
    stop("To include PA data, you *must* supply PAdat (argument) AND an entry in responseNames AND an entry in biasFormulas")
  #are the pieces available for PO?
  if( any( c( is.null( AAdat), is.na( responseNames["AA"]), is.null( biasFormulas[["AA"]]))))
    if( !all( c( is.null( AAdat), is.na( responseNames["AA"]), is.null( biasFormulas[["AA"]]))))
    stop("To include AA data, you *must* supply AAdat (argument) AND an entry in responseNames AND an entry in biasFormulas")


  #set up FULL control list (priors, nthreads, indexes etc)
  control <- makeControl( control)

  #define boundary if not already done
  if( is.null( mesh))
    stop("No mesh provided. Please create one using the function makeMesh() and pass it in using the mesh argument")
  if( is.null( boundary)){
    message("No boundary provided. Using the extent of covarBrick")
    boundary <- as( raster::extent( covarBrick), "SpatialPolygons")
  }

  #make the complete formula for the INLA call
  fullForm <- makeFormula( distributionFormula, biasFormulas)
  assign( "fullForm", fullForm, envir=environment())
  #make the 'full' formula containing all the raster terms
  rasterForm <- makeFormula( distributionFormula, list( PO=biasFormulas$PO))#, PA=NULL, AA=NULL))
  assign( "rasterForm", rasterForm, envir=environment())

  #Get variable names in formulas
  varNames <- getVarNames( distributionFormula, biasFormulas)
  rasterVarNames <- getVarNames( distributionFormula, list( PO=biasFormulas$PO))

  #create offset (for areas) if not already present
  if( is.null( sampleAreaNames)){
    sampleAreaNames <- rep("area1",2)
    names( sampleAreaNames) <- c("PA","AA")
    if( is.null( PAdat))
      sampleAreaNames <- sampleAreaNames[names( sampleAreaNames) != "PA"]
    else
      PAdat$area1 <- 1

    if( is.null( AAdat))
      sampleAreaNames <- sampleAreaNames[names( sampleAreaNames) != "AA"]
    else
      AAdat$area1 <- 1
  }
  #else assume that varnames are contained in the data.frames

  #create a mesh for inference and prediction
  #Make the mesh and get the node weights (produces warnings on my machine, but only depreciated warnings)
  FullMesh <- MakeSpatialRegion( data=NULL, coords = c("x", "y"), mesh=mesh, bdry=boundary, dataBrick=covarBrick, varNames=rasterVarNames, proj = raster::crs( covarBrick))

  #set up INLA data stacks
  stck <- makeCombinedStack( POdat, PAdat, AAdat, covarBrick, FullMesh, control, varNames, rasterVarNames, sampleAreaNames, responseNames)

  #make a family vector for the INLA call.  Also make the link list in the INLA call
  famLink <- makeFamLink( attr( stck, "ind"))

  #priors etc for fixed effects
  tmp <- list( mean=control$prior.mean)
  tmp$prec <- list( default=control$other.prec)
  if( "Intercept.PO" %in% colnames( stck$effects$data))
    tmp$prec$Intercept.PO <- control$int.prec
  if( "Intercept.PA" %in% colnames( stck$effects$data))
    tmp$prec$Intercept.PA <- control$int.prec
  if( "Intercept.AA" %in% colnames( stck$effects$data))
    tmp$prec$Intercept.AA <- control$int.prec
  my.control.fixed <- tmp

  ####  priors for the spatial model, and build the spatial effects
  my.spde <- INLA::inla.spde2.pcmatern(mesh = FullMesh$mesh,
                              # PC-prior on range: P(practic.range < 0.05) = 0.01
                              prior.range = control$prior.range,
                              # PC-prior on sigma: P(sigma > 1) = 0.01
                              prior.sigma = control$prior.space.sigma)

  #the INLA call...
  mod <- INLA::inla( fullForm, #the formula (combined distribution and bias formulas)
               family = famLink$fam,  #the families for the different data types.
               control.family = famLink$link,  #the link functions to use for each of the data types.
               data = INLA::inla.stack.data(stck),  #the data stack from all data types and the predictions
               verbose = control$verbose,  #don't print out much.
               #control.results = list(return.marginals.random = FALSE, return.marginals.predictor = FALSE),  #some control arguments
               control.predictor = list(A = INLA::inla.stack.A(stck), link=famLink$linkID, compute = TRUE),  #how to predict.  Note the link is generated in makeFamLink()
               control.fixed = my.control.fixed,  #priors etc for the fixed effects.
               Ntrials = INLA::inla.stack.data(stck)$Ntrials, E = INLA::inla.stack.data(stck)$e,
               offset= INLA::inla.stack.data(stck)$offy,
               num.threads = control$n.threads,
               control.compute = list(waic = control$calcICs, dic = control$calcICs, return.marginals = FALSE, return.marginals.predictor = FALSE))

  id <- INLA::inla.stack.index(stck, "pred")$data
  pred <- data.frame(mean = mod$summary.fitted.values$mean[id], stddev = mod$summary.fitted.values$sd[id])

  predRaster <- raster::rasterFromXYZ( cbind( attr( stck, "predcoords"), pred), crs=raster::crs(covarBrick)) ## tmp_predPoints does not exist - changed to predLocs - its only after the appropriate crs right? Matt.
  #Matt: predLocs is not of the right dimension (after NAs removed).  Would work only for regular shaped boundaries (not coastlines). Scott.

  return( list( mod=mod, preds=predRaster))
}
