
isdm <- function( obsertvationList=list( POdat=NULL, PAdat=NULL, AAdat=NULL, DCdat=NULL),
                                covarBrick=NULL,
                                mesh=NULL,
#                                boundary=NULL,
                                responseNames=NULL,
                                sampleAreaNames=NULL,  #a character vector of max length 3. Elements are named "PA" and/or "AA" and/or "DC". Refers to one of the cols in PAdat, AAdat, DCdat respectively
                                distributionFormula=NULL,
                                biasFormula=~1,
                                artefactFormula=list( PA=~1, AA=~1, DC=~1),  #an intercept for each data type.
                                DCobserverInfo=list( SurveyID="surveyID", Obs1="Obs1", Obs2="Obs2", Both="Both"), #info for observer stuff
                                control=list()) {

  #preliminary checks.  Very loose ATM
  flag <- checkInput( responseNames, biasFormula, artefactFormula, DCobserverInfo, observationList)
  
  #set up FULL control list (priors, nthreads, indexes etc)
  control <- makeControl( control)

  #define boundary if not already done
  if( is.null( mesh))
    stop("No mesh provided. Please create one using the function makeMesh() and pass it in using the mesh argument")

  #set up relevant information for DC data, if required.
  if( !is.null( DCdat)){
    responseNames['DC'] <- "DCcountDC"
    warning( "Observer probabilities (for double count data) are assumed to be *within* survey.")
    biasFormulas$DC <- update( biasFormulas$DC, paste0("~.+", DCobserverInfo$SurveyID,"/(alpha1Coef+alpha2Coef)"))
    tmpDCoffset <- rep( DCdat[,sampleAreaNames['DC']], times=3)
    DCdat <- prepareDCdata( DCdat, DCobserverInfo)
    DCdat <- cbind( DCdat, tmpDCoffset)
    colnames( DCdat)[ncol( DCdat)] <- sampleAreaNames['DC']
  }

  #make the complete formula for the INLA call
  fullForm <- makeFormula( distributionFormula, biasFormulas, control$addRandom)
  assign( "fullForm", fullForm, envir=environment())
  #make the 'full' formula containing all the raster terms
#  rasterForm <- makeFormula( distributionFormula, list( PO=biasFormulas$PO))#, PA=NULL, AA=NULL))
#  assign( "rasterForm", rasterForm, envir=environment())

  #Get variable names in formulas
  varNames <- getVarNames( distributionFormula, biasFormulas)
  rasterVarNames <- getVarNames( distributionFormula, list( PO=biasFormulas$PO))

  #create offset (for areas) if not already present
  if( is.null( sampleAreaNames)){
    sampleAreaNames <- rep("area1",3)
    names( sampleAreaNames) <- c("PA","AA","DC")
    if( is.null( PAdat))
      sampleAreaNames <- sampleAreaNames[names( sampleAreaNames) != "PA"]
    else
      PAdat$area1 <- 1

    if( is.null( AAdat))
      sampleAreaNames <- sampleAreaNames[names( sampleAreaNames) != "AA"]
    else
      AAdat$area1 <- 1
    
    if( is.null( DCdat))
      sampleAreaNames <- sampleAreaNames[names( sampleAreaNames) != "DC"]
    else
      DCdat$area1 <- 1
  }

  #add the offset for DC expansion 
  if( !is.null( DCdat))
    DCdat[,sampleAreaNames['DC']] <- exp( log( DCdat[,sampleAreaNames['DC']]) + DCdat$expansOffset)
  
  #create a mesh for inference and prediction
  #Make the mesh and get the node weights (produces warnings on my machine, but only depreciated warnings)
  FullMesh <- MakeSpatialRegion( data=NULL, coords = c("x", "y"), mesh=mesh, bdry=mesh$risdmBoundary$poly, dataBrick=covarBrick, varNames=rasterVarNames, proj = raster::crs( covarBrick))

  #set up INLA data stacks
  stck <- makeCombinedStack( POdat, PAdat, AAdat, DCdat, covarBrick, FullMesh, control, varNames, rasterVarNames, sampleAreaNames, responseNames)

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
  if( "Intercept.DC" %in% colnames( stck$effects$data))
    tmp$prec$Intercept.DC <- control$int.prec
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
               control.predictor = list(A = INLA::inla.stack.A(stck), link=famLink$linkID, compute = FALSE),  #how to predict.  Note the link is generated in makeFamLink()
               control.fixed = my.control.fixed,  #priors etc for the fixed effects.
               Ntrials = INLA::inla.stack.data(stck)$Ntrials, 
               E = INLA::inla.stack.data(stck)$e,
               offset = INLA::inla.stack.data(stck)$offy,
               num.threads = control$n.threads,
               control.compute = list(config=TRUE, waic = control$calcICs, dic = control$calcICs, return.marginals = FALSE, return.marginals.predictor = FALSE))

#  id <- INLA::inla.stack.index(stck, "pred")$data
#  pred <- data.frame(mean = mod$summary.fitted.values$mean[id], stddev = mod$summary.fitted.values$sd[id])

#  predRaster <- raster::rasterFromXYZ( cbind( attr( stck, "predcoords"), pred), crs=raster::crs(covarBrick)) ## tmp_predPoints does not exist - changed to predLocs - its only after the appropriate crs right? Matt.
  #Matt: predLocs is not of the right dimension (after NAs removed).  Would work only for regular shaped boundaries (not coastlines). Scott.

  res <- list( mod=mod, mesh=FullMesh$mesh, distFormula=distributionFormula)#, preds=predRaster)
  if( exists( "DCobserverInfo")){
    attr( res, "DCobserverInfo") <- DCobserverInfo
    attr( res, "DCSurveyIDLevels") <- unique( DCdat[,DCobserverInfo$SurveyID])
  }
  else
    attr( res, "DCobserverInfo") <- attr( res, "DCSurveyIDLevels") <- NULL

  class( res) <- "isdm"
  
  return( res)
}

