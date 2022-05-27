
isdm <- function( observationList=list( POdat=NULL, PAdat=NULL, AAdat=NULL, DCdat=NULL),
                                covarBrick=NULL,
                                mesh=NULL,
#                                boundary=NULL,
                                responseNames=NULL,
                                sampleAreaNames=NULL, 
                                distributionFormula=NULL,
                                biasFormula=NULL,
                                artefactFormulas=list( PO=NULL, PA=~1, AA=~1, DC=~1),  #an intercept for each data type.
                                DCobserverInfo=list( SurveyID="surveyID", Obs1="Obs1", Obs2="Obs2", Both="Both"), #info for observer stuff
                                control=list()) {

  #preliminary checks.  Very loose ATM
  flag <- checkInput( responseNames, biasFormula, artefactFormulas, DCobserverInfo, observationList)
  
  #make variable names in artefact models unique -- so that factor levels etc are not shared between data types
  tmp <- uniqueVarNames( observationList, artefactFormulas, DCobserverInfo$SurveyID)
  artefactFormulas <- tmp$arteForm
  observationList <- tmp$obsList
  DCobserverInfo$SurveyID <- tmp$DCsurvID
  
  #set up FULL control list (priors, nthreads, indexes etc)
  control <- makeControl( control)

  #define boundary if not already done
  if( is.null( mesh))
    stop("No mesh provided. Please create one using the function makeMesh() and pass it in using the mesh argument")

  #set up relevant information for DC data, if required.
  if( !is.null( observationList$DCdat)){
    responseNames['DC'] <- "DCcountDC"
##    warning( "Observer probabilities (for double count data) are assumed to be *within* survey.")
    observationList$DCdat <- prepareDCdata( DCdat=observationList$DCdat, DCobserverInfo=DCobserverInfo, sampArea=sampleAreaNames["DC"], control$DCmethod)
    if( control$DCmethod == "TaylorsLinApprox")
      artefactFormulas$DC <- update( artefactFormulas$DC, paste0("~.+", DCobserverInfo$SurveyID,":(alpha1Coef+alpha2Coef)"))
  }

  #make the complete formula for the INLA call
  fullForm <- makeFormula( distributionFormula, biasFormula, artefactFormulas, control$addRandom)
  assign( "fullForm", fullForm, envir=environment())

  #Get variable names in formulas
  varNames <- getVarNames( distributionFormula, biasFormula, artefactFormulas)
  rasterVarNames <- getVarNames( distributionFormula, biasFormula, list(PA=NULL, AA=NULL, DC=NULL))

  #create offset (for areas) if not already present
  tmp <- createOffsets( sampleAreaNames, observationList)
#  sampleAreaNames <- tmp$names
  observationList <- tmp$dat
    
  #create a mesh for inference and prediction
  #Make the mesh and get the node weights (produces warnings on my machine, but only depreciated warnings)
  FullMesh <- MakeSpatialRegion( mesh=mesh, dataBrick=covarBrick, varNames=rasterVarNames)

  #set up INLA data stacks
  stck <- makeCombinedStack( observationList, covarBrick, FullMesh, control, varNames, rasterVarNames, sampleAreaNames, responseNames)

  #wts for loglikelihood to attempt to exclude the influence on PO data likelihood of obs outside region (hull and expansion)
  loglWts <- makeLoglWts( stck)
  
  #make a family vector for the INLA call.  Also make the link list in the INLA call
  famLink <- makeFamLink( attr( stck, "ind"))

  #priors etc for fixed effects
  tmp <- setPriors( control, stck)
  my.control.fixed <- tmp

  ####  priors for the spatial model, and build the spatial effects
  my.spde <- INLA::inla.spde2.pcmatern(mesh = FullMesh$mesh,
                              # PC-prior on range: P(practic.range < 0.05) = 0.01
                              prior.range = control$prior.range,
                              # PC-prior on sigma: P(sigma > 1) = 0.01
                              prior.sigma = control$prior.space.sigma)

  print( class( my.spde))
  
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
               weights=loglWts,
               offset = INLA::inla.stack.data(stck)$offy,
               num.threads = control$n.threads,
               control.compute = list(config=TRUE, waic = control$calcICs, dic = control$calcICs, return.marginals = FALSE, return.marginals.predictor = FALSE))

  res <- list( mod=mod, distributionFormula=distributionFormula, biasFormula=biasFormula, artefactFormulas=artefactFormulas, mesh=FullMesh$mesh)
  if( control$returnStack)
    res$stack <- stck
  if( exists( "DCobserverInfo")){
    attr( res, "DCobserverInfo") <- DCobserverInfo
    attr( res, "DCSurveyIDLevels") <- unique( observationList$DCdat[,DCobserverInfo$SurveyID])
  }
  else
    attr( res, "DCobserverInfo") <- attr( res, "DCSurveyIDLevels") <- NULL

  class( res) <- "isdm"
  
  return( res)
}

