
###############################################################################################
###############################################################################################
####	SDM analysis from multiple sources of dispirate data
####	Takes different data sources, as a list of data.frames, and various formulas for 
####	different model components.
####	Will need to have area surveyed for each of the non-PO data soures
####
####	Returns an object of class "isdm" 
####			(a list whose main component is an object of class "INA")
####
####	Programmed by Scott in the first half of 2022
####	Genesis of code started with Isaac's et al. (2020) R function in fitModel.R, however
####		the code has expanded considerably from that to include more data types
####		and more flexibility and some different methods. Code is completely
####		reliant on INLA and spatial data pacakges
####
###############################################################################################
###############################################################################################

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

  #preliminary checks.  Aimed to catch obvious mis-specification only
  flag <- checkInput( responseNames, biasFormula, artefactFormulas, DCobserverInfo, observationList)
  #Check for mesh too
  if( is.null( mesh))
    stop("No mesh provided. Please create one using the function makeMesh() and pass it in using the mesh argument")
  
  #make variable names in artefact models unique -- so that factor levels etc are not shared between data types
  tmp <- uniqueVarNames( observationList, artefactFormulas, DCobserverInfo$SurveyID)
  artefactFormulas <- tmp$arteForm
  #same for observation list (the species data)
  observationList <- tmp$obsList
  #same for double count information
  DCobserverInfo$SurveyID <- tmp$DCsurvID
  
  #set up FULL control list (priors, nthreads, indexes etc)
  control <- makeControl( control)

  #set up relevant information for DC data, if required.
  if( !is.null( observationList$DCdat)){
    responseNames['DC'] <- "DCcountDC"
    observationList$DCdat <- prepareDCdata( DCdat=observationList$DCdat, DCobserverInfo=DCobserverInfo, sampAreaDC=sampleAreaNames["DC"], control$DCmethod)
    #add terms to the formula for the Taylor series approach.  But not if method is 'plugin'
    if( control$DCmethod == "TaylorsLinApprox")
      artefactFormulas$DC <- update( artefactFormulas$DC, paste0("~.+", DCobserverInfo$SurveyID,":(alpha1Coef+alpha2Coef)"))
  }

  #make the complete formula for the INLA call.  With all the right interactions etc.
  fullForm <- makeFormula( distributionFormula, biasFormula, artefactFormulas, control$addRandom)
  #this is an annoying scoping thing.  Took me ages to figure this out, and I still don't quite believe it.
  assign( "fullForm", fullForm, envir=environment())

  #Get variable names in formulas
  varNames <- getVarNames( distributionFormula, biasFormula, artefactFormulas)
  #and those that relate to rasters.
  rasterVarNames <- getVarNames( distributionFormula, biasFormula, list(PA=NULL, AA=NULL, DC=NULL))

  #create offset (for areas) if not already present
  tmp <- createOffsets( sampleAreaNames, observationList)
  #in an accident of history observationList gets rewritten again.
  observationList <- tmp$dat
    
  #make the mesh for approximating random effect over.
  FullMesh <- MakeSpatialRegion( mesh=mesh, dataBrick=covarBrick, varNames=rasterVarNames)

  #set up INLA data stacks
  stck <- makeCombinedStack( observationList, covarBrick, FullMesh, control, varNames, rasterVarNames, sampleAreaNames, responseNames)

  #wts for loglikelihood to attempt to exclude the influence on PO data likelihood of obs outside region (hull and expansion)
  #note that this currently does nothing and is therefore commented out.
#  loglWts <- makeLoglWts( stck)
  
  #make a family vector for the INLA call.  Also make the link list in the INLA call
  famLink <- makeFamLink( attr( stck, "ind"))

  #priors etc for fixed effects
  tmp <- setPriors( control, stck)
  my.control.fixed <- tmp

  #priors for the spatial model, and build the spatial effects
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
#               weights=loglWts,  #currently not used but could be if Berman-Turner implemented.
               offset = INLA::inla.stack.data(stck)$offy,
               num.threads = control$n.threads,
               control.compute = list(config=TRUE, waic = control$calcICs, dic = control$calcICs, return.marginals = FALSE, return.marginals.predictor = FALSE))

  #the return object, which contains some of the bits and pieces calcuated above.
  res <- list( mod=mod, distributionFormula=distributionFormula, biasFormula=biasFormula, artefactFormulas=artefactFormulas, mesh=FullMesh$mesh, control=control, responseNames=responseNames)
  #include the stack if requested
  if( control$returnStack){
    res$stack <- stck
    res$observationList <- observationList
  }
  #include DC information if presente
  if( exists( "DCobserverInfo")){
    attr( res, "DCobserverInfo") <- DCobserverInfo
    attr( res, "DCSurveyIDLevels") <- unique( observationList$DCdat[,DCobserverInfo$SurveyID])
  }
  else
    attr( res, "DCobserverInfo") <- attr( res, "DCSurveyIDLevels") <- NULL
  attr( res, "n.threads") <- control$n.threads
  attr( res, "coord.names") <- control$coord.names

  class( res) <- "isdm"
  
  return( res)
}
