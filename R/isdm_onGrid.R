
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
####	The first versions of this function had a structure that was based on 
####		Dambly et al (2019). https://doi.org/10.5281/zenodo.3363936 and their fucntion
####		code fitModel.R
####	The code of Damble et al (2019) also provided the pattern for dealing with INLA stacks 
####		and combining them, see makeCombinedStack().
####
####	This current code has evolved substantially from those early days, 
####		in terms of functionality and computational methods, user interface, 
####		return objects, associated functions etc
####
####	This code is completely reliant on INLA and spatial data pacakges
####
####	Re-jigged in March 2023 to make explicit for INLA the particular pattern
####		of constraints for the factor level effects
####
###############################################################################################
###############################################################################################

isdm <- function( observationList=list( POdat=NULL, PAdat=NULL, AAdat=NULL, DCdat=NULL),
                                covars=NULL,
				habitatArea=NULL,
                                mesh=NULL,
                                responseNames=NULL,
                                sampleAreaNames=NULL, 
                                distributionFormula=NULL,
                                biasFormula=NULL,
                                artefactFormulas=list( PO=NULL, PA=~1, AA=~1, DC=~1),  #an intercept for each data type.
                                DCobserverInfo=list( SurveyID="surveyID", Obs1="Obs1", Obs2="Obs2", Both="Both"), #info for observer stuff
                                control=list(), ...) {

  #preliminary checks.  Aimed to catch obvious mis-specification only
  flag <- checkInput( responseNames, biasFormula, artefactFormulas, DCobserverInfo, observationList)
  #covnert raster -> terra
  if( inherits( covars, "Raster"))
    covars <- terra::rast( covars)
  #Check for mesh too
  if( is.null( mesh))
    stop("No mesh provided. Please create one using the function makeMesh() and pass it in using the mesh argument")

  #set up FULL control list (priors, nthreads, indexes etc)
  control <- makeControl( control, covar.ext=terra::ext( covars))

  #set up relevant information for DC data, if required.
  if( !is.null( observationList$DCdat)){
    if( !control$DCmethod %in% c( "plugin", "TaylorsLinApprox"))
      stop( "Specified DC method not supported. Currently must be either 'TaylorsLinApprox' or 'plugin'.  Please see control$DCmethod argument.")
    responseNames['DC'] <- "DCcountDC"
    observationList$DCdat <- prepareDCdata( DCdat=observationList$DCdat, DCobserverInfo=DCobserverInfo, sampAreaDC=sampleAreaNames["DC"], DCmethod=control$DCmethod, coord.names=control$coord.names)
    #add terms to the formula for the Taylor series approach.  But not if method is 'plugin'
    if( control$DCmethod == "TaylorsLinApprox"){
      if( attr( observationList$DCdat, "nsurvey") > 1)
	artefactFormulas$DC <- update( artefactFormulas$DC, paste0("~.+", DCobserverInfo$SurveyID,":logDetectPi"))
      else
	artefactFormulas$DC <- update( artefactFormulas$DC, "~.+logDetectPi")
    }
  }

  #make variable names in artefact models unique -- so that factor levels etc are not shared between data types
  newInfo <- uniqueVarNames( obsList=observationList, covarBrick=covars, distForm=distributionFormula, biasForm=biasFormula, arteForm=artefactFormulas, habitatArea=habitatArea, DCsurvID=DCobserverInfo$SurveyID, coord.names=control$coord.names, responseNames=responseNames, sampleAreaNames=sampleAreaNames, stdCovs=control$standardiseCovariates, na.action=control$na.action)

  #same for double count information
  DCobserverInfo$SurveyID <- newInfo$DCsurvID

  #create offset (for areas) if not already present
  newInfo$offy <- createOffsets( sampleAreaNames, newInfo$obsList)
    
  #make the mesh for approximating random effect over.
  FullMesh <- MakeSpatialRegion( mesh=mesh, dataBrick=newInfo$covarBrick, varNames=names( newInfo$covarBrick))

  #set up INLA data stacks
  stck <- makeCombinedStack( newInfo$obsList, newInfo$covarBrick, habitatArea, newInfo$distributionFormula, newInfo$artefactFormulas, FullMesh, control, responseNames, ind=newInfo$ind, sampleAreaNames)#varNames, rasterVarNames)

  #make the complete formula for the INLA call.  With all the right interactions etc.
  fullForm <- combineFormulae( list( newInfo$distForm, newInfo$biasForm, newInfo$arteForm), addRE=control$addRandom)
#  fullForm <- makeFormula( distributionFormula, biasFormula, attr( stck, "newArtForms"), control$addRandom, interactArtefact=FALSE)
  environment( fullForm) <- environment()

  #wts for loglikelihood to attempt to exclude the influence on PO data likelihood of obs outside region (hull and expansion)
  #note that this currently does nothing and is therefore commented out.
#  loglWts <- makeLoglWts( stck)
  
  #make a family vector for the INLA call.  Also make the link list in the INLA call
  famLink <- makeFamLink( attr( stck, "ind"), attr( stck, "nObs"))

  #priors etc for fixed effects
  my.control.fixed <- setPriors( control, stck)

  #priors for the spatial model, and build the spatial effects
  my.spde <- INLA::inla.spde2.pcmatern(mesh = FullMesh$mesh, constr=control$re.constr,
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
               control.predictor = list(A = INLA::inla.stack.A(stck), link=1, compute = FALSE),  #how to predict.  Note the link is generated in makeFamLink()
               control.fixed = my.control.fixed,  #priors etc for the fixed effects.
               Ntrials = INLA::inla.stack.data(stck)$Ntrials, 
               E = INLA::inla.stack.data(stck)$e,
#               weights=loglWts,  #currently not used but could be if Berman-Turner implemented.
               offset = INLA::inla.stack.data(stck)$offy,
               num.threads = control$n.threads,
               control.compute = list(config=TRUE, waic = control$calcICs, dic = control$calcICs, return.marginals = FALSE, return.marginals.predictor = FALSE),
	       control.inla = list( tolerance=control$coverg.tol, control.vb=list(enable=control$vb.correction)),
	       safe=TRUE,
	       inla.mode=control$inla.mode,
	       ...)

  #the return object, which contains some of the bits and pieces calcuated above.
  res <- list( mod=mod, distributionFormula=distributionFormula, biasFormula=biasFormula, artefactFormulas=artefactFormulas, mesh=FullMesh$mesh, control=control, responseNames=responseNames, 
		data=list( covars=newInfo$covarBrick, obsList=newInfo$obsList))
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

