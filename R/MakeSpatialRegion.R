
###############################################################################################
###############################################################################################
####	
####	Currently this fucntion just calculates the covariates for the observations.
####	In future, and in the past, it has also provide information for the PPP approximation.
####
####	Returns a scalar double
####
####	Programmed by Scott in the first half of 2022.  
####	Based originally on code from Isaacs et al. (2020), Kainski et al (2019) and Flag & Hoegh (2022)
####
###############################################################################################
###############################################################################################


MakeSpatialRegion <- function ( mesh, dataBrick=NULL, varNames=NULL) {

  #find the covariate values under the mesh
  covars <- NULL
  if( !is.null( dataBrick)){
    #check if present
    if( ! all( varNames %in% names( dataBrick)))
      stop( "Covariates in formula not supplied in raster data")
    #extract them!
    tmp <- ExtractCovarsAtNodes( mesh=mesh, covars=dataBrick)
    #drop unwanted.
    covars <- tmp[,varNames,drop=FALSE]
  }

  return(list(mesh = mesh, covars=covars))
}
