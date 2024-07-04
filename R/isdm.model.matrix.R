

###############################################################################################
###############################################################################################
####	Yet Another version of model.matrix.
####	Does like model.matrix, except that NAs are removed before creating the bases and then
####	reinserted afterwared.
####	
####	Returns a design matrix
####
####	Programmed by Scott in March/April 2023
####	Touched up by Scott July 2024
####	refactored by Scott July 2024
####
###############################################################################################
###############################################################################################

isdm.model.matrix <- function( formmy, obsy, namy=NULL, includeNA=FALSE) {

  #standardise row names internally
  rownames( obsy) <- 1:nrow( obsy)

  #determining the data set for printing
  printnamy <- namy
  if( is.null( printnamy))	#if it is from the covariate raster
    printnamy <- "covariate brick"

  if( as.character( formmy)[2] == "1")
    XX <- stats::model.matrix( formmy, obsy)
  else{
    #prune down the data to just those variables required.
    varnames <- all.vars( formmy)
    obsy.pruned <- obsy[,varnames,drop=FALSE]  
  
    #remove any row with NAs -- basis expansion won't work with them
    allNAs.id <- apply( obsy.pruned, 1, function(zz) all( is.na(zz)))
    anyNAs.id <- apply( obsy.pruned, 1, function(zz) any( is.na(zz)))
    if( !printnamy %in% c("covariate brick","PO")){
      if( !all( allNAs.id == anyNAs.id))
	warning( paste0(printnamy,": Missing data (NA) present in a partial number of covariates for an observation (or cell). There are ",sum(anyNAs.id)-sum(allNAs.id)," such observations (cells). They are excluded from analysis."))
    }
    obsy.pruned.noNA <- obsy.pruned[!anyNAs.id,,drop=FALSE]
  
    #create the design matrix (no NAs in data
    XX <- stats::model.matrix( formmy, obsy.pruned.noNA)
  
    #put NAs back in, if requested
    if( includeNA){
      XXexpand <- matrix( NA, nrow=nrow( obsy.pruned), ncol=ncol( XX))
      rownames( XXexpand) <- 1:nrow( obsy.pruned)
      colnames( XXexpand) <- colnames( XX)
      XXexpand[rownames( XX),] <- XX
      XX <- XXexpand
    }
  }
  #attend to the names of the variables, by prefixing for data type
  if( !is.null( namy)){
    #rename the intercept (unique between data types)
    if( "(Intercept)" %in% colnames( XX))
      colnames( XX)["(Intercept)" == colnames( XX)] <- paste0( namy,"_Intercept")
    #make sure other variable names are also unique
    colnames( XX)[ !grepl( "_Intercept", colnames( XX))] <- paste0(namy,"_",colnames( XX)[ !grepl( "_Intercept",colnames( XX))])
  }

  #remove special characters for parsing in inla()
  colnames( XX) <- removeParsingChars( colnames( XX))


  return( XX)
}






