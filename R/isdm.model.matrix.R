

###############################################################################################
###############################################################################################
####	ANother version of model.matrix.
####	Does like model.matrix, except that NAs are removed before creating the bases and then
####	reinserted afterwared.
####	
####	Returns a design matrix
####
####	Programmed by Scott in March/April 2023
####
###############################################################################################
###############################################################################################


isdm.model.matrix <- function( formmy, obsy, namy=NULL) {

  #check NA patterns
  allNAs.id <- apply( obsy, 1, function(zz) all( is.na( zz)))
  anyNAs.id <- apply( obsy, 1, function(zz) any( is.na( zz)))
  if( !all( allNAs.id == anyNAs.id))
    stop( paste0("NA present in some, but not all, covariates: data name: ",namy," (if blank then in the raster)"))
  obsy.noNA <- obsy[!allNAs.id,,drop=FALSE]

  #make the model frame
  modFrame <- stats::model.frame( formmy, obsy.noNA, na.action='na.pass')
  
  #The model matrix
  XX_noNA <- stats::model.matrix( formmy, modFrame)  #should work fine with poly etc, NAs have already been dropped.
  
  #paste back together
  XX <- matrix( NA, nrow=length(allNAs.id), ncol=ncol( XX_noNA))
  XX[!allNAs.id,] <- XX_noNA
  colnames( XX) <- colnames( XX_noNA)
  
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
