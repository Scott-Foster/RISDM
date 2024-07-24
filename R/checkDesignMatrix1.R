
checkDesignMatrix <- function( covarBrick, distForm=NULL, biasForm=NULL, habitatArea=NULL, do.pairs=TRUE, ask=TRUE){
  
  #plotting up the residuals. RQR versus fitted and for PO data the raster.
  if( ask){
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit( grDevices::devAskNewPage(FALSE))
  }
  
  if( is.null( distForm) & is.null( biasForm))
    stop( "Must specify at least one of distForm or biasForm")
  
  # if( is.null( habitatArea)){
  #   habitatArea <- "cellSize"
  #   covarBrick <- c( covarBrick, terra::cellSize( covarBrick))
  #   names( covarBrick)[terra::nlyr( covarBrick)] <- habitatArea
  # }
  
  
  designMat <- RISDM:::uniqueVarNames( obsList=list(), 
                                     covarBrick=covarBrick,
                                     distForm=distForm,
                                     biasForm=biasForm,
                                     arteForm=list(), 
                                     habitatArea=habitatArea, stdCovs=TRUE, DCsurvID=NA)$covarBrick
  designMat <- as.data.frame( designMat)
  #Intercepts will cause trouble, remove them
  message( "Removing all intercepts")
  designMat <- designMat[,!grepl( "_Intercept", colnames( designMat))]
  
  corrplot::corrplot(
      stats::cor( designMat, method="pearson", use="complete.obs"),
      method = 'number',
      number.digits = 2,
      type = 'lower',
      diag = FALSE)
  
  tmp <- eigen( cov(designMat))
  cat( paste0( "Largest e-value: ", round( tmp$values[1], 4)," Smallest e-value: ", round( tmp$values[ncol(designMat)], 4),". Condition number: ", round( tmp$values[1]/tmp$values[ncol(designMat)],4),"\n"))

  if( do.pairs)
    pairs( designMat[sample.int(nrow( designMat), min( 5000, nrow( designMat))),], pch='.', lower.panel=function(x,y,...){points(x,y,...)}, upper.panel=NULL)
  
  return( invisible( TRUE))
}
