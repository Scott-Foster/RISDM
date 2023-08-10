

###############################################################################################
###############################################################################################
####	Make a boundary of the analysis space.  Default is to lower the resolution of the 
####	raster, as a form of smoothing, and base the boundary on that.
####	An attempt to remove holes in the boundary is made, successfull?
####
####	Returns a list with single object of raster. This format is an acident of history
####	
####	Programmed by Scott in the first half of 2022
####	
###############################################################################################
###############################################################################################

makeBoundary <- function( r, doPlot=TRUE){
  #r is raster to base poly off.  If rasterStack or Brick, then only first dim used.
  #doPlot is bool to indicate if plots of polygon should be produced.

  myDims <- dim( r)
  #just use the first raster in a rasterLayer/brick to find dimensions
  if( myDims[3]>1)
    r <- r[[1]]
  
  #make the lower res poly (for mesh creation)
  #lower the resoulation
  facs <- myDims[1:2] %/% 100
  if( any( facs >1 ))
    r <- terra::aggregate( r, fact=facs, expand=TRUE, na.rm=TRUE)  #expand means that the new raster contains the old one.

  terra::values( r) <- ifelse( !is.na( terra::values( r, na.rm=FALSE)), 1, NA)
  #covert to polygon for hole finding.
  #statePoly.low <- raster::rasterToPolygons( x=raster::raster( r), dissolve=TRUE, fun=function(xx){xx==1}) 
  statePoly.low <- terra::as.polygons( x=r, values=FALSE, aggregate=TRUE) 
  #weed out the interior holes.  They aren't so interesting here.  I HOPE!
  statePoly.low <- terra::fillHoles( statePoly.low)
  statePoly.low <- terra::union( statePoly.low)
  
#  tmp.low <- as( remove.holes( statePoly.low), "SpatialPolygons")

  res <- list( lower.res=statePoly.low, lower.res.ras=r)  
  
  return( res)
  
}

