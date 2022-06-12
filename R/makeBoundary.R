

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

  r.orig <- r
  myDims <- dim( r)
  #just use the first raster in a rasterLayer/brick to find dimensions
  if( myDims[3]>1)
    r <- r[[1]]
  
#  #so we know what is inside nd outside
#  raster::values( r) <- ifelse( !is.na( raster::values( r)), 1, 0)#ifelse( values( tmp)>-1, 1, NA)
  
  #make the lower res poly (for mesh creation)
  r <- r.orig
  #lower the resoulation
  facs <- myDims[1:2] %/% 100
  if( any( facs >1 ))
    r <- raster::aggregate( r, fact=facs, expand=TRUE, na.rm=TRUE)  #expand means that the new raster contains the old one.

  raster::values( r) <- ifelse( !is.na( raster::values( r)), 1, 0)
  #covert to polygon for hole finding.
  statePoly.low <- raster::rasterToPolygons( x=r, dissolve=TRUE, fun=function(xx){xx==1}) 
  #weed out the interior holes.  They aren't so interesting here.  I HOPE!
  tmp.low <- as( remove.holes( statePoly.low), "SpatialPolygons")

  res <- list( lower.res=tmp.low)  
  
  return( res)
  
}

