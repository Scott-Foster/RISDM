
###############################################################################################
###############################################################################################
####	Make a mesh for 2d analysis (especially for isdm analyses)
####	Based on rules of thumb etc 
####	Tries to simplify as far as possible the complexities of mesh generation by automation
####	Ambitious goal -- the user should be careful of outputs: is it a useful mesh?
####
####	Returns an object of class "inla.mesh" with a few additions to the list
####	
####	Programmed by Scott in the first half of 2022
####	Function is a wrapper to INLA::inla.mesh.2d. The wrapper aims to find some useful
####	arguments to call inla.mesh.2d.
####
###############################################################################################
###############################################################################################

makeMesh <- function( ras, max.n=NULL, dep.range=NULL, expandRegion=TRUE, expans.mult=NULL, hull.res=100, max.edge=NULL, cutoff=NULL, offset=NULL, barrier=list(useBarrier=FALSE, range.fraction=0.1), doPlot=TRUE, ...){
  
  #check input for important parameters.  Set (possibly dumb) defaults.
  if( is.null( max.n)){
    message( "No max.n given. ARBITRARILY setting the max number of mesh nodes to be c(500,200) for within the inner domain and outer domain respectively.\n")
    max.n <- c(500,200)
  }
  #check input for important parameters.  Set (possibly dumb) defaults.  
  if( is.null( dep.range)){
    message( "No range of dependence specified (dep.range argument). Assuming that this range is 1/5 of the extent of the raster. Maybe(?) this isn't a good value.  Please check.\n")
    tmp <- matrix( as.vector( terra::ext( ras)), ncol=2, byrow=TRUE)
    diagLen <- sqrt(sum(apply(tmp, 1, diff)^2))
    dep.range <- diagLen / 3
    rm(tmp, diagLen)
  }
  #check barrier input
  if( barrier$useBarrier){
    if( !is.double( barrier$range.fraction))
      stop( "Inappropriate value for barrier$range.fraction. Please specify a double. Default is 0.1.")
    if( barrier$range.fraction < 0)
      message( "You have specified a barrier$range.fraction < 0.  Continuing, but this is unlikely to be what you wanted. You will want this to be 0<range.franction<1, with the (possibly useful) default being 0.1.")
    if( barrier$range.fraction > 1)
      message( "You have specified a barrier$range.fraction > 1.  Continuing, but this is unlikely to be what you wanted. You will want this to be 0<range.franction<1, with the (possibly useful) default being 0.1.")
  }
  
  #boundary of the sampling area
  boundary <- list( poly=makeBoundary( ras))
  #as a raster
  boundary$ras <- list( lower.res=boundary$poly$lower.res.ras, specified.res=terra::rasterize( boundary$poly$lower.res, ras))
  #clean up poly list
  boundary$poly <- boundary$poly[1]
  
  #if the region is going to be expanded, then expanded. This is often a good idea.
#  if( expandRegion){
    #a default ('cause something has to be right?
    if( is.null( expans.mult)){
      message( "No value for (convex) expansion multiplier given.  Assuming 1.5 (so that convex.expansion is 1.5 * spatial dependence.\n)")
      expans.mult <- 1.5 #recommended, from HaakonBakkagit.github.io, about 1*range but allow 1.5 for safety (I hope)
    }
    convex.expansion <- expans.mult * dep.range  
  
    #the hull surrounding the spatial domain.
#    hully <- INLA::inla.nonconvex.hull( terra::crds( boundary$ras$lower.res, na.rm=TRUE), convex=convex.expansion, resolution=rep( hull.res,2))
    hully <- fmesher::fm_nonconvex_hull( x=terra::crds( boundary$ras$lower.res, na.rm=TRUE), convex=convex.expansion)

  #more checks for defaults
  if( is.null( max.edge)){
    message( "No max.edge given. Assuming that the inner max.edge is 1/5 of spatial dependence range (dep.range argument) and outer max.edge is 1/2 spatial dependence range.\n")
    #advice from (HaakonBakkagit.github.io), except for the outer max.edge
    max.edge <- c( 0.2, 0.5) * dep.range
  }

  #setting defaults
  if( is.null( cutoff)){
    message( "No cutoff given. Assuming that points less than min( max.edge) / 5 are considered to be the same.\n")
    #advice from (HaakonBakkagit.github.io) to get triangles that aren't too 'pointy'
    cutoff <- min( max.edge) / 5
  }
  
  #defaults
#  if( expandRegion){
    if( is.null( offset)){
      message( "No offset given. Assuming that the outer domain is approximately an expansion of the inner domain by the amount of a spatial dependence range.\n")
      offset <- c( -0.0001, 1*dep.range)
    }
    else
      offset <- c(-0.0001,1*offset)
#  }

  #make the mesh.  Note that sometimes the arguments seem to be treated as a 'guide' rather than gospel -- I guess
  #	that sometimes mesh solutions cannot be found for particular arguments.
#    meshy <- INLA::inla.mesh.2d(boundary=hully, max.edge = max.edge, cutoff=cutoff, max.n=max.n, offset=offset)
    meshy <- fmesher::fm_mesh_2d(boundary=hully, max.edge = max.edge, cutoff=cutoff, max.n=max.n, offset=offset)
  
  #Set the CRS
  fmesher::fm_crs( meshy) <- terra::crs( ras)
  
  if( !barrier$useBarrier){
    triBarrier.xy <- NULL
    range.fraction <- NULL
    poly.barrier <- NULL
  }
  else{    
  #define the barrier nodes
    tmpRas <- !is.na( ras)
    tmpRas[!tmpRas] <- NA
    my.bound.ALL.poly_agg <- terra::as.polygons( tmpRas)
    my.bound.ALL.poly_agg <- sf::st_as_sf( my.bound.ALL.poly_agg)

#    my.bound.ALL.poly <- sf::combineGeoms( boundaries[1,], boundaries[2,])
#    my.bound.ALL.poly <- sf::st_as_sf( my.bound.ALL.poly)

    #which nodes are NOT in the barrier region
    notTriBarrier <- unlist( fmesher::fm_contains(
      x = my.bound.ALL.poly_agg, 
      y = meshy, 
      type = "centroid"))
      
    #visual check of outside of barriers
    triCentres.xy <- cbind(
        meshy$loc[meshy$graph$tv[,1], 1:2] +
        meshy$loc[meshy$graph$tv[,2], 1:2] +
        meshy$loc[meshy$graph$tv[,3], 1:2])/3
  
    triBarrier.xy <- setdiff( 1:nrow( triCentres.xy), notTriBarrier)
    range.fraction = barrier$range.fraction
    poly.barrier = sf::st_as_sf( INLA::inla.barrier.polygon( meshy, triBarrier.xy)) #polygon for plotting etc
  }

  
  #is the result going to be plotted?
  if( doPlot){
#    INLA:::plot.inla.mesh( meshy, asp=1)
    fmesher:::plot.fm_mesh_2d( x=fmesher::fm_as_fm( meshy), xlim=range( meshy$loc[,1], na.rm=TRUE), ylim=range( meshy$loc[,2], na.rm=TRUE))
#    sp::plot( boundary$poly$lower.res, add=TRUE, border='green')
    if( !barrier$useBarrier){
      terra::plot( terra::buffer( boundary$poly$lower.res, width=1e-5), add=TRUE, border='green')
      legend("topright", col=c("blue","green"), legend=c("Nonconvex hull defining inner/outer domains","Boundary of raster values (low res)"), lty=c(1,1), lwd=c(2,2))
    }
    else{
#      plot( poly.barrier, add=TRUE, alpha=0.5, col=grey(0.7))
      points( triCentres.xy[triBarrier.xy,1:2], pch=20, col='gold', cex=0.3)
      legend("topright", col=c("blue","green",'gold'), legend=c("Nonconvex hull defining inner/outer domains","Boundary of raster values (low res)","Barrier"), lty=c(1,1,NA), lwd=c(2,2,NA), pch=c(NA, NA, 20))
    }
  }

  #add boundary to object
  meshy$risdmBoundary <- boundary
  meshy$hull <- hully
  meshy$dep.range <- dep.range
  meshy$barrier <- list( triBarrier.xy=triBarrier.xy, range.fraction=range.fraction, barrier.poly=poly.barrier)
  
  return( meshy)
} 
