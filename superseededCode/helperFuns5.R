makeMesh <- function( ras, max.n=NULL, dep.range=NULL, expandRegion=TRUE, expans.mult=NULL, hull.res=100, max.edge=NULL, cutoff=NULL, offset=NULL, doPlot=TRUE, ...){
  
  if( is.null( max.n)){
    message( "No max.n given. ARBITRARILY setting the max number of mesh nodes to be c(500,200) for within the inner domain and outer domain respectively.\n")
    max.n <- c(500,200)
  }
  
  if( is.null( dep.range)){
    message( "No range of dependence specified (dep.range argument). Assuming that this range is 1/5 of the extent of the raster. Maybe(?) this isn't a good value.  Please check.\n")
    tmp <- matrix( as.vector( raster::extent( ras)), ncol=2, byrow=TRUE)#as.matrix( raster::extent(ras))
    diagLen <- sqrt(sum(apply(tmp, 1, diff)^2))
    dep.range <- diagLen / 3
    rm(tmp, diagLen)
  }
  
  #boundary of the sampling area
  boundary <- list( poly=makeBoundary( ras, lowerRes=TRUE, doPlot=FALSE))
  #  boundary$ras <- list( higher.res=raster::rasterize( boundary$poly$higher.res, ras), lower.res=raster::rasterize( boundary$poly$lower.res, ras))
  boundary$ras <- list( lower.res=raster::rasterize( boundary$poly$lower.res, ras))  
  
  if( expandRegion){
    if( is.null( expans.mult)){
      message( "No value for (convex) expansion multiplier given.  Assuming 1.5 (so that convex.expansion is 1.5 * spatial dependence.\n)")
      expans.mult <- 1.5 #recommended, from HaakonBakkagit.github.io, about 1*range but allow 1.5 for safety (I hope)
    }
    convex.expansion <- expans.mult * dep.range  
  
    #the hull surrounding the spatial domain.
    hully <- INLA::inla.nonconvex.hull( raster::coordinates( ras)[!is.na( raster::values( ras)),], convex = convex.expansion, resolution = rep(hull.res,2))
  }
  else
    hully <- INLA::inla.mesh.segment( loc=boundary$poly)
  
  if( is.null( max.edge)){
    message( "No max.edge given. Assuming that the inner max.edge is 1/5 of spatial dependence range (dep.range argument) and outer max.edge is 1/2 spatial dependence range.\n")
    #advice from (HaakonBakkagit.github.io), except for the outer max.edge
    max.edge <- c( 0.2, 0.5) * dep.range
  }
  
  if( is.null( cutoff)){
    message( "No cutoff given. Assuming that points less than max.edge / 5 are considered to be the same.\n")
    #advice from (HaakonBakkagit.github.io) to get triangles that aren't too 'pointy'
    cutoff <- max.edge / 5
  }
  
  if( expandRegion){
    if( is.null( offset)){
      message( "No offset given. Assuming that the outer domain is approximately an expansion of the inner domain by the amount of a spatial dependence range.\n")
      offset <- c( -0.0001, 1*dep.range)
    }
    else
      offset <- c(-0.0001,1*offset)
  }

  meshy <- INLA::inla.mesh.2d(boundary=hully, max.edge = max.edge, cutoff=cutoff, max.n=max.n, offset=offset)
  
  if( doPlot){
    plot( meshy, asp=1)
    sp::plot( boundary$poly$lower.res, add=TRUE, border='green')
#    sp::plot( boundary$poly$higher.res, add=TRUE, border='red')
    legend("topright", col=c("blue","green",'red'), legend=c("Nonconvex hull defining inner/outer domains","Boundary of raster values (low res)","Boundary of raster values (high res)"), lty=c(1,1), lwd=c(2,2))
  }
  
  meshy$risdmBoundary <- boundary
  meshy$hull <- hully
  
  return( meshy)
} 

findBuffWidth <- function( pts, propOfExtent=1e-6){
  ext <- as.matrix( raster::extent( pts))
  tmp <- apply( ext, 1, diff)
  tmp <- tmp * propOfExtent  
  
  return( min( tmp))
}

cosRule <- function( s){
  tmpRad <- acos( ( s[1]^2+s[2]^2-s[3]^2) / (2*s[1]*s[2]))
  tmpDeg <- 180*( tmpRad / pi)
  return( tmpDeg)
}

find1TriInfo <- function( in1Tri, nl){
  pts <- coordinates( nl[in1Tri])
  sides <- as.vector( dist( pts))
  angles <- rep( NA, 3)
  angles[1] <- cosRule( sides[c(1,2,3)])
  angles[2] <- cosRule( sides[c(1,3,2)])
  angles[3] <- cosRule( sides[c(3,2,1)])
  
  #Heron's formula
  semiPerim <- sum( sides)/2
  area <- sqrt( semiPerim*prod( semiPerim-sides))
  
  res <- c( area, angles)
  return( res)
}

checkMesh <- function( mesh, hull){
  #check inner domain area
  #find inner domain
  nodeLocs <- sp::SpatialPoints( mesh$loc[,1:2])
  dom <- sp::SpatialPolygons( list( sp::Polygons( list( sp::Polygon( hull$loc)), 1)))
  buffWidth <- findBuffWidth(nodeLocs)
  domBuff <- raster::buffer( dom, width=0.1)
  polyNum <- sp::over( nodeLocs, domBuff)
  
  #find angles and areas  
  ouPts <- which( is.na( polyNum))
  
  inTris <- !apply( mesh$graph$tv, 1, function(xx) any( xx %in% ouPts))
  inTris <- mesh$graph$tv[inTris,] #the idx of the locations of each inside triangle's vertice.
  
  triSumms <- apply( inTris, 1, find1TriInfo, nl=nodeLocs)
  areas <- triSumms[1,]
  angles <- as.vector( triSumms[-1,])
  
#  par( mfrow=c(1,3))
  plot( mesh)
  points( nodeLocs[!is.na( polyNum),], col=polyNum[!is.na( polyNum)]+2, pch=20, cex=0.5)
  points( nodeLocs[is.na( polyNum),], col='red', pch=3, cex=0.5)
  
  legend("topright", legend=c("Inner Nodes","Outer Nodes"), pch=c(20,3), col=c(2+1:(length( unique(polyNum))-1),'red'))
  
  hist( areas, main="Inner domain triangle AREAS")
  hist( angles, main="Inner domain triangle ANGLES")
  abline( v=60, col='red')
  
  message( "If mesh is good then inner domain should consist of triangles that are: \n *Visually quite uniform (Left plot),\n *Be comprised of relatively homogeneously sized triangles (middle plot), and \n *Be comprised of triangles that are not weirdly shaped (most like equilateral triangles, right plot).")
  
  return( invisible( NULL))
}

makeBoundary <- function( r, lowerRes=TRUE, doPlot=TRUE){
  #r is raster to base poly off.  If rasterStack or Brick, then only first dim used.
  #lowerRes is Bool to indicate if a lower res version of the boundary is required.  Almost always TRUE, unless original raster is already low res
  #doPlot is bool to indicate if plots of polygon should be produced.
  #returns a SpatialPolygons object containing the outer boundary of the non-NA raster elements.
  r.orig <- r
  myDims <- dim( r)
  if( myDims[3]>1)
    r <- r[[1]]
  #higher res poly
  raster::values( r) <- ifelse( !is.na( raster::values( r)), 1, 0)#ifelse( values( tmp)>-1, 1, NA)
#  statePoly.hi <- raster::rasterToPolygons( x=r, dissolve=TRUE, fun=function(xx){xx==1})  #fun = function(xx) !is.na( xx),
#  #weed out the interior holes.  They aren't so interesting here.  I HOPE!
#  tmp.hi <- as( remove.holes( statePoly.hi), "SpatialPolygons")
  
  #now to make the lower res poly (for mesh creation)
  r <- r.orig
  if( lowerRes){
    facs <- myDims[1:2] %/% 100
    if( any( facs >1 ))
      r <- raster::aggregate( r, fact=facs, expand=TRUE, na.rm=TRUE)  #expand means that the new raster contains the old one.
  }
  raster::values( r) <- ifelse( !is.na( raster::values( r)), 1, 0)
  statePoly.low <- raster::rasterToPolygons( x=r, dissolve=TRUE, fun=function(xx){xx==1}) 
  #weed out the interior holes.  They aren't so interesting here.  I HOPE!
  tmp.low <- as( remove.holes( statePoly.low), "SpatialPolygons")

#  res <- list( higher.res=tmp.hi, lower.res=tmp.low)  
  res <- list( lower.res=tmp.low)  
  
  return( res)
  
}

remove.holes <- function(x) {
  #taken verbatim from package spatialEco https://rdrr.io/cran/spatialEco/src/R/remove.holes.R
  #taken May 6 2022
  
  if(!any(which(utils::installed.packages()[,1] %in% "maptools")))
    stop("please install maptools package before running this function")
  xp <- slot(x, "polygons")
  holes <- lapply(xp, function(x) sapply(methods::slot(x, "Polygons"), methods::slot, "hole"))
  res <- lapply(1:length(xp), function(i) methods::slot(xp[[i]], "Polygons")[!holes[[i]]])
  IDs <- row.names(x)
  x.fill <- sp::SpatialPolygons(lapply(1:length(res), function(i)
    sp::Polygons(res[[i]], ID=IDs[i])), 
    proj4string=sp::CRS(sp::proj4string(x)))
  methods::slot(x.fill, "polygons") <- lapply(methods::slot(x.fill, "polygons"), 
                                              maptools::checkPolygonsHoles)   
  methods::slot(x.fill, "polygons") <- lapply(methods::slot(x.fill, "polygons"), "comment<-", NULL)   
  pids <- sapply(methods::slot(x.fill, "polygons"), function(x) methods::slot(x, "ID"))
  x.fill <- sp::SpatialPolygonsDataFrame(x.fill, data.frame(row.names=pids, ID=1:length(pids)))	   
  return( x.fill )	   
}
