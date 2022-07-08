
###############################################################################################
###############################################################################################
####	
####	Check elementary validity of mesh for analysis (and helper functions).
####
####	Returns NULL (invisibly) but produces a graph as a side-effect
####
####	Programmed by Scott in the first half of 2022
####
###############################################################################################
###############################################################################################

checkMesh <- function( mesh, hull){
  #check inner domain area
  #find inner domain
  nodeLocs <- sp::SpatialPoints( mesh$loc[,1:2])
  #and the bounding hull.
  dom <- sp::SpatialPolygons( list( sp::Polygons( list( sp::Polygon( hull$loc)), 1)))
  #buffer the node locations
  domBuff <- raster::buffer( dom, width=0.1)
  polyNum <- sp::over( nodeLocs, domBuff)
  
  #find angles and areas outside area
  ouPts <- which( is.na( polyNum))
  
  inTris <- !apply( mesh$graph$tv, 1, function(xx) any( xx %in% ouPts))
  inTris <- mesh$graph$tv[inTris,] #the idx of the locations of each inside triangle's vertice.
  
  #area and angle properties of triangles inside area
  triSumms <- apply( inTris, 1, find1TriInfo, nl=nodeLocs)
  areas <- triSumms[1,]
  angles <- as.vector( triSumms[-1,])
  
  par( mfrow=c(1,3))
  
  #plotting mesh and inner/outer
  plot( mesh)
  points( nodeLocs[!is.na( polyNum),], col=polyNum[!is.na( polyNum)]+2, pch=20, cex=0.5)
  points( nodeLocs[is.na( polyNum),], col='red', pch=3, cex=0.5)
  
  legend("topright", legend=c("Inner Nodes","Outer Nodes"), pch=c(20,3), col=c(2+1:(length( unique(polyNum))-1),'red'))
  
  #distribution of triangles inside area.
  hist( areas, main="Inner domain triangle AREAS")
  hist( angles, main="Inner domain triangle ANGLES")
  abline( v=60, col='red')
  
  message( "If mesh is good then inner domain should consist of triangles that are: \n *Visually quite uniform (Left plot),\n *Be comprised of relatively homogeneously sized triangles (middle plot), and \n *Be comprised of triangles that are not weirdly shaped (most like equilateral triangles, right plot).")
  
  return( invisible( NULL))
}

cosRule <- function( s){
  #angle within triangles.
  tmpRad <- acos( ( s[1]^2+s[2]^2-s[3]^2) / (2*s[1]*s[2]))
  tmpDeg <- 180*( tmpRad / pi)
  return( tmpDeg)
}

find1TriInfo <- function( in1Tri, nl){
  #for each triangle, find the area and the angle.
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

