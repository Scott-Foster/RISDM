

maternfun <- function( dist1, nu, sigma, kappa){
  matty <- sigma^2 * (kappa * dist1)^nu * besselK( x=kappa * dist1, nu=nu) / ( 2^(nu-1) * gamma( nu))  #matern function, I think.
  return( matty)
}

materncov <- function( locats, nu=1, sigma, rho){
  #locats is matrix or something from which coordinates can be extracted
  #nu is smoothness of Matern
  #sigma is spatial variance
  #rho is effective range (see Lindgren and co)
  if( inherits( locats, "SpatialPoints"))
    locats <- sp::coordinates( locats)
  if( inherits( locats, "raster"))
    locats <- raster::coordinates( locats)
  
  kappa <- sqrt( 8*nu) / rho  #from effective range to scale
  disty <- as.matrix( stats::dist( locats, method="euclidean"))
  covvy <- maternfun( disty, nu, sigma, kappa)
  covvy[disty==0] <- sigma^2
  
  return( covvy)
}

simMaternField <- function( locats, nu, sigma, rho, nugget=sigma/100){
  if( inherits( locats, "SpatialPoints"))
    locats <- sp::coordinates( locats)
  if( inherits( locats, "raster") | inherits( locats, "RasterStack") | inherits( locats, "RasterBrick") | inherits( locats, "RasterLayer"))
    locats <- raster::coordinates( locats)
#  dupes <- duplicated( locats)
#  if( sum( dupes) >0 ){
#    warning( "Removing duplicated locations.  You'll have to reinsert them.")
#    locats
#  }
  
  covvy <- materncov( locats, nu, sigma, rho)
  diag( covvy) <- diag( covvy) + nugget^2
  decomp <- eigen( covvy)
  D <- diag( sqrt( decomp$values))
  mult <- decomp$vectors %*% D
  
  y_star <- rnorm( n=nrow( locats), mean=0, sd=1)
  y <- mult %*% y_star
  
  return( y)
  
}
