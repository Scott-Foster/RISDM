
###############################################################################################
###############################################################################################
####	
####	Simulate a GP with Matern covariance with nu =1/2 (exponential), 3/2 or 5/2.
####
####	Returns a matrix with three columns: x, y and random field
####
####	Programmed by Scott in the first half of 2022
####
###############################################################################################
###############################################################################################

myGrfSim <- function( xlim, ylim, nx, ny, sig2, rho, nu=1/2, nug.sig2=0){
  #Matern covariance function ONLY (but smooth=nu=1/2 or 3/2  or 5/2 ONLY -- nu=1/2 is exponential)
  #xlim (ylim) gives the range of x (y) values
  #nx (ny) gives the number of cells in the x and y directions
  #sig2 gives the marginal variance of the process
  #rho is effective range
  #nu is smoothness
  #nug.sig2 is nugget variance
  
  xStep <- diff( xlim) / (nx-1)  #-1 for *both* endpoints
  yStep <- diff( ylim) / (ny-1)  #-1 for *both* endpoints
  
  xseq <- seq( from=xlim[1], to=xlim[2], length=nx)
  yseq <- seq( from=ylim[1], to=ylim[2], length=ny)
  X <- expand.grid( xseq, yseq)

  #quicker(?) way to generate distances, exploiting gridded structure.
  A <- matrix( rep( 0:(nx-1), times=ny), nrow=nx, ncol=ny)
  B <- matrix( rep( 0:(ny-1), each=nx), nrow=nx, ncol=ny)
  
  s <- matrix( NA, nrow=nx*ny, ncol=nx*ny)
  for( xx in (1:nx)-1){
    A.star <- abs( A-xx)
    for( yy in (1:ny)-1){
      B.star <- abs( B-yy)
      distanceFromxx_yy <- sqrt( (A.star*xStep)^2 + (B.star*yStep)^2)
      id <- 1 + xx + nx*yy
      s[id,] <- distanceFromxx_yy
    }
  }
  #calculate the covariance
  if( nu == 1/2)
    covvy1 <- sig2 * exp( -s/rho)
  if( nu == 3/2)
    covvy1 <- sig2 * (1+sqrt(3)*s/rho) * exp( -sqrt(3)*s/rho)
  if( nu == 5/2)
    covvy1 <- sig2 * (1+sqrt(5)*s/rho+5*(s^2)/(3*(rho^2))) * exp( -sqrt(5)*s/rho)
  #take the decomposition, slow.  Is there some sort of structure that can be exploited?
  L <- chol(covvy1)
  
  #take the sample
  xstar <- rnorm( nrow( s), mean=0, sd=1)
  x <- t( L) %*% xstar + rnorm( nrow( s), mean=0, sd=sqrt( nug.sig2))
  
  #return
  simVal <- cbind( X, x)
  
  return( simVal)
}
