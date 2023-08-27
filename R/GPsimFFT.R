
#simulate a Matern GP based on a grid.  Uses the FFT to do computations much quicker than distance based stuff.
####  A product of Skip Woolley's spare time (~October 2022)
####  All credit goes to Skip.  This is excellent.


bcb <- function (x, y, sig2, rho, nu=1/2) {
  # get a the basis vector for a block-circulant matrix representing the
  # dispersal matrix between equally-spaced grid cells on a torus, as some
  # function `f` of euclidean distance. `x` and `y` are vectors containing
  # equally-spaced vectors for the x and y coordinates of the grid cells on a
  # plane (i.e. the real coordinates). These should have been extended in order
  # to approximate some centre portion as a 2D plane.

  # number and dimension of grid cells
  m <- length(x)
  n <- length(y)
  x_size <- x[2] - x[1]
  y_size <- y[2] - y[1]

  # create indices for x and y on the first row of the dispersal matrix
  xidx <- rep(1:m, n)
  yidx <- rep(1:n, each = m)

  # project onto the torus and get distances from the first cell, in each
  # dimension
  xdist <- abs(xidx - xidx[1])
  ydist <- abs(yidx - yidx[1])
  xdist <- pmin(xdist, m - xdist) * x_size
  ydist <- pmin(ydist, n - ydist) * y_size

  # flatten distances into Euclidean space, apply the dispersal function and
  # return
  d <- sqrt(xdist ^ 2 + ydist ^ 2)
  covFun(d, sig2, rho, nu)
}


covFun <- function(d, sig2, rho, nu=1/2){
  #Matern covariance function ONLY (but smooth=nu=1/2 or 3/2  or 5/2 ONLY -- nu=1/2 is exponential)
  #d is the euclidean distance - changed this to accept the euclidean distance as per bcb function
  #sig2 gives the marginal variance of the process
  #rho is effective range
  #nu is smoothness
  #nug.sig2 is nugget variance - no nugget here.

  #calculate the covariance
  if( nu == 1/2)
    covvy1 <- sig2 * exp( -d/rho)
  if( nu == 3/2)
    covvy1 <- sig2 * (1+sqrt(3)*d/rho) * exp( -sqrt(3)*d/rho)
  if( nu == 5/2)
    covvy1 <- sig2 * (1+sqrt(5)*d/rho+5*(d^2)/(3*(rho^2))) * exp( -sqrt(5)*d/rho)
  #take the decomposition, slow.  Is there some sort of structure that can be exploited?
  # L <- chol(covvy1)

  #take the sample
  # xstar <- stats::rnorm( nrow( s), mean=0, sd=1)
  # x <- t( L) %*% xstar + stats::rnorm( nrow( s), mean=0, sd=sqrt( nug.sig2))

  #return
  # simVal <- cbind( X, x)

  return(covvy1)
}

## Test the block circulant matrix
# image(matrix(bcb(xseq,yseq, sig2 = 1,rho = 1, nu = 1/2),nx,ny))


extend <- function (x, factor = 2) {
  # given an evenly-spaced vector `x` of cell center locations, extend it to the
  # shortest possible vector that is at least `factor` times longer, has a
  # length that is a power of 2 and nests the vector `x` approximately in the
  # middle. This is used to define each dimension of a grid mapped on a torus
  # for which the original vector x is approximately a plane.

  # get cell number and width
  n <- length(x)
  width <- x[2] - x[1]

  # the smallest integer greater than or equal to than n * factor and an
  # integer power of 2
  n2 <- 2 ^ ceiling(log2(factor * n))

  # find how much to pad each end of n
  pad <- n2 - n
  if (pad %% 2 == 0) {
    # if it's even then that's just tickety boo
    pad1 <- pad2 <- pad / 2
  } else {
    # otherwise put the spare cell on the right hand side
    pad1 <- (pad - 1) / 2
    pad2 <- (pad + 1) / 2
  }

  # define the padding vectors
  padx1 <- x[1] - rev(seq_len(pad1)) * width
  padx2 <- x[n] + seq_len(pad2) * width

  # combine these and add an attribute returning the start and end indices for
  # the true vector
  newx <- c(padx1, x, padx2)
  attr(newx, 'idx') <- pad1 + c(1, n)
  newx
}

seq_range <- function (range, by = 1) seq(range[1], range[2], by = by)

setupFFTgp <- function (x, y, factor = 2,  sig2 = 0.5, rho=0.1, nu=1/2) {

  # extend the vectors (to project our plane on <= 1/4 of a torus)
  xe <- extend(x, factor)
  ye <- extend(y, factor)

  # get indices to true vectors
  xidx <- seq_range(attr(xe, 'idx'))
  yidx <- seq_range(attr(ye, 'idx'))

  # get fft basis for covariance on a torus
  bcb_vec <- bcb(ye, xe, sig2, rho, nu)
  ## Sanity check
  # image(matrix(bcb_vec,length(ye),length(xe)))

  # create an empty covariance on all grid cells of the torus
  cov_torus <- matrix(0, length(ye), length(xe))

  # return as a named list for use in each iteration
  list(bcb_vec = bcb_vec,
       cov_torus = cov_torus,
       xidx = xidx,
       yidx = yidx)
}

# x <- y <- seq(-2,2,len=50)
# fs <- setupFFTgp(x,y,2)

fft_GP <- function(meanmat,fs){

# duplicate meanmat to create 'before' condition
meanmat_orig <- meanmat
missing <- is.na(meanmat)

# check for missing values and replace with zeros
if (any(missing)) {
  meanmat[missing] <- 0
}

# insert population matrix into the torus population
# if(!is.null(nugget)) {
  # fs$cov_torus[fs$yidx, fs$xidx] <- meanmat + stats::rnorm(length(fs$cov_torus[fs$yidx, fs$xidx]), mean=0, sd=sqrt(nugget))
# } else {
  fs$cov_torus[fs$yidx, fs$xidx] <- meanmat
# }


cov_fft <- stats::fft(t(fs$cov_torus))
bcb_fft <- stats::fft(fs$bcb_vec)
cov_new_torus_fft <- stats::fft(cov_fft * bcb_fft, inverse = TRUE)

# convert back to real domain, apply correction and transpose
cov_torus_new <- t(Re(cov_new_torus_fft) / length(fs$cov_torus))

cov_new <- cov_torus_new[fs$yidx, fs$xidx]

# return NA values to matrix
cov_new[missing] <- NA

return(cov_new)

}

fftGPsim <- function(x, y, sig2 = 1, rho = 0.5, nu = 1/2, nugget = NULL){

  fs <- setupFFTgp(x = x, y = y, factor = 2, sig2 = sig2, rho = rho, nu = nu)

  nx <- length(x)
  ny <- length(y)

  meanmat <- matrix( stats::rnorm(nx*ny,mean = 0,sd = 1),nx,ny)

  sim <- fft_GP(meanmat = meanmat, fs = fs)

  if(!is.null(nugget))
  sim <- sim + stats::rnorm(nx*ny,mean = 0,sd = sqrt(nugget))

  return(sim)

}


GPMaternSim <- function(x, y, sig2 = 1, rho = 0.5, nu = 1/2, 
                        nugget = NULL, n=1){
  # Input
  #  x: an evenly-spaced vector of cell center locations on the x-axis (e.g. easting)
  #  y: an evenly-spaced vector of cell center locations on the y-axis (e.g. northing)
  #  sig2: spatial variance
  #  rho: scale (i.e. range) parameter of Matérn correlation function
  #  nu: smooth parameter of Matérn correlation function. 
  #      Support nu=1/2, 3/2, and 5/2 only.
  #  nugget: variance of Gaussian white noise
  #  n: number of samples 
  #
  # Output
  #  A matrix of the first two columns representing x coordinate, y coordinate and 
  #  the rest columns the simulated values.
  #

  locs <- as.matrix(expand.grid( x=x, y=y)) # Use matrix to reduce memory size

  #dd <- dist(locs, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)

  # Create the covariacne matrix based on distance
  Sigma <- as.matrix(covFun(stats::dist(locs, method = "euclidean"), sig2, rho, nu))
  
  # Add spatial variance and/or nugget to the diagnoal 
  if (!is.null(nugget)) {
    diag(Sigma) <- sig2 + nugget
  } else {
    diag(Sigma) <- sig2
  }

  invisible(gc()) # clean memory
  
  # Cholesky decomposition on Sigma (Output an upper triangular matrix)
  if( nchar( system.file(package="Rfast")) == 0){
    LL <- chol(Sigma, pivot = FALSE) # pivot = TRUE messes up the order.
  }
  else
    LL <- Rfast::cholesky(Sigma, parallel = ifelse(nrow(Sigma)>1000,TRUE,FALSE))

  rm(Sigma)
  invisible(gc()) # clean memory
  
  # Generate a Gaussian random field
  ZZ <- matrix(stats::rnorm(nrow(locs)*n, mean=0, sd=1), 
         nrow=nrow(locs), ncol=n)

  sim <- crossprod(LL, ZZ) # Conduct t(LL) %*%  

  return(cbind( locs, sim))
 
}
