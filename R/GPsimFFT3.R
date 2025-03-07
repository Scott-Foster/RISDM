#simulate a Matern GP based on a grid.  Uses the FFT to do computations much quicker than distance based stuff.
####  A product of Skip Woolley's spare time (~October 2022)
####  All credit goes to Skip.  This is excellent.

# WH:
# Some functions were modified to simulate a 2D Matern GP on a grid using
# FFT. Those modified functions with suffix "2" are
#
#  * extend2()
#  * fft_GP2()
#  * fftGPsim2()
#  * setupFFTgp2()
#
# New functions:
#
#  * rho_criteria()


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

################################################################################
# extend2()
extend2 <- function (x, factor = 2, min_exponent2=12, padEndOnly = TRUE) {
  # given an evenly-spaced vector `x` of cell center locations, extend it to the
  # shortest possible vector that is at least `factor` times longer, has a
  # length that is a power of 2 and nests the vector `x` approximately in the
  # middle. This is used to define each dimension of a grid mapped on a torus
  # for which the original vector x is approximately a plane.

  # get cell number and width
  n <- length(x)
  width <- x[2] - x[1]

  # set the power of two has the exponent at least equal to min_exponent2;
  #  otherwise, FFT generates negative eigenvalues 
  minp <- max(ceiling(log2(factor * n)), min_exponent2)
  n2 <- 2^minp

  # find how much to pad each end of n
  pad <- n2 - n

  if (padEndOnly == TRUE) {
    # pad to the end
    # define the padding vectors
    padx2 <- x[n] + seq_len(pad) * width

    # combine these and add an attribute returning the start and end indices for
    # the true vector
    newx <- c(x, padx2)
    attr(newx, 'idx') <- c(1, n)
    #newx

  } else {
    # pad to the front and end

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
    #newx

  }

  return(newx) 

}

# fft_GP2()
fft_GP2 <- function(fs){
# Simulate a surface of the Gaussian process based on FFT. The function 
# uses the equation above Eq 2.44 in page 61 of Rue and Held, 2005) 
#
# Input
#  fs: a list created by setupFFTgp2
#    bcb_vec: the base of a block circulant covariance matrix
#    xidx: x index of the interesting matrix on a torus   
#    yidx: y index of the interesting matrix on a torus
#
#  References
#   * Rue, H. and Held, L., 2005. Gaussian Markov random fields: theory and applications. CRC press.
#

# Get eigenvalues using the circulant base
Lambda <- Re(fft(fs$bcb_vec))

# Check if there are negative eigenvalues
if (sum(Lambda<0)>0){

  # Calcuate the rho values of Wood and Chan (1994) 
  rho1 <- rho_criteria(Lambda, mbar=fs$dimmy$mbar)
  
  if (rho1$rho1 > 0.999) {  # SDF: should this be based on sigma 2?  
    warning("Negative eigenvalues appear but are negligible. Those values are set to zero.")
  } else {
    stop("Negative eigenvalues are significant. One solution is to increase the minimum exponent")
  }
  # SDF: moved this line outside of if statement.  The negatives will have to be dealt with irrespective of whether rho1 is large or not
  Lambda[Lambda<0] <- 0 
  # SDF: the adjustment in Wood and Chan
  Lambda <- (rho1$rho1^2) * Lambda
}

# Generate ZZ from N(0,1)
ZZ <- matrix(stats::rnorm(length(fs$bcb_vec)),
         nrow=nrow(fs$bcb_vec),ncol=ncol(fs$bcb_vec))

# IDFT(ZZ) 
FhZZ <- fft(ZZ,inverse=TRUE)/length(fs$bcb_vec)

# GP surface
GPall <- Re(fft(sqrt(Lambda)*FhZZ))

# get the right GP surface
cov_new <- GPall[fs$xidx, fs$yidx]

return(cov_new)

}

# setupFFTgp2()
setupFFTgp2 <- function (x, y, factor = 2,  sig2 = 0.5, rho=0.1, nu=1/2, 
    min_exponent2=12) {

  # extend the vectors on torus to have length at least 2^min_exponent2
  xe <- extend2(x, factor=factor, min_exponent2=min_exponent2, padEndOnly = TRUE)
  ye <- extend2(y, factor=factor, min_exponent2=min_exponent2, padEndOnly = TRUE)

  # get indices to true vectors
  xidx <- seq_range(attr(xe, 'idx'))
  yidx <- seq_range(attr(ye, 'idx'))

  # get the basis of the bloack circulant covariance matrix
  bcb_vec <- matrix(bcb(x=xe, y=ye, sig2=sig2, rho=rho, nu=nu),length(xe), length(ye))
  
  # SDF: get the basic dimensions of the grid
  dimmy <- list( m1=length( xidx), m2=length( yidx), mbar=length( xidx) * length( yidx))
  
  # return as a named list for use in each iteration
  list(bcb_vec = bcb_vec, 
       xidx = xidx,
       yidx = yidx,
       dimmy=dimmy)
}


# fftGPsim2()
fftGPsim2 <- function(x, y, sig2 = 1, rho = 0.5, nu = 1/2, nugget = NULL,
   min_exponent2=NULL) {

  if( is.null( min_exponent2))
    min_exponent2 <- ceiling(log2(2*(max(length(x),length(y))-1)))

  fs <- setupFFTgp2(x = x, y = y, factor = 2, sig2 = sig2, rho = rho, nu = nu,
         min_exponent2=min_exponent2)

  sim <- fft_GP2(fs = fs)

  if(!is.null(nugget))
  sim <- sim + stats::rnorm(length(sim),mean = 0,sd = sqrt(nugget))

  return(sim)

}

rho_criteria <- function (Lambda, mbar=NULL) {
  # Calculate the rho1 value of Wood and Chan (1994). This value 
  # is used to approxmate the embedding matrix when negative eigenvalues
  # occur and are set to zero (see Section 4). 
  #
  #  References
  #   * Wood, A. T., & Chan, G. (1994). Simulation of stationary Gaussian 
  #     processes in [0, 1] d. Journal of computational and graphical 
  #     statistics, 3(4), 409-432.
  #

  # indices of positive eigenvalues 
  aa <- (Lambda>0)
  
  # rho1 = tr(Gamma)/tr(Gamma+)
  rho1 <- sum(Lambda)/sum(Lambda[aa])

  # rho2 = sqrt(rho2)
  #rho2 <- sqrt(rho1)

  # SDF: added for calculation of sigma2
  negEigen <- -sum(Lambda[!aa]) 
  posEigen <- sum(Lambda[aa])
  allEigen <- sum(Lambda)
  
  if( !is.null( mbar))
    sig2 <- allEigen * negEigen / ( mbar*posEigen)
  else
    sig2 <- NA

  return(list(rho1=rho1, negEigen=negEigen, posEigen=posEigen, allEigen=allEigen, sig2=sig2)) 
}
