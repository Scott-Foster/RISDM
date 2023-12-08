
###############################################################################################
###############################################################################################
####	
####	Calcuate -logl for the very simple observation/detection model for DC data
####
####	Returns a scalar double
####
####	Programmed by Scott in the first half of 2022
####
###############################################################################################
###############################################################################################

DC_loglikFun_CommonPi <- function( par, dat){
  #estimating log (pi), so that is what is passed in.
  par <- exp( par)
  #adjusting exact zeros and ones
  par <- pmax( pmin( par, 1-1e-4), 1e-4)
  #the various combinations of obs1 (not obs2), obs2 (not obs1) and both
  lambda <- matrix( NA, nrow=nrow( dat), ncol=3)
  lambda[,1] <- par*(1-par)
  lambda[,2] <- (1-par)*par
  lambda[,3] <- par*par
 
  #the multinomial contributions from each observation and each category
  logl_ij <- dat * log( lambda)
  #the contribution from each observation
  logl_i <- rowSums( logl_ij)
  #the overal logl.
  logl <- sum( logl_i)
  
  return( -logl)
}
