
###############################################################################################
###############################################################################################
####	
####	Estimate pis from a single double count survey.
####
####	Returns a legnth2 vector of estimated probs.
####
####	Programmed by Scott in the first half of 2022
####
###############################################################################################
###############################################################################################

estimatePisDoubleCount <- function( dat){
  #dat is an nx3 matrix with colnames c("obs1","obs2","both")
  #data should already be ordered by the way that this function is called...
  start.vals <- cbind( (dat[,1]+dat[,3]), dat[,2]+dat[,3]) / rowSums( dat, na.rm=TRUE)
  start.vals <- colMeans( start.vals, na.rm=TRUE)
  start.vals <- log( start.vals + 0.01)
  
  tmp <- stats::optim( par=start.vals, fn=DC_loglikFun, dat=dat)
  if( tmp$convergence != 0)
    warning("Nelder-Mead optimisation for expansion points/probabilities has not worked properly (convergence code != 0)")
  expans.pts <- exp( tmp$par)  
  
  return( expans.pts)
}

