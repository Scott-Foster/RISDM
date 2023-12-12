
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
  
  if( sum( as.matrix( dat)) <= 0)
    stop( "At least one of the DC surveys hasn't seen a single individual in any transect by either observer.\n There is no way that detection probabilities can be estiamted.\n Please reconsider model (e.g. add these points as SC as an off-the-cuff suggestion/hack).")
  
  #note that the next line will produce NaNs when neither observer sees any animals at a site
  start.vals <- cbind( (dat[,1]+dat[,3]), dat[,2]+dat[,3]) / rowSums( dat, na.rm=TRUE)
#  #this line will not produce NaNs and will work for surveys where no koalas are seen (is that sensible?)
#  start.vals <- cbind( (dat[,1]+dat[,3]), dat[,2]+dat[,3]) / pmax( rowSums( dat, na.rm=TRUE), 1e-4)
  start.vals <- colMeans( start.vals, na.rm=TRUE)
  #if there are NaNs everywhere, define them
  start.vals <- ifelse( is.na( start.vals), 0, start.vals)
  start.vals <- log( pmax( 0.001, pmin( 0.999, start.vals)))
#  start.vals <- log( start.vals + 0.01)
  
  tmp <- stats::optim( par=start.vals, fn=DC_loglikFun, dat=dat)
  if( tmp$convergence != 0)
    warning("Nelder-Mead optimisation for expansion points/probabilities has not worked properly (convergence code != 0)")
  expans.pts <- exp( tmp$par)  
  
  return( expans.pts)
}

