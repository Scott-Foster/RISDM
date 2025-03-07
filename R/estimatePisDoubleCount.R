
###############################################################################################
###############################################################################################
####	
####	Estimate pis from a single double count survey.
####
####	Returns a legnth2 vector of estimated probs.
####
####	Programmed by Scott in the first half of 2022
####	Re-factored by Scott Feb 2024
#### 		to have only 1 prob, shared for both observers (sometimes the labels are
####		swapped and observers change and...)
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
#  start.vals <- colMeans( start.vals, na.rm=TRUE)
  start.vals <- mean( start.vals, na.rm=TRUE)
  #if there are NaNs everywhere, define them
#  start.vals <- ifelse( is.na( start.vals), 0.5, start.vals) #0.5 seems kind of reasonable as it could provide a place to avoid local maxima
  start.vals <- if( is.na( start.vals)) 0.5 else start.vals
#  start.vals <- log( pmax( 0.001, pmin( 0.999, start.vals)))
  start.vals <- log( max( 0.001, min( 0.999, start.vals))) #avoid boundaries
  
  tmp <- stats::optim( par=start.vals, fn=DC_loglikFun, dat=dat, method='Brent', upper=0, lower=log(0.01))  #revert to nelder-mead if there are two parameters.
  if( tmp$convergence != 0)
    warning("Optimisation for expansion points/probabilities has not worked properly (convergence code != 0)")
  expans.pts <- exp( tmp$par)  
  
  return( expans.pts)
}

