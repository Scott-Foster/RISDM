
###############################################################################################
###############################################################################################
####	
####	Estimate detection probability using the method of moments when the observed data 
####  ("obs1", "obs2", "both") follows a multinominal distribution with parameters 
####  (N, p_1 x (1-p_2), (1-p_1) x p_2, p_1 x p_2).   
####
####	Return 
####   * a vector of two probabilities when p_1 != p_2 
####   * one probability when p_1 = p_2 
####
####	Programmed by Wen-Hsi (26/5/2025)
####
###############################################################################################
###############################################################################################

DC_detProbMM <- function( dat, is_p1_p2_same = TRUE){
  #dat is an nx3 matrix with colnames c("obs1","obs2","both")
  #data should already be ordered by the way that this function is called...
  
  if( sum( as.matrix( dat)) <= 0)
    stop( "At least one of the DC surveys hasn't seen a single individual in any transect by either observer.\n There is no way that detection probabilities can be estiamted.\n Please reconsider model (e.g. add these points as SC as an off-the-cuff suggestion/hack).")

  y1_mean <- mean(dat[,1])
  y2_mean <- mean(dat[,2])
  y3_mean <- mean(dat[,3])

  if (is_p1_p2_same) {
    # when p_1 = p_2
    expans.pts <- (2*y3_mean)/(y1_mean + y2_mean + 2*y3_mean)
  } else {
    # when p_1 != p_2
    expans.pts <- c(
      # p_1 
      y3_mean/(y2_mean + y3_mean),
      # p_2 
      y3_mean/(y1_mean + y3_mean)
    )
  }
 
  return( expans.pts)
}

