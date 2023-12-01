
###############################################################################################
###############################################################################################
####	
####	Control for the simulation function
####
####	Returns a list
####
####	Programmed by Scott in the first half of 2022
####
###############################################################################################
###############################################################################################


set.sim.control <- function( contr){

  #should the random seed be set by user
  if( ! "set.random.seed" %in% names( contr))
    contr$set.random.seed <- FALSE
  #what random seed
  if( ! "random.seed" %in% names( contr))  
    contr$random.seed <- 787
  #what raster dimension to simulate over
  if( ! "raster.dim" %in% names( contr))
    contr$raster.dim <- rep(50,2)
  #plot the simulation data?
  if( ! "doPlot" %in% names( contr))
    contr$doPlot <- TRUE
  #should the random effect be added to the simulation
  if( ! "addRandom" %in% names( contr))
    contr$addRandom <- TRUE
  #what is the re std dev?
  if( ! "sd" %in% names( contr))
    contr$sd <- 0.5
  #what is the re effective range?  This is completely dependent on the spatial region
  if( ! "range" %in% names( contr))
    contr$range <- NA
  if( ! "useMBHdesign" %in% names( contr))
    contr$useMBHdesign <- FALSE
  if( ! "exact.n.PO" %in% names( contr))
    contr$exact.n.PO <- FALSE
  
  #add to as we go along
  if( !all( names( contr) %in% c("set.random.seed","random.seed","raster.dim","doPlot","addRandom", "sd", "range", "useMBHdesign","exact.n.PO")))
    warning( "There are control parameters specified that are not used.")
  return( contr)
  
}

