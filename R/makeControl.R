
###############################################################################################
###############################################################################################
####	
####	Organise the control parameters and set defaults.
####
####	Returns a list for controling isdm 
####
####	Programmed by Scott in the first half of 2022
####
###############################################################################################
###############################################################################################

makeControl <- function( contr) {

  #number of cores to use.  Default is greedy but not super-super greedy
  if( ! "n.threads" %in% names( contr))  
    contr$n.threads <- parallel::detectCores()-1
  #what are the coordinates called in the observation list data frames and the rasters.
  if( ! "coord.names" %in% names( contr))
    contr$coord.names <- c("Easting","Northing")
  #priors for the intercepts -- location
  if( ! "prior.mean" %in% names( contr))
    contr$prior.mean <- 0
  #priors for the intercepts -- precision
  if( ! "int.prec" %in% names( contr))
    contr$int.prec <- 0.0001
  #priors for the slopes -- precision
  if( ! "other.prec" %in% names( contr))
    contr$other.prec <- 0.0001
  #Should the info. crit be calculated.  Default is no.
  if( ! "calcICs" %in% names( contr))
    contr$calcICs <- FALSE
  #prior for the spatial range of the random effect
  #default is suitable for a unit square, possibly
  if( !"prior.range" %in% names( contr))
    contr$prior.range <- c(0.2, 0.1)
  #prior for the variance of the spatial rand effect
  if( !"prior.space.sigma" %in% names( contr))
    contr$prior.space.sigma <- c(1, 0.01)
  #should the INLA call be run as verbose?
  if( !"verbose" %in% names( contr))
    contr$verbose <- FALSE
  #should the spatial random effect be included in the model.
  #default is 'yes'.  If FALSE, then priors for random effect etc will be ignored.
  if( !"addRandom" %in% names( contr))
    contr$addRandom <- TRUE
  #should the INLA stack be returned with the object
  if( !"returnStack" %in% names( contr))
    contr$returnStack <- TRUE
  #which DC method should be used.  plugin is the other option (pre-estimate detection.
  if( !"DCmethod" %in% names( contr))
    contr$DCmethod <- "TaylorsLinApprox"  #the other option is "plugin"
  
  #add to as we go along
  #check remaining input -- user may have specified stuff that is not there accidentally.
  if( !all( names( contr) %in% c("n.threads","tag.pred","spat.index", "coord.names", "verbose",
                                 "prior.mean","int.prec","other.prec", "calcICs", "prior.range", 
                                 "prior.space.sigma", "addRandom", "returnStack", "DCmethod")))
    warning( "There are control parameters specified that are not used.")
  return( contr)
}


