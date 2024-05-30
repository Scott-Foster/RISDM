
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

makeControl <- function( contr, covar.ext) {

  #number of cores to use.  Default is greedy but not super-super greedy
  if( ! "n.threads" %in% names( contr))  
    contr$n.threads <- parallel::detectCores()-1
  #what are the coordinates called in the observation list data frames and the rasters.
  if( ! "coord.names" %in% names( contr))
    contr$coord.names <- c("Easting","Northing")
  #prior list
  if( ! "prior.list" %in% names( contr))
    contr$prior.list <- NULL
  if( is.null( contr$prior.list)){
    #priors for the intercepts -- location
    if( ! "prior.mean" %in% names( contr))
      contr$prior.mean <- 0
    #priors for the intercepts -- stdev and precision
    if( !"int.sd" %in% names( contr)){
      if( !"int.prec" %in% names( contr)){
	contr$int.sd <- 1000
	contr$int.prec <- 1 / (contr$int.sd^2)
      }
      if( "int.prec" %in% names( contr))
	contr$int.sd <- sqrt( 1 / contr$int.prec)
    }
    else  #if sd is specified then it overrides prec.
      contr$int.prec <- 1 / ( contr$int.sd^2)
    #priors for the slopes -- stdec & precision
    if( !"other.sd" %in% names( contr)){
      if( !"other.prec" %in% names( contr)){
	contr$other.sd <- 10
	contr$other.prec <- 1 / ( contr$other.sd^2)
      }
      if( "other.prec" %in% names( contr))
	contr$other.sd <- sqrt( 1 / contr$other.prec)
    }
    else
      contr$other.prec <- 1 / ( contr$other.sd^2)
  }
  else
    contr$other.prec <- contr$other.sd <- contr$int.prec <- contr$int.sd <- prior.mean <- NULL
  #Should the info. crit be calculated.  Default is no.
  if( ! "calcICs" %in% names( contr))
    contr$calcICs <- FALSE
  #prior for the spatial range of the random effect
  #default is suitable for a unit square, possibly
  if( !"prior.range" %in% names( contr)){
    maxDist <- sqrt( ( covar.ext[2] - covar.ext[1])^2 + ( covar.ext[4] - covar.ext[3])^2)
    contr$prior.range <- c( maxDist/25, 0.1)
  }
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
  if( !"standardiseCovariates" %in% names( contr))
    contr$standardiseCovariates <- TRUE
  if( !"na.action" %in% names( contr))
    contr$na.action <- na.omit
  if( !"inla.mode" %in% names( contr))
    contr$inla.mode <- "compact"  #other option is "classic"
  if( !"re.constr" %in% names( contr))
    contr$re.constr <- TRUE
  if( !"converg.tol" %in% names( contr))
    contr$converg.tol <- 1e-8	#tighter than usual (advice from H. Rue)
  if( !"vb.correction" %in% names( contr))
    contr$vb.correction <- FALSE
  
  
  #add to as we go along
  #check remaining input -- user may have specified stuff that is not there accidentally.
  if( !all( names( contr) %in% c("n.threads","tag.pred","spat.index", "coord.names", "verbose",
                                 "prior.list", "prior.mean","int.prec","other.prec", "int.sd","other.sd","calcICs", "prior.range", 
                                 "prior.space.sigma", "addRandom", "returnStack", "DCmethod", "standardiseCovariates", "na.action", 
				 "inla.mode", "re.constr", "converg.tol", "vb.correction")))
    warning( "There are control parameters specified that are not used.")
  return( contr)
}


