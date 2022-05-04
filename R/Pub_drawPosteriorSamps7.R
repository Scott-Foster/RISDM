
##############################################################
####  Function to draw posterior samples from an INLA object
####  Taken from online version, and put into function
####
####  Tidied (a very little bit) 8-July-2020
####
####  Scott Foster (scott.foster@data61.csiro.au)
##############################################################


#function to draw samples from posterior (hyperparams independet of fixed, unfortunately)
draw.posterior.samps <- function(inla.fm, B=100, what="params", field="iSpat") 
{
  #Arguemtn "what" is either "params", in which case function will return the fixed effects and hyper-params (all in one matrix)
  #or it is "effects" in which case it returns a matrix of fixed effects and a matrix of random effects.
  if( what=="params"){
    p.fixed <- nrow( inla.fm$summary.fixed)
    p.var <- 1
    p.dep <- 1
    p.smooth <- 1
    p.tot <- p.fixed + p.var + p.dep + p.smooth
  
    draws <- matrix( NA, nrow=B, ncol=p.tot)
    colnames( draws) <- c( rownames( inla.fm$summary.fixed), "Spatial.Variance", "kappa", "Smooth")
  
    tmp <- INLA::inla.posterior.sample( n=B, result=inla.fm, intern=FALSE, add.names=FALSE)
    latent.ids.wanted <- tail( rownames( tmp[[1]]$latent), p.fixed)
  
    for( ii in 1:p.fixed)
      draws[,ii] <- sapply( tmp, function(x) x$latent[rownames( tmp[[1]]$latent) == latent.ids.wanted[ii],])
    draws[,"Smooth"] <- 2 - 2/2 #alpha was 2 from above.
    ranges <- sapply( tmp, function(x) x$hyperpar[grep( "Range", names( tmp[[1]]$hyperpar))])
    sds <- sapply( tmp, function(x) x$hyperpar[grep( "Stdev", names( tmp[[1]]$hyperpar))])
    draws[,c("Spatial.Variance")] <- sds^2
    draws[,"kappa"] <- sqrt( 8*draws[,"Smooth"]) / ranges
    return( draws)
  }
  #else not needed as return already called.
  p.fixed <- nrow( inla.fm$summary.fixed)
  #take the sample
  tmp <- INLA::inla.posterior.sample( n=B, result=inla.fm, num.threads=parallel::detectCores(), add.names=FALSE)
  #table( Reduce( rbind, strsplit(rownames( tmp[[1]]$latent), ":"))[,1])
  iSpat.id <- grep( field, rownames( tmp[[1]]$latent))
  res  <- list()
  res$fieldAtNodes <- sapply( tmp, function(x) x$latent[iSpat.id])
  res$fixedEffects <- sapply( tmp, function(x) tail( x$latent, p.fixed))
  
  return( res)
  
}
