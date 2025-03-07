###############################################################################################
###############################################################################################
####	
####	Set control parameters for population estimation
####
####	Returns a list of predictions
####
####	Programmed by Scott in August 2024
####
###############################################################################################
###############################################################################################

makePopEstiControl <- function( cntr){

  if( !"winsor" %in% names( cntr))
    cntr$winsor <- TRUE
  if( !"tail" %in% names( cntr))
    cntr$tail <- "upper"
  if( !"percent" %in% names( cntr))
    cntr$percent <- 0.01
  if( !"probs" %in% names( cntr))
    cntr$probs<-c(0.025,0.975)
    
  if( !all( names( cntr) %in% c("probs","winsor","tail","percent")))
    warning( "There are control parameters specified that are not used.")
    
  return( cntr)
}
