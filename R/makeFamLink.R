
###############################################################################################
###############################################################################################
####	
####	Set the family and link functions etc for the various parts of the model.
####
####	Returns a list
####
####	Programmed by Scott in the first half of 2022
####	Based originally on Dambly et al. (2019)'s code. 
####
###############################################################################################
###############################################################################################

makeFamLink <- function( ind, nObs) {

  #the names (datatype) contained within the data
  nammy <- names( ind)

  #vectors for the families and links
  fammy <- rep( NA, length( nammy))
  names( fammy) <- nammy
  linkky <- list()

  #each of the data types in turn.
  if( ind["PO"]==1){
    fammy["PO"] <- "poisson"
    linkky[["PO"]] <- list(link="log")
  }
  if( ind["PA"]==1){
    fammy["PA"] <- "binomial"
    linkky[["PA"]] <- list(link="cloglog")
  }
  if( ind["AA"]==1){
    fammy["AA"] <- "poisson"
    linkky[["AA"]] <- list(link="log")
  }
  if( ind["DC"]==1){
    fammy["DC"] <- "poisson"
    linkky[["DC"]] <- list(link="log")
  }
  #remove the NA entries
  fammy <- fammy[!is.na( fammy)]

#  #the index of a log link, if PO AA data present.  Otherwise its gotta be cloglog.
#  linkID <- NA  
#  if( is.na( linkID) & ( "PO" %in% names( fammy)))
#    linkID <- which( names( fammy) == "PO")
#  if( is.na( linkID) & ( "AA" %in% names( fammy)))
#    linkID <- which( names( fammy) == "AA")
#  if( is.na( linkID) & ( "PA" %in% names( fammy)))
#    linkID <- which( names( fammy) == "PA")
#  if( is.na( linkID) & ( "DC" %in% names( fammy)))
#    linkID <- which( names( fammy) == "DC")

  linkID <- rep( 1:length( nObs), times=nObs)

  #to satisfy one of inla's internal checks.  Sigh...
  names( linkky) <- NULL  

  res <- list( fam=fammy, link=linkky, linkID=linkID)

  return( res)
}

