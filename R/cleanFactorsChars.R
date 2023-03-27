
###############################################################################################
###############################################################################################
####	Utility function to strip character values that could cause problems with parsing 
####	formulas etc...
####
####	Returns the vector with those terms stripped.
####
####	Programmed by Scott in Marfch 2023
####
###############################################################################################
###############################################################################################

cleanFactorsChars <- function( datList){
  
  for( ll in 1:length( datList)){
    for( ii in 1:ncol( datList[[ll]])){
      if( inherits( datList[[ll]][,ii], "character")){
	datList[[ll]][,ii] <- gsub( " ", "", datList[[ll]][,ii], fixed=TRUE)
	datList[[ll]][,ii] <- gsub( ":", "", datList[[ll]][,ii], fixed=TRUE)
	datList[[ll]][,ii] <- gsub( "*", "", datList[[ll]][,ii], fixed=TRUE)
	datList[[ll]][,ii] <- gsub( "/", "", datList[[ll]][,ii], fixed=TRUE)
	datList[[ll]][,ii] <- gsub( "-", "", datList[[ll]][,ii], fixed=TRUE)
	datList[[ll]][,ii] <- gsub( "+", "", datList[[ll]][,ii], fixed=TRUE)
      }
      if( inherits( datList[[ll]][,ii], "factor")){
	levels( datList[[ll]][,ii]) <- gsub( " ", "", levels( datList[[ll]][,ii]), fixed=TRUE)
	levels( datList[[ll]][,ii]) <- gsub( ":", "", levels( datList[[ll]][,ii]), fixed=TRUE)
	levels( datList[[ll]][,ii]) <- gsub( "*", "", levels( datList[[ll]][,ii]), fixed=TRUE)
	levels( datList[[ll]][,ii]) <- gsub( "/", "", levels( datList[[ll]][,ii]), fixed=TRUE)
	levels( datList[[ll]][,ii]) <- gsub( "-", "", levels( datList[[ll]][,ii]), fixed=TRUE)
	levels( datList[[ll]][,ii]) <- gsub( "+", "", levels( datList[[ll]][,ii]), fixed=TRUE)
      }
    }
  }
  return( datList)
}
#  
#  
#  removeFormTerms <- function(xx){
#    if( inherits( xx, "character")){
#      xx <- gsub( " ", "", xx, fixed=TRUE)
#      xx <- gsub( ":", "", xx, fixed=TRUE)
#      xx <- gsub( "*", "", xx, fixed=TRUE)
#      xx <- gsub( "/", "", xx, fixed=TRUE)
#      xx <- gsub( "-", "", xx, fixed=TRUE)
#      xx <- gsub( "+", "", xx, fixed=TRUE)
#    }
#    if( inherits( xx, "factor")){
#      levels( xx) <- gsub( " ", "", levels( xx), fixed=TRUE)
#      levels( xx) <- gsub( ":", "", levels( xx), fixed=TRUE)
#      levels( xx) <- gsub( "*", "", levels( xx), fixed=TRUE)
#      levels( xx) <- gsub( "/", "", levels( xx), fixed=TRUE)
#      levels( xx) <- gsub( "-", "", levels( xx), fixed=TRUE)
#      levels( xx) <- gsub( "+", "", levels( xx), fixed=TRUE)
#    }
#    return( xx)
#  }
#  
#  tmp <- lapply( datList, function(zz) lapply( as.data.frame( zz), removeFormTerms))
#  return( tmp)
#}
#  
#
