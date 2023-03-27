

removeParsingChars <- function( yy){
  
    yy <- gsub( ":", ".", yy, fixed=TRUE)
    yy <- gsub( "/", ".", yy, fixed=TRUE)
    yy <- gsub( "*", ".", yy, fixed=TRUE)
    yy <- gsub( "-", ".", yy, fixed=TRUE)
    yy <- gsub( "+", ".", yy, fixed=TRUE)
    yy <- gsub( "%in%", ".", yy, fixed=TRUE)
    yy <- gsub( "^", ".", yy, fixed=TRUE)
    yy <- gsub( "(", ".", yy, fixed=TRUE)
    yy <- gsub( ")", ".", yy, fixed=TRUE)
    yy <- gsub( ", ", ".", yy, fixed=TRUE)
    yy <- gsub( "=", ".", yy, fixed=TRUE)
    yy <- gsub( " ", "", yy, fixed=TRUE)
    
    return( yy)
}
