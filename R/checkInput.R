
###############################################################################################
###############################################################################################
####	Check the input arguments to an isdm call.
####	Checks for gross misspecifications, not detailed little things.
####
####	Returns a boolean if the input seems to be plausible (ie no gross errors).
####
####	Programmed by Scott in the first half of 2022
####
###############################################################################################
###############################################################################################

checkInput <- function( respNames, biasForm, arteForm, DCobsInfo, obsList){
  
  #checking the input for obvious inconsistencies
  
  #is the DC information supplied agree in structure to that which is wanted?  
  if( length( DCobsInfo) != 4)
    stop("Please specify DCobserverInfo as a named list of exactly length 4.  See help for guidance.")
  
  #has the user specified *any* data?
  if( all( is.null( obsList)))
    stop("No data specified.  Please do so through observationList argument.")
  
  #are the pieces available for PO?
  if( any( c( is.null( obsList$POdat), is.null( biasForm))))
    if ( !all( c( is.null( obsList$POdat), is.null( biasForm))))
      stop("To include PO data, you *must* supply observationList$POdat (argument) AND a biasFormula")
  #are the pieces available for PA?
  if( any( c( is.null( obsList$PAdat), is.na( respNames["PA"]), is.null( arteForm[["PA"]]))))
    if( !all( c( is.null( obsList$PAdat), is.na( respNames["PA"]), is.null( arteForm[["PA"]]))))
      stop("To include PA data, you *must* supply observationList$PAdat (argument) AND an entry in responseNames AND an entry in artefactFormulas")
  #are the pieces available for AA?
  if( any( c( is.null( obsList$AAdat), is.na( respNames["AA"]), is.null( arteForm[["AA"]]))))
    if( !all( c( is.null( obsList$AAdat), is.na( respNames["AA"]), is.null( arteForm[["AA"]]))))
      stop("To include AA data, you *must* supply observationList$AAdat (argument) AND an entry in responseNames AND an entry in artefactFormulas")
  #are the pieces available for DC?
  if( any( c( is.null( obsList$DCdat), is.null( arteForm[["DC"]])))){#, is.null( DCobsInfo)))){
    if( !all( c( is.null( obsList$DCdat), is.null( arteForm[["DC"]]))))#, is.null( DCobsInfo))))
      stop("To include DC data, you *must* supply observationList$DCdat (argument) AND an entry in artefactFormulas AND the information about observers (DCobserverInfo")
    else
      rm( DCobsInfo)
  }
  
  return( TRUE)  #if execution makes it this far
}
