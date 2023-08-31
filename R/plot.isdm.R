
###############################################################################################
###############################################################################################
####	
####	create diagnostic plots for isdm objects.
####
####	Returns NULL but does so invisibily
####
####	Programmed by Scott in the first half of 2022
####
###############################################################################################
###############################################################################################

plot.isdm <- function( x, nFigRow=1, ask=TRUE, ...){
 
  ressy <- residuals.isdm( x)
 
  #plotting up the residuals. RQR versus fitted and for PO data the raster.
  if( ask){
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit( grDevices::devAskNewPage(FALSE))
  }
  graphics::par( mfrow=c(nFigRow,3), oma=rep( 0,4), mar=c(2,2,2,4))
  
#  oask <- grDevices::devAskNewPage(FALSE) #don't ask for the first page of plots
  if( "DC_Intercept" %in% x$mod$names.fixed){
    plot.new()
    DCresids <- ressy$DC
    plot( DCresids$fitted, DCresids$residual, pch=20, ylab="DC residuals", xlab="DC fitted", main="Double Count")
    graphics::abline( h=0, col='green')
    stats::qqnorm( DCresids$residual, pch=20, ylab="DC quantile", main="Double Count")
    stats::qqline( DCresids$residual, col='green')
  }
  if( "AA_Intercept" %in% x$mod$names.fixed){
    plot.new()
    AAresids <- ressy$AA
    plot( AAresids$fitted, AAresids$residual, pch=20, ylab="AA residuals", xlab="AA fitted", main="Abundance-Absence")
    graphics::abline( h=0, col='green')
    stats::qqnorm( AAresids$residual, pch=20, ylab="AA quantile", main="Abundance-Absence")
    stats::qqline( AAresids$residual, col='green')
  }
  if( "PA_Intercept" %in% x$mod$names.fixed){
    plot.new()
    PAresids <- ressy$PA
    plot( PAresids$fitted, PAresids$residual, pch=20, ylab="PA residuals", xlab="PA fitted", main="Presence-Absence")
    graphics::abline( h=0, col='green')
    stats::qqnorm( PAresids$residual, pch=20, ylab="PA quantile", main="Presence-Absence")
    stats::qqline( PAresids$residual, col='green')
  }
  if( "PO_Intercept" %in% x$mod$names.fixed){
    POresids <- ressy$PO
    terra::plot( POresids$ras, main="Presence-Only", col=colorRamps::matlab.like2(25))
    plot( POresids$POresids$fitted, POresids$POresids$residual, pch=20, ylab="PO residuals", xlab="PO fitted", main="Presence-Only")
    graphics::abline( h=0, col='green')
    stats::qqnorm( POresids$POresids$residual, pch=20, ylab="PO quantile", main="Presence Only")
    stats::qqline( POresids$POresids$residual, col='green')
  }

  invisible( NULL)
  
}
  
