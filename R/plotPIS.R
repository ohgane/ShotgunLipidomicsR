#' Plot PIS data set (comparison of two spectra)
#' 
#' Plot PIS (precursor ion scanning) data
#' 
#' @param pisdata A data frame constructed from xcmsRaw object via \code{constructPIS()} function. 
#' @param col A vector of colors for lines. Default, blue and red.
#' @param lty Types of lines. Same as base graphics.
#' @param mzrange The m/z range to plot. A numeric vector of the length 2.
#' @param fill.col A vector of colors for area under the curve. Default, transparent colors.
#' 
#' @export
#' 
plotPIS=function(pisdata, col=c("blue", "red"), lty=c(1,1),
                 mzrange=NULL, 
                 fill.col=adjustcolor(c("skyblue","salmon"), alpha.f=0.5), ...){
  # Crop mzrange
  if (!is.null(mzrange)){
    pisdata=subset(pisdata, (mzrange[1]<mz)&(mz<mzrange[2]))
  }
  # Plot blank axis with appropriate xlim and ylim
  plot(pisdata$mz, pisdata$int, type="n",
       xlab="m/z", ylab="Relative intensity (%)", ...)
  # Extract group label (fragment m/z)
  frag=levels(pisdata$fragment)
  # Plot the first fragment PIS
  subdat=subset(pisdata, fragment==frag[1])
  polygon(subdat$mz, subdat$int, col=fill.col[1], border=NA)
  points(subdat$mz, subdat$int, type="l", col=col[1], lty=lty[1])
  # Plot the second fragment PIS
  subdat=subset(pisdata, fragment==frag[2])
  polygon(subdat$mz, subdat$int, col=fill.col[2], border=NA)
  points(subdat$mz, subdat$int, type="l", col=col[2], lty=lty[2])
}