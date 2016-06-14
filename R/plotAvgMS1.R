#' A function to plot averaged MS1 spectrum from a \code{xcmsRaw} object
#'
#' This function accepts \code{xcmsRaw} object with MS/MS data (direct infusion,
#'    data-dependent acquisition) and plot the spectrum and overlay colored lines over the precursors. Additionally, label some selected peaks.
#'
#' @param xcmsraw \code{xcmsRaw} object.
#' @param normalization Type of normalization for plotting.
#'     Default, constant. Either constant, relative, or TIC (relative to total ion current).
#' @param int.divide Intensity normalization factor (simple division). Default, 1000.
#' @param col Default, 1 (black).
#' @param type Default, l (line).
#' @param plot.precursor Logical (default, TRUE). If true, plot lines at the precursor m/z values.
#' @param custom.precursor Default, NULL.
#' @param col.precursor Default, light blue.
#' @param lwd.precursor Default, 2.
#' @param lty.precursor Default, 1.
#' @param baseline Logical (default TRUE)
#' @param col.baseline The color of baseline. Default, gray50
#' @param lty.baseline Line type of baseine. Default, 1.
#' @param scans A numeric vector of scan indexes. If provided, only the indicated scans are averaged and plotted.
#' @param mzrange m/z range to plot.
#' @param ylim Y axis range to plot.
#' @param ad.lab. mzrange m/z range to add additional peak labels.
#' @param ad.top The number of additional peak labels.
#' @param top The number of peak labels.
#' @param halfWindowSize See \code{labelPeaks} in \code{MALDIquant} for detail.
#' @param SNR See \code{labelPeaks} in \code{MALDIquant} for detail.
#' @param digits Digits for peak labels.
#'
#' @return Plot spectrrum. Plus, if assigned to a variable (invisibly),
#'    returns the averaged spectrum as a data frame.
#'
#' @export
plotAvgMS1=function(xcmsraw, normalization="constant", int.divide=1000,
                    col=1, type="l",
                    plot.precursor=TRUE,
                    custom.precursor=NULL,
                    col.precursor="light blue",
                    lwd.precursor=2, lty.precursor=1,
                    baseline=FALSE,
                    col.baseline="gray50", lty.baseline=1,
                    scans=NULL, mzrange=NULL, ylim=NULL,
                    ad.lab.mzrange=NA, ad.top=3,
                    top=12, halfWindowSize=NULL, SNR=50, digits=4, ...){
  require(MALDIquant)
  require(xcms)
  # Intensity normalization with a fixed value or TIC
  if (normalization %in% c("constant", "TIC")){
    if (normalization=="TIC"){
      if (is.null(scans)) {
        int.divide=mean(xcmsraw@tic)/100
      } else {
        int.divide=mean(xcmsraw@tic[xcmsraw@acquisitionNum %in% scans])/100
      }
    }
  } else if (normalization=="relative") {
    int.divide=1
  } else {
    warning("The argument normalization must be either constant or TIC.
            Default value (constant, int.divide=1000) is used for intensity normalization.")
  }

  # Extract averaged spectrum from a xcmsRaw object
  if (is.null(scans)){
    # getSpec returns matrix with mz and intensity columns
    spec=data.frame(getSpec_modified(xcmsraw, scanrange=NULL, mzrange=mzrange, ...))
    NumScans=length(xcmsraw@acquisitionNum)  # Number of scans to be used for plot title
  } else {
    # Extract averaged spectrum from the indicated scans (a vector of acquisitionNumber (= real scan))
    spec=data.frame(getSpec_modified(xcmsraw, mzrange=mzrange, scanrange=scans, ...))
    NumScans=length(scans)
  }

  spec[is.na(spec$intensity), "intensity"]=0 # Remove NA by fillin in 0
  # Extract mzrange from xcmsRaw object if mzrange is not supplied
  if(is.null(mzrange)) {
    mzrange=xcmsraw@mzrange
  } else {
    # If mzrange is supplied, crop the data
    spec=spec[mzrange[1] < spec$mz & spec$mz < mzrange[2], ]
  }
  # Copy averaged MS1 spectrum (to return upon exit)
  spec_original=spec
  # Normalization (for "relative")
  if (normalization=="relative"){
    spec$intensity=spec$intensity/max(spec$intensity)*100
  }
  # Scale intensity
  spec$intensity=spec$intensity/int.divide
  # Set plot range for intensity axis (range expanded by 5%), if ylim is not provided
  if (is.null(ylim)){
    intrange=range(spec$intensity)
    intrange.expanded=intrange+c(-0.01, 0.1)*diff(intrange)
  } else {
    intrange.expanded=ylim/int.divide
  }
  #
  if (normalization == "TIC") {
    ylabel="Intensity / Total Ion Current (%)"
  } else if (normalization == "relative") {
    ylabel = "Relative intensity"
  } else {
    ylabel=paste("Intensity/", int.divide, " (a.u.)", sep="")
  }
  # Plot axis first
  plot(NULL, xlim=mzrange, ylim=intrange.expanded,
       xlab="m/z", ylab=ylabel, cex.axis=0.9)
  # Color precursors
  if (is.null(custom.precursor)){
    precursor=xcmsraw@msnPrecursorMz
  } else {
    if (is.numeric(custom.precursor)) {
      precursor=custom.precursor
    } else {
      print("The argument custom.precursor should be a vector of m/z values. None plotted.")
    }
  }
  if (baseline) abline(h=0, col=col.baseline, lty=lty.baseline)
  if (plot.precursor){
    abline(v=precursor, col=col.precursor, lwd=lwd.precursor, lty=lty.precursor)
  }


  # Plot spectrum
  if (type=="l"){
    lines(spec$mz, spec$intensity, col=col, ...)
  } else {
    points(spec$mz, spec$intensity, col=col, type=type, ...)
  }

  # Convert data.frame to MALDIquant::MassSpectrum class
  spec=createMassSpectrum(mass=spec$mz, intensity=spec$intensity)
  # Perform peak picking
  if(is.null(halfWindowSize)) halfWindowSize=round(length(spec)/100) # set peak picking window
  peaks=detectPeaks(spec, method="MAD", halfWindowSize=halfWindowSize, SNR=SNR)

  if (length(peaks) == 0){
    warning("No peak was detected. Too high SNR (default 50) or too small mzrange?")
  } else {
    #MassPeaks object. Use mass() and intensity() to access the contents
    topN=MALDIquant::intensity(peaks) %in% sort(MALDIquant::intensity(peaks),
                                                decreasing=TRUE)[1:top]
    ### Additional label
    if(is.numeric(ad.lab.mzrange)){
      mz.ad=ad.lab.mzrange[1]<mass(peaks) & mass(peaks) < ad.lab.mzrange[2] & !topN
      topN.ad=intensity(peaks) %in% sort(intensity(peaks[mz.ad,]), decreasing=TRUE)[1:ad.top]
      topN=topN|topN.ad
    }
    # T/F vector of the length of "peaks"
    labelPeaks(peaks, index=topN, underline=FALSE, col="blue",
               avoidOverlap=TRUE,  digits=digits, arrowCol="blue")
  }

  if (normalization == "TIC"){
    tic=formatC(int.divide*100, digits=1, format="e")
    title(main=paste("Averaged MS1 (", NumScans, "scans, TIC:", tic, ")"),
          sub=xcmsraw@filepath@.Data,
          col.sub="gray50", col.main="gray50", cex.main=1, cex.sub=0.8, font.sub=1, font.main=1)
  } else {
    title(main=paste("Averaged MS1 Spectrum (", NumScans, "scans )"),
          sub=xcmsraw@filepath@.Data,
          col.sub="gray50", col.main="gray50", cex.main=1, cex.sub=0.8, font.sub=1, font.main=1)
  }
  invisible(spec_original)
}
