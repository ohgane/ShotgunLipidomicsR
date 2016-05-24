#' A function to pick peaks from an averaged MS1 spectrum from a \code{xcmsRaw} object
#'
#' This function accepts \code{xcmsRaw} object with MS/MS data (direct infusion,
#'    data-dependent acquisition) and return top N peaks. Correspond to \code{plotAvgMS1}.
#'
#' @param xcmsraw \code{xcmsRaw} object.
#' @param scans A numeric vector of scan indexes. If provided, only the indicated scans are averaged and plotted.
#' @param mzrange m/z range to plot.
#' @param ylim Y axis range to plot.
#' @param ad.lab. mzrange m/z range to add additional peak labels.
#' @param ad.top The number of additional peak labels.
#' @param top The number of peak labels.
#' @param halfWindowSize See \code{MALDIquant::labelPeaks()} for detail.
#' @param SNR See \code{MALDIquant::labelPeaks()} for detail.
#' @param digits Digits for peak labels.
#'
#' @return Top N peaks as a data frame.
#'
#' @export
getPeaksAvgMS1=function(xcmsraw,
                     scans=NULL, mzrange=NULL, ylim=NULL,
                     ad.lab.mzrange=NA, ad.top=3,
                     top=12, halfWindowSize=NULL, SNR=50, digits=4, ...){
  require(MALDIquant)
  require(xcms)


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
  }
  # Select Top N peaks from a MassPeaks object
  peaks=peaks[topN,]
  # Convert MassPeaks object into a dataframe
  res_df=data.frame(mz=peaks@mass, intensity=peaks@intensity, snr=peaks@snr)
  return(res_df)
}
