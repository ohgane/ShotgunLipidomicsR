#' Plot averaged MS2 spectrum for a selected m/z or scans
#'
#' From xcmsRaw object, this function allows plotting MS2 spectrum for a selected m/z
#'    or selected scans, with peak labels and PrecursorMz line. Additionally, common noise
#'    peaks can be excluded from the peak labels.
#'
#' @param xcmsraw xcmsRaw object
#' @param mz Precursor m/z value to plot (default NULL). Either one of the  \code{mz}
#'      or \code{scanNum} is used.
#' @param scanNum Scan numbers to plot (default NULL). A numeric vector.
#' @param mzabs Tolerance for scan selection in absolute m/z values.
#' @param digits Digits for peak label.
#' @param mzrange m/z range to plot.
#' @param halfWindowSize See \code{MALDIquant::labelPeaks()} for detail.
#' @param SNR See \code{MALDIquant::labelPeaks()} for detail.
#' @param top The number of peaks to label (default 10).
#' @param noisePeaks A numeric vector. Peaks that should not be labeled.
#'
#' @return Invisibly (only when assigned to a variable) returns averaged spectrum
#'      as a data frame. If no MS/MS spectrum is found for the precursor m/z or
#'      scan numbers, NULL is returned.
#'
#' @export
plotAvgMS2=function(xcmsraw, mz=NULL, scanNum=NULL, mzabs=0.1, digits=4, mzrange=NULL,
                     halfWindowSize=NULL, SNR=5, top=10, noisePeaks=NULL, ...){
  require(MALDIquant)
  require(xcms)
  # Require modified version of getSpec()
  # First, deepCopy xcmsRaw object, because msn2xcmsRaw() alter the original input data.
  xraw=deepCopy(xcmsraw)
  xms2=msn2xcmsRaw(xraw) # Copy MS2 data to xcmsRaw class's main slot

  # plot MS2 based on scan number or m/z values
  if (is.null(scanNum)) {
    scanNum=findScans(mz, xcmsraw, mzabs=mzabs)$ScanNum
  }
  if (is.null(scanNum)){
    return(NULL)
  }

  # extract index based on the provided "actual" scan number (= msnAcquisitionNum)
  MS2scan=which(xms2@msnAcquisitionNum %in% scanNum)
  # Average precursor m/z
  PrecursorMz=round(mean(xraw@msnPrecursorMz[MS2scan]), digits=digits)

  # Extract raw data and average over scans
  spec=data.frame(getSpec_modified(xms2, scanrange=MS2scan, ...))
  #spec=data.frame(getSpec(xms2, scanrange=MS2scan))
  #spec[is.na(spec$intensity),"intensity"]=0 # NA has been already removed in getSpec_modified
  # Crop data if mzrange is provided
  if (!is.null(mzrange)){
    spec=spec[mzrange[1] < spec$mz & spec$mz < mzrange[2], ]
  }
  # Convert to MALDIquant MassSpectrum class (S4 class object).
  #     Use mass() or intensity() to access to the m/z and intensity slots
  spec=createMassSpectrum(mass=spec$mz, intensity=spec$intensity/1000)
  # set peak picking window
  if(is.null(halfWindowSize)) {
    halfWindowSize=max(round(length(spec)/200), 3)
  }

  # Perform peak picking
  peaks=detectPeaks(spec, method="MAD", halfWindowSize=halfWindowSize, SNR=SNR)
  # Logical vector (TRUE for noise peaks)
  if (is.null(noisePeaks)){
    peakint=intensity(peaks)
    topLogical= peakint %in% sort(peakint, decreasing=TRUE)[1:top]
  } else {
    Noise=fuzzyMatch(query=mass(peaks), ref=noisePeaks, mzabs=0.1)
    peakint=intensity(peaks)
    topLogical=peakint %in% sort(peakint, decreasing=TRUE)[1:top]
    topLogical=topLogical & !Noise
  }
  inputfile=strsplit(xraw@filepath@.Data, split="/", fixed=TRUE)[[1]]
  inputfile=inputfile[length(inputfile)]

  ### Plot MS/MS spectrum
  # If mzrange is not manually provided, use precursor mz or max mz of extracted scan
  #    (in case of multiply charged ions, precursor mz < max mz of extracted scan)
  # Else, use the supplied mzrange
  if (is.null(mzrange)){
    mzmax=max(PrecursorMz*1.05, max(mass(spec))*1.05)
    mzmin=min(mass(spec))*0.95
  } else {
    mzmax=mzrange[2]
    mzmin=mzrange[1]
  }

  plot(spec, xlim=c(mzmin, mzmax),
       ylim=c(0, max(intensity(spec)*1.2)),
       xlab="m/z", ylab="Intensity (/1000)", ...)
  labelPeaks(peaks, index=topLogical, underline=FALSE, col="blue",
             avoidOverlap=TRUE, srt=0, digits=digits)
  abline(v=PrecursorMz, col="blue")
  title(main=paste("Precursor m/z: ", PrecursorMz, " Scan#: ", paste(scanNum, collapse="-")),
        sub=inputfile,
        col.main="gray50", cex.main=1, font.main=1, col.sub="gray50", cex.sub=0.8, font.sub=1)

  # Return a peaklist with precursor m/z, for LipidBlast search.
  peak.df=list(
    precursor=PrecursorMz,
    spec=data.frame(mz=peaks@mass, intensity=peaks@intensity)
  )
  invisible(peak.df)
}
