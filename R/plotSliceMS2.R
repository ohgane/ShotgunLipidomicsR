#' Plot "precursor ion scan"- or "neutral loss scan"-like spectrum from DDA xcmsRaw object.
#'
#' Note that multiple scans with the same precursor m/z values are averaged
#'      (after peak integration).
#'
#' @param xraw \code{xcmsRaw} object.
#' @param mz \code{xcmsRaw} object.
#' @param ppm Numeric (default 10 ppm). Integration range. If \code{mzabs} is provided,
#'      this parameter will not be ignored.
#' @param slicemode Character (default "FRAG").Select one from "FRAG" or "NL",
#'      to plot precursor ion spectrum with fragment or neutral loss intensity respectively.
#' @param scanType Character (default "profile"). Either "profile" or "centroid".
#' @param plotSpec Logical (default TRUE). Whether plot pseudo-precursor ion spectrum.
#'      If set to FALSE, only returns pseudo precursor ion spectrum as a matrix.
#' @param type A character vector containing \code{type} argments that are used in standard
#'      \code{plot()}. For example, \code{c("h", "p")} results in points over vertical lines.
#' @param col Colors used to plot multiple spectra. A character vector (default \code{c("blue", "red")})
#' @param pch Same as \code{plot()}. Type of points.
#' @param cex Same as \code{plot()}. Size of points.
#' @param labelPeaks Logical (default TRUE), that determines whether to plot peak labels.
#' @param top The number of peak labels (default 5).
#'      By default, top 5 peaks in each spectrum will be labelled.
#' @param SNR S/N ratio for peak picking (default 5).
#' @param halfWindowSize A parameter for peak picking (default NULL).
#'      By default, automatically determined based on data.
#'      Larger value prevents peak labeling of nearer peaks.
#' @param digits Digits of peak label (default 4).
#'
#' @return Plot spectrrum. Plus, if assigned to a variable, returns a data.frame
#'      with precursor m/z, and intensities (AUC for profile, sum for cetroid).
#'
#' @export

plotSliceMS2=function(xraw, mz, ppm=50, mzrange=NULL, ylim=NULL,
                      slicemode="FRAG", scanType="profile", plotSpec=TRUE,
                      type=c("h", "p"), col=c("blue", "red"), pch=1, cex=0.5,
                      labelPeaks=TRUE, top=5, SNR=5, halfWindowSize=NULL, digits=4, ...){
  require("caTools")
  ### Prepare precursor m/z - intensity data frame
  # Extract precursor m/z and scan number from xcmsRaw
  PrecursorMz=xraw@msnPrecursorMz
  scan=xraw@msnAcquisitionNum
  if (is.null(mzrange)){
    mzrange=range(PrecursorMz)
  }
  if (slicemode=="FRAG"){
    # Calculate Fragment intensity
    int.df=NULL
    for (i in 1:length(mz)){
      # Integration range
      mzabs=mz[i]*ppm/10^6
      # Calculate AUC (caTools::trapz() is used for integration of profiledata)
      #     For centroid data set, set type="centroid", which uses sum() instead of trapz().
      int=sapply(1:length(xraw@msnAcquisitionNum),
                 xcmsrawAuc, mzrange=mz[i]+c(-mzabs, mzabs), xraw=xraw,
                 type=scanType, msLevel=2)
      int.df=cbind(int.df, int)
    }
  } else if (slicemode=="NL") {
    # Calculate Neutral loss peak intensity
    int.df=NULL
    for (i in 1:length(mz)){
      # Integration range
      mzabs=mz[i]*ppm/10^6
      # Calculate AUC
      int=NULL
      for (k in 1:length(xraw@msnAcquisitionNum)){
        int=c(int,
              xcmsrawAuc(k, mzrange=PrecursorMz[k]-mz[i]+c(-mzabs, mzabs), xraw=xraw,
                         type=scanType, msLevel=2)
        )
      }
      int.df=cbind(int.df, int)
    }
  } else {
    stop("The slicemode argument must be either FRAG (to plot a slice with fragments)
         or NL (to plot a slice with neutral loss peaks)")
  }
  colnames(int.df)=paste("Int", round(mz), sep="_")


  ### Prepare spectrum data.frame
  # Plot blank axis with appropriate xlim and ylim
  spec=data.frame(PrecursorMz, int.df)
  spec=na.omit(spec)
  spec=spec[mzrange[1] < spec$PrecursorMz & spec$PrecursorMz < mzrange[2], ]

  # Average scans with the same precursor mz
  spec=summaryBy(.~PrecursorMz, data=spec, FUN=mean)

  ### Plot data
  if (plotSpec){
  if (is.null(ylim)){
    yrange=range(spec[,-1])/1000
    yrange=c(min(yrange)*0.9, max(yrange)*1.1)
  } else {
    yrange=ylim
  }

  # If no peak was found in the mz range, return NULL and exit with a warning.
  if (diff(yrange)==0){
    warning("No precursor ion peak was found in the specified region.")
    return(NULL)
  }


  # Plot spectrum
  plot(mzrange, yrange,
       type="n", xlab="m/z", ylab="Intensity/1000 (a.u.)", ...)
  for (type in type){
    for (i in 1:length(mz)){
      points(spec$PrecursorMz, spec[,i+1]/1000, type=type,
             col=col[i], cex=cex, pch=pch, ...)
      # Check if there is any peak for each mz
      PeakExist=TRUE
      if(diff(range(spec[,i+1]))==0){
        PeakExist=FALSE
      }
      if (labelPeaks==TRUE & PeakExist){
        # Convert to MALDIquant spectrum object (S4 class)
        specMALDI=createMassSpectrum(mass=spec$PrecursorMz, intensity=spec[,i+1]/1000)
        # set peak picking window
        if(is.null(halfWindowSize)) {
          halfWindowSize=max(round(length(spec)/200), 3)
        }
        # If no peak was picked, no need to label peaks.
        # If the number of peaks is less than 5, no need to pick peaks
        #   (to avoid error "halfWindowsize too large/small")
        if (length(specMALDI)==0){
          peaks=NULL
        } else if (length(specMALDI)<5){
          peaks=createMassPeaks(mass=mass(specMALDI), intensity=intensity(specMALDI))
        } else {
          # Perform peak picking
          peaks=detectPeaks(specMALDI, method="MAD", halfWindowSize=halfWindowSize, SNR=SNR)
          # If no peak was found, this returns MassPeaks class object of length==0
          #> head(peaks)
          #S4 class type            : MassPeaks
          #Number of m/z values     : 0
          #Range of m/z values      : NA
          #...
        }
        # Label peaks only when peaks were picked.
        if (length(peaks)!=0){
          topLogical=intensity(peaks) %in% sort(intensity(peaks), decreasing=TRUE)[1:top]
          labelPeaks(peaks, index=topLogical, underline=FALSE, col=col[i],
                     avoidOverlap=TRUE, srt=0, digits=digits)
        } else {
          warning("No peak was picked.")
        }
      }
    }
  }

  # Plot legend at outside of the plot area
  par(xpd=TRUE) # This (global) option allows to plot legend outside of the plot area
  legend(par()$usr[2], par()$usr[4], legend=mz, col=col, lty=1,
         bty="n", bg=NULL, xjust=1, yjust=0, cex=0.8)
  par(xpd=FALSE)

  # Extract file name from xcmsRaw object (the last element of path)
  inputfile=strsplit(xraw@filepath@.Data, split="/", fixed=TRUE)[[1]]
  inputfile=inputfile[length(inputfile)]
  # Plot title
  if (slicemode=="FRAG"){
    titlemode="Precursor ion spectrum (fragment)"
  } else {
    titlemode="Precursor ion spectrum (neutral loss)"
  }
  titlemode=paste(titlemode, " [", ppm, "ppm]", sep="")
  title(main=titlemode,
        sub=inputfile,
        col.main="gray50", cex.main=1, font.main=1, col.sub="gray50", cex.sub=0.8, font.sub=1)
}
  ### Return data.frame (invisibly)
  invisible(spec)
  }
