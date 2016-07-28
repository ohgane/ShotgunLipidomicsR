#' A function that retrieve averaged spectra from a xcmsRaw object as a hyperSpec object
#'
#' From a xcmsRaw object (with multiple scans), this function returns a hyperSpec object (averaged spectrum). As this function only accepts one object, later merging would be necessary to construct a hyperSpec object with multiple spectra (use rbind()).
#' @param xcmsraw An object of class \code{xcmsRaw} class
#' @param scans A numeric vector (default NULL returns average of all the scans). Scan index within the ms level, not an actual scan number.
#' @param msLevel MS level (MS1 or MS2). Either 1 or 2. Default, 1 (MS1).
#' @param step Numeric (default 0.01). The step size of x axis to which
#'      the raw data is interpolated.
#' @param mzrange Range of m/z to be included in the output \code{hyperSpec} object.
#'      By default, the min and max values of m/z in the data are used.
#'
#' @return Returns an object of class \code{hyperSpec}.
#'
#' @export

xcmsAvgHyperSpec=function(xcmsraw, scans=NULL, msLevel=1, step=0.01, mzrange=NULL){
  if (msLevel==2) {
    xcmsraw=deepCopy(xcmsraw)
    # Copy metadata for ms2 to main slot (in xcmsRaw object)
    xcmsraw=msn2xcmsRaw(xcmsraw)
  }
  if (is.null(scans)){
    scans=1:length(xcmsraw@acquisitionNum)
  }
  dat=getSpec_modified(xcmsraw, scan=scans) # Currently, mzrange can not be used here.
  if (is.null(mzrange)){
    xrange=range(dat[, "mz"])
  } else {
    xrange=mzrange
  }

  dat_approx=approx(dat[,1], dat[,2], xout=seq(xrange[1], xrange[2], by=step))
  spc_dat=new("hyperSpec", wavelength=dat_approx$x, spc=dat_approx$y,
              labels=list(.wavelength="m/z", spc="Intensity (a.u.)"))
  return(spc_dat)
}
