#' A function that retrieve averaged spectra from xcmsRaw object as a hyperSpec object
#'
#' @param xcmsraw An object of class \code{xcmsRaw} class
#' @param scans A numeric vector
#' @param step Numeric (default 0.01). The step size of x axis to which
#'      the raw data is interpolated.
#' @param mzrange Range of m/z to be included in the output \code{hyperSpec} object.
#'      By default, the min and max values of m/z in the data are used.
#'
#' @return Returns an object of class \code{hyperSpec}.
#'
#' @export

xcmsAvgHyperSpec=function(xcmsraw, scans, step=0.01, mzrange=NULL){
  xraw=deepCopy(xcmsraw)
  xraw=msn2xcmsRaw(xraw)
  dat=getSpec_modified(xraw, scan=scans)
  if (is.null(mzrange)){
    xrange=range(dat_all[, "mz"])
  } else {
    xrange=mzrange
  }

  dat_approx=approx(dat[,1], dat[,2], xout=seq(xrange[1], xrange[2], by=step))
  spc_dat=new("hyperSpec", wavelength=dat_approx$x, spc=dat_approx$y,
              labels=list(.wavelength="m/z", spc="Intensity (a.u.)"))
  return(spc_dat)
}
