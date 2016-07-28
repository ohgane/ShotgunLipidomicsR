#' A function that retrieve TIC from a list of xcmsRaw objects as a hyperSpec object
#'
#' From a list of xcmsRaw objects (LC-MS run), this function returns a hyperSpec object (TIC).
#' @param xraw_list A list of objects of class \code{xcmsRaw} class
#' @param step Numeric (default 1 sec). The step size of retention time axis to which
#'      the raw data is interpolated.
#' @param rtrange Range of m/z to be included in the output \code{hyperSpec} object.
#'      By default, the min and max values of m/z in the data are used.
#'
#' @return Returns an object of class \code{hyperSpec}.
#'
#' @export

xcmsTICHyperSpec=function(xraw_list, rtstep=1, rtrange=NULL, mzrange=NULL, mzstep=1){
  require(foreach)
  if (is.null(rtrange)){
    xrange=range(xraw_list[[1]]@scantime)
  } else {
    xrange=rtrange
  }
  ## If mzrange is not provided, return TIC, and if provided return EIC
  if (is.null(mzrange)){ # Total ion chromatogram from scan info
    spc_list=foreach (i = 1:length(xraw_list), .combine="rbind") %do% {
      dat=data.frame(scantime=xraw_list[[i]]@scantime, tic=xraw_list[[i]]@tic)
      dat_approx=approx(dat$scantime, dat$tic, xout=seq(xrange[1], xrange[2], by=rtstep))
      spc=new("hyperSpec", wavelength=dat_approx$x, spc=dat_approx$y,
              labels=list(.wavelength="Retention time (sec)", spc="Intensity (a.u.)"))
      return(spc)
    }
  } else { # Extracted ion chromatogram
    spc_list=foreach (i = 1:length(xraw_list), .combine="rbind") %do% {
      eic=getEIC(xraw_list[[i]], 
                 mzrange=matrix(mzrange, nrow=1, ncol=2), 
                 rtrange=matrix(xrange, nrow=1, ncol=2),
                 step=mzstep)
      dat=data.frame(scantime=eic@eic[[1]][[1]][,"rt"], eic=eic@eic[[1]][[1]][, "intensity"])
      dat_approx=approx(dat$scantime, dat$eic, xout=seq(xrange[1], xrange[2], by=rtstep))
      spc=new("hyperSpec", wavelength=dat_approx$x, spc=dat_approx$y,
              labels=list(.wavelength="Retention time (sec)", spc="Intensity (a.u.)"))
      return(spc)
    }
  }
  return(spc_list)
}
