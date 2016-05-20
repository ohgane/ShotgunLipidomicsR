#' Calculate area under the curve from \code{xcmsRaw} object
#'
#' @param mz Center m/z values
#' @param mzabs Width to integrate
#' @param xraw An object of class \code{xcmsRaw} class
#' @param scan An Numeric vector. Scan index within the indicated MS level..
#'      Note that this scan index is NOT the actual scan number (include both MS1 and MS2).
#'      This corresoinds to the index of the actual scan number in
#'      \code{xcmsRaw@msnAcquisitionNum} slot.
#' @param type "profile" or "centroid". Default, "profile".
#' @param msLevel MS level (default 1).
#'
#' @return
#'
#' @export

xcmsrawAuc=function(scan, mzrange=NULL, mz=NULL, mzabs=NULL,
                    xraw, type="profile", msLevel=1){
  require(xcms)
  ## integration range can be specified either mzrange or combination of mz and mzabs
  if (is.null(mz) & is.null(mzabs) & !is.null(mzrange)){
    mzrange=mzrange
  } else if (is.null(mzrange) & !is.null(mz) & !is.null(mzabs)) {
    mzrange=c(mz-0.5*mzabs, mz+0.5*mzabs)
  } else {
    stop("Range for integration must be specified either mzrange
         or combination of mz and mzabs.")
  }
  if (type=="profile"){
    if (msLevel==1){
      dat=getScan(xraw, scan=scan)
    } else if (msLevel == 2){
      dat=data.frame(getMsnScan(xraw, scan=scan)) # xcms::getMsnScan()
      # Note that the data is converted to data.frame,
      #   to avoid possible unwanted conversion to vector
      #   when nrow=1 ("dat[,1]" can cause error).
    } else {
      warning("The msLevel argment should be 1 or 2.")
    }
    # Extract region of interest
    dat=dat[mzrange[1]<dat$mz & dat$mz < mzrange[2], ]
    if (nrow(dat) != 0){
      auc=trapz(dat$mz, dat$intensity)
    } else {
      auc=0
    }
  } else if (type=="centroid") {
    # For centroided data, area under the curve can not be calculated.
    # Instead, sum of peak intensity within a range is used.
    if (msLevel==1){
      dat=getScan(xraw, scan=scan)
    } else if (msLevel == 2){
      dat=data.frame(getMsnScan(xraw, scan=scan)) # xcms::getMsnScan()
    } else {
      warning("The msLevel argment should be 1 or 2.")
    }
    dat=dat[mzrange[1]<dat$mz & dat$mz < mzrange[2], ]
    if (nrow(dat) != 0){
      auc=sum(dat$intensity)
    } else {
      auc=0
    }
  } else {
    stop("type must be either profile or centroid.")
  }
  return(auc)
}

# Older version: flux::auc() for integration, slow.
#xcmsrawAuc=function(mzrange, xraw, scan){
#  require(flux) # flux::auc()
#  require(xcms)
#  if (msLevel==1){
#    dat=getScan(xraw, scan=scan)
#  } else if (msLevel == 2){
#    dat=getMsnScan(xraw, scan=scan) # xcms::getMsnScan()
#  } else {
#    warning("The msLevel argment should be 1 or 2.")
#  }
#  dat=dat[mzrange[1]<dat[,1] & dat[,1] < mzrange[2], ]
#  if (nrow(dat) != 0){
#    auc=auc(dat[,1], dat[,2]) # flux::auc()
#  } else {
#    auc=0
#  }
#  return(auc)
#}

