#' Modification to xcms::getSpec()
#'
#' A function that extract scan(s) from a \code{xcmsRaw} object.
#'    A bug that prevented scan averaging of MS2 data was removed.
#'    Incontrast to the original getSpec, scanrange must be provided (default NULL).
#'
#' @param prmz A numeric value. Numeric vectors cannot be used
#'      (in that cases, only the first element will be used).
#' @param xraw An object of class \code{xcmsRaw}, which is read into by \code{xcmsRaw()}
#'      with \code{includeMSn = TRUE} option.
#' @param mzabs Search tolerance in m/z (not ppm).
#'
#' @return A data frame with columns of the names \code{mz} and \code{intensity}.
#'
#' @export
getSpec_modified=function (object, scanrange=NULL, ...)
{
  sel <- profRange(object, ...)
  if (is.null(scanrange)) scanrange=sel$scanidx
  scans <- list(length(scanrange))
  uniquemz <- numeric()
  for (i in seq(along = scanrange)) {
    scans[[i]] <- getScan(object, scan=scanrange[i]) #, mzrange=sel$mzrange
    uniquemz <- unique(c(uniquemz, scans[[i]][, "mz"]))
  }
  uniquemz <- sort(uniquemz)
  intmat <- matrix(nrow = length(uniquemz), ncol = length(scanrange))
  for (i in seq(along = scanrange)) {
    scan <- getScan(object, scan=scanrange[i]) #, mzrange=sel$mzrange
    intmat[, i] <- approx(scan, xout = uniquemz)$y
  }
  # na.rm=TRUE option in the following rowMeans() is important.
  points <- cbind(mz = uniquemz, intensity = rowMeans(intmat, na.rm=TRUE))
  invisible(points)
}

#> getMethod("getSpec", "xcmsRaw")
#Method Definition:
#
#  function (object, ...)
#  {
#     sel <- profRange(object, ...)
#     scans <- list(length(sel$scanidx))
#     uniquemz <- numeric()
#     for (i in seq(along = sel$scanidx)) {
#       scans[[i]] <- getScan(object, sel$scanidx[i], sel$mzrange)
#       uniquemz <- unique(c(uniquemz, scans[[i]][, "mz"]))
#     }
#     uniquemz <- sort(uniquemz)
#     intmat <- matrix(nrow = length(uniquemz), ncol = length(sel$scanidx))
#     for (i in seq(along = sel$scanidx)) {
#       scan <- getScan(object, sel$scanidx[i], sel$mzrange)
#       intmat[, i] <- approx(scan, xout = uniquemz)$y
#     }
#     points <- cbind(mz = uniquemz, intensity = rowMeans(intmat))
#     invisible(points)
#   }
# <environment: namespace:xcms>
#
#   Signatures:
#   object
# target  "xcmsRaw"
# defined "xcmsRaw"
