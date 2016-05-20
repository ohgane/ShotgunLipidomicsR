#' A function that returns scan indexes from a precursor m/z value and xcmsRaw object
#'
#' This function retrieves scan indexes ("real scan number" and scan number within MS2 scans)
#'      from precursor ion list (from meta-data).
#'
#' @param prmz A numeric value. Numeric vectors cannot be used
#'      (in that cases, only the first element will be used).
#' @param xraw An object of class \code{xcmsRaw}, which is read into by \code{xcmsRaw()}
#'      with \code{includeMSn = TRUE} option.
#' @param mzabs Search tolerance in m/z (not ppm).
#'
#' @return Returns a data frame, with columns \code{PrecursorMz}, \code{ScanNum},
#'      and \code{ScanNumMS2}.
#'
#' @export

findScans=function(prmz, xraw, mzabs=0.1){
  # depend on fuzzyMatch.r
  require(xcms)
  if (!is.numeric(prmz)){
    stop("The argment prmz must be a numeric.")
  }
  if (class(xraw) != "xcmsRaw") {
    stop("The argment xraw must be an object from xcmsRaw class.")
  }
  if (length(prmz)>1){
    warning("The argument prmz should be a numeric value. Only the first entry was used.")
    prmz=prmz[1]  # Use only the first element.
  }
  # Make index table (look up table for m/z and scan numbers)
  scanlist=data.frame(PrecursorMz=xraw@msnPrecursorMz,
                      ScanNum=xraw@msnAcquisitionNum,
                      ScanNumMS2=1:length(xraw@msnPrecursorMz))
  # Extract subset of the index table that fullfill the required matching criteria (tolerance).
  selected=scanlist[sapply(scanlist$PrecursorMz, fuzzyMatch, prmz, mzabs=mzabs), ]
  if (nrow(selected)==0){
    warning(paste("No MS/MS data was found for m/z", prmz, "within", mzabs, "tolerance."))
    return(NULL)
  } else {
    # Return as a data frame.
    return(selected)
  }
}
