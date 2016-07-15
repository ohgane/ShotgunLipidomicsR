#' Precursor m/z correction of xcmsRaw object imported from mzML (converted from Thermo RAW)
#'
#' This function is a kind of a bug fix for mzML (generated from ThermoRAW by
#'    msConvert) data acquired on LTQ-Orbitrap XL, with monoisotopic selection
#'    OFF under DDA settings. The setting results in the precursor mz entry in
#'    the mzML to be "monoisotopic mass" (even though "monoisotopic selection"
#'    is disabled). Therefore real precursor mz values are not read from the
#'    mzML file into xcmsRaw object. As the real precursor m/z values are still
#'    present in the mzML file (though only to 2 decimals), this function extract
#'    the real precursor m/z and pick peak near the m/z value.
#'
#' @param xraw An \code{xcmsRaw} object.
#' @param path A path to the original mzML file.
#' @param mzabs Window size to search real precursor m/z value (default 0.05).
#'
#' @return \code{xcmsRaw} object with corrected precursor m/z values
#'
#' @export


### A function that correct xcmsRaw's precursor m/z entries, using extractScanInfo().
correctPrecursorMz=function(xraw, path, mzabs=0.05){
  # Retrieve header information including precursor m/z (real precursor mz values) from mzML
  scaninfo=extractScanInfo(path)
  scaninfo=scaninfo[scaninfo$msLevel==2, ]
  precursor=scaninfo[, "precursorMz"]
  # Perform peak picking for MS1 data
  # As to apply matchedFilter peak picker for direct infusion data, fwhm is set to large value.
  peaks=findPeaks(xraw, method="matchedFilter", fwhm=10000, snthresh=1.5)

  # Perform precursor m/z correction, using peak list from MS1 data
  ND=0
  MD=0
  for (i in 1:nrow(scaninfo)){
    # Find peak nearest to the precursor m/z values extracted from mzML file (digits 2)
    peakext=peaks[fuzzyMatch(peaks[, "mz"], scaninfo[i, "precursorMz"], mzabs=mzabs), "mz"]
    if (length(peakext)==1) {
      precursor[i]=peakext
    } else if (length(peakext) > 1) {
      precursor[i]=median(peakext) # This may not be reasonable. Just selecting one value.
      MD=MD+1 # Count scans where median peak was used.
    } else if (length(peakext < 1)) {
      # If no peak is found in the indicated window, increase the window size 2-fold.
      peakext=peaks[fuzzyMatch(peaks[, "mz"], scaninfo[i, "precursorMz"], mzabs=2*mzabs), "mz"]
      if (length(peakext)==1) {
        precursor[i]=peakext
      } else if (length(peakext) > 1) {
        precursor[i]=median(peakext) # This may not be reasonable. Just selecting one value.
      } else {
      ND=ND+1 # Count scans where no peak can be assigned.
      }
    }
  }
  # Replace xcmsRaw@msnPrecursorMz with the retrieved values
  xraw@msnPrecursorMz=precursor
  # The number of MS/MS scans that could not be associated with MS1 peaks (peak picked values)
  #     is printed as a warning.
  if (ND > 0){
    warning(paste("Of the ", nrow(scaninfo), " MS/MS spectra, \n     ",
                  ND, " spectra could not be asssociated with precursor m/z values \n",
                  "     found in averaged MS1 spectrum.", sep=""))
  }
  if (ND > 0){
    warning(paste("Of the ", nrow(scaninfo), " MS/MS spectra, \n     ",
                  "for ", MD, " spectra, median of precursor m/z values found in averaged\n",
                  "    MS1 spectrum was associated with each scan.", sep=""))
  }
  return(xraw)
}


# Older version
#correctPrecursorMz=function(xraw, path){
#  # Retrieve "scan filter" that contains real precursor mz values
#  scanfilter=extractScanInfo(path)
#  precursor=scanfilter[scanfilter$msLevel==2, "precursorMz"]
#  # Replace xcmsRaw@msnPrecursorMz with the retrieved values
#  xraw@msnPrecursorMz=precursor
#  return(xraw)
#}



