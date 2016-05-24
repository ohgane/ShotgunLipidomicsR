#' Convert \code{xcmsFragments} object into a list (MGF format)
#'
#' @param xfrag \code{xcmsFragments} object
#'
#' @return Returns a nested list, each of which contains \code{precursor} (numeric)
#'      and \code{spec} (data.frame, with \code{mz} and \code{intensity} columns).
#'
#' @export

convertXFrag2Mgf=function(xfrag){
  MS1=xfrag@peaks[xfrag@peaks[, "msLevel"]==1,]
  MS2=xfrag@peaks[xfrag@peaks[, "msLevel"]==2,]
  mgf=list(1:nrow(MS1))
  for (i in 1:nrow(MS1)){
    # Limit the number of peaks to 20 at maximum for MS/MS search
    peakNum=min(sum(MS2[, "MSnParentPeakID"]==i), 20)
    spec=MS2[MS2[, "MSnParentPeakID"]==i, c("mz", "intensity")]
    mgf[[i]]=list(precursor=MS1[i, "mz"],
                  spec=spec[1:peakNum, ])
  }
  return(mgf)
}
