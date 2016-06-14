#' Construct precursor ion scan data from multiple xcmsRaw objects
#' 
#' This function constructs a summrized data frame from a set of xcmsRaw objects.
#' Currently, this is designed for use with LTQ-Orbitrap XL dataset acquired with
#' "Parent ion mapping" method. With this method, no MS1 data was acquired and only 
#' MS2 data set was acquired over a range of m/z values. The data should be mzML
#' converted from msConvert.
#' 
#' @param xrawlist A list of xcmsRaw objects to be used for construction of PIS 
#'    data.
#' @param mz A vector of m/z values of fragment ions.
#' @param mzabs Window size for peak area calculation. Default, 0.5.
#' @param ppm Window size for peak area calculation. If provided, this overrides 
#'    \code{mzabs} setting.
#' @param digits Digits used for fragment names.
#' @param centroid Default, TRUE. For centroided data, this should be TRUE.
#' @param summarizeData Default, TRUE. If set to FALSE, this function returns 
#'    PIS data without averaging.
#' 
#' @return A data frame containing PIS (precursor ion scanning) data set, with 
#'    column names of mz, fragment, and int.
#' 
#' @export
#' 
#' @examples 
#' # Think about a situation where 1_PIS.mzML etc are present in a directory "mzML".
#' # input=list.files("../mzML", pattern="PIS", full.names=TRUE)
#' # xraw_list=lapply(input, xcmsRaw, includeMSn=TRUE, profstep=1)
#' # pisdata=constructPIS(xraw_list[1:3], mz=c(184.1, 208.1))
#' # xyplot(int~mz, groups=fragment, data=pisdata, type="l")

constructPIS=function(xrawlist, mz, mzabs=0.5, ppm=NULL, digits=1,
                      centroid=TRUE, summarizeData=TRUE){
  require(reshape2)
  require(doBy)
  # If ppm was provided, calculate mzabs
  if (!is.null(ppm)){
    # Convert ppm to mzabs (vector)
    mzabs=mz*ppm/10^6
  } else {
    # mzabs as a vector
    mzabs=rep(mzabs, length(mz))
  }
  # Profile or centroid: used in peak area calculation
  if (centroid){
    datatype="centroid"
  } else {
    datatype="profile"
  }
  #
  pisdata=NULL
  for (i in 1:length(xrawlist)){
    temp=data.frame(mz=xrawlist[[i]]@msnPrecursorMz)
    for (k in 1:length(mz)){
      temp=cbind(temp,
                 sapply(1:length(xrawlist[[i]]@msnAcquisitionNum),
                        xcmsrawAuc, mzrange=c(mz[k]-mzabs[k]/2, mz[k]+mzabs[k]/2),
                        xraw=xrawlist[[i]], type=datatype, msLevel=2)
      ) # End of cbind
    } # End loop over mz
    colnames(temp)=c("mz", paste("", round(mz, digits=digits), sep=""))
    rownames(temp)=NULL
    pisdata=rbind(pisdata, temp)
  } # End loop over xrawlist
  pisdata=pisdata[order(pisdata$mz),]
  if (summarizeData){
    pislong=melt(pisdata, id.var="mz", variable.name="fragment", value.name="int")
    pis=summaryBy(int~mz+fragment, data=pislong, keep.names=TRUE,
                  FUN=function(x) sum(x)/length(xrawlist))
  } else {
    pis=melt(pisdata, id.var="mz", variable.name="fragment", value.name="int")
  }
  return(pis)
} # End of function