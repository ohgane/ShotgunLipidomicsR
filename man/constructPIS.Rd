% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constructPIS.R
\name{constructPIS}
\alias{constructPIS}
\title{Construct precursor ion scan data from multiple xcmsRaw objects}
\usage{
constructPIS(xrawlist, mz, mzabs = 0.5, ppm = NULL, digits = 1,
  centroid = TRUE, summarizeData = TRUE)
}
\arguments{
\item{xrawlist}{A list of xcmsRaw objects to be used for construction of PIS 
data.}

\item{mz}{A vector of m/z values of fragment ions.}

\item{mzabs}{Window size for peak area calculation. Default, 0.5.}

\item{ppm}{Window size for peak area calculation. If provided, this overrides 
\code{mzabs} setting.}

\item{digits}{Digits used for fragment names.}

\item{centroid}{Default, TRUE. For centroided data, this should be TRUE.}

\item{summarizeData}{Default, TRUE. If set to FALSE, this function returns 
PIS data without averaging.}
}
\value{
A data frame containing PIS (precursor ion scanning) data set, with 
   column names of mz, fragment, and int.
}
\description{
This function constructs a summrized data frame from a set of xcmsRaw objects.
Currently, this is designed for use with LTQ-Orbitrap XL dataset acquired with
"Parent ion mapping" method. With this method, no MS1 data was acquired and only 
MS2 data set was acquired over a range of m/z values. The data should be mzML
converted from msConvert.
}
\examples{
# Think about a situation where 1_PIS.mzML etc are present in a directory "mzML".
# input=list.files("../mzML", pattern="PIS", full.names=TRUE)
# xraw_list=lapply(input, xcmsRaw, includeMSn=TRUE, profstep=1)
# pisdata=constructPIS(xraw_list[1:3], mz=c(184.1, 208.1))
# xyplot(int~mz, groups=fragment, data=pisdata, type="l")
}

