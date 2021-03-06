% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotAvgMS1.R
\name{plotAvgMS1}
\alias{plotAvgMS1}
\title{A function to plot averaged MS1 spectrum from a \code{xcmsRaw} object}
\usage{
plotAvgMS1(xcmsraw, normalization = "constant", int.divide = 1000,
  col = 1, type = "l", plot.precursor = TRUE, custom.precursor = NULL,
  col.precursor = "light blue", lwd.precursor = 2, lty.precursor = 1,
  baseline = FALSE, col.baseline = "gray50", lty.baseline = 1,
  scans = NULL, mzrange = NULL, ylim = NULL, ad.lab.mzrange = NA,
  ad.top = 3, top = 12, halfWindowSize = NULL, SNR = 50, digits = 4,
  ...)
}
\arguments{
\item{xcmsraw}{\code{xcmsRaw} object.}

\item{normalization}{Type of normalization for plotting. Default, "constant".
Either "constant", "relative", or "TIC" (percent of total ion current).}

\item{int.divide}{Intensity normalization factor (simple division). Default, 1000.}

\item{col}{Default, 1 (=black).}

\item{type}{Default, "l" (line).}

\item{plot.precursor}{Logical (default, TRUE). If true, plot lines at the precursor m/z values.}

\item{custom.precursor}{Default, NULL.}

\item{col.precursor}{Default, "light blue".}

\item{lwd.precursor}{Default, 2.}

\item{lty.precursor}{Default, 1.}

\item{baseline}{Logical (default TRUE)}

\item{col.baseline}{The color of baseline. Default, "gray50"}

\item{lty.baseline}{Line type of baseine. Default, 1.}

\item{scans}{A numeric vector of scan indexes. If provided, only the indicated scans are averaged and plotted.}

\item{mzrange}{m/z range to plot.}

\item{ylim}{Y axis range to plot.}

\item{ad.top}{The number of additional peak labels.}

\item{top}{The number of peak labels.}

\item{halfWindowSize}{See \code{MALDIquant::labelPeaks()} for detail.}

\item{SNR}{See \code{MALDIquant::labelPeaks()} for detail.}

\item{digits}{Digits for peak labels.}

\item{ad.lab.}{mzrange m/z range to add additional peak labels.}
}
\value{
Plot spectrrum. Plus, if assigned to a variable (= invisibly),
   returns the averaged spectrum as a data frame.
}
\description{
This function accepts \code{xcmsRaw} object with MS/MS data (direct infusion,
   data-dependent acquisition) and plot the spectrum and overlay colored lines over the precursors. Additionally, label some selected peaks.
}

