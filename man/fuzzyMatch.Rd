% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fuzzyMatch.R
\name{fuzzyMatch}
\alias{fuzzyMatch}
\title{A function that performs matching with a fixed tolerance, and return T/F}
\usage{
sapply(vector, fuzzyMatch, ref, mzabs=1)
}
\arguments{
\item{query}{A numeric vector.}

\item{ref}{A numeric.}

\item{mzabs}{Tolerance in m/z unit. The larger of the mzabs or ppm is used.}

\item{ppm}{Tolerance in ppm. The larger of the mzabs or ppm is used.}
}
\value{
A logical vector of the same length as \code{query} vector.
}
\description{
This function accepts a query vector and a reference value to search,
    and return T/F vector of the same length as the query vector.
    For elements of query vector, it returns TRUE if the value is
    within the tolerance of reference value.
    The larger of the \code{mzabs} or \code{ppm} is used as tolerance.
}
\examples{
fuzzyMatch(1:5, ref=3, ppm=10)
}

