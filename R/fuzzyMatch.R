#' A function that performs matching with a fixed tolerance, and return T/F
#'
#' This function accepts a query vector and a reference value to search,
#'     and return T/F vector of the same length as the query vector.
#'     For elements of query vector, it returns TRUE if the value is
#'     within the tolerance of reference value.
#'     The larger of the \code{mzabs} or \code{ppm} is used as tolerance.
#'
#' @param query A numeric vector.
#' @param ref A numeric.
#' @param mzabs Tolerance in m/z unit. The larger of the mzabs or ppm is used.
#' @param ppm Tolerance in ppm. The larger of the mzabs or ppm is used.
#'
#' @return A logical vector of the same length as \code{query} vector.
#'
#' @usage sapply(vector, fuzzyMatch, ref, mzabs=1)
#'
#' @export

# Older version
#fuzzyMatch=function(query, ref, mzabs=0, ppm=5){
#  ppmabs=mzabs/ref*10^6
#  ppm=max(c(ppm, ppmabs))
#  refrange=ref+c(-ref*ppm/10^6, ref*ppm/10^6)
#  res=refrange[1] <= query & query <= refrange[2]
#  return(res)
#}

# Test version: vectorized
#fuzzyMatchVec=function(query, ref, mzabs=0, ppm=5){
#  res=apply(
#    sapply(query, fuzzyMatch, query=ref, mzabs=mzabs, ppm=ppm),
#    MARGIN=2, sum)
#  return(as.logical(res))
#}

fuzzyMatch=function(query, ref, mzabs=0, ppm=5){
  # Define an internal fuzzy matching function
  f=function(query, ref, mzabs, ppm){
    ppmabs=mzabs/ref*10^6
    ppm=max(c(ppm, ppmabs))
    refrange=ref+c(-ref*ppm/10^6, ref*ppm/10^6)
    res=refrange[1] <= query & query <= refrange[2]
    return(res)
  }
  if (length(ref)==1){
    res=f(query=query, ref=ref, mzabs=mzabs, ppm=ppm)
  } else {
    res=apply(
      sapply(query, f, query=ref, mzabs=mzabs, ppm=ppm),
      MARGIN=2, sum)
    res=as.logical(res)
  }

  return(res)
}

