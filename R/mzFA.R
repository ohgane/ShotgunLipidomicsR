#' Calculate exact mass for fatty acid
#'
#' This function calculate exact mass for fatty acid [n:m].
#'
#' @param n The number of carbon. A vector of length 2 (n,m) can also be used.
#'    In that case, \code{m} must be the default value \code{NULL}.
#' @param m The number of double bond. Default, \code{NULL}.
#' @param digits Decimal point (default 4)
#'
#' @export
#'
#' @examples
#' # Exact mass of FA[18:1] anion (281.2486).
#' mzFA(18, 1)

mzFA=function(n,m=NULL, digits=4){
  require(Rdisop)
  electron=0.000549
  if (length(n)==2 & is.null(m)){
    m=n[2]
    n=n[1]
  }
  fa=getMolecule(paste("C", n, "H", 2*n-1-2*m, "O2", sep=""), z=-1)
  mz=round(getMass(fa)+electron, digits=digits)
  return(mz)
}
