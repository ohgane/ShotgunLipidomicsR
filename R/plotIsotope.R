#' Plot isotope patterns
#' 
#' A wrapper function for Rdisop::getMolecule, getIsotopes etc., that allows easier plotting of the isotope patterns.
#' 
#' @param chemsc Chemical sum composition (elemental formula) as a character (e.g. "H2O").
#' @param z Charge
#' @param add Logical (default FALSE). If TRUE, overlay isotope pattern on an existing plot.
#' @param type Type argument for plotting (used within points() function). Default, "p" (point).
#' @param ymax Default=1. The height of the maximum peak.
#' @param max.iso Maximum number of isotopes to plot. Default, 5.
#' 
#' @export
#' 
plotIsotope=function(chemsc, z=0, add=FALSE, type="p", ymax=1, max.iso=5, ...){
  require(Rdisop)
  electron=0.000549
  mol=getMolecule(chemsc, z=z)
  # getIsotope returns a matrix with its first row mz and the second isotopic abundance
  iso=getIsotope(mol)[, 1:max.iso]
  if (!add){
    plot(iso[1,]+z*electron, ymax*iso[2,]/max(iso[2,]),
         type=type, xlab="m/z", ylab="Abundance", ...)
    title(paste(round(getMass(mol)+z*electron, digits=4), getFormula(mol), sep=" "))
  } else {
    points(iso[1,]+z*electron, ymax*iso[2,]/max(iso[2,]),
           type=type, ...)
  }
}