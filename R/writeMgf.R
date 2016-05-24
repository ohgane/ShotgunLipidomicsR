#' A function that export mgf file from filename, precursor mz, and a peak table.
#'
#' @param file A file name for MGF file to be exported.
#' @param precursor A precursor m/z value.
#' @param peak.table A data.frame with 2 column (mz, int).
#' @param charge Charge of the precursor. Character (default "+1").
#'
#' @export
writeMgf=function(file, precursor, peak.table, charge="+1"){
  out.file <- file(file, open = "w")
  writeLines("BEGIN IONS", out.file)
  writeLines(paste("PEPMASS=", precursor, sep=""), out.file)
  writeLines(paste("CHARGE=", charge, sep=""), out.file)
  for (i in 1:nrow(peak.table)) {
    writeLines(paste(peak.table[i,1], peak.table[i,2], sep=" "), out.file)
  }
  writeLines("END IONS", out.file)
  close(out.file)
}
