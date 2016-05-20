#' A function that extracts scan header information from a mzML file
#' 
#' This function extracts informations from a mzML file exported from 
#' ThermoRAW file with msConvert software. The information includes 
#' scan number, MS level, scanfilter message, total ion current, base peak, 
#' injection time etc. This function can be used to check injection time with 
#' the intention of optimizing the acquisition conditions.
#'
#' @param path The path to the mzML file.
#'
#' @return output A data frame 

#' @export



extractScanInfo=function(path){
  require(XML)
  doc=xmlTreeParse(path, useInternal=TRUE)
  top=xmlRoot(doc)
  specList=top[["mzML"]][["run"]][["spectrumList"]]

  # Helper function to retrieve scan number, msLevel, and scan filter text
  f=function(x){ # argument x each spectrum element, e.g. specList[[1]] etc
    scan=as.character(strsplit(xmlGetAttr(x, name="id"), split="scan=", fixed=TRUE)[[1]][2])
    msLevel=as.character(xmlGetAttr(x[[2]], name="value"))
    scanfilter=xmlGetAttr(x[["scanList"]][["scan"]][[2]], name="value")
    TIC=as.character(xmlGetAttr(x[[7]], name="value"))
    basePeakMz=as.character(xmlGetAttr(x[[5]], name="value"))
    basePeakInt=as.character(xmlGetAttr(x[[6]], name="value"))
    injectionTime=as.character(xmlGetAttr(x[["scanList"]][["scan"]][[4]], name="value"))
    scanStartTime=as.character(xmlGetAttr(x[["scanList"]][["scan"]][[1]], name="value"))
    if (msLevel != 1){
      precursorMz=xmlGetAttr(x[["precursorList"
                                ]][["precursor"]][["isolationWindow"]][[1]], name="value")
      precIsoWin=xmlGetAttr(x[["precursorList"
                               ]][["precursor"]][["isolationWindow"]][[2]], name="value")
      selectedMz=xmlGetAttr(x[["precursorList"
                               ]][["precursor"]][["selectedIonList"]][["selectedIon"]][[1]],
                            name="value")
      charge=xmlGetAttr(x[["precursorList"
                           ]][["precursor"]][["selectedIonList"]][["selectedIon"]][[2]],
                        name="value")
      activation=xmlGetAttr(x[["precursorList"]][["precursor"]][["activation"]][[1]],
                            name="name")
      nCE=xmlGetAttr(x[["precursorList"]][["precursor"]][["activation"]][[2]], name="value")
    } else {
      precursorMz=NA
      precIsoWin=NA
      selectedMz=NA
      charge=NA
      activation=NA
      nCE=NA
    }
    return(c(scan, msLevel, precursorMz, precIsoWin, selectedMz, charge, activation, nCE,
             TIC, basePeakMz, basePeakInt, injectionTime, scanStartTime, scanfilter))
  }

  temp=t(xmlSApply(specList, f))
  scanfilter=data.frame(scan=as.numeric(temp[,1]),
                        msLevel=as.numeric(temp[,2]),
                        precursorMz=as.numeric(temp[,3]),
                        isolationWindow=round(as.numeric(temp[,4])*2, digits=3),
                        selectedMz=as.numeric(temp[,5]),
                        charge=as.numeric(temp[,6]),
                        activation=temp[,7],
                        nCE=as.numeric(temp[,8]),
                        TIC=as.numeric(temp[,9]),
                        basePeakMz=as.numeric(temp[,10]),
                        basePeakInt=as.numeric(temp[,11]),
                        injectionTime=as.numeric(temp[,12]),
                        scanStartTime=as.numeric(temp[,13]),
                        scanfilter=temp[,14],
                        stringsAsFactors=FALSE)
  return(scanfilter)
}
