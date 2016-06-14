#' Species ratio correction for LipidXplorer identification result (negative mode)
#'
#' This function performs intensity correction based on intensities of precursor ion
#'      intensity and fatty acid anions. In a situation where there are several lipid
#'      species of the same m/z value but with different fatty acid compositions, this
#'      function allocates precursor ion intensity based on fatty acid fragment ion
#'      intensities.
#'
#' @param df A data frame. Result of LipidXplorer identification, without unnecessary lines.
#'      E.g., MFQL title lines that starts with ### should be removed.
#' @param anID Numeric vector. Column indexes for annotation part of input data frame.
#' @param precID Numeric vector. Column indexes for precursor ion intensities.
#' @param faID Numeric vector. Column indexes for fatty acid anion intensities.
#'
#' @return A data frame with annotation part, corrected precursor intensities,
#'      and original fatty acid anion intensities.
#'
#' @references Surma MA, Herzog R, Vasilj A, Klose C, Christinat N,
#'      Morin Rivron D, Simons K, Masoodi M & Sampaio JL (2015) An automated shotgun
#'      lipidomics platform for high throughput, comprehensive, and quantitative
#'      analysis of blood plasma intact lipids. European Journal of Lipid Science and
#'      Technology 117: 1540–1549
#'
#' @export
#'
correctLXneg=function(df, anID=NULL, precID=NULL, faID=NULL){
  lx_a=df[, anID] # Retrieve annotation columns
  lx_p=df[, precID] # Retrieve precursor intensities
  lx_f=df[, faID] # Retrieve fatty acid anion intensities
  lx_f[lx_f==" None"]=NA
  # Remove negative values
  lx_p[lx_p<0]=0
  lx_f[lx_f<0]=0
  # Convert to matrix (of numeric)
  #lx_p=matrix(as.numeric(as.matrix(lx_p)), nrow=dim(lx_p)[1], ncol=dim(lx_p)[2], byrow=T)
  #lx_f=matrix(as.numeric(as.matrix(lx_f)), nrow=dim(lx_f)[1], ncol=dim(lx_f)[2], byrow=T)
  # Perform correction, in a m/z-wise manner (correction within the same peak)
  lx_corrected=df[,c(anID, precID, faID)]
  masslist=unique(lx_a$MASS)
  for (i in seq_along(masslist)){
    index=which(lx_a$MASS==masslist[i])
    temp_a=lx_a[index,]
    temp_p=lx_p[index,]
    temp_f=lx_f[index,]
    # Define a helper function
    ratio=function(vec) {
      vec=as.numeric(vec)
      if (sum(is.na(vec))>0){ # If fatty acid data was missing, do nothing
        ratiovec=rep(1, length(vec))
      } else if (sum(vec)>0){
        ratiovec=vec/sum(vec)
      } else {
        ratiovec=rep(0, length(vec)) # return 0 vector of the same length
      }
      return(ratiovec)
    }
    if (is.null(dim(temp_f))){
      # In case there is only 1 entry, calculate ratio as a vector
      temp_fr=t(ratio(temp_f))
    } else {
      # Calculate ratio in a column-wise manner
      temp_fr=apply(temp_f, MARGIN=2, FUN=ratio)
    }

    # Convert vector to matrix
    if (is.null(dim(temp_fr))){
      temp_fr=t(temp_fr)
    }
    temp_p=temp_p*temp_fr
    lx_corrected[index,]=data.frame(cbind(temp_a, temp_p, temp_f))
  }
  return(lx_corrected)
}
