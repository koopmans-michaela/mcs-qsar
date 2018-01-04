#' Title
#'
#' @param SMIfile
#'
#' @return
#' @export
#'
#' @examples
handleSMIset <- function(SMIfile) {
  sampleSMI <- ChemmineR::read.SMIset(SMIfile) #Imports sample molecule in SMILES format
  sampleSDF <- ChemmineOB::convertFormat("SMI", "SDF", sampleSMI) #Converts user-input SMILES sample to SDF format
  fmcslookupfctn(sampleSDF)
}
