#' Title
#'
#' @param string
#'
#' @return
#' @export
#'
#' @examples
handleSMIstring <- function(string) {
  samplesmiles <- string #Stores SMILES string input
  sampleSDF <- ChemmineR::smiles2sdf(samplesmiles) #Converts user-input SMILES sample string to SDF format
  print(paste("Original Fragment SMILES: ", samplesmiles))
  fmcslookupfctn(sampleSDF)
}
