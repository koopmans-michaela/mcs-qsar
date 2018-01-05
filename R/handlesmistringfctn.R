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
  print(paste("Original Fragment SMILES: ", samplesmiles)) #Prints to console the SMILES string for the original sample fragment
  fmcslookupfctn(sampleSDF) #Performs MCS matching on sample fragment
}
