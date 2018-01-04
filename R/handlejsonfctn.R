#' Title
#'
#' @param JSON
#'
#' @return
#' @export
#'
#' @examples
#'
handleJSON <- function(samplejson) {
  samplesmiles <- samplejson[["data"]][["format"]][["smiles"]] #Parses JSON to get and store SMILES string
  sampleSDF <- ChemmineOB::convertFormat("SMI", "SDF", samplesmiles) #Converts user-input SMILES sample string to SDF format
  print(paste("Original Fragment SMILES: ", samplesmiles))
  fmcslookupfctn(sampleSDF)
}
