#' Title
#'
#' @param sampleSDF
#'
#' @return
#' @export
#'
#' @examples
fmcslookupfctn <- function(sampleSDF) {
  MCS = batchfull[1] #Loads comparison library
  Tanimoto = 0
  Index = 0 #Sets up variables Tanimoto, and Index
  samplesmi = ChemmineR::sdf2smiles(sampleSDF)
  for(i in 1:300)
  {

    compared <- fmcsR::fmcs(batchfull[i], sampleSDF[1], fast = TRUE) #Compares sample molecule to each molecule in comparison collection
    s <- readr::parse_number(compared[4]) #Removes words from output to isolate comparison value

    if (s > Tanimoto) #Checks if comparison value is the best match

    {
      Tanimoto = s
      MCS = batchfull[i]
      Index = i #Continues checking against the rest of the comparison collection
    }
  }

  if (Tanimoto > 0.7){
    smiles <- ChemmineR::sdf2smiles(MCS[1]) #Stores values for best match in SMILES format
    ChemmineR::write.SMI(smiles, file = "smiles.smi") #Creates output SMILES file
    sigma <- dbGetQuery(lookupdb2, 'SELECT "sigma" FROM testlookup2 WHERE "ID" == :z', params = list(z = Index)) #Searches for sigma values of substructure match
    sigma.meta <- dbGetQuery(lookupdb2, 'SELECT "sigma.meta" FROM testlookup2 WHERE "ID" == :z', params = list(z = Index))
    sigma.para <- dbGetQuery(lookupdb2, 'SELECT "sigma.para" FROM testlookup2 WHERE "ID" == :z', params = list(z = Index))
    char <- c(as.character(samplesmi), as.character(smiles), as.character(Tanimoto), as.character(sigma), as.character(sigma.meta), as.character(sigma.para)) #Creates list of output values
    outdf <- data.frame(matrix(unlist(char), nrow=1, byrow=T)) #Converts list into usable output format
    colnames(outdf) = c("Original Fragment SMILES","MCSS Match SMILES", "Tanimoto Index", "Sigma Value", "Sigma Meta Value", "Sigma Para Value")
    print(outdf) #Shows output values
    #print(paste("MCSS Match SMILES: ", smiles))
    #print(paste("Tanimoto Index: ", Tanimoto))
    #print(paste("Sigma Value: ", sigma))
    #print(paste("Sigma Meta Value: ", sigma.meta))
    #print(paste("Sigma Para Value: ", sigma.para)) #Prints relevant outputs
    #ChemmineR::plot(MCS[1], regenCoords = TRUE, print = FALSE) #Creates a visualization of output structure
    outmcs <- fmcsR::fmcs(MCS, sampleSDF[1])
    fmcsR::plotMCS(outmcs)
    jsonlite::write_json(outdf, "output-json.json", dataframe = "columns")
    }
  else {
    print("No similar matches were found.")
  }
}
