
getpara <- function(SMILES){
  
  originalSMILES <- SMILES #Stores SMILES string input
  originalSDF <- ChemmineR::smiles2sdf(originalSMILES) #Converts user-input SMILES sample string to SDF format
  MCS = PARAfrags[1] #Loads comparison library
  Tanimoto = 0
  Index = 0 #Sets up variables Tanimoto Match Value and Index Number
  for(i in 1:300) ###Should this be "nrow" instead of 300?
  {
    
    compared <- fmcsR::fmcs(PARAfrags[i], originalSDF[1], fast = TRUE) #Compares sample molecule to each molecule in comparison collection
    s <- readr::parse_number(compared[4]) #Removes words from output to isolate comparison value
    
    if (s > Tanimoto) #Checks if comparison value is the best match
      
    {
      Tanimoto = s
      MCS = PARAfrags[i]
      Index = i #Continues checking against the rest of the comparison collection
    }
  }
  
  matchSMILES <- ChemmineR::sdf2smiles(MCS[1]) #Stores the best match molecule in SMILES format
  sigma.para <- RSQLite::dbGetQuery(lookupparaDB, 'SELECT "sigma.para" FROM lookuppara WHERE "ID" == :z', params = list(z = Index))
  char <- c(as.character(originalSMILES), as.character(matchSMILES), as.character(Tanimoto), as.character(sigma.para)) #Creates list of output values
  output <- data.frame(matrix(unlist(char), nrow = 1, byrow = T)) #Converts list into usable output format as a dataframe
  colnames(output) = c("OriginalFragmentSMILES","MCSSMatchSMILES", "TanimotoMatchIndex", "SigmaParaValue") #Creates column labels for the dataframe
  return(output) #What format should the output be in for use in Python dataframe?
  
}