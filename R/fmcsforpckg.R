fmcsforpckg <- function(sampleSDF) {
  MCS = batchfull[1] #Loads comparison library
  Tanimoto = 0
  Index = 0 #Sets up variables Tanimoto, and Index
  samplesmi = ChemmineR::sdf2smiles(sampleSDF) #Converts sample fragment from SDF to SMILES to be used in the output dataframe
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
    smiles <- ChemmineR::sdf2smiles(MCS[1]) #Stores the best match molecule in SMILES format
    ChemmineR::write.SMI(smiles, file = "smiles.smi") #Creates output SMILES file for the match molecule
    sigma <- RSQLite::dbGetQuery(lookupdb2, 'SELECT "sigma" FROM testlookup2 WHERE "ID" == :z', params = list(z = Index)) #Retrieves sigma values of match molecule
    sigma.meta <- RSQLite::dbGetQuery(lookupdb2, 'SELECT "sigma.meta" FROM testlookup2 WHERE "ID" == :z', params = list(z = Index))
    sigma.para <- RSQLite::dbGetQuery(lookupdb2, 'SELECT "sigma.para" FROM testlookup2 WHERE "ID" == :z', params = list(z = Index))
    char <- c(as.character(samplesmi), as.character(smiles), as.character(Tanimoto), as.character(sigma), as.character(sigma.meta), as.character(sigma.para)) #Creates list of output values
    outdf <- data.frame(matrix(unlist(char), nrow=1, byrow=T)) #Converts list into usable output format as a dataframe
    colnames(outdf) = c("Original Fragment SMILES","MCSS Match SMILES", "Tanimoto Index", "Sigma Value", "Sigma Meta Value", "Sigma Para Value") #Creates column labels for the dataframe
    print(paste("Original Fragment SMILES: ", samplesmi))
    print(paste("MCSS Match SMILES: ", smiles))
    print(paste("Tanimoto Index: ", Tanimoto))
    print(paste("Sigma Value: ", sigma))
    print(paste("Sigma Meta Value: ", sigma.meta))
    print(paste("Sigma Para Value: ", sigma.para)) #Prints relevant outputs
    outmcs <- fmcsR::fmcs(MCS, sampleSDF[1]) #Runs MCS between match molecule and original fragment to get the output information in MCS format for use in the plotMCS visualization function
    fmcsR::plotMCS(outmcs) #Visualizes the original fragment and match molecules and highlights the similar substructure
    write.csv(outdf, "fmcs-output.csv") #Writes outdf to an output file
    }
  else {
    print("No similar matches were found.")
    print(paste("Original Fragment SMILES: ", samplesmi))
  }
}
