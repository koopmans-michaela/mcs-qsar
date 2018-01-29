fmcsforCTS<- function(fragSDF) {
  MCS = batchfull[1] #Loads comparison library
  Tanimoto = 0
  Index = 0 #Sets up variables Tanimoto Match Value and Index Number
  samplesmi = ChemmineR::sdf2smiles(fragSDF) #Converts sample fragment from SDF to SMILES to be used in the output dataframe
  for(i in 1:300)
  {
    
    compared <- fmcsR::fmcs(batchfull[i], fragSDF[1], fast = TRUE) #Compares sample molecule to each molecule in comparison collection
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
    smilesCHAR <- as.character(smiles) #Converts match molecule to character string for editing
    names(smilesCHAR) = "MCSS Match Structure" #Changes name of match molecule character vector to appropriate label for use in plotMCS
    matchSDF <- ChemmineR::smiles2sdf(smilesCHAR) #Converts renamed match molecule vector back to SDF for use in plotMCS
    fragSMI <- ChemmineR::sdf2smiles(fragSDF[1]) #Converts original fragment to SMIset format
    fragCHAR <- as.character(fragSMI) #Converts original fragment to character string for editing
    names(fragCHAR) = "Original Fragment Structure" #Changes name of original fragment character vector to appropriate label for use in plotMCS
    fragSDF2 <- ChemmineR::smiles2sdf(fragCHAR) #Converts renamed original fragment vector back to SDF for use in plotMCS
    outmcs <- fmcsR::fmcs(fragSDF2, matchSDF) #Runs MCS between match molecule and original fragment to get the output information in MCS format for use in the plotMCS visualization function
    fmcsR::plotMCS(outmcs) #Visualizes the original fragment and match molecules and highlights the similar substructure
    ChemmineR::write.SMI(smiles, file = "smiles.smi") #Creates output SMILES file for the match molecule
    sigma <- RSQLite::dbGetQuery(lookupdb2, 'SELECT "sigma" FROM testlookup2 WHERE "ID" == :z', params = list(z = Index)) #Retrieves sigma values of match molecule
    sigma.meta <- RSQLite::dbGetQuery(lookupdb2, 'SELECT "sigma.meta" FROM testlookup2 WHERE "ID" == :z', params = list(z = Index))
    sigma.para <- RSQLite::dbGetQuery(lookupdb2, 'SELECT "sigma.para" FROM testlookup2 WHERE "ID" == :z', params = list(z = Index))
    char <- c(as.character(samplesmi), as.character(smiles), as.character(Tanimoto), as.character(sigma), as.character(sigma.meta), as.character(sigma.para)) #Creates list of output values
    outdf <- data.frame(matrix(unlist(char), nrow = 1, byrow = T), stringsAsFactors = FALSE) #Converts list into usable output format as a dataframe
    colnames(outdf) = c("Original Fragment SMILES","MCSS Match SMILES", "Tanimoto Match Index", "Sigma Value", "Sigma Meta Value", "Sigma Para Value") #Creates column labels for the dataframe
    #print(outdf) #Shows output values in console
    #jsonlite::write_json(outdf, "output-json.json", dataframe = "columns") #Writes a JSON file containing the output values dataframe for use in CTS
    return(outdf)
  }
  else {
    Error = "Error"
    NoMatch = "No similar matches were found. The closest structure available is given." #Creates error message
    closest <- ChemmineR::sdf2smiles(MCS[1]) #Stores the closest MCSS match SMILES
    OutNoMatch <- c(Error, NoMatch, as.character(samplesmi), as.character(closest), as.character(Tanimoto)) #Starts creating dataframe to put error message into output JSON
    nomatchDF <- data.frame(matrix(unlist(OutNoMatch), nrow=1, byrow = T)) #Creates dataframe
    colnames(nomatchDF) = c("Error", "No Match Found", "Original Fragment SMILES", "Closest Match SMILES", "Tanimoto Match Index") #Gives dataframe column names
    #jsonlite::write_json(nomatchDF, "output-nomatch.json", dataframe = "columns") #Outputs error message in JSON format
    return(nomatchDF)
  }
}
