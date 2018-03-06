

dout = NULL #Create empty dataframe

for(t in 1:nrow(fraglisttest)) #iterate through all frags in list
{
  string <- as.character(fraglisttest[t,1]) #for each iteration, take the SMILES string as characters and store it as "string"
  samplesmiles <- string #Stores "string" as a SMILES file
  sampleSDF <- ChemmineR::smiles2sdf(samplesmiles) #Converts SMILES file to SDF file
  MCS = batchfull[1] #Loads comparison library
  Tanimoto = 0
  Index = 0 #Sets up variables Tanimoto Match Value and Index Number
  for(i in 1:300) #Runs MCS comparing sample frag to all frags in our library
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
  
  dout = rbind(dout, data.frame(string, Tanimoto)) #Creates output dataframe of the original fragment SMILES and the Tanimoto match value

}