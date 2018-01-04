
asSMIstringtest <- function(x) {
  samplesmiles <- x #Stores SMILES string input
  sampleSDF <- ChemmineOB::convertFormat("SMI", "SDF", samplesmiles) #Converts user-input SMILES sample string to SDF format
  MCS = batchSDF[1]
  Tanimoto = 0
  Index = 0 #Sets up variables MCS, Tanimoto, and Index
  for(i in 1:300)
  {

    test3 <- fmcsR::fmcs(batchSDF[i], sampleSDF[1], fast = TRUE) #Compares sample molecule to each molecule in comparison collection
    s <- readr::parse_number(test3[4]) #Removes words from output to isolate comparison value

    if (s > Tanimoto) #Checks if comparison value is the best match

    {
      Tanimoto = s
      MCS = batchSDF[i]
      Index = i #Continues checking against the rest of the comparison collection
    }
  }
  #Tanimoto
  #Index
  smiles <- ChemmineR::sdf2smiles(MCS[1]) #Stores values for best match
  ChemmineR::write.SMI(smiles, file = "smiles.smi") #Creates output file
  print(paste("MCSS Match SMILES: ", smiles))
  print(paste("Tanimoto Index: ", Tanimoto))
  print(paste("Library Index: ", Index))
  plot(MCS[1], regenCoords = TRUE, print = FALSE)
}

x <- "NC1=CC=CC=C1"
asSMIstringtest(x)
