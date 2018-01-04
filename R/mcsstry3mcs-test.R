samplejson <- jsonlite::read_json("samplejson.csv", simplifyVector = TRUE, flatten = FALSE) #Reads JSON input
samplesmiles <- samplejson[["data"]][["format"]][["smiles"]] #Parses JSON to get and store SMILES string
sampleSDF <- ChemmineOB::convertFormat("SMI", "SDF", samplesmiles) #Converts user-input SMILES sample string to SDF format

#batchSDF <- ChemmineR::read.SDFset("batch2.sdf") #Imports comparison collection in SDF format
#sampleSDF <- read.SDFset("sample.sdf") #Imports sample molecule in SDF format
MCS = batchSDF[1]
Tanimoto = 0
Index = 0 #Sets up variables MCS, Tanimoto, and Index
for(i in 1:300)
{

  test3 <- fmcsR::fmcs(batchSDF[i], sampleSDF[1], fast = TRUE) #Compares sample molecule to each molecule in comparison collection
  sigma <- readr::parse_number(test3[4]) #Removes words from output to isolate comparison value

  if (sigma > Tanimoto) #Checks if comparison value is the best match

  {
    Tanimoto = sigma
    MCS = batchSDF[i]
    Index = i          #Continues checking against the rest of the comparison collection
  }
}
#Tanimoto
#Index
smiles <- ChemmineR::sdf2smiles(MCS[1]) #Stores values for best match
ChemmineR::write.SMI(smiles, file = "smiles.smi") #Creates output file
print(paste("MCSS Match SMILES: ", smiles))
print(paste("Tanimoto Index: ", Tanimoto))
print(paste("MCSS Match Sigma: ", sigma))
print(paste("Library Index: ", Index))
