#CTS input from JSON

CTSin <- fromJSON("CTSinput2.json", simplifyVector = TRUE, simplifyDataFrame = simplifyVector, simplifyMatrix = simplifyVector, flatten = FALSE)

CTSdf <- as.data.frame(CTSin)
numfrags = nrow(CTSdf)
length = numfrags - 1
i = 1
finaldf <- data.frame("Original Fragment Smiles" = character(), "MCSS Match SMILES" = character(), "Tanimoto Match Index" = character(), "Sigma Value" = character(), "Sigma Meta Value" = character(), "Sigma Para Value" = character())
finallist <- list()

#while(i < numfrags){
#  frag <- CTSdf[1,1]
#}

for(i in 1:length){
  frag <- CTSdf[i,1]
  fragSDF <- smiles2sdf(as.character(frag))
  #put code/if else here to narrow down search list to aliphatic/meta/para/ortho
  fmcsout <- fmcsforCTS(fragSDF) #run fmcsforCTS function
  if (as.character(fmcsout[1,1]) == "Error"){
    
    break
  }
  #escape loop if return data.frame is nomatch/error
  #use return data.frames/lists from function to append onto finallist
  #
}

#check if for loop broke
#create final dataframe from finallist
#create output json from final data.frame


