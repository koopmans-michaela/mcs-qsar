stringfrag <- "NC1=CC=CC=C1*"
stringmatch <- "*C1=CC=CC=C1"
matchfrag <- smiles2sdf(stringmatch)
frag <- smiles2sdf(stringfrag)
outm <- fmcsR::fmcs(matchfrag, frag) #Runs MCS between match molecule and original fragment to get the output information in MCS format for use in the plotMCS visualization function
fmcsR::plotMCS(outm) 