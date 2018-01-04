
choosetype <- function(fpath) {
  samplefile <- as.character(fpath)
  if (file_ext(samplefile) == "json") {
    f <- read.csv(fpath)
    handleJSON(f)
  }
  else if (file_ext(samplefile) == "JSON") {
    f <- read.csv(fpath)
    handleJSON(f)
  }
  else if (file_ext(samplefile) == "sdf") {
    f <- read.SDFset(fpath)
    handleSDF(f)
  }
  else if (file_ext(samplefile) == "SDF") {
    f <- read.SDFset(fpath)
    handleSDF(f)
  }
  else if (file_ext(samplefile) == "smi") {
    f <- read.SMIset(fpath)
    handleSMI(f)
  }
  else if (file_ext(samplefile) == "SMI") {
    f <- read.SMIset(fpath)
    handleSMI(f)
  }
  else if (file_ext(samplefile) == "smiles") {
    f <- read.SMIset(fpath)
    handleSMI(f)
  }
  else handleSMIstring(fpath)
}
  #else if (class(fpath) == "character") {
   # handleSMIstring(fpath)
  #}
  #else if (file_ext(samplefile) == .csv) {

  #}
  #else print("Input is not an accepted file type.")
  #}

choosetype("batch2.sdf")
