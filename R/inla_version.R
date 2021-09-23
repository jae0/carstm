inla_version = function( test="21.09.13" ) {

  trim <- function(string) {
      string <- gsub("^[^:]*:", "", string)
      string <- gsub("^[ \t]+", "", string)
      string <- gsub("[ \t]+$", "", string)
      return(string)
  }
  desc <- readLines(system.file("DESCRIPTION", package = "INLA"))
  version <- trim(desc[grep("^Version:", desc)])

  v = tstrsplit(version, "[[:punct:]]" )
  w = tstrsplit(test,    "[[:punct:]]" )

  newer_than_test = FALSE
  for( i in 1:length(v) ) {
    if (as.numeric(v[[i]]) < as.numeric(w[[i]]) ) {
      newer_than_test = TRUE
      break()
    }  
  }
  return( list(version=version, test=test, current_version_newer_than_test=newer_than_test ))
}