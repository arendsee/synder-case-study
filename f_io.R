# io
#  * give prefixes to input files (from the Data Package)
#  * give prefixes to the output files (to the archive)

get_in <- function(f){
  fin <- file.path("arendsee", "synder-nfyc", "archive", f) 
  if(!file.exists(fin)){
    stop("Cannot find file ", fin)
  }
  fin
}

get_out <- function(f){
  if(!dir.exists("output")){
    dir.create("output")
  }
  file.path("output", f)
}
