SDI <- function(t, dfin){
  tmpind <- which(dfin$trait == t)
  fm <- unlist(dfin[tmpind, c("female", "male")])
  if(any(is.na(fm))){
    return(NA)
  } 
  if(fm[1] == fm[2]){
    return(0)
  }
  larger <- which.max(fm)
  smaller <- which.min(fm)
  index <- (fm[larger] / fm[smaller]) - 1
 if(larger == 1) return(index) else return(-1 * index)
}

