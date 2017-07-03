#' Sexual dimorphism index.
#'
#' Calculates the sexual dimorphism index of Lovich & Gibbons (1992):
#' \code{(larger sex / smaller sex) - 1} and arbitrarily made positive when
#' females are the larger sex or negative when males are the larger sex.
#'
#' @aliases SDI
#' @export
#' @param t Character object naming a trait in \code{dfin}.
#' @param dfin Data frame which must contain at least three columns: one
#'   containing trait names that match \code{t}, a column with female
#'   measurements that must match "female" exactly, and a column with male
#'   measurements that must match "male" exactly.
#' @return A numeric representing the calculated sexual dimorphism index for
#'   trait \code{t}.
#' @references
#' Lovich, JE & Gibbons, JW. 1992. A review of techniques for quantifying sexual
#'   size dimorphism. Growth, Development, and Ageing 56(4):269-281.
#' @author \email{matthewwolak@@gmail.com}
#' @examples
#' popDat <- data.frame(trait = c("length", "weight", "coloration", "intelligence"),
#'    female = c(4.2, 255.0, 0.57, 144),
#'    male = c(3.7, 255.0, 0.89, 143))
#' SDI("length", popDat)
#' sapply(popDat$trait, FUN = SDI, dfin = popDat)
#' stopifnot(SDI("weight", popDat) == 0.0)
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

