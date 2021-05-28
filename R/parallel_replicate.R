#' Parallel Implementation of replicate.
#'
#' Implements parallel version of \code{replicate} using
#' \code{parallel::mclapply()}.
#'
#' @aliases preplicate
#' @export
#' @import parallel
#' @param n integer: the number of replications.
#' @param expr the expression (a language object, usually a call) to evaluate
#'   repeatedly.
#' @param simplify logical or character string, see \code{replicate}.
#' @param ... other arguments to be passed to \code{mclapply}.
#' @return A list, vector, or matrix depending on arguments to \code{simplify}.
#'   See documentation for \code{replicate} for more details.
#' @author \email{matthewwolak@@gmail.com}
#' @examples
#' \dontrun{
#'   system.time(r <- replicate(10000, var(rnorm(500))))
#'   system.time(pr <- preplicate(10000, var(rnorm(500)), mc.cores = 3))
#' }
preplicate <- function(n, expr, simplify = "array", ...){
  require(parallel) 
  ans <- mclapply(integer(n), eval.parent(substitute(function(...) expr)), ...)
  if(!isFALSE(simplify) && length(ans))
    simplify2array(ans, higher = (simplify == "array"))
  else ans  
}  

