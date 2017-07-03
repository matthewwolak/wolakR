#' Posterior distribution summary.
#'
#' Reports a table containing the mode, mean, and highest posterior density of
#' the marginal posterior distribution.
#'
#' @aliases postTable
#' @export
#' @param mcpost Object of class \code{mcmc} containing marginal posterior
#'   distributions.
#' @param ind Column indices for the subset of \code{mcpost} for which the
#'   summary statistics are to be returned.
#' @param sigdig Integer value for the number of significant digits to report.
#' @param \dots Additional arguments.
#' @return A character \code{matrix} of two rows and \code{ncol(mcpost)} columns
#'   (\code{length(ind)} when \code{ind} is not \code{NULL}). The first row
#'   contains component names from \code{mcpost}. The second row contains a
#'   a string representing the posterior mode then posterior mean, separated by
#'   a comma, followed by the lower and upper limits of the highest posterior
#'   density, enclosed in parentheses.
#' @author \email{matthewwolak@@gmail.com}
#' @examples
#' # Simulate to example MCMC chains:
#' ## Both are standard normal distributions of only positive values
#' ## However, one is on the boundary (near zero) while the other is not
#' normMCMC <- as.mcmc(matrix(c(abs(rnorm(1000, 0.2, 1)), rnorm(1000, 10, 1)),
#'	ncol = 2))
#' postTable(normMCMC)
postTable <- function(mcpost, ind = NULL, sigdig = 3, ...){
  if(!is.null(ind)) mcpost <- mcpost[, ind]
  hpdint <- signif(HPDinterval(mcpost), sigdig)
  postmode <- signif(posterior.mode(mcpost), sigdig)
  postmean <- signif(if(length(postmode)==1) mean(mcpost) else apply(mcpost, MARGIN = 2, FUN = mean), sigdig)
  compnames <- if(is.null(colnames(mcpost))) names(mcpost) else colnames(mcpost)
    if(is.null(compnames)) compnames <- "V1"
  out <- t(data.frame(component = compnames, summary = paste(paste0(postmode, ", ", postmean), paste0("(", hpdint[, 1], ",", hpdint[, 2], ")"), sep = " ")))
  out[2, ] <- gsub(pattern = "e-", replacement = paste0("\U00D7", "10-"), x = out[2, ])
 out
}





