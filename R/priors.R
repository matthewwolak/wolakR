################################################################################
# 	Altered from the code written by Jarrod Hadfield
	################################################
#TODO add this documentation with `rpeIW()` so show up in single help file

#' Inverse Wishart (co)variance density.
#'
#' Density generating function for the inverse Wishart distribution.
#'
#' When \code{V = 1}, the inverse Wishart distribution is equivalent to an
#' inverse Gamma distribution with shape and scale parameters set to
#' \code{nu / 2}. In other words, the inverse Gamma is a special case of the
#' inverse Wishart (Hadfield 2015, p. 12), where the inverse Wishart is the
#' multivariate generalization of the inverse Gamma distribution.
#'
#' As \code{nu} goes to infinity, the point mass moves towards \code{V}. The
#' mode of the distributions is calculated by \code{(V * nu) / (nu + 2)}
#' (Hadfield 2015, p. 12).

#'
#' @aliases dIW
#' @export
#' @importFrom MCMCpack dinvgamma
#' 
#' @param x Vector of quantiles.
#' @param V Numeric matrix of the expected (co)variances
#' @param nu Numeric for the degree of belief parameter for the inverse-
#'   Wishart.
#' @param marginal Logical indicating whether the densities for a single
#'   variance (\code{FALSE}) or the marginal densities of the first variance
#'   in a covariance matrix \code{V[1, 1]} (\code{TRUE}) are to be returned.
#'
#' @return A \code{list} containing
#'   \describe{
#'     \item{x}{Numeric vector of quantiles originally passed to the function.}
#'     \item{y}{Numeric vector of densities.}
#'   }
#'
#' @author \email{matthewwolak@@gmail.com}
#' @references
#' Hadfield, J. 2015. MCMCglmm Course Notes. June 20, 2015.
#' @seealso \code{\link[MCMCglmm]{MCMCglmm}}, \code{\link[MCMCglmm]{rIW}},
#'   \code{rpeIW}
#' @family prior functions
#' @examples
#' xseq <- seq(from = 1e-16, to = 5, length = 1000)
#' # Plot range of prior distributions
#' ## start with inverse Gamma with small degree of belief and point mass ~ 0
#' IG0.002 <- dIW(xseq, V = 1, nu = 0.002)
#' IG0.02 <- dIW(xseq, V = 1, nu = 0.02)
#' IG0.2 <- dIW(xseq, V = 1, nu = 0.2)
#' ## end with point mass near V
#' IG1 <- dIW(xseq, V = 1, nu = 1)
#' plot(IG0.002, type = "n",
#'   main = "Inverse Gamma\nV = 1",
#'   xlab = "Variance", ylab = "Density",
#'   xlim = c(0, max(xseq)),
#'   ylim = c(0, max(c(IG0.002$y, IG0.02$y, IG0.2$y, IG1$y))))
#'  lines(IG0.02, lwd = 2, col = "red")
#'  lines(IG0.2, lwd = 2, col = "blue")
#'  lines(IG1, lwd = 2, col = "grey40")
#'  lines(IG0.002, lwd = 3, col = "black")
#'  legend("topright", lwd = 2, col = c("black", "red", "blue", "grey40"),
#' 	title = "nu", legend = as.character(c(0.002, 0.02, 0.2, 1)),
#' 	inset = 0.01)
#'
#' #######################
#' # Marginal variance
#' #######################
#' mar1 <- dIW(xseq, V = diag(2), nu = 1.002, marginal = TRUE)
#' # compare to IG0.002 above
#' plot(mar1, type = "n",
#'   main = "Marginal prior for a variance:\n IW(V = diag(2), nu = 1.002)",
#'   xlab = "Variance", ylab = "Density",
#'   xlim = c(0, max(xseq)), ylim = c(0, max(c(mar1$y, IG0.002$y))))
#'  lines(mar1, lwd = 2, col = "red")
#'  lines(IG0.002, lwd = 2, col = "black")
#'  legend("topright", col = c("red", "black"), lwd = 2,
#' 	legend = c("marginal prior",
#' 		"univariate prior\nIG(V=1, nu=0.002)"),
#' 	inset = 0.01)
dIW <- function(x, V = 1, nu = 1, marginal = FALSE){
  if(!marginal){
    if(is.matrix(V) && nrow(V) > 1){
      stop("'V' must be a scalar/single number when 'marginal = FALSE'")
    }
   return(list(x = x,
	y = MCMCpack:::dinvgamma(x, shape = nu / 2, scale = (nu * V) / 2)))
  } else{
      nu2 <- nu - dim(V)[1] + 1
      V2 <- (nu / nu2) * V[1, 1]
     return(list(x = x,
	y = MCMCpack:::dinvgamma(x, shape = nu2 / 2, scale =  (nu2 * V2) / 2)))
    }
}






################################################################################
# 	Altered from the code written by Jarrod Hadfield
	################################################
#' Parameter expanded inverse Wishart (co)variance.
#'
#' Simulate prior (co)variance matrix according to the inverse Wishart
#' distribution with parameter expansion on the hyperparameters.
#'
#' Details on the meaning of \code{V, nu, alpha.mu, alpha.V} can be found
#' in the \code{\link[MCMCglmm]{MCMCglmm}} documentation covering the way
#' to specify variance structure priors in the \code{prior} argument.
#'
#' @aliases rpeIW
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom MCMCglmm rIW
#' @param n Integer of how many samples to generate.
#' @param V Numeric matrix of the expected (co)variances
#' @param nu Numeric for the degree of belief parameter for the inverse-
#'   Wishart.
#' @param alpha.mu Numeric vector of means for the redundant working parameters.
#' @param alpha.V Numeric matrix of the covariance matrix for the redunant
#'   working parameters.
#'
#' @return A numeric vector (if the dimensions of \code{V} are 1) or matrix
#'   where each row designates a sampled (co)variance and the columns contain
#'   the matrix elements for the prior covariance matrix.
#'
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link[MCMCglmm]{MCMCglmm}}, \code{\link[MCMCglmm]{rIW}}
#' @family prior functions
#' @examples
#' set.seed(101)
#' (peIW <- rpeIW(n = 5, V = diag(2), nu = 3,
#'	alpha.mu = rep(0, 2), alpha.V = diag(2)*1000))
#' sapply(1:5, FUN = function(i){matrix(peIW[i, ], 2, 2)}, simplify = "array")
rpeIW <- function(n = 1, V, nu, alpha.mu, alpha.V){
  k <- nrow(V)
  alpha <- MASS::mvrnorm(n, alpha.mu, alpha.V)
  Valpha <- MCMCglmm::rIW(V, nu = nu, n = n)
  Vp <- sapply(1:n,
    FUN = function(i){amat <- diag(alpha[i,],k,k); amat %*%
	matrix(Valpha[i,],k,k) %*% amat})

  if(k == 1){
    return(Vp)
  } else{
      return(t(Vp))
    }
}



