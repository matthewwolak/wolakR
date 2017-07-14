################################################################################
# 	Altered from the code written by Jarrod Hadfield
	################################################
#' Parameter expanded inverse Wishart.
#'
#' Simulate inverse Wishart distribution with parameter expansion on the
#' hyperparameters.
#'
#' Details on the meaning of \code{V, nu, alpha.mu, alpha.V} can be found
#' in the \code{\link[MCMCglmm]{MCMCglmm}} documentation covering the way
#' to specify variance structure priors in the \code{prior} argument.
#'
#' @aliases rpeIW
#' @export
#' @param n Integer of how many samples to generate.
#' @param V Numeric matrix of the expected (co)variances
#' @param nu Numeric for the degree of belief parameter for the inverse-
#'   Wishart.
#' @param alpha.mu Numeric vector of means for the redundant working parameters.
#' @param alpha.V Numeric matrix of the covariance matrix for the redunant
#'   working parameters.
#'
#' @return A numeric vector (if the dimensions of \code{V} are 1) or matrix
#'   where each row designates a sample and the columns contain the matrix
#'   elements for a prior matrix.
#'
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link[MCMCglmm]{MCMCglmm}}, \code{\link[MCMCglmm]{rIW}}
#' @family prior functions
#' @examples
#' rpeIW(n = 10, V = diag(2), nu = 3,
#'	alpha.mu = rep(0, 2), alpha.V = diag(2)*1000)
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



