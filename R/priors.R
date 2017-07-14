################################################################################
# 	Altered from the code written by Jarrod Hadfield
	################################################
#' Parameter expanded inverse Wishart.
#'
#' Simulate inverse Wishart distribution with parameter expansion on the
#' hyperparameters.
#'
#'FIXME To plot a prior distribution, \code{prior} should either specify an object of
#'   class \code{\link[coda]{mcmc}} that contains a simulated prior distribution
#'
#'FIXME @section Warning:
#'
#' @aliases postPlot
#' @export
#' @param posterior Object containing a marginal posterior distribution from
#'   from which kernel density estimates can be calulated by the
#'   \code{coda::density} function.
#' @param plotHist Should a histogram of the density be plotted as well.
#' @param histbreaks If \code{plotHist = TRUE}, then the number of breaks. See
#'   \code{\link[graphics]{hist}} for details.
#'
#' @param prior Should a line be added representing the prior distribution. If
#'   \code{NULL} (the default) then no line is added. A prior is plotted if
#'   either a \code{\link[coda]{mcmc}} object or prior specification is supplied
#'   (see Details).
#'
#' @return A \code{list} returned invisibly (assign to keep):
#'   \describe{
#'     \item{call}{The original \code{call}.}
#'     \item{postDensity}{Object of \code{class} "density" for the posterior
#'       distribution.}
#'     \item{priorDensity}{Object of \code{class} "density" for the prior
#'       distribution.}
#'     \item{bandwith}{The density bandwith from the posterior used to set some
#'         the bandwith in the prior under some constraints.FIXME}
#'     \item{histOut}{Object of \code{class} "histogram" for the posterior
#'         distribution.}
#'     \item{constraint}{A \code{character} indicating if any densities were
#'         adjusted because of boundary constraints.}
#'   }
#'
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link[coda]{densplot}}, \code{\link[graphics]{hist}}
#' @family prior functions
#' @examples
#' # Simulate to example MCMC chains:
#' ## Both are standard normal distributions of only positive values
#' ## However, one is on the boundary (near zero) while the other is not
#' normMCMC <- coda::mcmc(matrix(c(abs(rnorm(1000, 0.2, 1)), rnorm(1000, 10, 1)),
#'	ncol = 2))
#' par(mfrow = c(2, 1))
#'   postPlot(normMCMC[, 1], ylim = c(0, 1))
#'   postPlot(normMCMC[, 2], ylim = c(0, 1))
rpx <- function(n = 1, V, nu, alpha.mu, alpha.V){
  k <- nrow(V)
  alpha <- MASS::mvrnorm(n, alpha.mu, alpha.V)
  Valpha <- MCMCglmm::rIW(V, nu=nu, n=n)
  Vp <- sapply(1:n, FUN=function(i){diag(alpha[i,],k,k)%*%matrix(Valpha[i,],k,k)%*%diag(alpha[i,],k,k)})

  if(k==1){
    return(Vp)
  }else{
    return(t(Vp))
  }
}

rpeIW_diagSet <- function(n = 1, V, nu, alpha.mu, alpha.V){
  k <- nrow(V)
  alpha <- MASS::mvrnorm(n, alpha.mu, alpha.V)
  Valpha <- MCMCglmm::rIW(V, nu=nu, n=n)
  Vp <- sapply(1:n, FUN=function(i){amat <- diag(alpha[i,],k,k); amat%*%matrix(Valpha[i,],k,k)%*%amat})

  if(k==1){
    return(Vp)
  }else{
    return(t(Vp))
  }
}


rpeIW_array <- function(n = 1, V, nu, alpha.mu, alpha.V){
  k <- nrow(V)
  alpha <- MASS::mvrnorm(n, alpha.mu, alpha.V)

  Valpha <- MCMCglmm::rIW(V, nu=nu, n=n)
  Valpha <- sapply(1:n, FUN = function(i){matrix(Valpha[i,], k, k)}, simplify = "array")

  Vp <- sapply(1:n, FUN=function(i){amat <- diag(alpha[i,]); amat %*% Valpha[,,i] %*% amat})

  if(k==1){
    return(Vp)
  }else{
    return(t(Vp))
  }
}

######################
# `rpeIW_diagSet()` faster than `rpx()` for kin=2 & 3
library(microbenchmark)
kin <- 3
system.time(mrk <- microbenchmark(rpx(10000, V=diag(kin), nu=kin+1, alpha.mu=rep(0,kin), alpha.V=diag(kin)*1000),
	rpeIW_diagSet(10000, V=diag(kin), nu=kin+1, alpha.mu=rep(0,kin), diag(kin)*1000),
		times = 100))
plot(mrk)
mrk

# `repIW_diagSet()` is faster than `rpeIW_array()` for kin=2, 3, and 5
kin <- 2
system.time(mrk <- microbenchmark(rpeIW_diagSet(10000, V=diag(kin), nu=kin+1, alpha.mu=rep(0,kin), diag(kin)*1000),
	rpeIW_array(10000, V=diag(kin), nu=kin+1, alpha.mu=rep(0,kin), diag(kin)*1000),
		times = 100))
plot(mrk)
mrk

