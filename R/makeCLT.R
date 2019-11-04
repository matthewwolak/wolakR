#' Cohort Life Table Constructor
#'
#' Given age-specific cohort size (n) and fecundity (m), construct a cohort life
#' table of class \sQuote{CLT}.
#' 
#' @aliases makeCLT is.CLT
#' @param n A vector of cohort-specific sizes/number of individuals. Should
#'   include a last cohort with n=0, else the function adds this.
#' @param m A vector of cohort-specific fecundities. Should be of length
#'   \code{length(n)-1} or else have \code{NA} at \code{n[length(n)]=0}.
#' @param x An object of class \sQuote{CLT}.
#' @param n_x A vector containing the number of individuals in each cohort
#' @param a_x A vector containing Annual survival for each cohort
#'
#' @return A \code{vector} or \code{data.frame} of class \sQuote{CLT} containing
#'   one or all of \code{names}:
#'   \describe{
#'     \item{x }{Cohort ages/classes}
#'     \item{n_x }{Number of individuals in each cohort}
#'     \item{l_x }{Survivorship}
#'     \item{a_x }{Annual survival}
#'     \item{d_x }{Number of individuals that die}
#'     \item{q_x }{Annual mortality}
#'     \item{m_x }{Fecundity}
#'   }
#' @author \email{matthewwolak@gmail.com}
#' @examples
#'   makeCLT(n = c(10, 5, 3, 2, 0), m = c(0, 1, 3, 2))
#' @export
makeCLT <- function(n, m){
  if(tail(n, 1) != 0){
    warning(cat("adding age class so ages span: 0-", length(n), " with n(x=", length(n), "=0\n"))
    n <- c(n, 0)
  }

  clt <- data.frame(x = seq(1, length(n), 1)-1,
	n_x = n,
	l_x = NA,
	a_x = NA,
	d_x = NA,
	q_x = NA,
	m_x = ifelse(length(n) == length(m), m, c(m, NA)))
 return(structure(clt, class = "CLT"))
}

#' @method is CLT
#' @rdname makeCLT
#' @export
is.CLT <- function(x) inherits(x, "CLT")

################################################################################
# makeCLT accesories
################################################################################
#####################################
#	l_x
#####################################
#' @rdname makeCLT
#' @export
lxFun <- function(x, ...){
  UseMethod("lxFun", x)
}

#' @rdname makeCLT
#' @method lxFun default
#' @export
  lxFun.default <- function(n_x){
    n_x / n_x[1]
  }
#' @rdname makeCLT
#' @method lxFun CLT
#' @export
  lxFun.CLT <- function(x){
    x$l_x <- x$n_x / x$n_x[1]
   return(x)
  }

#####################################
#	a_x
#####################################
#' @rdname makeCLT
#' @export
axFun <- function(x, ...){
  UseMethod("axFun", x)
}

#' @rdname makeCLT
#' @method axFun default
#' @export
  axFun.default <- function(n_vec){
    N <- length(n_vec)
   return(c(sapply(seq(1, N-1, 1), FUN = function(x){ n_vec[x+1] / n_vec[x]}), NA))
  }
#' @rdname makeCLT
#' @method axFun CLT
#' @export
  axFun.CLT <- function(x){
    n <- x$n_vec
    N <- length(n)
    x$a_x <- c(sapply(seq(1, N-1, 1), FUN = function(i){ n[i+1] / n[i]}), NA)
   return(x)
  }

#####################################
#	d_x
#####################################
#' @rdname makeCLT
#' @export
dxFun <- function(x, ...){
  UseMethod("dxFun", x)
}

#' @rdname makeCLT
#' @method dxFun default
#' @export
  dxFun.default <- function(n_vec){
    N <- length(n_vec)
   return(c(sapply(seq(1, N-1, 1), FUN = function(x){ n_vec[x] - n_vec[x+1]}), NA))
  }
#' @rdname makeCLT
#' @method dxFun CLT
#' @export
  dxFun.CLT <- function(x){
    n <- n_vec
    N <- length(n)
    x$d_x <- c(sapply(seq(1, N-1, 1), FUN = function(i){ n[i] - n[i+1]}), NA)
   return(x)
  }

#####################################
#	q_x
#####################################
#' @rdname makeCLT
#' @export
qxFun <- function(x, ...){
  UseMethod("qxFun", x)
}

#' @rdname makeCLT
#' @method qxFun default
#' @export
  qxFun.default <- function(a_vec){
    1 - a_vec
  }
#' @rdname makeCLT
#' @method qxFun CLT
#' @export
  qxFun.CLT <- function(x){
    x$q_x <- 1 - x$a_vec
   return(x)
  }
    

