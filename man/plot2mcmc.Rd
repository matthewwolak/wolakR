% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCMCaccessories.R
\name{plot2mcmc}
\alias{plot2mcmc}
\title{Two MCMC chain plot.}
\usage{
plot2mcmc(x1, x2 = NULL, smooth = FALSE, bwf, save = FALSE, ...)
}
\arguments{
\item{x1}{Object of class \code{mcmc} containing MCMC chain(s).}

\item{x2}{Optional object of class \code{mcmc} containing a second set of
MCMC chain(s).}

\item{smooth}{Logical. See \code{\link[coda]{traceplot}} for details.#FIXME}

\item{bwf}{Character indicating function to calculate the bandwith. See
\code{\link[coda]{traceplot}} for details.#TODO correspond with `postPlot()`}

\item{save}{Either a logical or character value. When \code{FALSE} (default)
plots are not saved to the hard drive. If \code{TRUE}, the plots are saved
as a pdf with a default filename to the current working directory. If a
\code{character} is specified, this will be used as the base filename to
save the plots as pdfs in the current working directory.}

\item{\dots}{Additional arguments passed to \code{densplot} and
\code{traceplot}.#FIXME update densplot and traceplot to own functions}
}
\value{
Nothing returned invisibly at the moment (assign to keep):#FIXME
}
\description{
Pretty plot of two MCMC chains side-by-side.
}
\section{Warning}{

The function attempts to match parameters by the column names, but may have
trouble where vastly different objects are in each of the two MCMC objects.
}

\examples{
# Simulate two sets of two example MCMC chains:
## Both are standard normal distributions of only positive values
## However, one is on the boundary (near zero) while the other is not
normMCMC <- coda::mcmc(matrix(c(abs(rnorm(1000, 0.2, 1)), rnorm(1000, 10, 1)),
ncol = 2))
normMCMC2 <- coda::mcmc(matrix(c(abs(rnorm(1000, 0.2, 1)), rnorm(1000, 10, 1)),
ncol = 2))
plot2mcmc(normMCMC, normMCMC2)
}
\seealso{
\code{postPlot} and \code{\link[coda]{densplot}}

Other MCMC posterior distribution helper functions: 
\code{\link{postPlot}()},
\code{\link{postTable}()}
}
\author{
\email{matthewwolak@gmail.com}
}
\concept{MCMC posterior distribution helper functions}
