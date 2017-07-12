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
#' @family MCMC posterior distribution helper functions
#' @examples
#' # Simulate to example MCMC chains:
#' ## Both are standard normal distributions of only positive values
#' ## However, one is on the boundary (near zero) while the other is not
#' normMCMC <- coda::mcmc(matrix(c(abs(rnorm(1000, 0.2, 1)), rnorm(1000, 10, 1)),
#'	ncol = 2))
#' postTable(normMCMC)
postTable <- function(mcpost, ind = NULL, sigdig = 3, ...){
  if(!is.null(ind)) mcpost <- mcpost[, ind]
  hpdint <- signif(coda::HPDinterval(mcpost), sigdig)
  postmode <- signif(MCMCglmm::posterior.mode(mcpost), sigdig)
  postmean <- signif(if(length(postmode)==1) mean(mcpost) else apply(mcpost, MARGIN = 2, FUN = mean), sigdig)
  compnames <- if(is.null(colnames(mcpost))) names(mcpost) else colnames(mcpost)
    if(is.null(compnames)) compnames <- "V1"
  out <- t(data.frame(component = compnames, summary = paste(paste0(postmode, ", ", postmean), paste0("(", hpdint[, 1], ",", hpdint[, 2], ")"), sep = " ")))
  out[2, ] <- gsub(pattern = "e-", replacement = paste0("\U00D7", "10-"), x = out[2, ])
 out
}





################################################################################
#' Posterior distribution plot.
#'
#' Pretty plot of a posterior distribution, including summary statistics such as the
#' posterior mean and credible interval.
#'
#' To plot a prior distribution, \code{prior} should either specify an object of
#'   class \code{\link[coda]{mcmc}} that contains a simulated prior distribution
#'   or it should contain a \code{list} with the prior specification according
#'   to the instructions for specifying G-structure priors in
#'   \code{\link[MCMCglmm]{MCMCglmm}}.
#'
#' Add more about prior plotting, specifically how parameter expanded priors are
#'   simulated versus inverse-Wishart/Gamma are derived straight from their
#'   distribution functions. Further, add how to specify which marginal component
#'   to plot, how to specify `n` samples, fix in residual covariance matrices,
#'   and what about fixed effects? Also, note how density is calculated
#'   after 'trimming' to range of the posterior samples.
#'
#' @section Warning:
#' The use of the \dots argument is untested. Arguments of the same name that
#' are shared among \code{densplot}, \code{hist}, or \code{abline} will affect
#' each of these aspects of the plot if supplied to \code{postPlot}. Proceed
#' with caution!
#'
#' @aliases postPlot
#' @export
#' @param posterior Object containing a marginal posterior distribution from
#'   from which kernel density estimates can be calulated by the
#'   \code{coda::density} function.
#' @param plotHist Should a histogram of the density be plotted as well.
#' @param histbreaks If \code{plotHist = TRUE}, then the number of breaks. See
#'   \code{\link[graphics]{hist}} for details.
#' @param prior Should a line be added representing the prior distribution. If
#'   \code{NULL} (the default) then no line is added. A prior is plotted if
#'   either a \code{\link[coda]{mcmc}} object or prior specification is supplied
#'   (see Details).
#' @param denslwd,hpdlwd,meanlwd,priorlwd The line widths for the lines to be plotted
#'   representing the kernel density estimate (\code{denslwd}), credible
#'   interval limits (\code{hpdlwd}), or posterior mean (\code{meanlwd}). See
#'   \code{\link[graphics]{par}} for details. 
#' @param denslty,hpdlty,meanlty,priorlty The line type for the lines to be plotted
#'   representing the kernel density estimate (\code{denslty}), credible
#'   interval limits (\code{hpdlty}), or posterior mean (\code{meanlty}). See
#'   \code{\link[graphics]{par}} for details. 
#' @param denscol,hpdcol,meancol,priorcol The line colors for the lines to be plotted
#'   representing the kernel density estimate (\code{denscol}), credible
#'   interval limits (\code{hpdcol}), or posterior mean (\code{meancol}). See
#'   \code{\link[graphics]{par}} for details.
#' @param ylim Y-axis limits to the plotting region, passed to
#'   \code{\link[coda]{densplot}}. See \code{\link[graphics]{plot.window}}.
#' @param at1,at2 The points at which tick-marks are to be drawn on the
#'   x-axis (\code{at1}) and y-axis (\code{at2}). See
#'   \code{\link[graphics]{axis}} for details.
#' @param labels1,labels2 Labels to be placed at the tick-marks of the
#'   x-axis (\code{labels1}) and and y-axis (\code{labels2}). See
#'   \code{\link[graphics]{axis}} for details.
#' @param \dots Additional arguments passed to \code{densplot}, \code{hist},
#'   \emph{and} \code{abline}.
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link[coda]{densplot}}, \code{\link[graphics]{hist}}
#' @family MCMC posterior distribution helper functions
#' @examples
#' # Simulate to example MCMC chains:
#' ## Both are standard normal distributions of only positive values
#' ## However, one is on the boundary (near zero) while the other is not
#' normMCMC <- coda::mcmc(matrix(c(abs(rnorm(1000, 0.2, 1)), rnorm(1000, 10, 1)),
#'	ncol = 2))
#' par(mfrow = c(2, 1))
#'   postPlot(normMCMC[, 1], ylim = c(0, 1))
#'   postPlot(normMCMC[, 2], ylim = c(0, 1))
postPlot <- function(posterior, plotHist = TRUE, histbreaks = 100,
	prior = NULL,
	denslwd = 6, denslty = "solid", denscol = "black",
	hpdlwd = 6, hpdlty = "dashed", hpdcol = "black",
	meanlwd = 7, meanlty = "12", meancol = "red",
	priorlwd = 3, priorlty = "solid", priorcol = "green",
	ylim,
	at1 = NULL, at2 = NULL,
	labels1 = NULL, labels2 = NULL, ...){

  if(missing(ylim)){
    pdens <- stats::density(posterior)$y
    histout <- if(plotHist){
      graphics::hist(posterior, breaks = histbreaks, plot = FALSE)$density
    } else 0
    ylimit <- c(0, max(c(pdens, histout)))
  } else ylimit <- ylim

  coda::densplot(posterior, show.obs = FALSE, axes = FALSE,
	ylim = ylimit, lwd = denslwd, ...)
  if(plotHist) graphics::hist(posterior, breaks = histbreaks, freq = FALSE,
	add = TRUE, ...)
  graphics::abline(v = coda::HPDinterval(posterior),
	lty = hpdlty, lwd = hpdlwd, col = hpdcol, ...)
  graphics::abline(v = mean(posterior),
	col = meancol, lty = meanlty, lwd = meanlwd, ...)
  axis(1, at = at1, labels = labels1)
  axis(2, at = at2, labels = labels2)

  #######
  # Prior
  #######
  if(!is.null(prior)){
    if(coda::is.mcmc(prior)){
      #TODO add some sort of bandwith more to max(posterior)
      ## but stop density there. This will keep it from dropping down to zero
      ## Also,
      pr <- range(posterior)
      if(pr[1L] > 0 && pr[1L] < 1){
        prDens <- stats::density(prior[which(prior >= pr[1L] & prior <= pr[2L])],
	  from = 0)
      } else{
          prDens <- stats::density(prior[which(prior >= pr[1L] & prior <= pr[2L])])
        }
    }
    if(is.list(prior)){
      stop("prior as a list is not yet implemented. Supply a `mcmc` object") #TODO
      if(!is.matrix(V)) V <- as.matrix(V)
      rc <- nrow(V)
      if(nrow(V) != ncol(V)) stop("V must be a square symmetric matrix")
      if(!all(diag(V) > 0)){
        stop("V is not positive definite (V has 0 or negative diagonal values)")
      } 
      #TODO add check of V for positive definiteness
      ## could use `MCMCglmm:::is.positive.definite()` or `eigen()`

      if(is.null(alpha.mu) && is.null(alpha.V)){
        stop("Need to add how to do non-parameter expanded priors")  #TODO
      }
    }
    if(!coda::is.mcmc(prior) && !is.list(prior)){
      warning("prior is not a `list` or `mcmc` object - no prior added to plot")
    }

    # Add prior density to plot
    graphics::lines(prDens,
	col = priorcol, lty = priorlty, lwd = priorlwd, ...)
  }
 

}



