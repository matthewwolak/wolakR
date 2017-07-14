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
#' Pretty plot of a posterior distribution that can include summary statistics,
#' such as the posterior mean and credible interval as well as a prior density.
#'
#' To plot a prior distribution, \code{prior} should either specify an object of
#'   class \code{\link[coda]{mcmc}} that contains a simulated prior distribution
#'   or it should contain a \code{list} which contains the prior density (where
#'   the first two elements must be named \code{x} and \code{y} and contain the
#'   quantiles and density, respectively.). Note that a function to generate the
#'   prior distribution may be specified. However, this must return either a
#'   \code{mcmc} or \code{list} object.
#'
#' Add more about prior plotting, specifically how parameter expanded priors are
#'   simulated versus inverse-Wishart/Gamma are derived straight from their
#'   distribution functions. Also, note how density is calculated
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
#'
#' @param prior The prior distribution for a line be added to the plot. If
#'   \code{NULL} (the default) then no line is added. A prior is plotted if
#'   either a \code{\link[coda]{mcmc}} object or \code{list} with the prior 
#'   density is supplied. See Details.
#' @param prange Character argument specifying which range
#'   \code{c("prior", "posterior")} the prior distribution should span. Currently,
#'   only works when class \code{mcmc} provided to \code{prior}.#FIXME
#'
#' @param main Overall title for the plot. See \code{\link[graphics]{plot}} and
#'   \code{\link[graphics]{title}}.
#' @param sub Sub-title for the plot. See \code{\link[graphics]{plot}} and
#'   \code{\link[graphics]{title}}.
#' @param ylim Y-axis limits to the plotting region, passed to
#'   \code{\link[coda]{densplot}}. See \code{\link[graphics]{plot.window}}.
#'
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
#' @param at1,at2 The points at which tick-marks are to be drawn on the
#'   x-axis (\code{at1}) and y-axis (\code{at2}). See
#'   \code{\link[graphics]{axis}} for details.
#' @param labels1,labels2 Labels to be placed at the tick-marks of the
#'   x-axis (\code{labels1}) and and y-axis (\code{labels2}). See
#'   \code{\link[graphics]{axis}} for details.
#' @param \dots Additional arguments passed to \code{densplot}, \code{hist},
#'   \emph{and} \code{abline}.
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
	prior = NULL, prange = c("prior", "posterior"),
	main, sub, ylim,
	denslwd = 6, denslty = "solid", denscol = "black",
	hpdlwd = 6, hpdlty = "dashed", hpdcol = "black",
	meanlwd = 7, meanlty = "12", meancol = "red",
	priorlwd = 3, priorlty = "solid", priorcol = "green",
	at1 = NULL, at2 = NULL,
	labels1 = NULL, labels2 = NULL, ...){

  cl <- match.call()
  # bw.nrd() is like `bwf()` in `coda::densplot()` 
  poDens <- stats::density(posterior, bw = "nrd", kernel = "gaussian", n = 2^13)
  bw <- poDens$bw
  main.title <- if(missing(main)) "" else main
  sub.title <- if(missing(sub)) "" else sub
  histout <- if(plotHist){
               graphics::hist(posterior, breaks = histbreaks, plot = FALSE)
             } else 0
  if(missing(ylim)){
    ylimit <- c(0, max(c(poDens$y, histout$density)))
  } else ylimit <- ylim

  # Adjust density at boundaries
  ## reflect density estimates back across a boundary
  constraint <- "unbounded"
  po <- posterior
  if(max(po) <=1 && 1 - max(po) < 2 * bw){
    if(min(po) >= 0 && min(po) < 2 * bw){
      constraint <- "zero-one"
      po <- c(po, -po, 2 - po)
    } 
    if(min(po) < 0 && 1 + min(po) < 2 * bw){
      constraint <- "correlation"
      po <- c(po, -2 - po, 2 - po)
    } else{
        constraint <- "upbound1"  #<-- most likely correlation converge on 1
        po <- c(po, 2-po)
      }
  } else{
      if(min(po) >= 0 && min(po) < 2 * bw){
        constraint <- "positive"
        po <- c(po, -po)
      }
      if(min(po) >= -1 && 1 + min(po) < 2 * bw){
        constraint <- "loboundNeg1"  #<-- most likely correlation converge on -1
        po <- c(po, -2 - po)
      }
    }
  poDens <- stats::density(po, kernel = "gaussian", width = 4 * bw, n = 2^13)
  if(constraint == "zero-one"){
    poDens$y <- 3 * poDens$y[poDens$x >= 0 & poDens$x <=1]
    poDens$x <- poDens$x[poDens$x >= 0 & poDens$x <= 1]
  }
  if(constraint == "correlation"){
    poDens$y <- 3 * poDens$y[poDens$x >= -1 & poDens$x <=1]
    poDens$x <- poDens$x[poDens$x >= -1 & poDens$x <= 1]
  }
  if(constraint == "upbound1"){  #<-- most likely correlation converge on 1
    poDens$y <- 2 * poDens$y[poDens$x <=1]
    poDens$x <- poDens$x[poDens$x <= 1]
  }
  if(constraint == "loboundNeg1"){  #<-- most likely correlation converge on -1
    poDens$y <- 2 * poDens$y[poDens$x >= -1]
    poDens$x <- poDens$x[poDens$x >= -1]
  }
  if(constraint == "positive"){
    poDens$y <- 2 * poDens$y[poDens$x >= 0]
    poDens$x <- poDens$x[poDens$x >= 0]
  }
  plot(poDens, type = "n", axes = FALSE,
	main = main.title, sub = sub.title, ylim = ylimit, ...)
  graphics::lines(poDens, lwd = denslwd, ...)
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
    if(is.list(prior)){
      if(names(prior)[1L] != "x" || names(prior)[2L] != "y"){
        stop("prior must be either an object of class `mcmc`, a `list` containing a density, or a function that returns either of these")
      }
      prDens <- prior
    }

    if(coda::is.mcmc(prior)){
      ## range of prior/posterior
      #TODO argument matching for prange
      prra <- if(prange == "posterior") range(posterior) else range(prior)
      # General formula for finding bounds
      ## Lower: 2*bound - prior (then select all >= bound)
      ## Upper: 2*bound - prior (then select all <= bound)
      if(constraint == "zero-one"){
        pr <- prior[prior >= 0 & prior <= 1]
        pr <- c(pr, -pr, 2 - pr)
        prDens <- stats::density(pr, kernel = "gaussian", width = 4 * bw, n = 2^13)
        prDens$y <- 3 * prDens$y[prDens$x >= 0 & prDens$x <=1]
        prDens$x <- prDens$x[prDens$x >= 0 & prDens$x <= 1]
      }
      if(constraint == "upbound1" | constraint == "loboundNeg1"){
        if(min(prior) >= -1 && 1 + min(prior) < 2 * bw &&
          max(prior) <= 1 && 1 - max(prior) < 2 * bw) constraint <- "correlation"
      }
      if(constraint == "correlation"){
        pr <- prior[prior >= -1 & prior <= 1]
        pr <- c(pr, -2 - pr, 2 - pr)
        prDens <- stats::density(pr, kernel = "gaussian", width = 4 * bw, n = 2^13)
        prDens$y <- 3 * prDens$y[prDens$x >= -1 & prDens$x <=1]
        prDens$x <- prDens$x[prDens$x >= -1 & prDens$x <= 1]
      }
      if(constraint == "upbound1"){  #<-- most likely correlation converge on 1
        pr <- prior[prior <= 1]
        pr <- c(pr, 2 - pr)
        prDens <- stats::density(pr, kernel = "gaussian", width = 4 * bw, n = 2^13)
        prDens$y <- 2 * prDens$y[prDens$x <=1]
        prDens$x <- prDens$x[prDens$x <= 1]
      }
      if(constraint == "loboundNeg1"){  #<-- most likely correlation converge on -1
        pr <- prior[prior >= -1]
        pr <- c(pr, -2 - pr)
        prDens <- stats::density(pr, kernel = "gaussian", width = 4 * bw, n = 2^13)
        prDens$y <- 2 * prDens$y[prDens$x >= -1]
        prDens$x <- prDens$x[prDens$x >= -1]
      }
      if(constraint == "positive"){
#FIXME issues when prior range is used on very flat parameter expanded prior
## Density is high near zero, so usually above posterior plot and off figure panel
        #XXX NOTE calculates own `bw` and does not use `width` with posterior bw
        pr <- prior[prior >= 0 & prior <= prra[2L]]
        pr <- c(pr, -pr, 2*prra[2L] - pr)
        prDens <- stats::density(pr, bw = "nrd", kernel = "gaussian", n = 2^13)
        prDens$y <- 3 * prDens$y[prDens$x >= 0 & prDens$x <= prra[2L]]
        prDens$x <- prDens$x[prDens$x >= 0 & prDens$x <= prra[2L]]
#        pr <- prior[prior >= 0]
#        pr <- c(pr, -pr)
#        prDens <- stats::density(pr, bw = "nrd", kernel = "gaussian", n = 2^13)
#        prDens$y <- 2 * prDens$y[prDens$x >= 0]
#        prDens$x <- prDens$x[prDens$x >= 0]
      }
      if(constraint == "unbounded"){
        #XXX NOTE calculates own `bw` and does not use `width` with posterior bw
        pr <- prior[prior >= prra[1L] & prior <= prra[2L]]
        pr <- c(pr, 2*prra[1L] - pr, 2*prra[2L] - pr)
        prDens <- stats::density(pr, bw = "nrd", kernel = "gaussian", n = 2^13)
        prDens$y <- 3 * prDens$y[prDens$x >= prra[1L] & prDens$x <= prra[2L]]
        prDens$x <- prDens$x[prDens$x >= prra[1L] & prDens$x <= prra[2L]]
      }
    }

    if(!coda::is.mcmc(prior) && !is.list(prior)){
      warning("prior is not a `list` or `mcmc` object - no prior added to plot")
    } else{
    # Add prior density to plot
        graphics::lines(prDens,
	  col = priorcol, lty = priorlty, lwd = priorlwd, ...)
      }
  }
 
 return(invisible(list(call = cl,
	postDensity = poDens, priorDensity = prDens,
	bandwith = bw,
	histOut = histout,
	constraint = constraint)))
	
}


