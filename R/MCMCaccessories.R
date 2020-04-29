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
#' # Simulate two example MCMC chains:
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
#'   quantiles and densities, respectively.). Note that a function to generate the
#'   prior distribution may be specified. However, this must return either a
#'   \code{mcmc} or \code{list} object.
#'
#' TODO: Add more about prior plotting, specifically how parameter expanded
#'   priors are simulated versus inverse-Wishart/Gamma are derived straight from
#'   their distribution functions. Also, note how density is calculated#TODO
#'
#' If the prior is of class \code{mcmc} then the argument \code{prange} can be
#'   used to either extend or restrict the range over which the prior density
#'   estimation is conducted, in relation to the posterior distribution.
#'   The range alteration is done by specifying \code{posteriorN}, or appending
#'   the \code{"posterior"} character string with a numeric value indicating
#'   the desired adjustment of the range. Prior density will then be calculated
#'   after taking a subset of the current samples that fall in the range. Note,
#'   this process only makes sense when used in conjunction with
#'   \code{posterior}. By default, a \code{1} is used if there is no numeric
#'   suffix; indicating that a 1-to-1 mapping of the range should be used.
#'   Values of the suffix >1 will extend the range, while those <1 will contract
#'   it.
#'
#' @section Warning:
#' The use of the \dots argument is untested. Arguments of the same name that
#' are shared among \code{densplot}, \code{hist}, or \code{lines} will affect
#' each of these aspects of the plot if supplied to \code{postPlot}. Proceed
#' with caution!
#'
#' @aliases postPlot
#' @export
#' @param posterior Object containing a marginal posterior distribution from
#'   from which kernel density estimates can be calulated by the
#'   \code{coda::density} function.
#' @param bw Either a \code{numeric} or \code{character} for calculating the
#'   smoothing bandwidth used in the kernel density estimation. If
#'   \code{numeric}, the value is used as the standard deviation of the
#'   smoothing kernel. If \code{character}, the string is used as a standard
#'   rule to choose the bandwith. See details in \code{\link[stats]{density}}.
#'   Note that the default value \code{"nrd"} is not the default used in
#'   \code{\link[stats]{density}}.
#'
#' @param plotHist Should a histogram of the density be plotted as well.
#' @param histbreaks If \code{plotHist = TRUE}, then the number of breaks. See
#'   \code{\link[graphics]{hist}} for details.
#'
#' @param prior The prior distribution for a line to be added to the plot. If
#'   \code{NULL} (the default) then no line is added. A prior is plotted if
#'   either a \code{\link[coda]{mcmc}} object or \code{list} with the prior 
#'   density is supplied. See Details.
#' @param prange Character argument specifying the relative range
#'   \code{c("prior", "posterior")} for the density estimation of the prior.
#'   The supplied argument must begin with either \code{prior} or
#'   \code{posterior} exactly, but can end with a numeric value indicating the
#'   relative range of the prior posterior to be used. See Details. Range
#'   adjustment currently only works when class \code{mcmc} provided to
#'   \code{prior}.#FIXME
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
#' @param denscol,hpdcol,meancol,priorcol,histcol The line colors for the lines
#'   to be plotted representing the kernel density estimate (\code{denscol}),
#'   credible interval limits (\code{hpdcol}), posterior mean (\code{meancol}),
#'   prior density (\code{priorcol}), or histogram fill (\code{histcol}). See
#'   \code{\link[graphics]{par}} for details.
#' @param at1,at2 The points at which tick-marks are to be drawn on the
#'   x-axis (\code{at1}) and y-axis (\code{at2}). See
#'   \code{\link[graphics]{axis}} for details.
#' @param labels1,labels2 Labels to be placed at the tick-marks of the
#'   x-axis (\code{labels1}) and and y-axis (\code{labels2}). See
#'   \code{\link[graphics]{axis}} for details.
#' @param plot Logical. If \code{TRUE} (default) the plotting is performed,
#'   otherwise the graphical output is suppressed. 
#' @param \dots Additional arguments passed to \code{densplot}, \code{hist},
#'   \emph{and} \code{lines}.
#'
#' @return A \code{list} returned invisibly (assign to an object to view):
#'   \describe{
#'     \item{call}{The original \code{call}.}
#'     \item{postDensity}{List containing objects describing the posterior
#'	 distribution. Items include a \code{character} indicating if the
#'	 density was adjusted because of boundary constraints, an object of
#'	 \code{class} "density", and an object of \code{class} "histogram".}
#'     \item{priorDensity}{List containing objects describing the prior
#'	 distribution. Items include a \code{character} indicating if the
#'	 density was adjusted because of boundary constraints and an object of
#' 	 \code{class} "density".}
#'   }
#'
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link[coda]{densplot}}, \code{\link[graphics]{hist}}
#' @family MCMC posterior distribution helper functions
#' @examples
#' # Simulate two example MCMC chains:
#' ## Both are standard normal distributions of only positive values
#' ## However, one is on the boundary (near zero) while the other is not
#' normMCMC <- coda::mcmc(matrix(c(abs(rnorm(1000, 0.2, 1)), rnorm(1000, 10, 1)),
#'	ncol = 2))
#' par(mfrow = c(2, 1))
#'   postPlot(normMCMC[, 1], ylim = c(0, 1))
#'   postPlot(normMCMC[, 2], ylim = c(0, 1))
#'
#' # Now include an inverse Gamma prior, using the density function
#' postPlot(normMCMC[, 1], ylim = c(0,1),
#'   prior = dIW(seq(1e-16, 5, length = 1000), V = diag(1), nu = 1))
#'
#' # Example with parameter expanded prior:
#' pepr <- coda::mcmc(rpeIW(1000, V = diag(2), nu = 2, alpha.mu = rep(0, 2),
#'   alpha.V = diag(2)*1000))
#' par(mfrow = c(2, 1))
#'   postPlot(normMCMC[, 1], ylim = c(0, 1),
#'     prior = pepr[,1], prange = "posterior6")
#'   postPlot(normMCMC[, 2], ylim = c(0, 1),
#'     prior = pepr[,4], prange = "posterior1.5") 
#' # different scalar to get same y-intercept
#'
postPlot <- function(posterior, bw = "nrd", #TODO make separate prior/posterior bws?
	plotHist = TRUE, histbreaks = 100, histcol = NULL,
	prior = NULL, prange = c("prior", "posterior"),
	main, sub, ylim,
	denslwd = 6, denslty = "solid", denscol = "black",
	hpdlwd = 6, hpdlty = "dashed", hpdcol = "black",
	meanlwd = 7, meanlty = "12", meancol = "red",
	priorlwd = 3, priorlty = "solid", priorcol = "green",
	at1 = NULL, at2 = NULL,
	labels1 = NULL, labels2 = NULL,
	plot = TRUE, ...){

  cl <- match.call()
  # bw.nrd() is like `bwf()` in `coda::densplot()` 
  poDens <- stats::density(posterior, bw = bw, kernel = "gaussian", n = 2^13)
  bw <- poDens$bw
  main.title <- if(missing(main)) "" else main
  sub.title <- if(missing(sub)) "" else sub
  histout <- if(plotHist){
               graphics::hist(posterior, breaks = histbreaks, plot = FALSE)
             } else list(breaks = NULL, counts = NULL, density = NULL, mids = NULL)
  maxYhistpoDens <- max(c(poDens$y, histout$density))
  if(missing(ylim)){
    ylimit <- c(0, maxYhistpoDens)
  } else ylimit <- ylim

  # Adjust density at boundaries
  ## reflect density estimates back across a boundary
  constraint <- "unbounded"
  po <- posterior
  if(max(po) <=1 && 1 - max(po) < 2 * bw){
    if(min(po) >= 0 && min(po) < 2 * bw){
      constraint <- "zero-one"
      po <- c(po, -po, 2 - po)
    } else{
        if(min(po) < 0 && 1 + min(po) < 2 * bw){
          if(min(po) >= -1){  #<-- otherwise 'unbounded' happens to have max near 1
            constraint <- "correlation"
            po <- c(po, -2 - po, 2 - po)
          }
        } else{
            constraint <- "upbound1"  #<-- most likely correlation converge on 1
            po <- c(po, 2-po)
          }
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

  if(plot){
    plot(poDens, type = "n", axes = FALSE,
	main = main.title, sub = sub.title, ylim = ylimit, ...)
    graphics::lines(poDens, lwd = denslwd, col = denscol, ...)
    if(plotHist) graphics::hist(posterior, breaks = histbreaks, col = histcol,
	freq = FALSE,
	add = TRUE, ...)
    yaxmax <- max(c(at2, par("yaxp")[2]))
    lineylim <- ifelse(yaxmax >= maxYhistpoDens, yaxmax, maxYhistpoDens)
    nout <- sapply(coda::HPDinterval(posterior), FUN = function(x){(
	graphics::lines(x = rep(x, 2), y = c(ylimit[1], lineylim),
	  lty = hpdlty, lwd = hpdlwd, col = hpdcol, ...)
      )})
    graphics::lines(x = rep(mean(posterior), 2), y = c(ylimit[1], lineylim),
	col = meancol, lty = meanlty, lwd = meanlwd, ...)
    axis(1, at = at1, labels = labels1)
    axis(2, at = at2, labels = labels2)
  }

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
      if(missing(prange)) prange <- "prior"
      if(!startsWith(prange, "prior") && !startsWith(prange, "posterior")){
        stop("argument to 'prange' must start with the full word of either 'prior' or 'posterior'")
      }
      if(startsWith(prange, "prior")){
        prraN <- as.numeric(strsplit(prange, "prior")[[1L]][2L])
        prra <- range(prior)
      }
      if(startsWith(prange, "posterior")){
        prraN <- as.numeric(strsplit(prange, "posterior")[[1L]][2L])
        prra <- range(posterior)
      }
      if(is.na(prraN)) prraN <- 1
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
## Density is high near zero, so usually above posterior plot and off figure
        #XXX NOTE if problem occurs, consider calculate own `bw` and not use
        ## `width` with posterior bw (see commented code below)
        if(max(prior) > prra[2L]*prraN){
          pr <- prior[prior >= 0 & prior <= prra[2L]*prraN]
          pr <- c(pr, -pr, 2*prra[2L]*prraN - pr)
#          prDens <- stats::density(pr, bw = "nrd", kernel = "gaussian", n = 2^13)
          prDens <- stats::density(pr, kernel = "gaussian", width = 4 * bw, n = 2^13)
          prDens$y <- 3 * prDens$y[prDens$x >= 0 & prDens$x <= prra[2L]*prraN]
          prDens$x <- prDens$x[prDens$x >= 0 & prDens$x <= prra[2L]*prraN]
        } else{
            pr <- prior[prior >= 0]
            pr <- c(pr, -pr)
#            prDens <- stats::density(pr, bw = "nrd", kernel = "gaussian", n = 2^13)
            prDens <- stats::density(pr, kernel = "gaussian", width = 4 * bw, n = 2^13)
            prDens$y <- 2 * prDens$y[prDens$x >= 0]
            prDens$x <- prDens$x[prDens$x >= 0]
          }
      }
      if(constraint == "unbounded"){
        #XXX NOTE if problem occurs, consider calculate own `bw` and not use
        ## `width` with posterior bw (see commented code below)
        if(prra[1L] < 0) prra[1L] <- prra[1L] * prraN
        if(prra[1L] >= 0) prra[1L] <- prra[1L] * (1/prraN)
        if(prra[2L] < 0) prra[2L] <- prra[2L] * (1/prraN)
        if(prra[2L] >= 0) prra[2L] <- prra[2L] * prraN
        pr <- prior[prior >= prra[1L] & prior <= prra[2L]]
        pr <- c(pr, 2*prra[1L] - pr, 2*prra[2L] - pr)
#        prDens <- stats::density(pr, bw = "nrd", kernel = "gaussian", n = 2^13)
        prDens <- stats::density(pr, kernel = "gaussian", width = 4 * bw, n = 2^13)
        prDens$y <- 3 * prDens$y[prDens$x >= prra[1L] & prDens$x <= prra[2L]]
        prDens$x <- prDens$x[prDens$x >= prra[1L] & prDens$x <= prra[2L]]
      }
    }

    if(!coda::is.mcmc(prior) && !is.list(prior)){
      warning("prior is not a `list` or `mcmc` object - no prior added to plot")
    } else{
    # Add prior density to plot
        if(plot){
	  # to keep asymptotically increasing priors from 'running away'
          ##TODO check the following restriction for lots of different priors
          prDensSub <- which(prDens$y <= lineylim)
          graphics::lines(prDens$y[prDensSub] ~ prDens$x[prDensSub],
	    col = priorcol, lty = priorlty, lwd = priorlwd, ...)
        }
      }
  } else prDens <- NULL

 return(invisible(list(call = cl,
	postDensity = list(constraint = constraint,
		density = poDens,
		histogram = histout),
	priorDensity = list(constraint = constraint,
		density = prDens)
	)))
	
}











################################################################################
#' Two MCMC chain plot.
#'
#' Pretty plot of two MCMC chains side-by-side.
#'
#' @section Warning:
#' The function attempts to match parameters by the column names, but may have
#' trouble where vastly different objects are in each of the two MCMC objects.

#' @aliases plot2mcmc
#' @export
#' @param x1 Object of class \code{mcmc} containing MCMC chain(s).
#' @param x2 Optional object of class \code{mcmc} containing a second set of
#'   MCMC chain(s).
#'
#' @param smooth Logical. See \code{\link[coda]{traceplot}} for details.#FIXME
#' @param bwf Character indicating function to calculate the bandwith. See
#'   \code{\link[coda]{traceplot}} for details.#TODO correspond with `postPlot()`
#'
#' @param save Either a logical or character value. When \code{FALSE} (default)
#'   plots are not saved to the hard drive. If \code{TRUE}, the plots are saved
#'   as a pdf with a default filename to the current working directory. If a
#'   \code{character} is specified, this will be used as the base filename to
#'   save the plots as pdfs in the current working directory.
#' @param \dots Additional arguments passed to \code{densplot} and
#'   \code{traceplot}.#FIXME update densplot and traceplot to own functions
#'
#' @return Nothing returned invisibly at the moment (assign to keep):#FIXME
#'
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{postPlot} and \code{\link[coda]{densplot}}
#' @family MCMC posterior distribution helper functions
#' @examples
#' # Simulate two sets of two example MCMC chains:
#' ## Both are standard normal distributions of only positive values
#' ## However, one is on the boundary (near zero) while the other is not
#' normMCMC <- coda::mcmc(matrix(c(abs(rnorm(1000, 0.2, 1)), rnorm(1000, 10, 1)),
#'	ncol = 2))
#' normMCMC2 <- coda::mcmc(matrix(c(abs(rnorm(1000, 0.2, 1)), rnorm(1000, 10, 1)),
#'	ncol = 2))
#' plot2mcmc(normMCMC, normMCMC2)
plot2mcmc <- function(x1, x2 = NULL, smooth = FALSE, bwf, save = FALSE, ...){
    oldpar <- NULL
    on.exit(par(oldpar))
    if(!is.logical(save) && is.character(save)){
      fname <- save
      save <- TRUE
    } else fname <- "set"
    if(is.null(x2)){
      mfrow <- coda:::set.mfrow(Nchains = nchain(x1), Nparms = coda::nvar(x1), 
          nplots = 2)
      for (i in 1:coda::nvar(x1)) {
        if(i %in% seq(from = 1, to = coda::nvar(x1), by = mfrow[1])){
          x11(w = 7, h = 9)
          par(mfrow = mfrow)
          cnt <- 0
        }

        y1 <- coda::mcmc(as.matrix(x1)[, i, drop = FALSE], start(x1), 
          end(x1), thin(x1))
#TODO add own version of `traceplot()` (`plot` + `rug`)
        coda::traceplot(y1, smooth = smooth, ...)
        if(length(unique(y1)) == 1){
          plot.new()
        } else{
          if (missing(bwf)) coda::densplot(y1, ...) else coda::densplot(y1, bwf = bwf, ...)
          }
        cnt <- cnt + 1
        if(save && (cnt == mfrow[1] | i == coda::nvar(x1))){
          dev.copy(pdf, paste0(fname, "_", i-(cnt-1), "to", i, ".pdf"), w = 7, h = 9)
          dev.off()
          cat(paste0("Plots saved to: ", getwd(), "/", fname, "_", i-(cnt-1), "to", i, ".pdf"), "\n")
        } 
      }

    } else{   
        nvar_x1x2 <- max(c(coda::nvar(x1), coda::nvar(x2)))
        largerX <- which(nvar_x1x2 == c(coda::nvar(x1), coda::nvar(x2)))[1]
        mfrow <- coda:::set.mfrow(Nchains = nchain(x1), Nparms = nvar_x1x2, 
            nplots = 2)
        mfrow[2] <- 4

        for(i in 1:nvar_x1x2){
          if(i %in% seq(from = 1, to = nvar_x1x2, by = mfrow[1])){
             x11(w = 13, h = 9)
              par(mfrow = mfrow)
              cnt <- 0
          }

          if(coda::nvar(x1) == coda::nvar(x2)){
            y1 <- coda::mcmc(as.matrix(x1)[, i, drop = FALSE], start(x1), 
              end(x1), thin(x1))
            y2 <- coda::mcmc(as.matrix(x2)[, i, drop = FALSE], start(x2), end(x2), thin(x2))
            if(length(unique(y1)) != 1 & length(unique(y2)) != 1){
              denslist_y1y2 <- postPlot(y1, plotHist = FALSE,
		prior = y2, prange = "prior", plot = FALSE)
	      xlimit <- with(denslist_y1y2, range(c(postDensity$density$x,
		priorDensity$density$x)))
#	      ylimit <- with(denslist_y1y2, range(c(postDensity$density$y,
#		priorDensity$density$y)))
            }
            coda::traceplot(y1, smooth = smooth, ...)
            if(length(unique(y1)) == 1){
              plot.new()
            } else{
                if(missing(bwf)){
	          coda::densplot(y1, xlim = xlimit, ...)
                } else coda::densplot(y1, bwf = bwf, xlim = xlimit, ...)
              }
            coda::traceplot(y2, smooth = smooth, ...)
            if(length(unique(y1)) == 1){
              plot.new()
            } else{
                if(missing(bwf)){
	          coda::densplot(y2, xlim = xlimit, ...)
	        } else coda::densplot(y2, bwf = bwf, xlim = xlimit, ...)
              }
          } else{
              if(largerX == 1){
                y1 <- coda::mcmc(as.matrix(x1)[, i, drop = FALSE],
	          start(x1), end(x1), thin(x1))
                if(colnames(x1)[i] %in% colnames(x2)){
                  j <- match(colnames(x1)[i], colnames(x2))
                  y2 <- coda::mcmc(as.matrix(x2)[, j, drop = FALSE],
		    start(x2), end(x2), thin(x2))

                  if(length(unique(y1)) != 1 & length(unique(y2)) != 1){
                    denslist_y1y2 <- postPlot(y1, plotHist = FALSE,
		      prior = y2, prange = "prior", plot = FALSE)
	            xlimit <- with(denslist_y1y2, range(c(postDensity$density$x,
		      priorDensity$density$x)))
#	            ylimit <- with(denslist_y1y2, range(c(postDensity$density$y,
#		      priorDensity$density$y)))
                  }
		  #TODO use own version
                  coda::traceplot(y1, smooth = smooth, ...)
                  if(length(unique(y1)) == 1){
                    plot.new()
                  } else{
                      if(missing(bwf)){
	                coda::densplot(y1, xlim = xlimit, ...)
                      } else coda::densplot(y1, bwf = bwf, xlim = xlimit, ...)
                    }
                  coda::traceplot(y2, smooth = smooth, ...)
                  if(length(unique(y2)) == 1){
                    plot.new()
                  } else{
                      if(missing(bwf)){
	                coda::densplot(y2, xlim = xlimit, ...)
	              } else coda::densplot(y2, bwf = bwf, xlim = xlimit, ...)
                    }
                } else{
                    coda::traceplot(y1, smooth = smooth, ...)
                    if(length(unique(y1)) == 1){
                      plot.new()
                    } else{
                        if(missing(bwf)){
	                  coda::densplot(y1, ...)
                        } else coda::densplot(y1, bwf = bwf, ...)
                      }
                    plot.new()
                    plot.new()
                  }

              } else{
                  y2 <- coda::mcmc(as.matrix(x2)[, i, drop = FALSE],
		    start(x2), end(x2), thin(x2))
                  if(colnames(x2)[i] %in% colnames(x1)){
                    j <- match(colnames(x2)[i], colnames(x1))
                    y1 <- coda::mcmc(as.matrix(x1)[, j, drop = FALSE],
		      start(x1), end(x1), thin(x1))
 
                    if(length(unique(y1)) != 1 & length(unique(y2)) != 1){
                      denslist_y1y2 <- postPlot(y1, plotHist = FALSE,
		        prior = y2, prange = "prior", plot = FALSE)
	              xlimit <- with(denslist_y1y2, range(c(postDensity$density$x,
		        priorDensity$density$x)))
#	              ylimit <- with(denslist_y1y2, range(c(postDensity$density$y,
#		        priorDensity$density$y)))
                    }
                    coda::traceplot(y1, smooth = smooth, ...)
                    if(length(unique(y1)) == 1){
                      plot.new()
                    } else{
                        if(missing(bwf)){
		          coda::densplot(y1, xlim = xlimit, ...)
		        } else coda::densplot(y1, bwf = bwf, xlim = xlimit, ...)
                      }
                    coda::traceplot(y2, smooth = smooth, ...)
                    if(length(unique(y2)) == 1){
                      plot.new()
                    } else{
                        if(missing(bwf)){
		          coda::densplot(y2, xlim = xlimit, ...)
		        } else coda::densplot(y2, bwf = bwf, xlim = xlimit, ...)
                      }
                  } else{
                      plot.new()
                      plot.new()
                      coda::traceplot(y2, smooth = smooth, ...)
                      if(length(unique(y2)) == 1){
                        plot.new()
                      } else{
                          if(missing(bwf)){
		            coda::densplot(y2, ...)
		          } else coda::densplot(y2, bwf = bwf, ...)
                        }
                    }  
                }
            }
          cnt <- cnt + 1
          if(save && (cnt == mfrow[1] | i == coda::nvar(x1))){
            dev.copy(pdf, paste0(fname, "_", i-(cnt-1), "to", i, ".pdf"), w = 13, h = 9)
            dev.off()
            cat(paste0("Plots saved to: ", getwd(), "/", fname, "_", i-(cnt-1), "to", i, ".pdf"), "\n")
          }

        }
      }
}

