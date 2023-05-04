#' List of all color palettes
#'
#' Run the function \code{\link{wolak_palette}} to create the palette objects.
#'
#' \dQuote{aucfriend} is a color palette that is accessible across all types of
#'   color vision and was checked with the package \code{colorblindcheck} that
#'   implements a model of human color vision.
#' @author \email{matthewwolak@@gmail.com}
#' @export
wolak_palettes <- list(
  aucfriend = c(blue = "#03244d", cyan = "#00A08A", orange = "#e86823")
)


#' Create wolakR package palette
#'
#' Sets up the R object of class palette that contains one of the pre-specified
#'   color palettes.
#'
#' Function is based off the \code{wesanderson::wes_palette} function.
#' @param name Character specifying name of palette to use. See
#'   \code{\link{wolak_palettes}} for list of choices.
#' @param n Integer number of colors to return. 
#' @param type Character indicating either \dQuote{discrete} for colors already
#'   listed in the paletter or \dQuote{continuous} to interpolate between colors.
#'
#' @return A vector of colors of the class \dQuote{palette}.
#' @author \email{matthewwolak@@gmail.com}
#' @export
#' @examples
#' clr <- wolak_palette("aucfriend")
#' \dontrun{
#'   # check:
#'     colorblindcheck::palette_check(clr, plot = TRUE)
#' }
wolak_palette <- function(name = "aucfriend", n,
	type = c("discrete", "continuous")){
  type <- match.arg(type)
  pal <- wolak_palettes[[name]]
  if(is.null(pal)){
    stop("Palette name not found in wolak_palettes")
  }
  
  if(missing(n)) n <- length(pal)
  if(type == "discrete" & n > length(pal)){
    stop("Number of colors `n` is greater than total number in requested palette. Switch palettes or use `type = continuous`")  
  }
  
  pout <- switch(type,
    discrete = pal[1:n],
    continuous = grDevices::colorRampPalette(pal)(n))
  
 structure(pout, class = "palette", name = name)
}


#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link{wolak_palette}}
#' @export
#' @method plot palette
plot.palette <- function(pal, ...){
  n <- length(pal)
  bplotOut <- barplot(rep(10, n), col = pal, axes = FALSE)
    axis(1, at = bplotOut, labels = names(pal), lty = "blank")   	
}
