\name{SDI}
\alias{SDI}
\title{
Sexual Dimorphism Index
}
\description{
Calculates the sexual dimorphism index of Lovich & Gibbons (1992): \code{(larger sex / smaller sex) - 1} and arbitrarily made positive when females are the larger sex or negative when males are the larger sex.
}
\usage{
SDI(t, dfin)
}
\arguments{
  \item{t }{
Character object naming a trait in \code{dfin}.
  }
  \item{dfin }{
Data frame which must contain at least three columns: one containing trait names that match \code{t}, a column with female measurements that must match "female" exactly, and a column with male measurements that must match "male" exactly.
  }
}
\value{
A numeric representing the calculated sexual dimorphism index for trait \code{t}.
}
\references{
Lovich, JE & Gibbons, JW. 1992. A review of techniques for quantifying sexual size dimorphism. Growth, Development, and Ageing 56(4):269-281.
}
\author{\email{matthewwolak@generic.email.com}}
\examples{
popDat <- data.frame(trait = c("length", "weight", "coloration", "intelligence"),
    female = c(4.2, 255.0, 0.57, 144),
    male = c(3.7, 255.0, 0.89, 143))
SDI("length", popDat)
sapply(popDat$trait, FUN = SDI, dfin = popDat)
stopifnot(SDI("weight", popDat) == 0.0)
}
