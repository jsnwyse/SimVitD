\name{SimVitD-package}
\alias{SimVitD-package}
\docType{package}
\title{Simulation tools for vitamin D studies.}
\description{
A suite of tools for statistical planning of vitamin D trials. 
}

\details{ 

\tabular{ll}{
Package: \tab SimVitD\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2020-12-22\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}

The \code{SimVitD} package uses simulation to aid in statistical planning of vitamin D trials. The core elements of the simulation based tools involve 
  
  (i) simulation of an individual vitamin D staus profiles taking into account natural cyclic variation in summer and winter
  
  (ii) simulation of exposures to potential infections in an individual's envrionment
  
  (iii) simulation of the possibility of developing infection conditional on vitamin D status at exposure
  
  (iv) Monte Carlo estimation of power of detecting supplementation effects for a selection of different supplementation scheme and two classes of hypothesis test.
  
The fuctions \code{\link{vitd.curve}} to generate individual status profiles and \code{\link{exposure.levels}} to simulate exposures. The accompanying vignette provides a gateway into the simulation tools. 
}

\examples{
# simulate individual Vitamin D profiles
indprofiles <- vitd.curve( n=10, type="placebo" )

# make a plot of all these profiles
plot( indprofiles )

# intensity function for exposures to infection
intensfun <- intensity.function( summer.rate=0.1, winter.rate=0.9 )

# exposure times
expos <- exposure.levels( indprofiles, rate=2, intensfun, end=2 )

# plot of exposures on top of Vitamin D profiles
plot( expos )

# disease calculation
infect <- infection.count( expos )

# plot disease points on top of exposure points
infection.count.plot( expos, infect )
}

\author{Rebecca Mangan, Jason Wyse, Lina Zgaga}
