\name{rr.profile.plot}
\alias{rr.profile.plot}
\title{Plotting a vitamin D profile and relative risk curve.}
\description{
  Plot a vitamin D status profile for a single participant and the relative risk curve for vitamin D (with exposure times and whether a participant was infected at that exposure time) side by side.}
\usage{
rr.profile.plot( x, expos, infect, idx = 1, \dots )
}
\arguments{
  \item{x}{An object of class \code{vitd.curve}.}
  \item{expos}{An object of class \code{exposure.levels}.}
  \item{infect}{An object of class \code{infection.count}.}
  \item{idx}{Index of curve to plot.}
  \item{\dots}{Additional arguments to \code{plot}.}
}

\author{Rebecca Mangan and Jason Wyse}

\examples{
# individual profiles
indprofiles <- vitd.curve( n=10, type="placebo" )

# intensity function for exposures to infection
intensfun <- intensity.function( summer.rate=0.1, winter.rate=0.9 )

# exposure times
expos <- exposure.levels( indprofiles, rate=2, intensfun, end=2 )

# disease calculation
infect <- infection.count( expos )

# plot RR profile
rr.profile.plot( indprofiles, expos, infect )
}
