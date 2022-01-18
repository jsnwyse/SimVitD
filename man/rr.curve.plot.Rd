\name{rr.curve.plot}
\alias{rr.curve.plot}
\title{Plotting the relative risk curve for vitamin D status.}
\description{
  Plot the relative risk curve for vitamin D showing times of exposure and whether a participant developed infection from exposure.}
  \usage{
    rr.curve.plot( expos, infect, idx = 1, main = NULL, xlab = "25-hydroxyvitamin D",
        ylab = "Risk scaling", col = "blue", \dots )
  }
  \arguments{
    \item{expos}{An object of class \code{exposure.levels}.}
    \item{infect}{An object of class \code{infection.count}.}
    \item{idx}{A vector of indexes of specific exposures to plot.}
    \item{main}{Main title for the plot.}
    \item{xlab}{A title for the x-axis.}
    \item{ylab}{A title for the y-axis.}
    \item{col}{Character, value for the colour of the points.}
    \item{\dots}{Additional arguments to \code{plot}.}
  }

\author{Rebecca Mangan and Jason Wyse}
