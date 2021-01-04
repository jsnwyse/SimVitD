\name{infection.count.plot}
\alias{infection.count.plot}
\title{Plotting the times of disease on the vitamin D curves.}
\description{
  Plot whether a participant was disease or not as points on top of the vitamin D curves plotted from \code{plot.vitd.curve}.
  }
  \usage{
    infection.count.plot( expos, infect, pch = 20, cex = 1.5, col = "red" )
  }
  \arguments{
    \item{expos}{An object of class \code{exposure.levels}.}
    \item{infect}{An object of class \code{infection.count}.}
    \item{pch}{Plotting symbol for the points.}
    \item{cex}{Standard graphical parameter.}
    \item{col}{Colour of the plotted points.}
  }

\author{Rebecca Mangan and Jason Wyse}