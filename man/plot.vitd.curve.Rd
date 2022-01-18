\name{plot.vitd.curve}
\alias{plot.vitd.curve}
\title{Plotting a vitd.curve object}
\description{
  Plot vitamin D status curves.
  }
\usage{
  \method{plot}{vitd.curve}( x, main = " ", xlab = " ", ylab = "25-hydroxyvitamin D", 
    col=1:6, add = FALSE, ylim = NULL, \dots )
}
\arguments{
  \item{x}{An object of class \code{vitd.curve}.}
  \item{main}{Main title for the plot.}
  \item{xlab}{A title for the x-axis.}
  \item{ylab}{A title for the y-axis.}
  \item{col}{A vector of colours for the plotted status curves.}
  \item{add}{If \code{TRUE}, add to current plotting window}
  \item{ylim}{Gives the plotting limits on tye y-axis.}
  \item{\dots}{Additional arguments to the \code{plot} function.}
}

\author{Rebecca Mangan and Jason Wyse}
