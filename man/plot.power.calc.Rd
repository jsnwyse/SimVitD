\name{plot.power.calc}
\alias{plot.power.calc}
\title{Plotting a power.calc object}
\description{
  Plot the power between two groups as the number of participants and the relative risk vary.
  }
  \usage{
    \method{plot}{power.calc}( x, col = "hotpink", lwd = 1.5, lty = 1, ylab = NULL, 
          x.legend = NULL, y.legend = NULL, main.legend = "Risk scaling", legend.size = 1, 
          target.power = NA, which = 1L, \dots )
  }
  \arguments{
    \item{x}{An object of class \code{power.calc}.}
    \item{col}{Colour of the line.}
    \item{lwd}{Line width.}
    \item{lty}{Line type.}
    \item{ylab}{Title for the y-axis.}
    \item{x.legend}{The x position of the legend. If not given, tries to default to a sensible value.}
    \item{y.legend}{The y position of the legend. If not given, tries to default to a sensible value.}
    \item{main.legend}{Title for the legend.}
    \item{legend.size}{Size of the legend.}
    \item{target.power}{Add a horizontal dotted line at a target power. A value of \code{NA} suppresses plotting.}
    \item{which}{Type of plot. Value 1 gives a plot of power, and value 2 plot of estimated effect size.}
    \item{\dots}{Additional arguments to be passed to the \code{plot} function.}
  }

\author{Rebecca Mangan and Jason Wyse}
