\name{exposure.levels}

\alias{exposure.levels}

\title{Generate exposures to infectious agent and vitamin D status}

\description{
  \code{\link{exposure.levels}} creates a list of the exposures to an infectious agent and the status of vitamin D at the time of each exposure for each participant. The exposure times are simulated from a non homogeneous poisson process.
}

\usage{
  exposure.levels( x, rate, intensity.func = intensity.function(), end = 1 )
}

\arguments{
  \item{x}{An object of class \code{vitd.curve}.}
  \item{rate}{Expected number of exposures per week, the rate at which exposures occur in an equivalent homogeneous Poisson process.}
  \item{intensity.func}{A function that takes time as sole argument and returns value between 0 and 1, scaling the rate argument. For example, exposures may be higher at certain times for seasonal infections.}
  \item{end}{End time of study in units of years.}
  }

\value{

  \code{\link{exposure.levels}} returns an object of class \code{exposure.levels} that is a list.
  The list has the following slots.

  \item{exposures}{The exposure times to infection for each participant. }
  \item{levels}{The vitamin D staus levels for each participant at the time of exposure to infection.}
  }
  
\author{Rebecca Mangan and Jason Wyse}
