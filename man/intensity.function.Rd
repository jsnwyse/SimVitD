\name{intensity.function}
\alias{intensity.function}
\title{Generate an intensity function}
\description{
  Generates intensity function, function with time as sole argument.
}
\usage{
  intensity.function( summer.rate = 0, winter.rate = 1, flu = TRUE )
}
\arguments{
  \item{summer.rate}{Value bwteeen 0 and 1, the rate at which exposures occur in summer months.}

  \item{winter.rate}{Value between 0 and 1, the rate at which exposures occur in winter months.}

  \item{flu}{If \code{TRUE} define winter months as flu season (Sep-May in Northern Hemisphere, Mar-Nov in Southern Hemisphere). If \code{FALSE} have constant exposure rate all year.}

}

\value{

  \code{\link{intensity.function}} returns a function that takes
          time as sole argument and returns value between 0 and 1,
          the returned function acts as
          \code{intensity.func} for input to \code{\link{exposure.levels}}.
}


\author{Rebecca Mangan and Jason Wyse}
