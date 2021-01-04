\name{infection.count}
\alias{infection.count}
\title{Generate whether a participant contracted a infection at exposure times.}
\description{
  Generates probabilities of developing infection, the relative risk and whether a participant becomes infected or not at exposure times.
}
\usage{
  infection.count( expos, baseline = 0.03, RR = 3, holding.time = 2, lohi.vit = c(10,70) )
}
\arguments{
  \item{expos}{An object of class \code{exposure.levels}.}

  \item{baseline}{Numeric, baseline prevalence of developing infection at any exposure time. Probability of developing infection when a participant has maximum (fully replete) vitamin D levels.}

  \item{RR}{Numeric, the relative risk of the study. The difference between becoming infected at highest and lowest vitamin D levels.}

  \item{holding.time}{Numeric, the mean number of weeks for after becoming infected when the participant can not become re-infected. Times are simulated independently from an exponential distribution with this mean.}
  
  \item{lohi.vit}{A vector of length 2 giving the level at which a participant is considered to have insufficient and sufficient vitamin D status levels.}
}

\value{

  \code{\link{infection.count}} returns an object of class \code{infection.count} that is a list.
  The list has the following slots.

  \item{baseline}{The baseline prevalence of disease.}

  \item{RR}{The relative risk of the study.}
  
  \item{inflection}{The inflection points of the relative risk curve.}

  \item{probs}{Matrix, the probability of developing infection at each exposure time for each participant.}

  \item{relativerisk}{Matrix, the relative risk of infection at each exposure time for each participant.}

  \item{infection}{Matrix, 1 if participant was infected at the corresponding exposure time, 0 if participant was not infected at exposure time.}

  \item{count}{Number of infections for each participant over the course of the study.}

  \item{mean}{The mean number of infections for the group over the course of the study.}

}

\author{Rebecca Mangan and Jason Wyse}

