\name{vitd.curve}
\alias{vitd.curve}
\title{Generate Vitamin D curves}
\description{
  Generates a vitamin D status profile curve for each individual in a group
}
\usage{
vitd.curve( n = 1, type = c("placebo","fixed-dose","dynamic-dose"), start = 0, end = 1,  
    mu = 45, amplitude = 35, dyn.dose.thresh = 50, sd.mu = 5, sd.amplitude = 5, 
    sd.dyn.dose.thresh = 5, supp.dose = 20, supp.dose.rate = Inf, weight = 1, 
    sd.weight = 0, min.thresh = 10, north.hemi = TRUE, res = 40 )
}
\arguments{
  \item{n}{Number of curves to simulate.}

  \item{type}{One of "placebo", "fixed-dose", "dynamic-dose".}

  \item{start}{Time in units of years when trial started.}

  \item{end}{Time in units of years when trial ended.}

  \item{mu}{The mean level of 25OHD in the trial arm around which there is cosine variation.}

  \item{amplitude}{Amplitude of cosine function describing variation around \code{mu}}

  \item{dyn.dose.thresh}{Threshold for the concentration-controlled scheme.}

  \item{sd.mu}{Standard deviation levels around mean; this is the standard deviation of H in vignette.}

  \item{sd.amplitude}{Standard deviation of the amplitude.}

  \item{sd.dyn.dose.thresh}{Standard deviation of the participant retained concentration in the concentration-controlled trial.}

  \item{supp.dose}{The 25OHD nmol/l equivalent for dosage in fixed-dose supplementation.}
  
  \item{supp.dose.rate}{Concentration parameter for fixed-dose scheme uptake. Large values imply that all participants derive the same equivalent.}
  
  \item{weight}{For fixed-dose supplementation, this is the expected value of the proportion of the dose which is always utilized.}
  
  \item{sd.weight}{Standard deviation of weight}
  
  \item{min.thresh}{The minimum detectable threshold of 25OHD. Defaults to 10 nmol/l.}
  
  \item{north.hemi}{Summer/winter months as in the Northern Hemisphere if \code{TRUE} and Southern Hemisphere if \code{FALSE}.}
  
  \item{res}{Resolution parameter for plotting of vitamin D curves.}
}

\value{
  \code{\link{vitd.curve}} returns an object of class \code{vitd.curve}. Curve parameter settings are returned. The slot \code{curves} give the curve specific parameters for each of the \code{n} generated curves. Additionally, \code{time} used for plotting gives time values passed to \code{plot.vitd.curve}.
}

\author{Rebecca Mangan and Jason Wyse}


