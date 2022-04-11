\name{power.calc}
\alias{power.calc}
\title{Generate the power to detect the difference between two groups}
\description{
  generates a value for the power between two groups
}
\usage{
  power.calc( n, ratio=1, N = 500, test.type, sig.level = 0.05,
    vitdcurves.placebo = NULL, vitdcurves.treatment = NULL,
  	baseline = 0.03, RR = 3, rate = 1, intensity.func = intensity.function(), 
  	holding.time = 2, lohi.vit = c(10,70), clt = NULL, mc.error = 1, boot.rep = 500, 
  	parallel = FALSE, num.cores = NULL, verbose=FALSE )
}
\arguments{
  \item{n}{Numeric, the number of participants in the control group.}
  
  \item{ratio}{Ratio greater than or equal to 1 giving size of treatment group as a multiple of \code{n}. Defaults to groups of equal size.}

  \item{N}{Number of simulations of the study to run.}

  \item{test.type}{Type of test to calculate the power, one of "count" or "proportions". }

  \item{sig.level}{Significance level used to test for a statistically significant difference between the groups.}
  
  \item{vitdcurves.placebo}{An object of class \code{vitd.curve}.}

  \item{vitdcurves.treatment}{An object of class \code{vitd.curve}.}

  \item{baseline}{Baseline prevalence of getting diseased at any exposure time. Probability of getting diseased when a participant has sufficient vitamin D levels.}

  \item{RR}{Fold risk difference between getting infection between the most deficient and most sufficient vitamin D levels.}

  \item{rate}{Expected number of exposures per week, the rate at which exposures occur in the equivalent homogeneous Poisson process.}

  \item{intensity.func}{Function taking time as sole argument and returns value between 0 and 1, input to nhpp function see poisson package.}

  \item{holding.time}{Expected number of weeks for the holding time.}
  
  \item{lohi.vit}{Inflection points of the relative risk curve used in \code{infection.count}.}
  
  \item{clt}{Vector or single logical variable indicating whether Central Limit Theorem should be used to approximate the p-value.}
  
  \item{mc.error}{Number of times to repeat the experiment at each \code{n} value to explore Monte Carlo error.}
  
  \item{boot.rep}{Number of bootstrap samples to carry out non-parametric tests of hypotheses.}
  
  \item{parallel}{Use parallel processing to carry out the simulations. This will parallelise over mc.error.}
  
  \item{num.cores}{Number of cores to exploit in parallel mode. Defaults to (cores available) - 1. }
  
  \item{verbose}{If \code{TRUE} print out regular updates on progression. Not available when parallel is set to \code{TRUE}.}
}

\value{

   \code{\link{power.calc}} returns an object of class \code{power.calc} that is a list.
  The list has the following slots.

  \item{test.type}{The type of study the power has been calculated on.}

  \item{baseline}{The baseline prevalence for disease.}

  \item{RR}{The relative risk of the study.}

  \item{npergroup}{The number of participants per group in the study.}
  
  \item{mc.error}{Number of repetitions of experiment to approximate Monte Carlo error.}

  \item{power}{A \code{length(RR)} x \code{mc.error} x \code{length(npergroup)} array of estimated power values. }
  
  \item{eff.size}{A \code{length(RR)} x \code{mc.error} x \code{length(npergroup)} array of estimated treatment effects (differences between groups).}

}

\author{Rebecca Mangan and Jason Wyse}


\examples{
# simulate placebo group
placebo <- vitd.curve( n = 10, type = "placebo" )

# simulate treatment group
treatment <- vitd.curve( n = 10, type = "dynamic-dose" )

# intensity function for exposures to infection
intensfun <- intensity.function( summer.rate = 0, winter.rate = 1 )
  
# calculate power: example only- run for much larger num.sims
pow <- power.calc( n = c(10,20,30), 
        N = 10,  test.type = 'count', 
        vitdcurves.placebo = placebo, vitdcurves.treatment = treatment,
        baseline = 0.03, RR = c(2,3,4), rate = 1, intensity.func = intensfun )
            
### NOT RUN ###            
# approximate the Monte Carlo error in estimation of the power-- takes longer to run
#pow <- power.calc( n = c(10,20,30), 
#         N = 100, test.type = 'count', 
#         vitdcurves.placebo = placebo, vitdcurves.treatment = treatment,
#         baseline = 0.03, RR = c(2,3,4), rate = 1, 
#         intensity.func = intensfun, mc.error = 10 )

# plot power curves
plot( pow, xlab = "n", x.legend = 10, y.legend = 1, main.legend = "Relative Risk" )
}
