########################################################################
# x is output from exposure.levels function x[[2]] is matrix of vitamin D levels at time of exposures
infection.count <- function( expos, baseline = 0.03, RR = 3, holding.time = 2, lohi.vit = c(10,70) ){
  
  x <- expos
  
  if( class(x) != "exposure.levels" ) stop("Argument 'expos' not of class 'exposure.levels'")
  
  if( RR < 1 ) stop( "RR must be > 1" )
  if( length(lohi.vit)  != 2 ) stop( "lohi.vit must be of form (low,high)" ) 
  if( lohi.vit[1] > lohi.vit[2] ) stop( "lohi.vit must be of form (low,high)" )
  
  ht.rate <- ( 1/holding.time ) * ( 52/pi ) # holding time rate 

  lower <- baseline 
  upper <- baseline*RR # upper asymptote for generalised logisitic function
  
  if( upper > 1 ) stop("Argument 'RR' scales 'baseline' to probability greater than 1.")
  
  lo <- lohi.vit[1]
  hi <- lohi.vit[2]
  
  tau <- 0.045*(hi-lo) # tau to control inflection points
  
  slope <- log( (hi-lo-tau)^2/tau^2 ) / (hi-lo)
  intercept <- log( (hi-lo-tau)/tau ) - slope * hi
  
  ## where is lower defined in this?
  
  n <- nrow( x[[2]] ) # number of participants in the study
  num.exp <- ncol( x[[2]] ) # greatest number of possible exposures
  
  prob <- glf( x[[2]], lower, upper, intercept, slope )
  
  risk <- prob / lower
  
  u <- matrix( runif( n*num.exp, 0, 1 ), nrow=n, ncol=num.exp )
  
  disease <- u < prob # 1 if diseased, 0 if not diseased
  if(holding.time != 0) disease <- holding.time( disease, x[[1]], ht.rate ) # function that adds holding time if infection
  
  count <- rowSums( disease , na.rm = TRUE )
  mean <- mean( count )

  y <- list( baseline, RR, lohi.vit, prob, risk, disease, count, mean )
  names( y ) <- c("baseline", "RR", "inflection", "probs", "relativerisk", "infection", "count", "mean")
  class( y ) <- "infection.count"
  return( y )
}

