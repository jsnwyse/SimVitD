vitd.curve <- function( n = 1, type = c("placebo","fixed-dose","dynamic-dose"), start = 0, end = 1,  
                        mu = 45, amplitude = 35, dyn.dose.thresh = 50, sd.mu = 5, sd.amplitude = 5, sd.dyn.dose.thresh = 5, 
                        supp.dose = 20, supp.dose.rate = Inf, weight = 1, sd.weight = 0, min.thresh = 10, north.hemi = TRUE, res = 40  )
{
  
  type <- type[1]
  if( !( type %in% c("placebo","fixed-dose","dynamic-dose") ) ) stop("Argument 'type' is not valid")

  if( any( c(start < 0, end < 0, end <= start ) ) ) stop("Arguments 'start' and 'end' must be positive with 'end' > 'start'")
  if( any( c( mu < 10, mu  < amplitude ) ) ) stop("Please check arguments mu and amplitude are sensible")
  if( any( c(sd.mu < 0, sd.amplitude < 0, sd.dyn.dose.thresh < 0 ) ) ) stop( "Arguments 'sd' must be positive " ) 
  if( supp.dose < 0 ) stop("Argument 'supp.dose' must be positive")
  
  cross <- .5*(start + end) # not used
  
  years <- end - start
  time1 <- seq( start*pi, end*pi, length.out=years*(12*res) + 1 ) # pi is one year
  cross1 <- cross*pi
  time <- ( time1 * (12/pi) ) # gives time in terms of months
  
  if( type == "placebo" ) 
  {
    supp.dose <- 0
    delta <- rep(supp.dose,n)
    shape1 <- 1 # irrelevant- does not contribute here
    shape2 <- 1
    sd.dyn.dose.thresh <- 0
  }else if( type == "fixed-dose"){
    delta <- rtrexp( n, supp.dose.rate, supp.dose)   # allow for variable uptake in supplementation
    #delta <- rep(supp.dose,n) 
    sd.weight <- ( (weight == 0 || weight == 1) & sd.weight == 0 ) * 0 + ( (weight > 0 & weight < 1) & sd.weight == 0 ) *  1e-4
    var.weight <- sd.weight^2
    c0 <- weight*(1-weight)/var.weight - 1
    shape1 <- weight*c0
    shape2 <- (1-weight)*c0
    if( sd.weight == 0 )
    {
      if( weight == 1 ){ shape1 <- Inf; shape2 <- 1 }
      if( weight == 0 ){ shape1 <- 1; shape2 <- Inf }
    }
    sd.dyn.dose.thresh <- 0
  }else if( type == "dynamic-dose"){
    supp.dose <- 0
    delta <- rep(supp.dose,n)
    shape1 <- 1 # irrelevant- does not contribute here
    shape2 <- 1
  }
  
  if( sd.amplitude == 0 ) sd.amplitude <- 1e-4
  
  y <- list( time, random.cosine( n , k.A = sd.amplitude, k.H = sd.mu, k.T = sd.dyn.dose.thresh, H.0 = mu,
                                  A.0 = amplitude, T.0 = dyn.dose.thresh, time = time1, delta=delta, tau=min.thresh, shape1=shape1, shape2=shape2 ), 
             mu, amplitude, dyn.dose.thresh,
             sd.mu, sd.amplitude, sd.dyn.dose.thresh, type, supp.dose, supp.dose.rate, weight, sd.weight, min.thresh, north.hemi, res, start, end, res )
  names( y ) <- c( "time", "curves", "mu", "amplitude", "dyn.dose.thresh", "sd.mu", "sd.amplitude", "sd.dyn.dose.thresh",
                   "type", "supp.dose", "supp.dose.rate", "weight", "sd.weight", "min.thresh", "north.hemi", "res", "start", "end", "res" )
  
  class( y ) <- "vitd.curve"
  return( y )
  
}

