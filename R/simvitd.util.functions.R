#################################################################

# Functions to be used for Vitamin D Simulation

#################################################################
# Rate Function for nhpp that has zero in summer and 1 in winter
# allows for many years

intensity.extreme <- function(t)
{
  
  intens <- numeric(length(t))
  
  intens[ t < pi/6 ] <- 1
  
  a <- ( t - pi / 6 ) / ( pi / 12 )
  
  b <- floor(a) %% 12
  
  intens[ b < 4 & t >= pi/6 ] <- 0
  
  intens[ b >= 4 & t >= pi/6 ] <- 1
  
  return( intens )

}

#################################################################
# Function that returns n random cosine curves, amplitudes and heights
random.cosine <- function( n=1, k.H=1, k.A=1, k.T=1, H.0=45, A.0=35, T.0=50, time=NULL, cross=-1, delta=20, tau=10, shape1=1, shape2=1 )
{
  
 timemat <- matrix( rep(time, n), ncol=length(time),  byrow=TRUE )

 # convert k.A and A.0 into shape and rate
 A <- rgamma( n, shape = (A.0/k.A)^2, rate = A.0/k.A^2 )
 
 H <- rnorm( n, mean = H.0, sd = k.H )
 
 omega <- rbeta( n, shape1 = shape1, shape2 = shape2)
 
 thresh <- rgamma( n, shape=(T.0/k.T)^2, rate=T.0/k.T^2 )
 
 supp <- delta * ( omega  + 0.5 * (1-omega) * (1 + sin(2*timemat - 3*pi/2) ) )
 
 curves <- H + A * cos( 2 * timemat - pi ) + supp 
 
 U <- curves < tau
 
 curves <- (1-U) * curves + U * tau
 
 rm(U)
 
 Z <- curves - thresh
 
 Z <- Z >= 0
 
 curves <- Z * curves + ( 1 - Z ) * thresh 
 
 return( list( curves.eval=curves, mu=H, amplitude=A, weights = omega, thresh = thresh ) )  
 
}

#################################################################
# Input event times from nhpp and returns vitamin d levels
# for flat sine curves when there are na event times
event.height.na <- function( thresh, H, A, t, tau, dose, weights)
{
	# evaluate curves
  output <- H + A * cos( 2*t - pi ) + dose * ( weights + 0.5 * (1-weights) * ( 1 + sin( 2*t - 3*pi/2) ))  
  
  # threhold this at 10 minimum
  U <- output < tau
  
  # correct for threshold
  output <- (1-U) * output + U * tau
  
  # threshold
  Z <- output - thresh
  Z <- Z > 0

  output <- Z * output + (1 - Z) * thresh
  
  return( output )
}



#################################################################
# Generate generalised logisitic regression curve

glf <- function( x, l.asym, u.asym, intercept, slope )
{
  return( l.asym + (u.asym-l.asym)/(1 + exp(intercept+slope*x)) )
}

#################################################################
# Input vector with binary disease (0, 1), vector of event times,
# rate - ouput vector with binary disease including holding time

holding.time <- function( infect, t, rate)
{
  
  n <- nrow( infect )
  m <- ncol( infect )
  
  hold.times <- matrix( rexp( n*m, rate ), nrow=n, ncol=m )
  
  tphold <- t + hold.times
  if( n == 1 ){ 
    tphold <- matrix( c( t[,1], tphold[,1:m-1]), ncol=m, nrow=1 ) 
  }else{ 
    tphold <- cbind( t[,1], tphold[,1:(m-1)] ) ## the issue is here with scalar
  }
  
  tphold[ is.na(tphold) ] <- 0 # to make sure no times following NA are knocked out
  
  Z <- ( t < tphold ) * infect # gives TRUE where event is before next time
  
  infect <- infect - Z # knocks out bad times
  
  return( infect )

}

#################################################################
# Truncated exponential for fixed increase
# 
rtrexp <- function(n, lambda, ub)
{
  u <- runif(n)
  ub + log( exp(- lambda*ub) + (1-exp(- lambda* ub) ) * u ) / lambda
}

dtrexp <- function(x, lambda, ub)
{
  lambda * exp(-lambda*(ub-x)) / (1-exp(- lambda* ub) )
}

