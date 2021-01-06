exposure.levels <- function( x, rate, intensity.func = intensity.function(), end = 1 )
{
  if( class(x) != "vitd.curve" )
    stop("Argument 'x' is not of class 'vitd.curve'")
  
  n <- nrow(x[[2]][[1]]) # gives number of participants in the study
  L <- rate * ( 52 / pi ) # 52/pi is one exposure per week 
  start <- x$time[1] * (pi/12)
  if( end < ( x$time[1] / 12 ) ) stop( "End time must be greater than start time" )

  num.exposures <- ceiling( 1.5 * L * end * pi )  # simulate a conservative amount to buffer at end

  eventtimes <- matrix( rexp( n * num.exposures, rate=L  ), ncol=n )
  eventtimes <- apply( eventtimes, 2, cumsum )
  eventtimes[ eventtimes > end*pi ] <- NA
  
  prob.acc <- intensity.func( eventtimes )
  U <- runif( n * num.exposures )
  Z <- U < prob.acc
  eventtimes[ !Z ] <- NA
  
  eventtimes <- t(eventtimes)

  mu <- x$curves$mu
  amplitude <- x$curves$amplitude
  weights <- x$curves$weights
  
  if( x$type == "dynamic-dose" ) thresh <- x$curves$thresh else thresh <- rep(0,n)
  
  eventheights <- event.height.na( thresh, mu, amplitude, eventtimes, x$min.thresh, x$supp.dose, weights )

  y <- list( eventtimes, eventheights, x$type )
  names( y ) <- c("exposures", "levels", "type")

  class( y ) <- "exposure.levels"
  return( y )
}
