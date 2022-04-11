intensity.function <- function( summer.rate = 0, winter.rate = 1, flu = TRUE ){

  if( summer.rate < 0 || summer.rate > 1 ) stop( "summer.rate must be value between 0 and 1" )
  if( winter.rate < 0 || winter.rate > 1 ) stop( "winter.rate must be value between 0 and 1" )

  if( flu ){
    
    intensity <- function(t)
    {
      t <- 6 * ( t - floor(t/pi) * pi )
      l <- (t <= pi ) * ( t >= 4*pi)
      t[l] <- winter.rate
      t[!l] <- summer.rate
      return(t)
    }
    
  }else{
    
    intensity <- function(t)
    {
      return( rep( winter.rate, length(t) ) )
    }
  }
  
  
  return( intensity )
}
