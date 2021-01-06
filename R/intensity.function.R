intensity.function <- function( summer.rate = 0, winter.rate = 1, flu = TRUE ){

  if( summer.rate < 0 || summer.rate > 1 ) stop( "summer.rate must be value between 0 and 1" )
  if( winter.rate < 0 || winter.rate > 1 ) stop( "winter.rate must be value between 0 and 1" )

  if( flu ){
    
    intensity <- function(t)
    {
      intens <- numeric( length(t) )
      
      intens[ t <= pi/6 ] <- winter.rate
      
      a <- ( t - pi/6 )/( pi/12 )
      
      b <- floor(a) %% 12
      
      intens[ ( b < 4 ) & t > pi/6 ] <- summer.rate
      
      intens[ ( b >= 4 ) & t > pi/6  ] <- winter.rate
      
      return(intens)
      
    }
    
  }else{
    
    intensity <- function(t)
    {
      return( rep( winter.rate, length(t) ) )
    }
  }
  return( intensity )
}
