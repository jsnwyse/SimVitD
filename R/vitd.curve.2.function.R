# create a vitd.function from parameters
vitd.fn <- function( offset=0, scale=1, mu, a)
{
  force(offset); force(scale); force(mu); force(a)
  function(t){offset + scale * pmax( mu + a * cos( 2*pi*t - pi), 10 ) }
}

vitd.curve.2.function <- function( x, offset=0, scale=1 )
{
  if( class(x) != "vitd.curve" ) stop("Not a vitd.curve object")
  n <- length(x$curves$mu)
  fn <- vector(n, mode="list") # list to store functions
  for( k in 1:n ) fn[[k]] <- vitd.fn( offset, scale, x$curves$mu[k], x$curves$amplitude[k] )
  return( fn )
}
