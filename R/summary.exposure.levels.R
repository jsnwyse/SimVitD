summary.exposure.levels <- function( object, ... ){
  x <- object
  # number of exposures per participants
  U <- apply( x$exposures, 1, function(z) !is.na(z) )
  l <- colSums(U)
  if( length(l) <= 10 )
    cat("\tNumber of exposures per participant is: ",l )
  else
    cat("\tNumber of exposures for first 10 participants is: ",l[1:10])
}

