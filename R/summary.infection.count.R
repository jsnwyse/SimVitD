summary.infection.count <- function( object, ... ){
  x <- object
  cat( "\tBaseline prevalence is", x$baseline )
  cat( "\n\tRelative Risk is", x$RR )
  n <- length( x$count )
  if( n > 10 ) 
    cat( "\n\tNumber of infections for first 10 participants is: ", x$count[1:10] )
  else 
    cat( "\n\tNumber of infections per participant is: ", x$count )
  cat( "\n\tProportion of infected in group is ",  length( which(x$count >= 1) ) / length(x$count) )
  cat( "\n\tMean number of infections in group is ", sum(x$count) / length(x$count) )
}

