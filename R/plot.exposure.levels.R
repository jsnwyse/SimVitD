plot.exposure.levels <- function( x, col = "blue", ... ){
  if( .Device == "null device") stop("plot.vitd.curve must be called before plot.exposure.levels")
  n <- nrow(x[[1]])
  time <- x[[1]]
  for( i in 1:n ) points( time[i,], x[[2]][i,], col=col, ... )
}

