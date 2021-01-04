# expos output from exposure.levels
# infect output from infection.count
infection.count.plot <- function( expos, infect, pch = 20, cex = 1.5, col = "red" ){
  x <- expos
  y <- infect
  if( class(x) != "exposure.levels" ) stop("Argument x is not of class 'exposure.levels'")
  if( class(y) != "infection.count" ) stop("Argument y is not of class 'infection.count'") 
  n <- nrow(x[[1]]) # number of participants
  time <- x[[1]] #/ (12/pi)
  for( i in 1:n ){
    index <- which( y$infection[i,] == 1 ) # vector of infection times
    points( time[i,index], x[[2]][i,index], pch = pch, cex = cex, col = col )
  }
}

