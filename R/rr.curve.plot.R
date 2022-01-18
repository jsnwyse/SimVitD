# expos is output from exposure levels
# infect is output from infection.count
rr.curve.plot <- function( expos, infect, idx=1, main = NULL, xlab = "25-hydroxyvitamin D",
                           ylab = "Risk scaling", col = "blue", ... ){
   
  x <- expos
  y <- infect
                          
  if( class(x) != "exposure.levels" ) stop("Argument 'expos' not of class 'exposure.levels'")
  if( class(y) != "infection.count" ) stop("Argument 'infect' not of class 'infection.count'")
  
  n <- nrow(x[[1]])
  
  if( length(idx) > 1 | floor(idx) != idx ) stop("Set idx to a single index from 1:n")
  
  
  lvl <- seq( 0, 120, 1 )
  
  lower <- y$baseline
  upper <- lower * y$RR
  
  lo <- y$inflection[1]
  hi <- y$inflection[2]
  
  tau <- 0.045*(hi-lo) # tau to control inflection points
  
  slope <- log( (hi-lo-tau)^2/tau^2 ) / (hi-lo)
  intercept <- log( (hi-lo-tau)/tau ) - slope * hi
  
  OR.curve <- glf( lvl, lower, upper, intercept, slope )
  
  plot( lvl, OR.curve, type = "l", xlab = xlab, ylab = ylab, axes = FALSE, ... )
  axis( 2, at = seq( lower, upper, length.out = 6 ), labels = seq( 1, y$RR, length.out = 6 ) )
  axis( 1, at = seq( 0, 120, 20 ), labels = seq( 0, 120, 20 ) )
  if( is.null(main) ) main <- paste0( expos$type, " group" )
  if( length(which) == 1 ) main <- paste0( main, ": participant ", idx)
  title( main = main, ... )
  
  for( i in idx ){
    points( x$levels[i,], y$probs[i,], col = col, ... )
  }
  
  for( i in idx ){
    index <- which( y$infect[i,] == 1 )
    points( x$levels[i,index], y$probs[i,index], pch = 20, col = "red", ... )
  }
  
}
