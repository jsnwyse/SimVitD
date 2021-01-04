# vitdcurves is output from vitd.curve
# expos is output from exposure levels
# infect is output from infection.count

rr.profile.plot <- function( x, expos, infect, idx = 1, ... )
{
  
  vitdcurves <- x
  
  if( class(vitdcurves) != "vitd.curve" ) stop("Argument 'vitdcurves' not of class 'vitd.curve'")
  if( class(expos) != "exposure.levels" ) stop("Argument 'expos' not of class 'exposure.levels'")
  if( class(infect) != "infection.count") stop("Argument 'infect' not of class 'infection.count'")
  if( length(idx) > 1 ) warning("Length of idx is greater than 1: only first element used")
  idx <- idx[1]
  
  par( mfrow = c(1,2) )
  m <- 12*vitdcurves$res + 1 # take one yr
  time <- ( vitdcurves$time[1:m] ) / ( 12 / pi )
  h <- max( vitdcurves$curve$curves.eval[idx,] ) + 10
  
  # plot vit d profile for single participant
  plot( time, vitdcurves$curve$curves.eval[ idx, 1:m] , type = "l", xlab = "Time",
        ylab = "25 Hydroxy Vitamin D", ylim = c(0, h), axes = FALSE )
  
  axis( 2, ... )
  months <- c( month.abb[3:12], month.abb[1:3] )
  axis( 1, at = seq( time[1], time[length(time)], length.out = 13 ), labels = months, ... )
  
  year <- max( which(expos$exposures[ idx, ]  < pi) )
  points( expos$exposures[idx, 1:year], expos$levels[idx,1:year], col = "blue", cex=1.3, ...  )
  
  index <- which( infect$infection[idx,1:year] == 1 )
  points( expos$exposures[idx, index], expos$levels[idx,index], pch = 20, cex = 1.5, col = "red" )
  
# plot relative risk curve
  lvl <- seq( 0, h - 10 , 1 )
  
  lower <- infect$baseline
  upper <- lower * infect$RR
  
  lo <- infect$inflection[1]
  hi <- infect$inflection[2]
  
  tau <- 0.045*(hi-lo) # tau to control inflection points
  
  slope <- log( (hi-lo-tau)^2/tau^2 ) / (hi-lo)
  intercept <- log( (hi-lo-tau)/tau ) - slope * hi
  
  OR.curve <- glf( lvl, lower, upper, intercept, slope )
  
  plot( OR.curve, lvl, type = "l", 
        xlab="Relative Risk", ylab="", ylim=c(0, h), axes = FALSE, ... )
  axis( 1, at = seq( lower, upper, length.out = 6 ), labels = seq( 1, infect$RR, length.out = 6 ) )
  
  points( infect$probs[idx,1:year], expos$levels[idx,1:year], cex = 1.3, col = "blue" )
  points( infect$probs[idx,index], expos$levels[idx,index], type = "p", pch = 20, cex = 1.5, col = "red" )
  
  par( mfrow = c(1,1) )
  
}
