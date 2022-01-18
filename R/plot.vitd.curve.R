plot.vitd.curve <- function( x, main = " ", xlab = " ", ylab = "25-hydroxyvitamin D", col=1:6, add = FALSE, ylim = NULL, ... ){
  time <- ( x[[1]] ) / ( 12 / pi )
  if( add ) 
  {
    matplot( time, t(x[[2]][[1]]), col=col, type="l", axes=FALSE, add=TRUE, ... )
    return()
  }
  if( is.null(ylim) ) ylim <- c(0,min(125,max(x[[2]][[1]]+10)))
  matplot( time, t(x[[2]][[1]]), type="l", col=col, main=main, xlab=xlab, ylab=ylab, axes=FALSE, ylim = ylim, ... )
  axis( side = 2,   ... )
  year <- ( length(x[[1]]) - 1 ) / (x$res * 12)
  
  if( year < 1 )
  {
    # this is the case of less than 1 year
    month.start <- floor( x$start * 12 )
    month.end <- floor( x$end * 12 )
    if( x$north.hemi ){
      if( 3 + month.end > 12 ){
        months <- c( month.abb[(3+month.start):12], month.abb[1:((3+month.end)%%12)] )
      }else{
        months <- month.abb[(3 + month.start):(3 + month.end)]
      }
    }else{
      if( 9 + month.start > 12 ){
        st <- (9 + month.start)%%12 
        months <- month.abb[st:((9+month.end)%%12)]
      }else{
        st <- 9 + month.start
        if( month.end > 4){ 
          months <- c( month.abb[st:12], month.abb[1:((9+month.end)%%12)] )
        }else{
          months <- c( month.abb[st:(9+month.end)])
        }
      }
    }
    axis( 1, at = seq(time[1], time[length(time)], length.out = floor(year*12) + 1 ),
          labels = months, ... )
  }else{
    if( x$north.hemi ){
      months <- c( month.abb[3:12], rep(month.abb, year - 1), month.abb[1:3] )
    }else{
      months <- c( month.abb[9:12], rep(month.abb, year - 1), month.abb[1:9] )
    }
    axis( 1, at = seq(time[1], time[length(time)], length.out = year*12 + 1 ),
          labels = months, ... )
  }
  
}

