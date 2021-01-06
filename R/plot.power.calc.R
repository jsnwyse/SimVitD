plot.power.calc <- function( x, col = "hotpink", lwd = 1.5, lty = 1,
                             ylab = NULL, x.legend = NULL, y.legend = NULL,
                             main.legend = "Risk scaling", legend.size = 1, target.power = NA, which=1L, ... ){
  
  if( x$ratio == 1 ) xlab <- "n" else xlab <- paste0("n [ n supplement is ",x$ratio,"*n]")
  mainti <- paste( x$curve.type[1], "vs", x$curve.type[2], ":", x$test.type)

  if( is.null(ylab) & which==1L ) ylab <- "power" else ylab <- "effect size" 
  
  if( is.null(x.legend) ) x.legend <- min(x$npergroup)
  if( is.null(y.legend) ) y.legend <- 1
  
  lnrr <- length( x$RR )
  
  summaryPOWER <- vector( lnrr, mode="list" )
  
  if( which == 1L ) tarr <- x$power else tarr <- x$eff.size
  
  for( i in 1:lnrr ) 
  {
    if( x$mc.error > 1 )
    {
      if( length(x$npergroup) > 1 )
      { 
        summaryPOWER[[i]] <- apply( tarr[i, , ], 2, summary )
      }else{ 
        summaryPOWER[[i]] <- matrix(nrow=6,ncol=length(x$npergroup))
        summaryPOWER[[i]][,1] <- as.vector(summary( tarr[i,,1]))
      }
    }else{
      summaryPOWER[[i]] <- matrix(nrow=6,ncol=length(x$npergroup))
      summaryPOWER[[i]][3,] <- tarr[i,1,]
    }
    
  }  
  
  if( lnrr > 1 & length( x$npergroup ) == 1 ){ 
    
    if( x$mc.error == 1 )
    {
      
      plot( x$RR, tarr[,1,1], type = "l", col = col, lwd = lwd, lty = lty,
            main = mainti, ylim = c(0,1), xlab = "Risk scaling", ylab = ylab )
      
    }else{
      
      # make a vector of the medians
      meds <- numeric(lnrr)
      for( i in 1:lnrr ) meds[i] <- summaryPOWER[[i]][3,1]
      
      plot( x$RR, meds, type = "l", col = col, lwd = lwd, lty = lty,
            main = mainti, ylim = c(0,1), xlab = "Risk scaling", ylab = ylab, xlim=c(min(x$RR)-.5, max(x$RR)+.5) )
      
      for( i in 1:lnrr ) boxplot( tarr[i, ,1], add=TRUE, at=x$RR[i], xaxt="n", yaxt="n", boxfill="white", boxcol=col, outline=FALSE, medcol=col, whisklty="solid", whisklwd=1.5, whiskcol=col, staplecol=col )
      
    }
    
  }else if( lnrr > 0 & length( x$npergroup ) > 1 ){
    
    colour <- rainbow( lnrr )
    
    plot( x$npergroup, summaryPOWER[[1]][3,], type = "l", col = col, lwd = lwd, lty = lty,
          main = mainti, ylim = c(0,1), xlab = xlab, ylab = ylab, xaxt="n", xlim=c( min(x$npergroup)-1, max(x$npergroup)+1 ) )
    axis( 1, at=x$npergroup, labels=x$npergroup)
    
    if( x$mc.error > 1 ) boxplot( x$power[1, ,], use.cols=TRUE, add=TRUE, at=x$npergroup, xaxt="n", yaxt="n", boxfill="white", boxcol=col, outline=FALSE, medcol=col, whisklty="solid", whisklwd=1.5, whiskcol=col, staplecol=col )
    
    if( lnrr > 1 )
    {
      for( i in 2:lnrr ){
        lines( x$npergroup, summaryPOWER[[i]][3,], col = colour[i], lwd = lwd, lty = lty )
        if( x$mc.error > 1 ) boxplot( x$power[i, ,], use.cols=TRUE, add=TRUE, at=x$npergroup, xaxt="n", yaxt="n", boxfill="white", boxcol=colour[i], outline=FALSE, medcol=colour[i], whisklty="solid", whisklwd=1.5, whiskcol=colour[i], staplecol=colour[i] )
      }
    }
    
    if( !is.na(target.power) ) abline(h = target.power, lty="dotted", lwd=1.5, col="ForestGreen")
    legend(x.legend, y.legend, legend = x$RR, col = c(col, colour[2:lnrr]), title = main.legend, lwd = lwd, cex = legend.size )
    
  }else if( length( x$RR ) == 1 & length( x$npergroup ) == 1 ){
    
    if( x$mc.error == 1 )
    {
      stop("Only one relative risk and n value- plotting is not sensible.")
    }else{
      boxplot( tarr[1, ,], use.cols=TRUE, add=FALSE, boxfill="white", boxcol=col, outline=FALSE, medcol=col, whisklty="solid", whisklwd=1.5, whiskcol=col, staplecol=col, xlab=xlab, ylab=ylab, xaxt="n", main = mainti )
      axis(1, at=1, labels=x$npergroup )
    }
  }
}
