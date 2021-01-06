print.power.calc <- function( x, ... ){
  if( x$mc.error == 1 )
  {
    if( length(x$npergroup) > 1 & length(x$RR) > 1 )
    {
      A <- x$power[,1,]
      colnames( A ) <- paste0("n=",x$npergroup)
      rownames( A ) <- paste0("RR=",x$RR)
    }else{
      A <- x$power[,1,]
      if( length(x$npergroup) > 1 ) 
        names(A) <- paste0("n=",x$npergroup)
      else
        names(A) <- paste0("RR=",x$RR)
    }
  }else{
    pr <- vector( length(x$RR), mode="list" )
    names(pr) <- paste0("RR=",x$RR)
    for( i in 1:length(x$RR) )
    {
      if( length(x$npergroup) > 1){
        pr[[i]] <- apply( x$power[i,,], 2, summary ) 
        colnames( pr[[i]] ) <- paste0("n=",x$npergroup)
      }else{
        pr[[1]] <- x$power[i,,]
        names( pr[[i]] ) <- rep( paste0("n=",x$npergroup), length(pr[[i]]) )
      }
    }
  }
  cat( "\tPower" )
  if( x$mc.error == 1 )
  {
    if( length(x$npergroup) > 1 & length(x$RR) > 1 )
    {
      cat("\n")
      print( round(A, 3) )
    }else{
      if( length(x$npergroup) > 1 ) 
      {
        cat( paste0(" @ RR=",x$RR, "\n" ) )
        print( round(A,3) )
      }else{
        cat( paste0(" @ n=",x$npergroup, "\n")  )
        print( round(A,3) )
      }
    }
  }else{
    cat("\n\t")
    print( pr )
  }
  
  if( x$mc.error > 1 ) cat( "\n\tMonte Carlo reps: ", x$mc.error ) 
}
