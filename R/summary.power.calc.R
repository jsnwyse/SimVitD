summary.power.calc <- function( object, ... ){
  x <- object
  cat( "\tGroups tested are", x$curve.type[1],"and",x$curve.type[2])
  cat( "\n\tTest type is", x$test.type )
  cat( "\n\tBaseline prevalence is", x$baseline )
  cat( "\n\tRelative Risk is", x$RR )
}

