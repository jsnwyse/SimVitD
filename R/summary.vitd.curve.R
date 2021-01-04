summary.vitd.curve <- function( object, ... ){ # have to put in if loop if treatment / placebo
  
  x <- object
  
  type <- x$type # test for type of vitamin D curve 3 - placebo/trad 4 - treatment

  cat("\tLength of study:",(length( x$time )-1)/(12*x$res),"years")
  cat("\n\tNumber of participants:",nrow(x$curves$curves.eval))
  cat("\n\tExpected 25OHD level:",round(x$mu),"nmol/L")
  
  if( type == "placebo" | type ==  "fixed-dose"){
    if( type == "fixed-dose" ) cat("\n\tSupplement 25OHD contribution:",round(x$supp.dose),"nmol/L")
  }else{
    cat("\n\tThreshold 25OHD level:",round(x$dyn.dose.thresh),"nmol/L")
  }
}


