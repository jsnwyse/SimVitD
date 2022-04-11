power.calc.x <- function( idx, n, ratio, N, test.type, sig.level,
                          pl, tr,
                          baseline, rel.risk,
                          rate, intensity.func, 
                          holding.time, lohi.vit, boot.rep, verbose, clt ){
  
  # this uses .External to call C implementation
  lower <- baseline
  upper <- lower * rel.risk
  intercept <- -4
  slope <- 0.1
  glfcall <- function(x) glf( x, lower, upper, intercept, slope ) 
  glf.pars <- c(lower, upper, intercept, slope )
  
  start <- tr$time[1] / 12
  end <- tr$time[length(tr$time)] / 12
  
  n.treat <- floor( ratio * n )
  exper <- length( n )
  
  # placebo has no supp
  pars.pl <- c( pl$sd.amplitude, pl$sd.mu, pl$sd.dyn.dose.thresh,
                pl$mu, pl$amplitude, pl$dyn.dose.thresh,
                pl$supp.dose, pl$supp.dose.rate, 
                pl$min.thresh, 
                1.0, 1.0 ) # shapes set to 1.0 as not used
  
  if( tr$type == "fixed-dose")
  {
    sd.weight <- if( (tr$weight == 0 || tr$weight == 1) & tr$sd.weight == 0 )
      0 else if( (tr$weight > 0 & tr$weight < 1) & tr$sd.weight == 0 )  1e-4
    var.weight <- sd.weight^2
    c0 <- tr$weight*(1-tr$weight)/var.weight - 1
    shape1 <- tr$weight*c0
    shape2 <- (1-tr$weight)*c0
    if( sd.weight == 0 )
    {
      if( weight == 1 ){ shape1 <- Inf; shape2 <- 1 }
      if( weight == 0 ){ shape1 <- 1; shape2 <- Inf }
    }
  }else{ 
    shape1 <- 1 
    shape2 <- 1 
  } 

  pars.tr <- c( tr$sd.amplitude, tr$sd.mu, tr$sd.dyn.dose.thresh,
                tr$mu, tr$amplitude, tr$dyn.dose.thresh,
                tr$supp.dose, tr$supp.dose.rate, 
                tr$min.thresh, 
                shape1, shape2 ) 
  
  pars <- c( pars.pl, pars.tr )
  
  res <- .External2( "SimVitD_powersim", intensity.func, pars, tr$type, c(start,end), 
                       rate, 1, holding.time, exper, c(n, n.treat), N, test.type, 
                     sig.level, boot.rep, glf.pars, as.integer(clt), as.integer(verbose), PACKAGE = "SimVitD" )

  return( list( power=res$power, eff.size=res$effsize, group.1=res$group1, group.2=res$group2 ) )
  
}
