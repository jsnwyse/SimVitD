power.calc <- function( n, ratio=1, N = 500, test.type, sig.level = 0.05,
                        vitdcurves.placebo = NULL, vitdcurves.treatment = NULL,
                        baseline = 0.03, RR = 3,
                        rate = 1, intensity.func = intensity.function(), 
                        holding.time = 2, lohi.vit = c(10,70), clt = NULL, mc.error = 1, 
                        boot.rep = 9999, parallel = FALSE, num.cores = NULL, verbose=FALSE ){
  
  placebo.group <- vitdcurves.placebo
  treatment.group <- vitdcurves.treatment
  if( ratio < 1 ) stop("Argument ratio should not be less than 1.")
  if( sig.level < 0 | sig.level > 1 ) stop("Argument 'sig.level' is a significance level.")
  if( baseline < 0 ) stop("Argument 'baseline' is a probability and must be positive.")
  if( any(RR < 0 ) ) stop("Some entries of 'RR' are negative.")
  if( rate < 0 ) stop("Argument 'rate' must be positive.")
  if( holding.time < 0 ) stop("Argument 'holding.time' must be positive." )
  if( mc.error < 0 | mc.error != floor(mc.error) ) stop("Argument 'mc.error' must be a positive integer.")
  if( any( baseline*RR > 1 ) ) stop("Argument 'RR' scales 'baseline' to probabilities greater than 1.")
  if( !( test.type %in% c('proportions','count') )) stop("Argument 'test.type' must be one of 'proportions','count'.")
  if( (!is.null(placebo.group) & class(placebo.group) != 'vitd.curve') | (!is.null(treatment.group) & class(treatment.group)!= 'vitd.curve' ) ) stop("Arguments 'vitdcurves.placebo' and 'vitdcurves.treatment' must be of class 'vitd.curve'.")
  if( mc.error == 1 & parallel ){ message("Setting parallel equal to true will only impact scenarios where mc.error > 1."); parallel <- FALSE }
  
  if( verbose & parallel ) message("Printing progress updates is not possible with parallel set to TRUE.")
  
  # put in a check here for sizes of n and print warning for clt approximation
  clt <- if( !is.null(clt) )
  {
    if( length(clt) > 1 & length(clt) < length(n) ) stop("If argument clt is a vector, it must have the same length as n.")
    if( length(clt) == 1 )
    {
      if( clt==FALSE & any(n >= 35) ) message("Some n values seem large: setting clt to TRUE may give a faster approximation than bootstrapping for these values.")
      if( clt==TRUE & any(n < 35) ) message("Some n values seem small: a clt approximation may not give reliable results")
      rep(clt,length(n))
    }
  }else{
    # determine automatically
    message("Power approximation (CLT or bootstrap) has been determined for each instance in n: these can be set using argument clt")
    ap <- rep( TRUE, length(n) )
    ap[ n < 35 ] <- FALSE
    ap
  }
  
  num.participants <- n
  num.sims <- N
  
  if( parallel & !requireNamespace("parallel", quietly=TRUE) )  stop("Package \"parallel\" needed for this parallel functionality to work. Please install it.", call. = FALSE)
  
  if( parallel & !is.null(num.cores) )
  {
    # make sure parallelisation is sensible...
    if( num.cores >  mc.error ) num.cores <- mc.error
  }
  
  if( parallel & mc.error > 1 & is.null(num.cores) ) num.cores <- parallel::detectCores() - 1
  
  # make the clusters for parallel simulation
  if( parallel )
  {
    cl <- parallel::makePSOCKcluster( num.cores )
    parallel::clusterSetRNGStream( cl )
  }
  
  curve.type <- c(placebo.group$type, treatment.group$type)
  
  pow.sim <- array( dim=c( length(RR), mc.error, length(num.participants) ) )
  eff.size <- array( dim=c( length(RR), mc.error, length(num.participants) ) )
  group1.mean <- array( dim=c( length(RR), mc.error, length(num.participants) ) )
  group2.mean <- array( dim=c( length(RR), mc.error, length(num.participants) ) )
  
  for( iteration in 1:length(RR))
  {
    rel.risk <- RR[ iteration ]
    if( parallel )
    {
      powpara <- parallel::parLapplyLB( cl, 1:mc.error, power.calc.x, num.participants, ratio, num.sims, test.type, sig.level,
                                        placebo.group, treatment.group,
                                        baseline, rel.risk,
                                        rate, intensity.func, 
                                        holding.time, lohi.vit, boot.rep, verbose, clt )
      for( r in seq_along(powpara) )
      {
        pow.sim[ iteration, r, ] <- powpara[[ r ]]$power
        eff.size[ iteration, r, ] <- powpara[[ r ]]$eff.size
        group1.mean[ iteration, r, ] <- powpara[[ r ]]$group.1
        group2.mean[ iteration, r, ] <- powpara[[ r ]]$group.2
      }
    }else{
      for( r in 1:mc.error )
      {
        # need to incorporate the ratio into here...
        powpara <- power.calc.x( 0, num.participants, ratio, num.sims, test.type, sig.level,
                                 placebo.group, treatment.group,
                                 baseline, rel.risk,
                                 rate, intensity.func, 
                                 holding.time, lohi.vit, boot.rep, verbose, clt )
        pow.sim[ iteration, r,  ] <- powpara$power
        eff.size[ iteration, r, ] <- powpara$eff.size
        group1.mean[ iteration, r, ] <- powpara$group.1
        group2.mean[ iteration, r, ] <- powpara$group.2
      }
    }
  }
  
  if( verbose & !parallel ) cat("All runs complete.")
  
  if( parallel ) parallel::stopCluster(cl)
  
  
  y <- list( curve.type, test.type, baseline, RR, num.participants, mc.error, pow.sim, eff.size, group1.mean, group2.mean, ratio )
  names( y ) <- c( "curve.type", "test.type", "baseline", "RR", "npergroup", "mc.error", "power", "eff.size", "group1.mean", "group2.mean", "ratio" )
  class( y ) <- "power.calc"
  return( y )
  
}
