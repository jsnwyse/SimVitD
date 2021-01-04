power.calc <- function( n, ratio=1, N = 500, test.type, sig.level = 0.05,
                        vitdcurves.placebo = NULL, vitdcurves.treatment = NULL,
                        baseline = 0.03, RR = 3,
                        rate = 1, intensity.func = intensity.function(), 
                        holding.time = 2, lohi.vit = c(10,70), mc.error = 1, boot.rep = 500, parallel = FALSE, num.cores = NULL, verbose=FALSE ){
  
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
  if( mc.error == 1 & parallel ) warning("Setting parallel equal to true will only impact scenarios where mc.error > 1.")
  
  if( verbose & parallel ) cat("Printing progress updates is not possible with parallel set to TRUE.")
  
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

  for( iter in 1:length(num.participants) )
  {
    ni <- num.participants[iter]
    
    if( verbose & !parallel ) cat("Beginning run for",ni,"participants: ")
    
    for( iteration in 1:length(RR) )
    {
      
      rel.risk <- RR[ iteration ]
      
      if( verbose & !parallel ) cat("\n\t... with RR =", rel.risk )
      
      if( parallel ){
        # parallel implementation 
        
        powpara <- parallel::parLapplyLB( cl, 1:mc.error, power.calc.0, ni, ratio, num.sims, test.type, sig.level,
                                          placebo.group, treatment.group,
                                          baseline, rel.risk,
                                          rate, intensity.func, 
                                          holding.time, lohi.vit, boot.rep )
        
        for( r in seq_along(powpara) ) pow.sim[ iteration, r, iter ] <- powpara[[ r ]]$power
        for( r in seq_along(powpara) ) eff.size[ iteration, r, iter ] <- powpara[[ r ]]$eff.size
        
      }else{
        # not parallel- execute in sequence
        for( r in 1:mc.error )
        {
          # need to incorporate the ratio into here...
          powpara <- power.calc.0( 0, ni, ratio, num.sims, test.type, sig.level,
                                   placebo.group, treatment.group,
                                   baseline, rel.risk,
                                   rate, intensity.func, 
                                   holding.time, lohi.vit, boot.rep )
          pow.sim[ iteration, r, iter ] <- powpara$power
          eff.size[ iteration, r, iter ] <- powpara$eff.size
        }
        
      }
      
      
    }
    if( verbose & !parallel ) cat("\n")
  }
  
  if( verbose & !parallel ) cat("All runs complete.")
  
  if( parallel ) parallel::stopCluster(cl)
  
  
  y <- list( curve.type, test.type, baseline, RR, num.participants, mc.error, pow.sim, eff.size, ratio )
  names( y ) <- c( "curve.type", "test.type", "baseline", "RR", "npergroup", "mc.error", "power", "eff.size", "ratio" )
  class( y ) <- "power.calc"
  return( y )
  
}
