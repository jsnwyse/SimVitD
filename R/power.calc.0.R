power.calc.0 <- function( idx, n, ratio, N, test.type, sig.level,
                          pl, tr,
                          baseline, rel.risk,
                          rate, intensity.func, 
                          holding.time, lohi.vit, boot.rep ){
  
  num.sims <- N
  
  p <- numeric( num.sims )
  eff.size <- numeric( num.sims ) 
  
  lower <- baseline
  upper <- lower * rel.risk
  intercept <- -4
  slope <- 0.1
  
  start <- tr$time[1] / 12
  end <- tr$time[length(tr$time)] / 12
  cross  <- tr$cross/pi
  years <- end - start
  time <- seq( start * pi, end * pi, length.out = years * 24 + 1 )
  
  n.treat <- floor( ratio * n )
  
  x.tr <- numeric( n * num.sims ) # store all of the treatment counts/indicators
  x.pl <- numeric( n * num.sims )
  
  # how many blocks
  L <- min( 1e+05, n * num.sims )
  B <- n * num.sims / L
  
  # correction at end if needed
  corr <- (n * num.sims) %% L
  corr <- c( rep(0,B-1), corr )
  
  for( b in 1:B )
  {
    
    placebo <- vitd.curve( n = L+corr[b], type = pl$type, start = start, end = end, mu = pl$mu, amplitude = pl$amplitude,
                           dyn.dose.thresh = pl$dyn.dose.thresh, sd.mu = pl$sd.mu, sd.amplitude = pl$sd.amplitude,
                           sd.dyn.dose.thresh = pl$sd.dyn.dose.thresh, supp.dose = pl$supp.dose, supp.dose.rate = pl$supp.dose.rate,
                           weight = pl$weight, sd.weight = pl$sd.weight, min.thresh = pl$min.thresh, north.hemi = pl$north.hemi)
    
    placebo.levels <- exposure.levels( placebo, rate, intensity.func = intensity.func, start, end )
    
    rm(placebo)
    
    x.pl[((b-1)*L+1):(b*L+corr[b])] <- infection.count( placebo.levels, baseline, rel.risk, holding.time, lohi.vit=lohi.vit )$count
    
    rm(placebo.levels)
    
    treat <- vitd.curve( n = L+corr[b], type = tr$type, start = start, end = end, mu = tr$mu, amplitude = tr$amplitude,
                         dyn.dose.thresh = tr$dyn.dose.thresh, sd.mu = tr$sd.mu, sd.amplitude = tr$sd.amplitude,
                         sd.dyn.dose.thresh = tr$sd.dyn.dose.thresh, supp.dose = tr$supp.dose, supp.dose.rate = tr$supp.dose.rate,
                         weight = tr$weight, sd.weight = tr$sd.weight, min.thresh = tr$min.thresh, north.hemi = tr$north.hemi)
    
    treat.levels <- exposure.levels( treat, rate,intensity.func = intensity.func, start, end )
    
    rm(treat)
    
    x.tr[((b-1)*L+1):(b*L+corr[b])] <- infection.count( treat.levels, baseline, rel.risk, holding.time, lohi.vit=lohi.vit )$count
    
    rm(treat.levels)
    
  }
  
  if( test.type == 'proportions' )
  {
    x.pl <- as.integer( x.pl > 0)
    x.tr <- as.integer( x.tr > 0)
  }
  
  CLT.approx <- FALSE
  
  #boot.test <- boot.two.bca( x.pl, x.tr, stacked=FALSE, mean, null.hyp=0, alternative="greater", R=boot.rep )
  for( k in 1:num.sims )
  {
    idx <- ((k-1)*n + 1):(k*n)
    if( CLT.approx )
    {
      z <- (mean(x.pl[idx]) - mean(x.tr[idx])) / sqrt( (var(x.pl[idx]) + var(x.tr[idx])) / n )
      p[k] <- 1 - pnorm(z)
    } else p[k] <- boot.two.SimVitD( x.pl[idx], x.tr[idx], B=boot.rep, alternative="greater" ) #boot.test$p.value
    eff.size[k] <- mean(x.pl[idx]) - mean(x.tr[idx])
  }
  
  power <- sum( p < sig.level ) / num.sims
  
  return( list( power=power, eff.size= mean(eff.size)) )
  
}
