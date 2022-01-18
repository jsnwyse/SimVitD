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
  
  for( dataset in 1:num.sims )
  {
    
    placebo <- vitd.curve( n = n, type = pl$type, start = start, end = end, mu = pl$mu, amplitude = pl$amplitude,
                           dyn.dose.thresh = pl$dyn.dose.thresh, sd.mu = pl$sd.mu, sd.amplitude = pl$sd.amplitude,
                           sd.dyn.dose.thresh = pl$sd.dyn.dose.thresh, supp.dose = pl$supp.dose, supp.dose.rate = pl$supp.dose.rate,
                           weight = pl$weight, sd.weight = pl$sd.weight, min.thresh = pl$min.thresh, north.hemi = pl$north.hemi)
    
    placebo.levels <- exposure.levels( placebo, rate, intensity.func = intensity.func, start, end )
    
    x.pl <- infection.count( placebo.levels, baseline, rel.risk, holding.time, lohi.vit=lohi.vit )$count
    
    treat <- vitd.curve( n = n, type = tr$type, start = start, end = end, mu = tr$mu, amplitude = tr$amplitude,
                         dyn.dose.thresh = tr$dyn.dose.thresh, sd.mu = tr$sd.mu, sd.amplitude = tr$sd.amplitude,
                         sd.dyn.dose.thresh = tr$sd.dyn.dose.thresh, supp.dose = tr$supp.dose, supp.dose.rate = tr$supp.dose.rate,
                         weight = tr$weight, sd.weight = tr$sd.weight, min.thresh = tr$min.thresh, north.hemi = tr$north.hemi)
    
    treat.levels <- exposure.levels( treat, rate,intensity.func = intensity.func, start, end )
    
    x.tr <- infection.count( treat.levels, baseline, rel.risk, holding.time, lohi.vit=lohi.vit )$count
    
    if( test.type == 'proportions' )
    {
      x.pl <- as.integer( x.pl > 0)
      x.tr <- as.integer( x.tr > 0)
    }
    
    boot.test <- boot.two.bca( x.pl, x.tr, stacked=FALSE, mean, null.hyp=0, alternative="greater", R=boot.rep )
    
    p[dataset] <- boot.test$p.value
    eff.size[dataset] <- mean(x.pl) - mean(x.tr)
  }
  
  power <- sum( p < sig.level ) / num.sims
  
  return( list( power=power, eff.size= mean(eff.size)) )
  
}
