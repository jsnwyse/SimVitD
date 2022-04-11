
boot.two.SimVitD <- function( x, y, B, alternative = c("two.sided", "less", "greater") )
{
  # implement for p-value only : altertative is population x in relation to y e.g. x less than y etc.
  type <- match( alternative, c("two.sided", "less", "greater") ) - 1
  pval <- 0
  r <- .C("TWO_SAMP_BOOT", as.double(x), as.double(y), 
          as.integer(length(x)), as.integer(length(y)),
          as.integer(B), as.integer(type), pval = as.double(pval), 
          PACKAGE = "SimVitD" )
  return(r$pval[1])
}

test.fn.call <- function(fn,par)
{
  .External2("testfunc", fn, par, PACKAGE = "SimVitD")
}