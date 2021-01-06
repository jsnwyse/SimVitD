# SimVitD 1.0 

## Major changes

- Generative model of vitamin D status curves from a squared sine to cosine with threshold.
- Rather than maximum and minimum levels of vitamin D a more flexible mean level is now used.
- Changes to arguments passed to `vitd.curves()`.
- Cross-over designs have been disabled in current version.
- Parallel computation is now available for Monte Carlo error quantification in `pow.calc()`.
- Addition of internal function `power.calc.0()`.
- `prop.test()` and `wilcox.test()` used in version 0.1.2 replaced by non-parametric bootstrap from `wBoot` package.
- Dependency on `poisson` package removed as it was possible to vectorize our own sampling.
- Updated print and summary function for `exposure.levels` and `infection.count` objects.
- Changes to package vignette to document user facing changes.


