# SimVitD 2.0

- New C implementation of simulation included (faster)
- Bootstrap BCa method from wBoot translated to C code to be accessed directly during simulation
- Added new internal function `power.calc.x()` to access the C implementation

# SimVitD 1.0.2

- Added functions from wBoot by Neil A Weiss into package due to archival
- Due to wBoot incorporation dependency changed to simpleboot
- Included CLT approximation to be used when larger samples
- Determination of CLT/bootstrap is carried out automatically- clt argument to power.calc can be used to set
- Change default number of bootstrap re-samples to 9999

# SimVitD 1.0.1

## Minor changes

- Added extra argument `start` to `exposure.levels()` to pass the time unit in years of study start
- Corrected bug in time offset of exposure levels when start > 0
- Changed default plotting label for 25OHD from '25 Hydroxy Vitamin D' to '25-hydroxyvitamin D'
- Added new function 'vitd.curve.2.function' which will convert the output of vitd.curve into a list of callable functions (not available at user level)

# SimVitD 1.0.0 

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


