# R package `adaptBayes`

This package contains `R` functions implementing the adaptive priors described in Boonstra and Barbaro (2018). To install and load this package, run the following code in `R`

```r
library(devtools)
# may take some time:
install_github('umich-biostatistics/adaptBayes') 

library(adaptBayes)
```

See the companion repository at https://www.github.com/psboonstra/AdaptiveBayesianUpdates (specifically, the `exemplar.R` script) to see how these functions can be used

### Compilation warning

STAN is smart enough to recognize the need for the normalizing constant and so, upon initial installation of the package, the compiler will throw the following warning (once for each of the stan files):

```r
DIAGNOSTIC(S) FROM PARSER:
Warning (non-fatal):
Left-hand side of sampling statement (~) may contain a
non-linear transform of a parameter or local variable.
If it does, you need to include a target += statement with
the log absolute determinant of the Jacobian of the
transform.
Left-hand-side of sampling statement:
    normalized_beta ~ normal(...)
```

This warning can be safely ignored because we do, in fact, calculate the normalizing constant

#### Note, 10-Jul-2018:

After updating to version 3.5.0, <samp>R</samp> occasionally throws the following 'error':

`Error in x$.self$finalize() : attempt to apply non-function`

Error is used in quotes because it does not interrupt any processes and does not seem to affect any results. Searching online, this has been asked about by others and seems to be related to garbage collection:

http://discourse.mc-stan.org/t/very-mysterious-debug-error-when-running-rstanarm-rstan-chains-error-in-x-self-finalize-attempt-to-apply-non-function/4746


### Current Suggested Citation

Boonstra, Philip S. and Barbaro, Ryan P., "Incorporating Historical Models
with Adaptive Bayesian Updates" (2018) *Biostatistics* 
https://doi.org/10.1093/biostatistics/kxy053

<a href="https://biostats.bepress.com/umichbiostat/paper124">Authors' Copy </a>
