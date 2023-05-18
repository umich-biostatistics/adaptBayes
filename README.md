# R package `adaptBayes`

The adaptBayes R packages provides code to accompany the methodology developed
in Boonstra and Barbaro (2020). The sensible adaptive bayesian update is
implemented in `glm_sab()` and the naive adaptive bayesian update is implemented
in `glm_nab()`. A non-adaptive version of the sensible bayesian update is 
provided in `glm_sb()`. An extended version of the sensible and sensible adaptive
priors are provided in `glm_sb2()` and `glm_sab2()`. 

```r
library(devtools)
# may take some time:
install_github('umich-biostatistics/adaptBayes') 

library(adaptBayes)
```

A companion repository for this package exists at 
https://www.github.com/psboonstra/AdaptiveBayesianUpdates, 
which contains a vignette (`vignette.pdf`) on using the adaptive priors in this 
package as well as code for running the simulation studies in Boonstra and 
Barbaro (2020). 

### Current Suggested Citation

Boonstra, Philip S. and Barbaro, Ryan P., "Incorporating Historical Models
with Adaptive Bayesian Updates" (2020) *Biostatistics* 21, e47--e64
https://doi.org/10.1093/biostatistics/kxy053

Boonstra, Philip S. and Kleinsasser, Michael, "adaptBayes: R package for adaptive Bayesian updates" (2022) R package version 2.0.0.
https://github.com/umich-biostatistics/adaptBayes

<a href="https://biostats.bepress.com/umichbiostat/paper124">Authors' Copy </a>
