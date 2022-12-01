## Assessing the adequacy of binomial models and their predicted uncertainty for the monitoring of the L50 from sampled fish
by Mainguy et al.

This repository contains code and data that can be used to reproduce the figures and simulation studies in the paper.

The file `confint_L.R` contains the function `confint_L`, which can be use to calculate confidence intervals of different types for the L50 (and any other percentile). Methods included are

- Delta
- Fieller
- Profile likelihood
- Parametric bootstrap
- Nonparametric bootstrap
- Monte Carlo
- Bayesian

For the resampling methods (parametric and nonparametric bootstrap and monte carlo), three types of intervals can be calculated, namely, equal-tailed intervals (ETI; or percentile), bias corrected accelerated intervals (BCa) or highest density intervals (HDI).

Usage of the `confint_L` function involves the following arguments:

- `object`: a `glm` object fitted using the `binomial` family
- `p`: the percentile to be used to calculate the L (defaults to 0.5, i.e. the L50)
- `cf`: a vector of length two indicating the position of the intercept and slope of the fitted ogive in the vector of model coefficients (defaults to 1 and 2)
- `level`: confidence level of the interval (defaults to 0.95)
- `nboot`: number of bootstrap/monte carlo samples (defualts to 10000)
- `method`: method to be used to compute the confidence interval. Options are `delta`, `fieller`, `proflik`, `parboot`, `nonparboot`, `montecarlo`, and `bayesian`
- `interval_type`: for the resampling-type intervals, options are `eti` (for equal-tailed interval; or percentile), `bca` (for bias-corrected and accelerated interval), `hdi` (for highest density interval), or `all` (for all previous three types in one function run)
- `...`: arguments passed to internal functions
