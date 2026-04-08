# Table of transmission parameter estimates

Returns a clean data frame of posterior parameter estimates with
credible intervals and effective sample sizes. Community and household
parameters are reported on the probability scale.

## Usage

``` r
table_parameters(fit, probs = c(0.025, 0.975), show_ess = FALSE)
```

## Arguments

- fit:

  An object of class `hhdynamics_fit`.

- probs:

  Numeric vector of length 2 for credible interval bounds. Default:
  `c(0.025, 0.975)` (95% CrI).

- show_ess:

  Logical. If `TRUE`, include an ESS (effective sample size) column.
  Default: `FALSE`.

## Value

A data frame with columns: Parameter, Mean, Median, Lower, Upper,
Acceptance (and ESS if `show_ess = TRUE`).

## Examples

``` r
# \donttest{
data(inputdata)
fit <- household_dynamics(inputdata, ~sex, ~age,
  n_iteration = 1000, burnin = 500, thinning = 1)
#> The running time is 2 seconds
table_parameters(fit)
#>                      Parameter        Mean      Median       Lower       Upper
#> 1 Daily P(community infection) 0.004151039 0.003977196 0.002269415 0.006883804
#> 2    P(household transmission) 0.058053480 0.057670222 0.038660396 0.084865105
#> 3                       sex1.0 0.915157650 0.891361769 0.440160108 1.519210989
#> 4                       age1.0 0.998497397 0.968224149 0.546231978 1.585079549
#> 5                       age2.0 0.734305166 0.720502991 0.472992437 1.048756551
#>   Acceptance
#> 1  0.2512513
#> 2  0.3663664
#> 3  0.7137137
#> 4  0.7107107
#> 5  0.6936937
# }
```
