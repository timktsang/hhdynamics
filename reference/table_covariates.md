# Table of covariate effects (odds ratios)

Returns a data frame of covariate effects on infectivity and
susceptibility, with exponentiated estimates interpretable as relative
risks.

## Usage

``` r
table_covariates(fit, probs = c(0.025, 0.975))
```

## Arguments

- fit:

  An object of class `hhdynamics_fit`.

- probs:

  Numeric vector of length 2 for credible interval bounds.

## Value

A data frame with columns: Covariate, Type, Estimate, Lower, Upper,
exp_Estimate, exp_Lower, exp_Upper. Returns an empty data frame if no
covariates were used.

## Examples

``` r
# \donttest{
data(inputdata)
fit <- household_dynamics(inputdata, ~sex, ~age,
  n_iteration = 15000, burnin = 5000, thinning = 1)
#> Iteration: 1000
#> Iteration: 2000
#> Iteration: 3000
#> Iteration: 4000
#> Iteration: 5000
#> Iteration: 6000
#> Iteration: 7000
#> Iteration: 8000
#> Iteration: 9000
#> Iteration: 10000
#> Iteration: 11000
#> Iteration: 12000
#> Iteration: 13000
#> Iteration: 14000
#> The running time is 41 seconds
table_covariates(fit)
#>   Covariate           Type    Estimate      Lower     Upper exp_Estimate
#> 1    sex1.0    Infectivity -0.07614257 -0.6968124 0.5288521    0.9266841
#> 2    age1.0 Susceptibility -0.07725659 -0.5576413 0.3855814    0.9256523
#> 3    age2.0 Susceptibility -0.33159280 -0.8451054 0.1434486    0.7177795
#>   exp_Lower exp_Upper
#> 1 0.4981708  1.696983
#> 2 0.5725580  1.470469
#> 3 0.4295121  1.154248
# }
```
