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
#> 1    sex1.0    Infectivity -0.09794125 -0.7581836 0.5257743    0.9067022
#> 2    age1.0 Susceptibility -0.08529700 -0.5542926 0.4009306    0.9182395
#> 3    age2.0 Susceptibility -0.34225861 -0.8363658 0.1881606    0.7101645
#>   exp_Lower exp_Upper
#> 1 0.4685167  1.691768
#> 2 0.5744785  1.493214
#> 3 0.4332823  1.207027
# }
```
