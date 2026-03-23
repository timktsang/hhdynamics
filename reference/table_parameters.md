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
#> The running time is 42 seconds
table_parameters(fit)
#>                      Parameter        Mean      Median       Lower      Upper
#> 1 Daily P(community infection) 0.004206964 0.004104985 0.002169363 0.00688956
#> 2    P(household transmission) 0.058300558 0.057707305 0.035072077 0.08707642
#> 3                       sex1.0 0.946666185 0.906852737 0.499818086 1.66506970
#> 4                       age1.0 0.950137619 0.928354839 0.572672622 1.45903547
#> 5                       age2.0 0.749904453 0.725732014 0.437018622 1.19733720
#>   Acceptance
#> 1  0.5491033
#> 2  0.4890993
#> 3  0.5835056
#> 4  0.5960397
#> 5  0.5738383
# }
```
