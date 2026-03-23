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
#>                      Parameter       Mean      Median       Lower       Upper
#> 1 Daily P(community infection) 0.00425829 0.004020923 0.002240888 0.007595111
#> 2    P(household transmission) 0.05608197 0.056129193 0.027726669 0.083782978
#> 3                       sex1.0 0.99584741 0.972346707 0.528895175 1.693404894
#> 4                       age1.0 0.96174618 0.917322140 0.627487307 1.510067852
#> 5                       age2.0 0.74201833 0.720782831 0.428461544 1.182334871
#>   Acceptance
#> 1  0.2372372
#> 2  0.3063063
#> 3  0.6846847
#> 4  0.6656657
#> 5  0.6386386
# }
```
