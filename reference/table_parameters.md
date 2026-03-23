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
#> The running time is 41 seconds
table_parameters(fit)
#>                      Parameter        Mean      Median       Lower       Upper
#> 1 Daily P(community infection) 0.004297105 0.004180552 0.002157186 0.007241747
#> 2    P(household transmission) 0.057799646 0.056802601 0.034247973 0.085981164
#> 3                       sex1.0 0.959288634 0.916211880 0.492720579 1.670063592
#> 4                       age1.0 0.944429198 0.917964994 0.565029244 1.481706049
#> 5                       age2.0 0.743496508 0.724802870 0.442918118 1.183207883
#>   Acceptance
#> 1  0.4902994
#> 2  0.4898993
#> 3  0.5807721
#> 4  0.5812387
#> 5  0.5842389
# }
```
