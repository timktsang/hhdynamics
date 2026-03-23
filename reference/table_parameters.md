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
#>                      Parameter       Mean      Median      Lower       Upper
#> 1 Daily P(community infection) 0.00421830 0.003939608 0.00224782 0.007248609
#> 2    P(household transmission) 0.06005274 0.059974402 0.03793963 0.083718831
#> 3                       sex1.0 0.96779049 0.941273491 0.46800034 1.625625900
#> 4                       age1.0 0.90526250 0.889619477 0.56989564 1.303183909
#> 5                       age2.0 0.69595512 0.680699868 0.41513144 1.057859868
#>   Acceptance
#> 1  0.2762763
#> 2  0.3393393
#> 3  0.6826827
#> 4  0.6926927
#> 5  0.6866867
# }
```
