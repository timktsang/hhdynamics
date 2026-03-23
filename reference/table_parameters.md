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
#> 1 Daily P(community infection) 0.004286945 0.004123752 0.002127015 0.006686184
#> 2    P(household transmission) 0.054171342 0.053289169 0.033269222 0.078908730
#> 3                       sex1.0 0.992834907 0.956390609 0.480493147 1.761999638
#> 4                       age1.0 0.988432797 0.948722639 0.544689260 1.601053411
#> 5                       age2.0 0.809427995 0.754455129 0.437896738 1.595635701
#>   Acceptance
#> 1  0.2822823
#> 2  0.3283283
#> 3  0.7197197
#> 4  0.7067067
#> 5  0.6666667
# }
```
