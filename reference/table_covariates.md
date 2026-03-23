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
  n_iteration = 1000, burnin = 500, thinning = 1)
#> The running time is 2 seconds
table_covariates(fit)
#>   Covariate           Type    Estimate      Lower     Upper exp_Estimate
#> 1    sex1.0    Infectivity -0.16113348 -0.7449109 0.3976943    0.8511785
#> 2    age1.0 Susceptibility -0.08072098 -0.5388523 0.3652897    0.9224510
#> 3    age2.0 Susceptibility -0.31688556 -0.7011758 0.1266222    0.7284141
#>   exp_Lower exp_Upper
#> 1 0.4747766  1.488389
#> 2 0.5834174  1.440931
#> 3 0.4960017  1.134988
# }
```
