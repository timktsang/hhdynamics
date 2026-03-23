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
#> 1    sex1.0    Infectivity -0.11676775 -0.7188051 0.4812972    0.8897918
#> 2    age1.0 Susceptibility -0.07181166 -0.4580998 0.3765685    0.9307062
#> 3    age2.0 Susceptibility -0.33587788 -0.7930456 0.1503764    0.7147104
#>   exp_Lower exp_Upper
#> 1 0.4873342  1.618172
#> 2 0.6324843  1.457275
#> 3 0.4524647  1.162272
# }
```
