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
#>   Covariate           Type     Estimate      Lower     Upper exp_Estimate
#> 1    sex1.0    Infectivity -0.116264228 -0.7299445 0.5710986    0.8902400
#> 2    age1.0 Susceptibility  0.002280579 -0.4057281 0.4345640    1.0022832
#> 3    age2.0 Susceptibility -0.257778187 -0.6069611 0.1736721    0.7727666
#>   exp_Lower exp_Upper
#> 1 0.4819357  1.770211
#> 2 0.6664914  1.544290
#> 3 0.5450046  1.189665
# }
```
