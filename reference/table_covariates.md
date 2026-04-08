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
#>   Covariate           Type    Estimate      Lower      Upper exp_Estimate
#> 1    sex1.0    Infectivity -0.09092722 -0.6445150 0.55148445    0.9130842
#> 2    age1.0 Susceptibility -0.10014816 -0.6377861 0.40254734    0.9047034
#> 3    age2.0 Susceptibility -0.35661037 -0.8775240 0.05750059    0.7000452
#>   exp_Lower exp_Upper
#> 1 0.5249171  1.735828
#> 2 0.5284611  1.495630
#> 3 0.4158112  1.059186
# }
```
