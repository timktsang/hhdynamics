# Summary method for hhdynamics_fit

Computes posterior summaries (mean, 2.5%, 97.5% credible intervals) for
all model parameters. For community and household parameters, reports
the daily probability (via 1-exp(-x) transform). For covariate effects,
additionally reports exponentiated estimates.

## Usage

``` r
# S3 method for class 'hhdynamics_fit'
summary(object, ...)
```

## Arguments

- object:

  An object of class `hhdynamics_fit`.

- ...:

  Additional arguments (unused).

## Value

An object of class `summary.hhdynamics_fit`, which is a data frame with
columns: Variable, Point.estimate, Lower.bound, Upper.bound, and
optionally exp columns for covariates.

## Examples

``` r
# \donttest{
data(inputdata)
fit <- household_dynamics(inputdata, ~sex, ~age,
  n_iteration = 1000, burnin = 500, thinning = 1)
#> The running time is 2 seconds
summary(fit)
#>                                                    Variable Point estimate
#>               Daily probability of infection from community          0.004
#>  Probability of person-to-person transmission in households          0.053
#>                                                      sex1.0         -0.046
#>                                                      age1.0         -0.023
#>                                                      age2.0         -0.275
#>  Lower bound Upper bound exp(Point estimate) exp(Lower bound) exp(Upper bound)
#>        0.002       0.007                  NA               NA               NA
#>        0.032       0.079                  NA               NA               NA
#>       -0.682       0.601               0.955            0.505            1.825
#>       -0.477       0.444               0.977            0.621            1.558
#>       -0.673       0.169               0.759            0.510            1.184
# }
```
