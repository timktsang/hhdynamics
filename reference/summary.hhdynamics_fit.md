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
#>  Probability of person-to-person transmission in households          0.060
#>                                                      sex1.0         -0.140
#>                                                      age1.0         -0.091
#>                                                      age2.0         -0.353
#>  Lower bound Upper bound exp(Point estimate) exp(Lower bound) exp(Upper bound)
#>        0.002       0.007                  NA               NA               NA
#>        0.037       0.092                  NA               NA               NA
#>       -0.638       0.324               0.869            0.528            1.382
#>       -0.581       0.322               0.913            0.560            1.380
#>       -0.850       0.075               0.702            0.428            1.078
# }
```
