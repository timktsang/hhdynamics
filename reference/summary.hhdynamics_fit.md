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
summary(fit)
#>                                                    Variable Point estimate
#>               Daily probability of infection from community          0.004
#>  Probability of person-to-person transmission in households          0.058
#>                                                      sex1.0         -0.091
#>                                                      age1.0         -0.077
#>                                                      age2.0         -0.325
#>  Lower bound Upper bound exp(Point estimate) exp(Lower bound) exp(Upper bound)
#>        0.002       0.007                  NA               NA               NA
#>        0.034       0.087                  NA               NA               NA
#>       -0.732       0.512               0.913            0.481            1.668
#>       -0.566       0.390               0.926            0.568            1.477
#>       -0.802       0.169               0.723            0.448            1.184
# }
```
