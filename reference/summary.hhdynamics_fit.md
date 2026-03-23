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
#>  Probability of person-to-person transmission in households          0.059
#>                                                      sex1.0         -0.040
#>                                                      age1.0         -0.127
#>                                                      age2.0         -0.396
#>  Lower bound Upper bound exp(Point estimate) exp(Lower bound) exp(Upper bound)
#>        0.002       0.007                  NA               NA               NA
#>        0.040       0.082                  NA               NA               NA
#>       -0.673       0.512               0.960            0.510            1.669
#>       -0.655       0.382               0.881            0.519            1.465
#>       -0.835       0.199               0.673            0.434            1.220
# }
```
