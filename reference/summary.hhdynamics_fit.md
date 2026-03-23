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
#>               Daily probability of infection from community          0.005
#>  Probability of person-to-person transmission in households          0.056
#>                                                      sex1.0         -0.067
#>                                                      age1.0         -0.078
#>                                                      age2.0         -0.321
#>  Lower bound Upper bound exp(Point estimate) exp(Lower bound) exp(Upper bound)
#>        0.002       0.009                  NA               NA               NA
#>        0.032       0.082                  NA               NA               NA
#>       -0.705       0.427               0.935            0.494            1.533
#>       -0.577       0.452               0.925            0.562            1.571
#>       -0.950       0.294               0.726            0.387            1.342
# }
```
