# Print method for summary.hhdynamics_fit

Print method for summary.hhdynamics_fit

## Usage

``` r
# S3 method for class 'summary.hhdynamics_fit'
print(x, digits = 3, ...)
```

## Arguments

- x:

  An object of class `summary.hhdynamics_fit`.

- digits:

  Number of significant digits for printing. Default is 3.

- ...:

  Additional arguments (unused).

## Value

Invisibly returns `x`.

## Examples

``` r
# \donttest{
data(inputdata)
fit <- household_dynamics(inputdata, ~sex, ~age,
  n_iteration = 1000, burnin = 500, thinning = 1)
#> The running time is 2 seconds
s <- summary(fit)
print(s, digits = 4)
#>                                                    Variable Point estimate
#>               Daily probability of infection from community         0.0045
#>  Probability of person-to-person transmission in households         0.0567
#>                                                      sex1.0        -0.0797
#>                                                      age1.0        -0.1081
#>                                                      age2.0        -0.3307
#>  Lower bound Upper bound exp(Point estimate) exp(Lower bound) exp(Upper bound)
#>       0.0023      0.0073                  NA               NA               NA
#>       0.0340      0.0776                  NA               NA               NA
#>      -0.7384      0.3912              0.9234           0.4779           1.4788
#>      -0.5352      0.2809              0.8975           0.5856           1.3243
#>      -0.7729      0.1274              0.7184           0.4617           1.1358
# }
```
