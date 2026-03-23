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
#>  Probability of person-to-person transmission in households         0.0602
#>                                                      sex1.0        -0.1683
#>                                                      age1.0        -0.1072
#>                                                      age2.0        -0.3420
#>  Lower bound Upper bound exp(Point estimate) exp(Lower bound) exp(Upper bound)
#>       0.0023      0.0076                  NA               NA               NA
#>       0.0352      0.0906                  NA               NA               NA
#>      -0.8876      0.3880              0.8451           0.4116           1.4740
#>      -0.5981      0.3497              0.8984           0.5498           1.4186
#>      -0.7991      0.1765              0.7103           0.4498           1.1931
# }
```
