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
#>               Daily probability of infection from community         0.0040
#>  Probability of person-to-person transmission in households         0.0554
#>                                                      sex1.0        -0.0473
#>                                                      age1.0        -0.0134
#>                                                      age2.0        -0.2793
#>  Lower bound Upper bound exp(Point estimate) exp(Lower bound) exp(Upper bound)
#>       0.0021      0.0068                  NA               NA               NA
#>       0.0299      0.0835                  NA               NA               NA
#>      -0.7286      0.7509              0.9538           0.4826           2.1190
#>      -0.4303      0.4099              0.9867           0.6503           1.5067
#>      -0.8305      0.1741              0.7563           0.4358           1.1902
# }
```
