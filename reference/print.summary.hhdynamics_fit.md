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
s <- summary(fit)
print(s, digits = 4)
#>                                                    Variable Point estimate
#>               Daily probability of infection from community         0.0042
#>  Probability of person-to-person transmission in households         0.0576
#>                                                      sex1.0        -0.0825
#>                                                      age1.0        -0.0826
#>                                                      age2.0        -0.3285
#>  Lower bound Upper bound exp(Point estimate) exp(Lower bound) exp(Upper bound)
#>       0.0022      0.0071                  NA               NA               NA
#>       0.0364      0.0852                  NA               NA               NA
#>      -0.7131      0.4985              0.9208           0.4901           1.6463
#>      -0.5566      0.3981              0.9207           0.5732           1.4890
#>      -0.8159      0.1572              0.7200           0.4422           1.1702
# }
```
